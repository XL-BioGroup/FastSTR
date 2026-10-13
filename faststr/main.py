import argparse
import multiprocessing
import logging
import os
import time
import csv

from . import utils
from . import get_subread_trs
from . import make_consensus_TRs
from .parallel_windows import detect_windows
from .runtime_log import run_logging, log_input_file


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="FastSTR: A tool for identifying tandem repeats (STRs) from DNA sequences."
    )

    # --- 模式选择（必选，三选一） ---
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument("--strict", action="store_true", help="Strict mode (match=2, mismatch=7, gap_open=7, gap_extend=7)")
    mode_group.add_argument("--normal", action="store_true", help="Normal mode (match=2, mismatch=5, gap_open=7, gap_extend=3)")
    mode_group.add_argument("--loose", action="store_true", help="Loose mode (match=2, mismatch=3, gap_open=5, gap_extend=1)")

    # --- 模型选择（目前只有一个 default）---
    model_group = parser.add_mutually_exclusive_group(required=True)
    model_group.add_argument("--default", action="store_true", help="Default thresholds (p_indel=15%%, p_match=80%%, score=50); explicit values override them")

    # --- 输入序列（位置参数）---
    parser.add_argument("sequence", type=str, help="Path to the DNA sequence FASTA file")

    # --- 可选参数 ---
    parser.add_argument("-f", "--out_dir", type=str, default=".", help="Output directory (default: current directory)")
    parser.add_argument("-s", "--start", type=int, default=1, help="Start index (default: 1)")
    parser.add_argument("-e", "--end", type=int, default=0, help="End index (0 means full length)")
    parser.add_argument("-l", "--read_length", type=int, default=15000, help="Sub-read length (default: 15000)")
    parser.add_argument("-o", "--overlap", type=int, default=1000, help="Overlap length (default: 1000)")
    parser.add_argument("-p", "--processes", type=int, default=1, help="Number of CPU cores to use (default: 1)")
    parser.add_argument("-b", "--beta", type=float, default=0.045, help="Motif coverage threshold for alignment (default: 0.045)")
    parser.add_argument("--insertion-share", type=float, default=0.5, metavar="THETA",
                        help="Insertion share of total indels in final overlap gain (0-1; default: 0.5)")

    # Explicit values override the selected mode and --default model.
    override_group = parser.add_argument_group(
        "preset overrides", "Only specified values replace preset values; option order does not matter."
    )
    override_group.add_argument("--match", type=int, metavar="N", help="Match reward (preset: 2)")
    override_group.add_argument("--mismatch", type=int, metavar="N", help="Mismatch penalty magnitude")
    override_group.add_argument("--gap-open", "--gap_open", type=int, metavar="N", help="Gap-open penalty magnitude")
    override_group.add_argument("--gap-extend", "--gap_extend", type=int, metavar="N", help="Gap-extension penalty magnitude")
    override_group.add_argument("--p-indel", "--p_indel", type=float, metavar="PCT", help="Maximum indel percentage (0-100; preset: 15)")
    override_group.add_argument("--p-match", "--p_match", type=float, metavar="PCT", help="Minimum match percentage (0-100; preset: 80)")
    override_group.add_argument("--score", type=int, metavar="N", help="Minimum alignment score (preset: 50)")

    return parser.parse_args(argv)


def resolve_parameters(args):
    """Apply required presets, then explicit overrides, independently of CLI order."""
    modes = {
        'strict': (2, 7, 7, 7),
        'normal': (2, 5, 7, 3),
        'loose': (2, 3, 5, 1),
    }
    mode_name = next(name for name in modes if getattr(args, name))
    if not args.default:
        raise ValueError("--default is required")
    values = dict(zip(('match', 'mismatch', 'gap_open', 'gap_extend'), modes[mode_name]))
    values.update(p_indel=15.0, p_match=80.0, score=50)
    for name in values:
        override = getattr(args, name)
        if override is not None:
            values[name] = override
    for name in ('match', 'mismatch', 'gap_open', 'gap_extend'):
        if values[name] <= 0:
            raise ValueError(f"--{name.replace('_', '-')} must be greater than 0")
    for name in ('p_indel', 'p_match'):
        if not 0 <= values[name] <= 100:
            raise ValueError(f"--{name.replace('_', '-')} must be between 0 and 100 percent")
    if values['score'] < 0:
        raise ValueError("--score must be nonnegative")
    if not 0 <= args.insertion_share <= 1:
        raise ValueError("--insertion-share must be between 0 and 1")
    return mode_name, values



def Fast_TR(m, r, g, e, p_indel, p_match, score, input_path, out_path='', start_index=1,
            end_index=0, read_length=15000, overlap_length=1000, process=1, beta=0.045,
            insertion_share=0.5):
    if not 0 <= insertion_share <= 1:
        raise ValueError("insertion_share must be between 0 and 1")
    # 预处理
    current_dir = os.path.dirname(os.path.abspath(__file__))
    gene_name = os.path.splitext(os.path.basename(input_path))[0]
    parameter_name = str(m) + '_' + str(r) + '_' + str(g) + '_' + str(e) + '_' + str(p_indel) + '_' + str(
        p_match) + '_' + str(score)
    if insertion_share != 0.5:
        parameter_name += f'_theta_{insertion_share:g}'
    out_file_name = gene_name + '.' + parameter_name
    # logfile_path = os.path.join(current_dir, out_file_name + '.log')
    if out_path == '' or out_path == 'None':
        logfile_path = os.path.join(current_dir, out_file_name + '.log')
        detected_tr_detailed_path = os.path.join(current_dir, out_file_name + '_detailed.dat')
        detected_tr_aligned_path = os.path.join(current_dir, out_file_name + '_aligned.dat')
        detected_tr_summary_path = os.path.join(current_dir, out_file_name + '_summary.csv')
    else:
        logfile_path = os.path.join(out_path, out_file_name + '.log')
        detected_tr_detailed_path = os.path.join(out_path, out_file_name + '_detailed.dat')
        detected_tr_aligned_path = os.path.join(out_path, out_file_name + '_aligned.dat')
        detected_tr_summary_path = os.path.join(out_path, out_file_name + '_summary.csv')
    with run_logging(logfile_path, process) as run:
        r = -r
        g = -g
        e = -e
        p_indel = p_indel / 100
        p_match = p_match / 100
        logging.info(
            "Effective parameters: match=%s mismatch=-%s gap_open=-%s gap_extend=-%s "
            "p_indel=%s%% p_match=%s%% score=%s beta=%s start=%s end=%s "
            "read_length=%s overlap=%s processes=%s insertion_share=%s input=%s output=%s",
            m, -r, -g, -e, p_indel * 100, p_match * 100, score, beta,
            start_index, end_index, read_length, overlap_length, process, insertion_share,
            input_path, out_path,
        )

        # 获取基因序列
        run.stage('reading input')
        log_input_file(input_path)
        print('Start reading gene sequence')
        start_time_part = time.perf_counter()
        sequences_list = utils.read_fasta(input_path, start_index, end_index)
        end_time_part = time.perf_counter()
        run_time_part = end_time_part - start_time_part
        logging.info(f"read gene sequence running time: {run_time_part} seconds")

        # 划分各条染色体
        run.stage('dividing sequences')
        sequence_count = len(sequences_list)
        processed_bases = 0
        print('Start dividing gene sequences')
        num_subreads = []
        sub_reads_list = []
        sequences_name_list = []
        start_time_part = time.perf_counter()
        for seqs in sequences_list:
            processed_bases += len(seqs)
            sequences_name_list.append(seqs.description)
            sub_reads_list.extend(utils.make_sub_reads(seqs, read_length, overlap_length))
            num_subreads.append(len(sub_reads_list))
        end_time_part = time.perf_counter()
        run_time_part = end_time_part - start_time_part
        logging.info(f"divide gene sequence running time: {run_time_part} seconds")
        logging.info("Input workload: sequences=%s processed_bases=%s windows=%s (bases after start/end slicing; before window overlap)",
                     sequence_count, processed_bases, len(sub_reads_list))
        del sequences_list

        # 并行处理sub_reads
        run.stage('detecting window TRs')
        print('Start detecting TRs in parallel')
        start_time_part = time.perf_counter()
        Reads_TRs = detect_windows(sub_reads_list, (m, r, g, e, p_indel, p_match, score, beta), process)
        end_time_part = time.perf_counter()
        run_time_part = end_time_part - start_time_part
        logging.info(f"detect TRs in parallel running time: {run_time_part} seconds")

        # 合并各个sub_read上的trs,构建共识trs
        run.stage('cross-window candidate construction')
        print('Start constructing consensus TRs')
        start_time_part = time.perf_counter()
        substage_started = start_time_part
        before_handling_compatibility_TRs_dict = {}
        cross_sub_reads_TRs_list = []
        cro_sub_read_params = []
        num_cro_sub_read = [0]
        time_it = []
        for index, s_n in enumerate(sequences_name_list):
            before_handling_compatibility_TRs_dict[s_n] = []
            cross_sub_reads_TRs_list.clear()
            last_seq_subreads = 0 if index == 0 else num_subreads[index - 1]
            for pos, _ in enumerate(Reads_TRs[last_seq_subreads:num_subreads[index]]):
                sub_read = Reads_TRs[last_seq_subreads + pos]
                if pos + 1 < num_subreads[index] - last_seq_subreads:
                    sequence = sub_reads_list[last_seq_subreads + pos][:read_length - overlap_length] + sub_reads_list[
                        last_seq_subreads + pos + 1]
                    next_sub_read = Reads_TRs[last_seq_subreads + pos + 1]
                else:
                    sequence = sub_reads_list[num_subreads[index] - 1]
                    next_sub_read = []
                up_trs, mid_trs, if_change_motif = make_consensus_TRs.make_two_subreads_consensus(sub_read, next_sub_read,
                                                                                                  sequence, pos,
                                                                                                  read_length,
                                                                                                  overlap_length,
                                                                                                  start_index,
                                                                                                  cross_sub_reads_TRs_list,
                                                                                                  p_indel, p_match, m, r, g,
                                                                                                  e, score, beta)
                before_handling_compatibility_TRs_dict[s_n].extend(up_trs)
                if pos + 1 < num_subreads[index] - last_seq_subreads:
                    Reads_TRs[last_seq_subreads + pos + 1] = mid_trs

                if if_change_motif == 1:
                    cross_tr_seq = []
                    start_piece = (cross_sub_reads_TRs_list[-1][1] - start_index) // (read_length - overlap_length)
                    end_picee = (cross_sub_reads_TRs_list[-1][2] - start_index) // (read_length - overlap_length)
                    cross_tr_seq.append(str(
                        sub_reads_list[last_seq_subreads + start_piece][
                        cross_sub_reads_TRs_list[-1][1] - start_index - ((read_length - overlap_length)) * start_piece:]))
                    for i in range(start_piece + 1, end_picee):
                        cross_tr_seq.append(str(sub_reads_list[last_seq_subreads + i]))
                    cross_tr_seq.append(str(
                        sub_reads_list[last_seq_subreads + end_picee][
                        :cross_sub_reads_TRs_list[-1][2] - start_index - ((read_length - overlap_length)) * end_picee + 1]))
                    cross_tr_seq = ''.join(cross_tr_seq)
                    consensus_motif = utils.tri_gram_model(cross_tr_seq, len(cross_sub_reads_TRs_list[-1][0]))
                    cross_sub_reads_TRs_list[-1][0] = consensus_motif[0]


            if cross_sub_reads_TRs_list == []:
                num_cro_sub_read.append(len(cro_sub_read_params))
                continue
            for c_s_r in cross_sub_reads_TRs_list:
                cross_tr_seq = []
                start_piece = (c_s_r[1] - start_index) // (read_length - overlap_length)
                end_picee = (c_s_r[2] - start_index) // (read_length - overlap_length)
                cross_tr_seq.append(str(
                    sub_reads_list[last_seq_subreads + start_piece][
                    c_s_r[1] - start_index - ((read_length - overlap_length)) * start_piece:]))
                for i in range(start_piece + 1, end_picee):
                    cross_tr_seq.append(str(sub_reads_list[last_seq_subreads + i]))
                cross_tr_seq.append(str(
                    sub_reads_list[last_seq_subreads + end_picee][
                    :c_s_r[2] - start_index - ((read_length - overlap_length)) * end_picee + 1]))
                cross_tr_seq = ''.join(cross_tr_seq)
                cro_sub_read_params.append((c_s_r, cross_tr_seq, p_indel, p_match, m, r, g, e, score, beta))
            num_cro_sub_read.append(len(cro_sub_read_params))

        del Reads_TRs
        del sub_reads_list
        del cross_sub_reads_TRs_list

        logging.info("cross-window candidate construction running time: %.6f seconds; tasks=%s",
                     time.perf_counter() - substage_started, len(cro_sub_read_params))
        run.stage("cross-window realignment")
        substage_started = time.perf_counter()
        Cross_TRs = []
        if len(cro_sub_read_params) > 0:
            from .parallel_cross import calculate_cross_tasks
            Cross_TRs = calculate_cross_tasks(cro_sub_read_params, process)

        logging.info("cross-window realignment running time: %.6f seconds",
                     time.perf_counter() - substage_started)
        run.stage("final compatibility filtering and region merging")
        substage_started = time.perf_counter()
        del cro_sub_read_params
        Final_TRs_dict = {}
        Final_TRs_Region_dict = {}
        for index, s_n in enumerate(sequences_name_list):
            i = 0
            j = 0
            Final_TRs_dict[s_n] = []
            Final_TRs_Region_dict[s_n] = []
            merged_crosubtrs = []
            for sublist in Cross_TRs[num_cro_sub_read[index]:num_cro_sub_read[index + 1]]:
                merged_crosubtrs.extend(sublist)
            if len(before_handling_compatibility_TRs_dict[s_n]) > 0 and len(Cross_TRs) > 0:
                while i < len(before_handling_compatibility_TRs_dict[s_n]) and j < len(merged_crosubtrs):
                    if before_handling_compatibility_TRs_dict[s_n][i][1] <= merged_crosubtrs[j][1]:
                        Final_TRs_dict[s_n].append(before_handling_compatibility_TRs_dict[s_n][i])
                        i += 1
                    else:
                        Final_TRs_dict[s_n].append(merged_crosubtrs[j])
                        j += 1
            if before_handling_compatibility_TRs_dict[s_n] and i < len(before_handling_compatibility_TRs_dict[s_n]):
                Final_TRs_dict[s_n].extend(before_handling_compatibility_TRs_dict[s_n][i:])
            if Cross_TRs and j < len(merged_crosubtrs):
                Final_TRs_dict[s_n].extend(merged_crosubtrs[j:])

            Final_TRs_dict[s_n], Final_TRs_Region_dict[s_n] = make_consensus_TRs.handling_compatibility(Final_TRs_dict[s_n],
                                                                                                        p_match, p_indel, insertion_share)

        end_time_part = time.perf_counter()
        run_time_part = end_time_part - start_time_part
        logging.info("final compatibility filtering and region merging running time: %.6f seconds",
                     end_time_part - substage_started)
        logging.info(f"merge sub_reads and construct consensus TRs running time: {run_time_part} seconds")
        del before_handling_compatibility_TRs_dict

        # 写入识别结果
        run.stage('saving results')
        print('Start saving the detected TRs')
        start_time_part = time.perf_counter()
        '''
        写入摘要信息(dat)
        '''
        with open(detected_tr_summary_path, mode='w', newline='', encoding='utf-8', buffering=1048576) as file:
            writer = csv.writer(file)
            # 写入表头（手动创建）
            writer.writerow(
                ["Seq Name", "Start", "End", "Primary Motif", "STR Gain Score", "Second Motif", "STR Gain Score",
                 "Third Motif", "STR Gain Score"])
            for seq_name, final_trs_region in Final_TRs_Region_dict.items():
                for f_t_r in final_trs_region:
                    if f_t_r[0] == 1:
                        writer.writerow([seq_name, f_t_r[1], f_t_r[2], f_t_r[3], f_t_r[4]])
                    elif f_t_r[0] == 2:
                        writer.writerow([seq_name, f_t_r[1], f_t_r[2], f_t_r[3], f_t_r[4], f_t_r[5], f_t_r[6]])
                    else:
                        writer.writerow(
                            [seq_name, f_t_r[1], f_t_r[2], f_t_r[3], f_t_r[4], f_t_r[5], f_t_r[6], f_t_r[7], f_t_r[8]])

        '''
        写入检测报告（dat）
        '''
        with open(detected_tr_detailed_path, 'w', buffering=1048576) as file1, open(detected_tr_aligned_path, 'w', buffering=1048576) as file2:
            # 写入注释
            file1.write(
                'The report on FastSTR detected of STRs provides a detailed list of the distribution, quality, and structure of all STRs, with the following content template:\n\n\n')
            file1.write('**********************************************************************\n')
            file1.write('Gene sequence name\n')
            file1.write('----------------------------------------\n')
            file1.write(
                'start\t\tend\t\tregion length\t\tmotif length\t\tcopy number\t\tmotif\t\tindel percentage\t\tmatch percentage\t\talign score\t\talign uuid\n')
            file1.write('----------------------------------------\n')
            file1.write("The total number of detected STRs is: X\n")
            file1.write('**********************************************************************\n\n\n')
            # 写入tr
            for seq_name, final_trs in Final_TRs_dict.items():
                file1.write('**********************************************************************\n')
                file1.write(seq_name + '\n')
                file1.write('----------------------------------------\n')
                s_n = seq_name.replace(' ', '')
                uuid_prefix = s_n[:min(len(s_n), 20)]
                for f_t in final_trs:
                    uuid_line = '%'.join([uuid_prefix, str(f_t[1]), str(f_t[2]), f_t[0]]) + '\n'
                    file1.write(
                        f"{str(f_t[1]):<12}\t{str(f_t[2]):<12}\t{str(f_t[2] - f_t[1]):<10}\t{str(len(f_t[0])):<5}"
                        f"\t{f'{f_t[7]:.2f}':<10}\t{f_t[0]:<12}\t{f'{f_t[3]:.4f}':<10}\t{f'{f_t[4]:.4f}':<10}"
                        f"\t{str(f_t[5]):<10}\t" + uuid_line
                    )
                    file2.write(uuid_line)
                    file2.write(f_t[6] + '\n\n')
                # 写入tr总个数
                if len(final_trs) > 0:
                    file1.write('----------------------------------------\n')
                file1.write("The total number of detected STRs is: " + str(len(final_trs)) + '\n')
                file1.write('**********************************************************************\n\n\n')
        end_time_part = time.perf_counter()
        run_time_part = end_time_part - start_time_part

        logging.info(f"save the detected TRs running time: {run_time_part} seconds")
        logging.info("Results saved: detailed_calls=%s summary_regions=%s",
                     sum(len(calls) for calls in Final_TRs_dict.values()),
                     sum(len(regions) for regions in Final_TRs_Region_dict.values()))
        logging.info("Output files: summary=%s detailed=%s aligned=%s",
                     os.path.abspath(detected_tr_summary_path),
                     os.path.abspath(detected_tr_detailed_path),
                     os.path.abspath(detected_tr_aligned_path))


def main():
    args = parse_args()
    try:
        mode_name, values = resolve_parameters(args)
    except ValueError as exc:
        raise SystemExit(f"faststr: error: {exc}") from None

    print(f"Running FastSTR in {mode_name.upper()} mode with DEFAULT model")
    print(f"Input sequence: {args.sequence}")
    print(f"Output directory: {args.out_dir}")
    print(
        "Effective parameters: "
        f"match={values['match']}, mismatch=-{values['mismatch']}, "
        f"gap_open=-{values['gap_open']}, gap_extend=-{values['gap_extend']}, "
        f"p_indel={values['p_indel']:g}%, p_match={values['p_match']:g}%, "
        f"score={values['score']:g}, beta={args.beta:g}, "
        f"start={args.start}, end={args.end}, read_length={args.read_length}, "
        f"overlap={args.overlap}, processes={args.processes}, "
        f"insertion_share={args.insertion_share:g}"
    )

    Fast_TR(
        values['match'],
        values['mismatch'],
        values['gap_open'],
        values['gap_extend'],
        values['p_indel'],
        values['p_match'],
        values['score'],
        args.sequence,
        args.out_dir,
        args.start,
        args.end,
        args.read_length,
        args.overlap,
        args.processes,
        args.beta,
        args.insertion_share,
    )


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
