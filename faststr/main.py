import argparse
import multiprocessing
import logging
import os
import time
import csv

from . import utils
from . import get_subread_trs
from . import make_consensus_TRs


def parse_args():
    parser = argparse.ArgumentParser(
        description="FastSTR: A tool for identifying tandem repeats (STRs) from DNA sequences."
    )

    # --- 模式选择（必选，三选一） ---
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument("--strict", action="store_true", help="Strict mode (match=2, mismatch=3, gap_open=5, gap_extend=1)")
    mode_group.add_argument("--normal", action="store_true", help="Normal mode (match=2, mismatch=5, gap_open=7, gap_extend=3)")
    mode_group.add_argument("--loose", action="store_true", help="Loose mode (match=2, mismatch=7, gap_open=7, gap_extend=7)")

    # --- 模型选择（目前只有一个 default）---
    model_group = parser.add_mutually_exclusive_group(required=True)
    model_group.add_argument("--default", action="store_true", help="Use the default STR identification model")

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

    return parser.parse_args()



def Fast_TR(m, r, g, e, p_indel, p_match, score, input_path, out_path='', start_index=1,
            end_index=0, read_length=15000, overlap_length=1000, process=1, beta=0.045):
    # 预处理
    current_dir = os.path.dirname(os.path.abspath(__file__))
    gene_name = os.path.splitext(os.path.basename(input_path))[0]
    parameter_name = str(m) + '_' + str(r) + '_' + str(g) + '_' + str(e) + '_' + str(p_indel) + '_' + str(
        p_match) + '_' + str(score)
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
    logging.basicConfig(filename=logfile_path, level=logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', filemode='w')
    r = -r
    g = -g
    e = -e
    p_indel = p_indel / 100
    p_match = p_match / 100
    logging.info("Program started")

    # 获取基因序列
    print('Start reading gene sequence')
    start_time_part = time.time()
    sequences_list = utils.read_fasta(input_path, start_index, end_index)
    end_time_part = time.time()
    run_time_part = end_time_part - start_time_part
    logging.info(f"read gene sequence running time: {run_time_part} seconds")

    # 划分各条染色体
    print('Start dividing gene sequences')
    num_subreads = []
    sub_reads_list = []
    sequences_name_list = []
    start_time_part = time.time()
    for seqs in sequences_list:
        sequences_name_list.append(seqs.description)
        sub_reads_list.extend(utils.make_sub_reads(seqs, read_length, overlap_length))
        num_subreads.append(len(sub_reads_list))
    end_time_part = time.time()
    run_time_part = end_time_part - start_time_part
    logging.info(f"divide gene sequence running time: {run_time_part} seconds")
    del sequences_list

    # 并行处理sub_reads
    print('Start detecting TRs in parallel')
    start_time_part = time.time()
    params = []
    for sub_read in sub_reads_list:
        params.append((sub_read, m, r, g, e, p_indel, p_match, score, beta))
    pool = multiprocessing.Pool(process)
    Reads_TRs = pool.starmap(get_subread_trs.get_subread_trs, params)
    pool.close()
    pool.join()
    end_time_part = time.time()
    run_time_part = end_time_part - start_time_part
    logging.info(f"detect TRs in parallel running time: {run_time_part} seconds")
    del params

    # 合并各个sub_read上的trs,构建共识trs
    print('Start constructing consensus TRs')
    start_time_part = time.time()
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

    Cross_TRs = []
    if len(cro_sub_read_params) > 0:
        pool = multiprocessing.Pool(process)
        Cross_TRs = pool.starmap(make_consensus_TRs.calculate_cross_subread_tr, cro_sub_read_params)
        pool.close()
        pool.join()

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
                                                                                                    p_match, p_indel)

    end_time_part = time.time()
    run_time_part = end_time_part - start_time_part
    logging.info(f"merge sub_reads and construct consensus TRs running time: {run_time_part} seconds")
    del before_handling_compatibility_TRs_dict

    # 写入识别结果
    print('Start saving the detected TRs')
    start_time_part = time.time()
    '''
    写入摘要信息(dat)
    '''
    with open(detected_tr_summary_path, mode='w', newline='', encoding='utf-8') as file:
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
    with open(detected_tr_detailed_path, 'w') as file1, open(detected_tr_aligned_path, 'w') as file2:
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
            for f_t in final_trs:
                file1.write(
                    f"{str(f_t[1]):<12}\t{str(f_t[2]):<12}\t{str(f_t[2] - f_t[1]):<10}\t{str(len(f_t[0])):<5}"
                    f"\t{f'{f_t[7]:.2f}':<10}\t{f_t[0]:<12}\t{f'{f_t[3]:.4f}':<10}\t{f'{f_t[4]:.4f}':<10}"
                    f"\t{str(f_t[5]):<10}\t"
                )
                s_n = seq_name.replace(' ', '')
                file1.write('%'.join([s_n[:min(len(s_n), 20)], str(f_t[1]), str(f_t[2]), f_t[0]]) + '\n')
                file2.write('%'.join([s_n[:min(len(s_n), 20)], str(f_t[1]), str(f_t[2]), f_t[0]]) + '\n')
                file2.write(f_t[6] + '\n\n')
            # 写入tr总个数
            if len(final_trs) > 0:
                file1.write('----------------------------------------\n')
            file1.write("The total number of detected STRs is: " + str(len(final_trs)) + '\n')
            file1.write('**********************************************************************\n\n\n')
    end_time_part = time.time()
    run_time_part = end_time_part - start_time_part

    logging.info(f"save the detected TRs running time: {run_time_part} seconds")
    logging.info("Program completed")


def main():
    args = parse_args()

    # --- 模式选择 ---
    if args.strict:
        match, mismatch, gap_open, gap_extend = 2, 3, 5, 1
        mode_name = "strict"
    elif args.normal:
        match, mismatch, gap_open, gap_extend = 2, 5, 7, 3
        mode_name = "normal"
    elif args.loose:
        match, mismatch, gap_open, gap_extend = 2, 7, 7, 7
        mode_name = "loose"
    else:
        raise ValueError("You must specify one of --strict, --normal, or --loose")

    # --- 模型选择 ---
    if args.default:
        p_indel, p_match, score, quality_control = 0.15, 0.80, 50, False
    else:
        raise ValueError("You must specify one of --default")


    # --- 输出信息 ---
    print(f"Running FastSTR in {mode_name.upper()} mode with DEFAULT model")
    print(f"Input sequence: {args.sequence}")
    print(f"Output directory: {args.out_dir}")
    print(f"Using {args.processes} CPU core(s)")


    # --- 调用核心算法 ---
    Fast_TR(
        match,
        mismatch,
        gap_open,
        gap_extend,
        p_indel * 100,
        p_match * 100,
        score,
        args.sequence,
        args.out_dir,
        args.start,
        args.end,
        args.read_length,
        args.overlap,
        args.processes,
        args.beta
    )


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()