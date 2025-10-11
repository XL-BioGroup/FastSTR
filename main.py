import utils
import time
import get_subread_trs
import make_consensus_TRs
import logging
import multiprocessing
import argparse
import os
import csv


def parse_args():
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser(description='STR Identify Algorithm.')
    # 添加参数
    parser.add_argument('match', type=int, default=2, help='Match score')
    parser.add_argument('mismatch', type=int, default=5, help='Mismatch score')
    parser.add_argument('gap_open', type=int, default=7, help='Gap open score')
    parser.add_argument('gap_extend', type=int, default=3, help='Gap extend score')
    parser.add_argument('p_indel', type=int, default=15, help='Indel percentage threshold')
    parser.add_argument('p_match', type=int, default=80, help='Match percentage threshold')
    parser.add_argument('score', type=int, default=50, help='Alignment score threshold')
    parser.add_argument('DNA_file', type=str, help='The path of the DNA sequence file to be identified')
    parser.add_argument('-f', type=str, help='The directory path of the file to be output')
    parser.add_argument('-s', type=int, default=1, help='Start index')
    parser.add_argument('-e', type=int, default=0, help='End index')
    parser.add_argument('-l', type=int, default=15000, help='Sub_read length')
    parser.add_argument('-o', type=int, default=1000, help='Overlap length')
    parser.add_argument('-p', type=int, default=1, help='Number of CPU cores used')
    parser.add_argument('-b', type=float, default=0.045, help='Motif coverage of global sequence alignment algorithm')
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

        # t1 = time.time()
        # print(t1 - start_time_part)

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
    # print(end_time_part - t2)
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
                # file1.write('\t\t'.join(map(str, f_t[1:3])) + "\t\t" + str(f_t[2] - f_t[1]) + '\t\t' + str(len(
                #     f_t[0])) + '\t\t' + str(f"{f_t[7]:.2f}") + '\t\t' + f_t[0] + '\t\t' + str(
                #     f"{f_t[3]:.2f}") + '\t\t' + str(f"{f_t[4]:.2f}") + '\t\t' + str(f_t[5]) + '\t\t')
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
    Fast_TR(int(args.match), int(args.mismatch), int(args.gap_open), int(args.gap_extend), int(args.p_indel),
            int(args.p_match), int(args.score), str(args.DNA_file), str(args.f), int(args.s), int(args.e),
            int(args.l), int(args.o), int(args.p), float(args.b))
    # Fast_TR(2, 5, 7, 3, 15, 80, 50,
    #         r"/home/wenlong/project/20240330_TR_Detection/genome/maize/NC_050096.1_Zea_mays_cultivar_B73_chromosome_1_Zm-B73-REFERENCE-NAM-5.0_whole_genome_shotgun_sequence.fna",
    #         start_index=336000, end_index=351000, process=72,
    #         out_path='/home/wenlong/project/20240330_TR_Detection/FastSTR/temp_maize')


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
