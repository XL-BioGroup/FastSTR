import heapq
import regex

from . import trs_align
from . import utils


def judge_overlap_type(motif1, motif2, rear_overlap_pro, front_overlap_pro, start, end, overlap_length):
    if len(motif1) != len(motif2):
        return -1
    pattern = f"({motif1}){{s<=1}}"
    matches = regex.finditer(pattern, str(motif2 + motif2))
    if not any(matches):
        return -1
    if start <= overlap_length and end <= overlap_length and motif1 in motif2 + motif2 and rear_overlap_pro >= 0.9 and front_overlap_pro >= 0.9:
        return 0
    if start > overlap_length and end <= overlap_length and front_overlap_pro >= 0.8:
        return 1
    if start <= overlap_length and end > overlap_length and rear_overlap_pro >= 0.8:
        return 2
    if start > overlap_length and end > overlap_length:
        return 3
    return -1


def motif_sim(m1, m2):
    m1 = str(m1)
    m2 = str(m2)
    if len(m1) != len(m2):
        return -1
    if m1 in m2 + m2:
        return 0
    pattern = f"({m1}){{s<=1}}"
    matches = regex.finditer(pattern, str(m2 + m2))
    if any(matches):
        return 1
    return -1


def construct_consensus(rear_overlap, front_overlap, read_length, overlap_length, all_seq, global_index,
                        cross_sub_reads_TRs_list, p_indel, p_match, mat, mis, gap, ext, align_score, beta):
    consensus_trs = []
    j = 0
    global_index_p = global_index + read_length - overlap_length
    if_chage_motif = -1
    while rear_overlap and j < len(front_overlap) or rear_overlap == [] and j < len(front_overlap):
        # 下一个sub_read有全覆盖tr，可能是跨tr的起始或者中间
        if (front_overlap[j][2] - overlap_length) / (read_length - overlap_length) > 0.98:
            if_start = -1
            # 找出前一个重叠的tr，确定起始位置
            for index, r_o in enumerate(sorted(rear_overlap, key=lambda x: x[2], reverse=True)):
                if_chage_motif = motif_sim(front_overlap[j][0], r_o[0])
                if if_chage_motif >= 0:
                    tr_overlap_length = r_o[2] - max(front_overlap[j][1] + read_length - overlap_length, r_o[1]) + 1
                    rear_overlap_pro = tr_overlap_length / (overlap_length - front_overlap[j][1])
                    if rear_overlap_pro >= 0.8:
                        if_start = len(rear_overlap) - 1 - index
                        break
            # 如果找到，则该tr必为起始段
            if if_start > -1:
                cross_sub_reads_TRs_list.append(
                    [rear_overlap[if_start][0], rear_overlap[if_start][1] + global_index,
                     front_overlap[j][2] + global_index_p])
                del rear_overlap[if_start]
            # 如果没找到，则tr可能为起始段，可能为中间段，取决于跨sub_read字典是否存在前序tr重叠
            elif if_start == -1:
                # 如果该tr全覆盖sub_read，且与前序tr有重叠，motif相似，则为中间段
                if (front_overlap[j][2] - front_overlap[j][
                    1] + 1) / read_length >= 0.96 and cross_sub_reads_TRs_list and front_overlap[j][
                    1] + global_index_p < cross_sub_reads_TRs_list[-1][-1]:
                    if_chage_motif = motif_sim(front_overlap[j][0], cross_sub_reads_TRs_list[-1][0])
                    if if_chage_motif >= 0:
                        cross_sub_reads_TRs_list[-1][-1] = front_overlap[j][2] + global_index_p
                # 否则为起始段（无需全覆盖sub_read）
                else:
                    cross_sub_reads_TRs_list.append([front_overlap[j][0], front_overlap[j][1] + global_index_p,
                                                     front_overlap[j][2] + global_index_p])
            del front_overlap[j]
            j = 0
            continue
        # 到了跨tr的结束端
        if cross_sub_reads_TRs_list:
            if_chage_motif = motif_sim(front_overlap[j][0], cross_sub_reads_TRs_list[-1][0])
            if front_overlap[j][1] + global_index_p < cross_sub_reads_TRs_list[-1][2] and front_overlap[j][
                2] + global_index_p >= cross_sub_reads_TRs_list[-1][2] and if_chage_motif >= 0:
                cross_sub_reads_TRs_list[-1][2] = front_overlap[j][2] + global_index_p
                del front_overlap[j]
                j = 0
                continue
            else:
                if_chage_motif = -1
        if rear_overlap == []:
            break
        tr_overlap_length = min(rear_overlap[0][2], front_overlap[j][2] + read_length - overlap_length) - max(
            rear_overlap[0][1], front_overlap[j][1] + read_length - overlap_length) + 1
        rear_overlap_pro = tr_overlap_length / (rear_overlap[0][2] - rear_overlap[0][1] + 1)
        front_overlap_pro = tr_overlap_length / (front_overlap[j][2] - front_overlap[j][1] + 1)
        rear_tr = (rear_overlap[0][0], rear_overlap[0][1] + global_index, rear_overlap[0][2] + global_index,
                   rear_overlap[0][3], rear_overlap[0][4], rear_overlap[0][5], rear_overlap[0][6], rear_overlap[0][7])
        front_tr = (front_overlap[j][0], front_overlap[j][1] + global_index_p, front_overlap[j][2] + global_index_p,
                    front_overlap[j][3], front_overlap[j][4], front_overlap[j][5], front_overlap[j][6],
                    front_overlap[j][7])
        overlap_type = judge_overlap_type(rear_overlap[0][0], front_overlap[j][0], rear_overlap_pro, front_overlap_pro,
                                          read_length - rear_overlap[0][1], front_overlap[j][2] + 1, overlap_length)
        # 未重叠
        if front_overlap[j][1] + read_length - overlap_length >= rear_overlap[0][2] or j + 1 == len(front_overlap):
            consensus_trs.append(rear_tr)
            del rear_overlap[0]
            j = 0
            continue
        # 重叠但未匹配
        if overlap_type == -1:
            j += 1
            continue
        # 全重叠匹配
        if overlap_type == 0:
            if rear_tr[1] <= front_tr[1]:
                consensus_trs.append(rear_tr)
            else:
                consensus_trs.append(front_tr)
            del rear_overlap[0]
            del front_overlap[j]
            j = 0
            continue
        # 左重叠匹配
        if overlap_type == 1:
            consensus_trs.append(rear_tr)
            del rear_overlap[0]
            del front_overlap[j]
            j = 0
            continue
        # 右重叠匹配
        if overlap_type == 2:
            consensus_trs.append(front_tr)
            del rear_overlap[0]
            del front_overlap[j]
            j = 0
            continue
        # 嵌入重叠匹配
        if overlap_type == 3:
            to_scan_sequence = utils.get_realscan_sequence(
                all_seq[rear_overlap[0][1]:front_overlap[j][2] + 1 + read_length - overlap_length])
            motifs = utils.tri_gram_model(to_scan_sequence, len(rear_overlap[0][0]))
            if motifs == []:
                consensus_trs.append(rear_tr)
                consensus_trs.append(front_tr)
                del rear_overlap[0]
                del front_overlap[j]
                j = 0
                continue
            maybe_consensus = []
            for mot in motifs:
                maybe_consensus.extend(trs_align.trs_align_algorithm(
                    [(mot, rear_overlap[0][1], front_overlap[j][2] + 1 + read_length - overlap_length)], all_seq,
                    p_indel, p_match, mat, mis, gap, ext, align_score, beta))
            if maybe_consensus == []:
                consensus_trs.append(rear_tr)
                consensus_trs.append(front_tr)
                del rear_overlap[0]
                del front_overlap[j]
                j = 0
                continue
            for m_c in maybe_consensus:
                consensus_trs.append(
                    (m_c[0], m_c[1] + global_index, m_c[2] + global_index, m_c[3], m_c[4], m_c[5], m_c[6], m_c[7]))
            del rear_overlap[0]
            del front_overlap[j]
            j = 0
            continue

    # 上一个后重叠区还有tr，全部并入
    if rear_overlap:
        for reminder_tr in rear_overlap:
            rear_tr = (reminder_tr[0], reminder_tr[1] + global_index, reminder_tr[2] + global_index, reminder_tr[3],
                       reminder_tr[4], reminder_tr[5], reminder_tr[6], reminder_tr[7])
            consensus_trs.append(rear_tr)
    elif front_overlap:
        for reminder_tr in front_overlap:
            front_tr = (
                reminder_tr[0], reminder_tr[1] + global_index_p, reminder_tr[2] + global_index_p, reminder_tr[3],
                reminder_tr[4], reminder_tr[5], reminder_tr[6], reminder_tr[7])
            consensus_trs.append(front_tr)

    return sorted(consensus_trs, key=lambda x: x[1]), if_chage_motif


def make_two_subreads_consensus(up_read, mid_read, all_seq, subread_index, read_length, overlap_length, start_index,
                                cross_sub_reads_TRs_list, p_indel, p_match, mat, mis, gap, ext, align_score, beta):
    with_nosorted_trs = []
    to_be_merged_trs = []
    rear_repeating_area = read_length - overlap_length
    global_index = rear_repeating_area * subread_index + start_index

    if len(all_seq) <= read_length:
        if len(up_read) == 0:
            return [], [], -1
        for u_r in up_read:
            with_nosorted_trs.append(
                (u_r[0], u_r[1] + global_index, u_r[2] + global_index, u_r[3], u_r[4], u_r[5], u_r[6], u_r[7]))
        return sorted(with_nosorted_trs, key=lambda x: x[1]), [], -1
    if len(up_read) > 0:
        for u_r in up_read:
            if u_r[2] - rear_repeating_area < 24:
                with_nosorted_trs.append(
                    (u_r[0], u_r[1] + global_index, u_r[2] + global_index, u_r[3], u_r[4], u_r[5], u_r[6], u_r[7]))
            else:
                to_be_merged_trs.append(u_r)

    after_front_repeat = next((ind for ind, m_r in enumerate(mid_read) if overlap_length - m_r[1] <= 24), len(mid_read))
    if after_front_repeat == 0:
        for t_b_m in to_be_merged_trs:
            with_nosorted_trs.append(
                (t_b_m[0], t_b_m[1] + global_index, t_b_m[2] + global_index, t_b_m[3], t_b_m[4], t_b_m[5], t_b_m[6],
                 t_b_m[7]))
        return sorted(with_nosorted_trs, key=lambda x: x[1]), mid_read, -1
    after_consensus, if_chage_motif = construct_consensus(to_be_merged_trs, mid_read[:after_front_repeat], read_length,
                                                          overlap_length, all_seq, global_index,
                                                          cross_sub_reads_TRs_list, p_indel, p_match, mat, mis, gap,
                                                          ext, align_score, beta)
    with_nosorted_trs.extend(after_consensus)
    return sorted(with_nosorted_trs, key=lambda x: x[1]), mid_read[after_front_repeat:], if_chage_motif


def calculate_cross_subread_tr(cross_tr_tuple, cross_tr_seq, p_indel, p_match, mat, mis, gap, ext, align_score, beta):
    cross_tr = trs_align.trs_align_algorithm(
        [(cross_tr_tuple[0], cross_tr_tuple[1] - cross_tr_tuple[1], cross_tr_tuple[2] - cross_tr_tuple[1])],
        cross_tr_seq, p_indel, p_match, mat, mis, gap, ext, align_score, beta)
    if cross_tr == []:
        return []
    Cross_Tr = []
    for c_t in cross_tr:
        Cross_Tr.append(
            (c_t[0], c_t[1] + cross_tr_tuple[1], c_t[2] + cross_tr_tuple[1], c_t[3], c_t[4], c_t[5], c_t[6], c_t[7]))

    return Cross_Tr


def handling_compatibility(to_final_trs_list, p_match, p_indel):
    if len(to_final_trs_list) < 1:
        return to_final_trs_list, []
    if len(to_final_trs_list) < 2:
        return to_final_trs_list, [(1, to_final_trs_list[0][1], to_final_trs_list[0][2], to_final_trs_list[0][0],
                                    2 * (1 + to_final_trs_list[0][4]) * (1 - to_final_trs_list[0][3]) / (
                                            2 - to_final_trs_list[0][3]) - 1)]
    P_th_profit = 2 * (1 + p_match) * (1 - p_indel) / (2 - p_indel)
    region_minprop = 1 / P_th_profit
    final_trs_list = []
    can_retain = [1] * len(to_final_trs_list)
    cur_repeat_region = []
    cur_repeat_region_score = []
    repeat_region_mark = []
    max_repeat_region_end = 0
    for index, t_f in enumerate(to_final_trs_list):
        if can_retain[index] == 0:
            continue
        if index + 1 == len(to_final_trs_list):
            final_trs_list.append(t_f)
            if cur_repeat_region == []:
                cur_repeat_region.append(t_f)
                max_repeat_region_end = t_f[2]
                cur_repeat_region_score.append(2 * (1 + t_f[4]) * (1 - t_f[3]) / (2 - t_f[3]) - 1)
            else:
                one_region_legnth = max(max_repeat_region_end, t_f[2]) - cur_repeat_region[0][1]
                one_region_score = []
                for c_r_r in cur_repeat_region:
                    cache_score = 2 * (c_r_r[2] - c_r_r[1] + 1) * (1 + c_r_r[4]) * (1 - c_r_r[3]) / (
                            2 - c_r_r[3]) - one_region_legnth
                    if cache_score >= 0:
                        one_region_score.append(cache_score / one_region_legnth)
                    else:
                        top_three_tr = heapq.nlargest(3, zip(cur_repeat_region_score, cur_repeat_region))
                        if len(top_three_tr) == 1:
                            repeat_region_mark.append((1, cur_repeat_region[0][1], max_repeat_region_end,
                                                       cur_repeat_region[0][0], cur_repeat_region_score[0]))
                        elif len(top_three_tr) == 2:
                            repeat_region_mark.append((2, cur_repeat_region[0][1], max_repeat_region_end,
                                                       top_three_tr[0][1][0], top_three_tr[0][0], top_three_tr[1][1][0],
                                                       top_three_tr[1][0]))
                        else:
                            repeat_region_mark.append((3, cur_repeat_region[0][1], max_repeat_region_end,
                                                       top_three_tr[0][1][0], top_three_tr[0][0], top_three_tr[1][1][0],
                                                       top_three_tr[1][0], top_three_tr[2][1][0], top_three_tr[2][0]))
                        cur_repeat_region = [t_f]
                        max_repeat_region_end = t_f[2]
                        cur_repeat_region_score = [2 * (1 + t_f[4]) * (1 - t_f[3]) / (2 - t_f[3]) - 1]
                        break
                if len(one_region_score) == len(cur_repeat_region):
                    cur_repeat_region.append(t_f)
                    cur_repeat_region_score = one_region_score
                    cache_score = 2 * (t_f[2] - t_f[1] + 1) * (1 + t_f[4]) * (1 - t_f[3]) / (
                            2 - t_f[3]) - one_region_legnth
                    cur_repeat_region_score.append(cache_score / one_region_legnth)
                    max_repeat_region_end = max(max_repeat_region_end, t_f[2])
            continue
        for pos in range(index + 1, len(to_final_trs_list)):
            if to_final_trs_list[pos][1] >= t_f[2]:
                break
            region_length = max(t_f[2], to_final_trs_list[pos][2]) - min(t_f[1], to_final_trs_list[pos][1]) + 1
            tr_length_pro1 = (t_f[2] - t_f[1] + 1) / region_length
            tr_length_pro2 = (to_final_trs_list[pos][2] - to_final_trs_list[pos][1] + 1) / region_length
            if len(to_final_trs_list[pos][0]) == len(t_f[0]) and to_final_trs_list[pos][0] in t_f[0] + t_f[0]:
                if tr_length_pro1 < tr_length_pro2:
                    can_retain[pos] = 0
                else:
                    can_retain[index] = 0
                continue
            if tr_length_pro1 >= 0.8 and tr_length_pro2 <= region_minprop or tr_length_pro1 <= region_minprop and tr_length_pro2 >= 0.8:
                gain1 = 2 * (t_f[2] - t_f[1] + 1) * (1 + t_f[4]) * (1 - t_f[3]) / (2 - t_f[3]) - region_length
                gain2 = 2 * (to_final_trs_list[pos][2] - to_final_trs_list[pos][1] + 1) * (
                        1 + to_final_trs_list[pos][4]) * (1 - to_final_trs_list[pos][3]) / (
                                2 - to_final_trs_list[pos][3]) - region_length
                if gain2 < 0:
                    can_retain[pos] = 0
                if gain1 < 0:
                    can_retain[index] = 0
                    break
        if can_retain[index] == 1:
            final_trs_list.append(t_f)
            if cur_repeat_region == []:
                cur_repeat_region.append(t_f)
                max_repeat_region_end = t_f[2]
                cur_repeat_region_score.append(2 * (1 + t_f[4]) * (1 - t_f[3]) / (2 - t_f[3]) - 1)
            else:
                one_region_legnth = max(max_repeat_region_end, t_f[2]) - cur_repeat_region[0][1]
                one_region_score = []
                for c_r_r in cur_repeat_region + [t_f]:
                    cache_score = 2 * (c_r_r[2] - c_r_r[1] + 1) * (1 + c_r_r[4]) * (1 - c_r_r[3]) / (
                            2 - c_r_r[3]) - one_region_legnth
                    if cache_score >= 0:
                        one_region_score.append(cache_score / one_region_legnth)
                    else:
                        top_three_tr = heapq.nlargest(3, zip(cur_repeat_region_score, cur_repeat_region))
                        if len(top_three_tr) == 1:
                            repeat_region_mark.append((1, cur_repeat_region[0][1], max_repeat_region_end,
                                                       cur_repeat_region[0][0], cur_repeat_region_score[0]))
                        elif len(top_three_tr) == 2:
                            repeat_region_mark.append((2, cur_repeat_region[0][1], max_repeat_region_end,
                                                       top_three_tr[0][1][0], top_three_tr[0][0], top_three_tr[1][1][0],
                                                       top_three_tr[1][0]))
                        else:
                            repeat_region_mark.append((3, cur_repeat_region[0][1], max_repeat_region_end,
                                                       top_three_tr[0][1][0], top_three_tr[0][0], top_three_tr[1][1][0],
                                                       top_three_tr[1][0], top_three_tr[2][1][0], top_three_tr[2][0]))
                        cur_repeat_region = [t_f]
                        max_repeat_region_end = t_f[2]
                        cur_repeat_region_score = [2 * (1 + t_f[4]) * (1 - t_f[3]) / (2 - t_f[3]) - 1]
                        break
                if len(one_region_score) == len(cur_repeat_region) + 1:
                    cur_repeat_region.append(t_f)
                    cur_repeat_region_score = one_region_score
                    max_repeat_region_end = max(max_repeat_region_end, t_f[2])

    if cur_repeat_region and len(cur_repeat_region) == 1:
        repeat_region_mark.append((1, cur_repeat_region[0][1], max_repeat_region_end,
                                   cur_repeat_region[0][0], cur_repeat_region_score[0]))
    elif len(cur_repeat_region) > 1:
        top_three_tr = heapq.nlargest(3, zip(cur_repeat_region_score, cur_repeat_region))
        if len(top_three_tr) == 2:
            repeat_region_mark.append((2, cur_repeat_region[0][1], max_repeat_region_end,
                                       top_three_tr[0][1][0], top_three_tr[0][0], top_three_tr[1][1][0],
                                       top_three_tr[1][0]))
        else:
            repeat_region_mark.append((3, cur_repeat_region[0][1], max_repeat_region_end,
                                       top_three_tr[0][1][0], top_three_tr[0][0], top_three_tr[1][1][0],
                                       top_three_tr[1][0], top_three_tr[2][1][0], top_three_tr[2][0]))

    return final_trs_list, repeat_region_mark
