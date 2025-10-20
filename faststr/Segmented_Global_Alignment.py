import math
import re
from collections import deque
import numpy as np
import pandas as pd



# 寻找最佳解释：原则是如果有整motif长度，首选整motif长度，如果有多个，选取错误率最低的；如果没有整motif长度，则选取他周边是整motif长度，且编辑距
# 离最低的，如果也有多个，则选取错误率最低的
# RPTA
def realtime_path_tracking_alignment(standard_seq, seq, mat, mis, gap, logo=0):
    dp = [[0] * (len(standard_seq) + 1) for _ in range(len(seq) + 1)]
    upper_left = [[set() for _ in range(len(standard_seq) + 1)] for _ in range(len(seq) + 1)]
    priority_path = [[-1] * (len(standard_seq) + 1) for _ in range(len(seq) + 1)]

    # 初始化
    for i in range(len(seq) + 1):
        if i == 0:
            dp[i][0] = 0
        else:
            dp[i][0] = i * gap
    for j in range(len(standard_seq) + 1):
        if j == 0:
            dp[0][j] = 0
        else:
            dp[0][j] = j * gap
    for i in range(len(seq) + 1):
        upper_left[i][0].add(0)
    for j in range(len(standard_seq) + 1):
        upper_left[0][j].add(0)
    for i in range(len(seq) + 1):
        if i == 0:
            priority_path[i][0] = -1
        else:
            priority_path[i][0] = (i - 1, 0, -1)
    for j in range(len(standard_seq) + 1):
        if j == 0:
            priority_path[0][j] = -1
        else:
            priority_path[0][j] = (0, j - 1, -1)

    # 遍历序列，进行比对
    for i in range(1, len(seq) + 1):
        for j in range(1, len(standard_seq) + 1):
            # 选择状态基于插入、删除、替换，用三位二进制来记录
            chose_path = 0
            delete_cost = dp[i][j - 1] + gap
            insert_cost = dp[i - 1][j] + gap
            replace_cost = dp[i - 1][j - 1] + (mat if standard_seq[j - 1] == seq[i - 1] else mis)
            dp[i][j] = max(insert_cost, delete_cost, replace_cost)
            # 更新最佳比对的替换数量
            if dp[i][j] == delete_cost:
                upper_left[i][j] = upper_left[i][j - 1]
                chose_path += 2
            if dp[i][j] == insert_cost:
                # if not upper_left[i][j]:
                if chose_path == 0:
                    upper_left[i][j] = upper_left[i - 1][j]
                else:
                    upper_left[i][j] = upper_left[i][j - 1] | upper_left[i - 1][j]
                chose_path += 4
            if dp[i][j] == replace_cost:
                updated_set = {x + 1 for x in upper_left[i - 1][j - 1]}
                # if not upper_left[i][j]:
                if chose_path == 0:
                    upper_left[i][j] = updated_set
                else:
                    upper_left[i][j] = upper_left[i][j] | updated_set
                chose_path += 1
            # 更新比对路径,0表示不匹配、1表示匹配、-1表示indel
            if chose_path == 1:
                priority_path[i][j] = (i - 1, j - 1, 1) if standard_seq[j - 1] == seq[i - 1] else (i - 1, j - 1, 0)
            elif chose_path == 2:
                priority_path[i][j] = (i, j - 1, -1)
            elif chose_path == 4:
                priority_path[i][j] = (i - 1, j, -1)
            elif chose_path == 3:
                if min(upper_left[i - 1][j - 1]) + 1 <= min(upper_left[i][j - 1]):
                    priority_path[i][j] = (i - 1, j - 1, 1) if standard_seq[j - 1] == seq[i - 1] else (i - 1, j - 1, 0)
                else:
                    priority_path[i][j] = (i, j - 1, -1)
            elif chose_path == 5:
                if min(upper_left[i - 1][j - 1]) + 1 <= min(upper_left[i - 1][j]):
                    priority_path[i][j] = (i - 1, j - 1, 1) if standard_seq[j - 1] == seq[i - 1] else (i - 1, j - 1, 0)
                else:
                    priority_path[i][j] = (i - 1, j, -1)
            elif chose_path == 6:
                if min(upper_left[i][j - 1]) <= min(upper_left[i - 1][j]):
                    priority_path[i][j] = (i, j - 1, -1)
                else:
                    priority_path[i][j] = (i - 1, j, -1)
            else:
                min_ul = min(min(upper_left[i][j - 1]), min(upper_left[i - 1][j]), min(upper_left[i - 1][j - 1]) + 1)
                if min_ul == min(upper_left[i - 1][j - 1]) + 1:
                    priority_path[i][j] = (i - 1, j - 1, 1) if standard_seq[j - 1] == seq[i - 1] else (i - 1, j - 1, 0)
                elif min_ul == min(upper_left[i][j - 1]):
                    priority_path[i][j] = (i, j - 1, -1)
                else:
                    priority_path[i][j] = (i - 1, j, -1)

    if logo == 1:
        return dp, upper_left, priority_path
    return dp[len(seq)], upper_left[len(seq)], priority_path


# MCLA
def motif_constrained_local_align(standard_seq, seq, mat, mis, gap, n):
    dp = [[0] * (len(standard_seq) + 1) for _ in range(len(seq) + 1)]
    upper_left = [[{(0, -1, -1)} for _ in range(len(standard_seq) + 1)] for _ in range(len(seq) + 1)]
    priority_path = [[-1] * (len(standard_seq) + 1) for _ in range(len(seq) + 1)]
    # new_tr_start = [-1] * len(seq)
    max_score = 0

    for i in range(1, len(seq) + 1):
        for j in range(1, len(standard_seq) + 1):
            if j > 1 and priority_path[i - 1][j] == -1 and priority_path[i][j - 1] == -1 and priority_path[
                i - 1][j - 1] == -1:  # 确保B从第一个碱基开始比对
                continue
            match = mis
            if j == 1 or priority_path[i - 1][j - 1] != -1:
                match = dp[i - 1][j - 1] + (mat if seq[i - 1] == standard_seq[j - 1] else mis)
            delete = dp[i][j - 1] + gap
            insert = dp[i - 1][j] + gap
            dp[i][j] = max(0, match, delete, insert)
            if match < 0 and delete < 0 and insert < 0:
                continue
            chose_path = 0
            # 更新最佳比对的替换数量
            if dp[i][j] == delete:
                upper_left[i][j] = upper_left[i][j - 1]
                chose_path += 2
            if dp[i][j] == insert:
                if chose_path == 0:
                    upper_left[i][j] = upper_left[i - 1][j]
                else:
                    upper_left[i][j] = upper_left[i][j - 1] | upper_left[i - 1][j]
                chose_path += 4
            if dp[i][j] == match:
                updated_set = set()
                for x in upper_left[i - 1][j - 1]:
                    if priority_path[i - 1][j - 1] == -1:
                        updated_set.add((x[0] + 1, i, j))
                    else:
                        updated_set.add((x[0] + 1, x[1], x[2]))
                # updated_set = {(x[0] + 1, x[1], x[2]) for x in upper_left[i - 1][j - 1]}
                if chose_path == 0:
                    upper_left[i][j] = updated_set
                else:
                    upper_left[i][j] = upper_left[i][j] | updated_set
                chose_path += 1
            # 更新比对路径,0表示不匹配、1表示匹配、-1表示indel
            if chose_path == 1:
                priority_path[i][j] = (i - 1, j - 1, 1) if standard_seq[j - 1] == seq[i - 1] else (i - 1, j - 1, 0)
            elif chose_path == 2:
                priority_path[i][j] = (i, j - 1, -1)
            elif chose_path == 4:
                priority_path[i][j] = (i - 1, j, -1)
            elif chose_path == 3:
                if min(t[0] for t in upper_left[i - 1][j - 1]) + 1 <= min(t[0] for t in upper_left[i][j - 1]):
                    priority_path[i][j] = (i - 1, j - 1, 1) if standard_seq[j - 1] == seq[i - 1] else (
                        i - 1, j - 1, 0)
                else:
                    priority_path[i][j] = (i, j - 1, -1)
            elif chose_path == 5:
                if min(t[0] for t in upper_left[i - 1][j - 1]) + 1 <= min(t[0] for t in upper_left[i - 1][j]):
                    priority_path[i][j] = (i - 1, j - 1, 1) if standard_seq[j - 1] == seq[i - 1] else (
                        i - 1, j - 1, 0)
                else:
                    priority_path[i][j] = (i - 1, j, -1)
            elif chose_path == 6:
                if min(t[0] for t in upper_left[i][j - 1]) <= min(t[0] for t in upper_left[i - 1][j]):
                    priority_path[i][j] = (i, j - 1, -1)
                else:
                    priority_path[i][j] = (i - 1, j, -1)
            else:
                min_ul = min(min(t[0] for t in upper_left[i][j - 1]), min(t[0] for t in upper_left[i - 1][j]),
                             min(t[0] for t in upper_left[i - 1][j - 1]) + 1)
                if min_ul == min(t[0] for t in upper_left[i - 1][j]) + 1:
                    priority_path[i][j] = (i - 1, j - 1, 1) if standard_seq[j - 1] == seq[i - 1] else (
                        i - 1, j - 1, 0)
                elif min_ul == min(t[0] for t in upper_left[i][j - 1]):
                    priority_path[i][j] = (i, j - 1, -1)
                else:
                    priority_path[i][j] = (i - 1, j, -1)

            if i == len(seq) and j % n == 0 and dp[i][j] > max_score:
                max_score = dp[i][j]

    return dp[len(seq)], upper_left[len(seq)], priority_path, max_score


def format_alignment(align_path, seq, motif, stab_l):
    target = deque()
    query = deque()
    alignment = deque()
    i = len(seq)
    j = stab_l
    length = 0
    seq_length = len(seq) - 1
    motif_length = len(motif) - 1
    while True:
        if align_path[i][j] == -1:
            break
        length += 1
        if align_path[i][j][2] == 0:
            query.appendleft(seq[seq_length])
            seq_length -= 1
            target.appendleft(motif[motif_length])
            motif_length -= 1
            if motif_length == -1:
                motif_length += len(motif)
            alignment.appendleft('.')
            ii = align_path[i][j][0]
            j = align_path[i][j][1]
            i = ii
        elif align_path[i][j][2] == 1:
            query.appendleft(seq[seq_length])
            seq_length -= 1
            target.appendleft(motif[motif_length])
            motif_length -= 1
            if motif_length == -1:
                motif_length += len(motif)
            alignment.appendleft('|')
            ii = align_path[i][j][0]
            j = align_path[i][j][1]
            i = ii
        else:
            if align_path[i][j][1] == j:
                query.appendleft(seq[seq_length])
                seq_length -= 1
                target.appendleft('-')
            else:
                query.appendleft('-')
                target.appendleft(motif[motif_length])
                motif_length -= 1
                if motif_length == -1:
                    motif_length += len(motif)
            alignment.appendleft('-')
            ii = align_path[i][j][0]
            j = align_path[i][j][1]
            i = ii

    format_align = [''.join(target), ''.join(alignment), ''.join(query), stab_l, length, len(seq) - 1 - seq_length]

    return format_align


def global_motif_multiple_align(seq, motif, logo, indel, mat, mis, ope, ext, marks):
    standard_seq_length, real_seq_length = math.ceil(len(seq) * (1 + indel + 0.005)), len(seq)
    n = len(motif)
    standard_seq = motif[marks % n:] + motif * (standard_seq_length // n)
    # 如果是开头
    if logo == 1:
        standard_seq = motif + motif * (standard_seq_length // n)
        align_score, diagonal_num, align_path, max_score = motif_constrained_local_align(standard_seq, seq, mat, mis,
                                                                                         ope, n)
        candidate_length = []
        indel = real_seq_length + 1
        mismatch = real_seq_length + 1
        insert_num = 0
        delete_num = 0
        for tr_l, a_s in enumerate(align_score):
            if tr_l % n == 0 and a_s == max_score:
                candidate_length.append(tr_l)
        stab_l = candidate_length[0]
        for c_l in candidate_length:
            min_up_left = min(diagonal_num[c_l], key=lambda x: x[0])
            indel_num = c_l + (real_seq_length - min_up_left[1] + 1) - 2 * min_up_left[0]
            mismatch_num = (align_score[c_l] - indel_num * ope - min_up_left[0] * mat) / (mis - mat)
            mismatch = mismatch_num if indel_num < indel else mismatch
            insert_num = (real_seq_length - min_up_left[1] + 1) - min_up_left[0] if indel_num < indel else insert_num
            delete_num = indel_num - insert_num if indel_num < indel else delete_num
            stab_l = c_l if indel_num < indel else stab_l
            indel = indel_num if indel_num < indel else indel
        if max_score == 0 and align_path[len(seq)][stab_l] == -1:
            return 0, 0, 0, '', len(seq)
        format_align = format_alignment(align_path, seq, motif, stab_l)
        left = real_seq_length - format_align[5]
        return insert_num, delete_num, mismatch, format_align, left
    # 如果是结尾，则不需要整motif长度,选取最佳匹配长度
    if logo == 2:
        align_score, diagonal_num, align_path = realtime_path_tracking_alignment(standard_seq, seq, mat, mis, ope, 1)
        max_score = np.max(np.array(align_score))
        if max_score < 0:
            return 0, 0, 0, '', len(seq)
        candidate_length = []
        indel = real_seq_length + 1
        mismatch = real_seq_length + 1
        insert_num = 0
        delete_num = 0
        row = 0
        for ro in range(len(seq), 0, -1):
            for tr_l, a_s in enumerate(align_score[ro]):
                if a_s == max_score:
                    candidate_length.append(tr_l)
            if candidate_length != []:
                row = ro
                break
        if candidate_length == []:
            return 0, 0, 0, '', len(seq)
        # 选取indel错误率最低的
        stab_l = candidate_length[0]
        for c_l in candidate_length:
            indel_num = c_l + row - 2 * min(diagonal_num[row][c_l])
            mismatch_num = (align_score[row][c_l] - indel_num * ope - min(diagonal_num[row][c_l]) * mat) / (mis - mat)
            mismatch = mismatch_num if indel_num < indel else mismatch
            insert_num = row - min(diagonal_num[row][c_l]) if indel_num < indel else insert_num
            delete_num = indel_num - insert_num if indel_num < indel else delete_num
            stab_l = c_l if indel_num < indel else stab_l
            indel = indel_num if indel_num < indel else indel
        format_align = format_alignment(align_path, seq[:row], motif, stab_l)
        right = real_seq_length - row
        return insert_num, delete_num, mismatch, format_align, right
    # 如果不是结尾，那么必须是整motif长度
    align_score, diagonal_num, align_path = realtime_path_tracking_alignment(standard_seq, seq, mat, mis, ope)
    candidate_length = []
    indel = real_seq_length + 1
    mismatch = real_seq_length + 1
    insert_num = 0
    delete_num = 0
    for tr_l, a_s in enumerate(align_score):
        if (tr_l - marks % n) % n == 0 and a_s == max(align_score):
            candidate_length.append(tr_l)
    # 如果有整motif长度,选取indel错误率最低的
    if candidate_length:
        stab_l = candidate_length[0]
        for c_l in candidate_length:
            indel_num = c_l + real_seq_length - 2 * min(diagonal_num[c_l])
            mismatch_num = (align_score[c_l] - indel_num * ope - min(diagonal_num[c_l]) * mat) / (mis - mat)
            mismatch = mismatch_num if indel_num < indel else mismatch
            insert_num = real_seq_length - min(diagonal_num[c_l]) if indel_num < indel else insert_num
            delete_num = indel_num - insert_num if indel_num < indel else delete_num
            stab_l = c_l if indel_num < indel else stab_l
            indel = indel_num if indel_num < indel else indel
        format_align = format_alignment(align_path, seq, motif, stab_l)
        return insert_num, delete_num, mismatch, format_align
    # 如果没有整motif长度，那么选取比对分数最高的整motif长度，并从中选取indel最低的
    integralmultiples_motif = []
    for tr_l, _ in enumerate(align_score):
        if (tr_l - marks % n) % n == 0:
            integralmultiples_motif.append((tr_l, align_score[tr_l]))
    integralmultiples_motif.sort(key=lambda x: x[1])
    for i_m in integralmultiples_motif:
        if i_m[1] == integralmultiples_motif[-1][1]:
            candidate_length.append(i_m[0])
    stab_l = candidate_length[0]
    for c_l in candidate_length:
        indel_num = c_l + real_seq_length - 2 * min(diagonal_num[c_l])
        mismatch_num = (align_score[c_l] - indel_num * ope - min(diagonal_num[c_l]) * mat) / (mis - mat)
        mismatch = mismatch_num if indel_num < indel else mismatch
        insert_num = real_seq_length - min(diagonal_num[c_l]) if indel_num < indel else insert_num
        delete_num = indel_num - insert_num if indel_num < indel else delete_num
        stab_l = c_l if indel_num < indel else stab_l
        indel = indel_num if indel_num < indel else indel
    format_align = format_alignment(align_path, seq, motif, stab_l)
    return insert_num, delete_num, mismatch, format_align


def segmented_global_align_algorithm(sequence, motif, indel, mat, mis, ope, ext, beta=0.045):
    sequence = str(sequence)
    marks = len(motif)
    if len(motif) == 2:
        marks = 5
    elif len(motif) == 3:
        marks = 6
    probe = str(marks // len(motif) * motif + motif[:marks % len(motif)])
    pattern = re.compile('(?={})'.format(re.escape(str(probe))))
    com_match = [match.start() for match in pattern.finditer(sequence)]

    if len(com_match) * len(probe) / len(sequence) < beta:
        return None, None, None, None, None, None, None, None

    distances = [com_match[i + 1] - com_match[i] for i in range(len(com_match) - 1)]
    dis_series = pd.Series(distances)
    most_common_dis = dis_series.value_counts().idxmax()
    if most_common_dis != len(motif):
        return None, None, None, None, None, None, None, None

    Insert = 0
    Delete = 0
    Mismatch = 0
    Target = []
    Align = []
    Query = []
    target_length = 0
    align_length = 0
    pre_end = 0
    chaos_seq_start = -1
    left = 0
    for index, cm in enumerate(com_match):
        # 未出现混沌序列
        if cm - pre_end == 0 and chaos_seq_start < 0:
            pre_end = pre_end + marks
            if cm > 0 and marks % len(motif) > 0:
                Delete += (len(motif) - marks % len(motif))
                Target.append(motif[marks % len(motif):])
                Align.append('-' * (len(motif) - marks % len(motif)))
                Query.append('-' * (len(motif) - marks % len(motif)))
                align_length += (len(motif) - marks % len(motif))
                target_length += (len(motif) - marks % len(motif))
            if index + 1 == len(com_match) or com_match[index + 1] >= pre_end:
                Target.append(probe)
                Align.append('|' * marks)
                Query.append(probe)
                align_length += marks
                target_length += marks
            continue
        if cm - pre_end > 0 and chaos_seq_start < 0:
            '''
            起始端单独考虑
            '''
            if pre_end == 0:
                insert_num, delete_num, mismatch_num, format_align, left = global_motif_multiple_align(
                    sequence[pre_end:cm], motif, 1, indel, mat, mis, mis, mis, marks)
                Insert += insert_num
                Delete += delete_num
                Mismatch += mismatch_num
                pre_end = cm + marks
                if format_align == '':
                    if index + 1 == len(com_match) or com_match[index + 1] >= pre_end:
                        Target.append(probe)
                        Align.append('|' * marks)
                        Query.append(probe)
                        target_length += marks
                        align_length += marks
                    continue
                Target.append(format_align[0])
                Align.append(format_align[1])
                Query.append(format_align[2])
                target_length += format_align[3]
                align_length += format_align[4]
                if index + 1 == len(com_match) or com_match[index + 1] >= pre_end:
                    Target.append(probe)
                    Align.append('|' * marks)
                    Query.append(probe)
                    target_length += marks
                    align_length += marks
                continue
            insert_num, delete_num, mismatch_num, format_align = global_motif_multiple_align(sequence[pre_end:cm],
                                                                                             motif, 0, indel,
                                                                                             mat, mis, mis, mis, marks)
            Insert += insert_num
            Delete += delete_num
            Mismatch += mismatch_num
            pre_end = cm + marks
            Target.append(format_align[0])
            Align.append(format_align[1])
            Query.append(format_align[2])
            target_length += format_align[3]
            align_length += format_align[4]
            if index + 1 == len(com_match) or com_match[index + 1] >= pre_end:
                Target.append(probe)
                Align.append('|' * marks)
                Query.append(probe)
                target_length += marks
                align_length += marks
            continue
        # 发现混沌序列起始部分
        if chaos_seq_start == -1:
            chaos_seq_start = pre_end - marks
            pre_end = cm + marks
            continue
        # 继续发现混沌序列部分
        if chaos_seq_start >= 0 and cm - pre_end < 0:
            pre_end = cm + marks
            continue
        # 发现混沌序列结尾部分
        # motif为2和3对应的混沌序列包含一个大的重构标记序列，不可能是混沌序列
        if len(motif) in [2, 3]:
            Target.append(sequence[chaos_seq_start:pre_end])
            Align.append('|' * len(sequence[chaos_seq_start:pre_end]))
            Query.append(sequence[chaos_seq_start:pre_end])
            target_length += len(sequence[chaos_seq_start:pre_end])
            align_length += len(sequence[chaos_seq_start:pre_end])
            chaos_seq_start = pre_end
        insert_num, delete_num, mismatch_num, format_align = global_motif_multiple_align(sequence[chaos_seq_start:cm],
                                                                                         motif, 0, indel, mat, mis, mis,
                                                                                         mis, marks)
        Insert += insert_num
        Delete += delete_num
        Mismatch += mismatch_num
        chaos_seq_start = -1
        pre_end = cm + marks
        Target.append(format_align[0])
        Align.append(format_align[1])
        Query.append(format_align[2])
        target_length += format_align[3]
        align_length += format_align[4]
        # 如果该新的标记序列处于最后一个或者不是下一个混沌序列的开始，那么可以直接比对，否则纳入到下一次比对中
        if index == len(com_match) - 1 or com_match[index + 1] - pre_end >= 0:
            Target.append(probe)
            Align.append('|' * marks)
            Query.append(probe)
            target_length += marks
            align_length += marks

    # 处理剩余序列
    # 如果不是处于混沌序列部分,直接比对
    right = 0
    if chaos_seq_start == -1 and pre_end < len(sequence):
        insert_num, delete_num, mismatch_num, format_align, right = global_motif_multiple_align(sequence[pre_end:],
                                                                                                motif, 2, indel, mat,
                                                                                                mis, mis, mis, marks)
        if format_align != '':
            Insert += insert_num
            Delete += delete_num
            Mismatch += mismatch_num
            Target.append(format_align[0])
            Align.append(format_align[1])
            Query.append(format_align[2])
            target_length += format_align[3]
            align_length += format_align[4]
        # 否则如果还处于混沌序列之中，
    elif chaos_seq_start >= 0:
        insert_num, delete_num, mismatch_num, format_align, right = global_motif_multiple_align(
            sequence[chaos_seq_start:], motif, 2, indel, mat, mis, mis, mis, marks)
        if format_align != '':
            Insert += insert_num
            Delete += delete_num
            Mismatch += mismatch_num
            Target.append(format_align[0])
            Align.append(format_align[1])
            Query.append(format_align[2])
            target_length += format_align[3]
            align_length += format_align[4]

    score = (Insert + Delete) * ope + Mismatch * mis + ((len(sequence) - right) - Mismatch - Insert) * mat
    Target = ''.join(Target)
    Align = ''.join(Align)
    Query = ''.join(Query)
    align = 'target            0 ' + Target + ' ' + str(
        target_length) + '\n' + '                  0 ' + Align + ' ' + str(
        align_length) + '\n' + 'query             0 ' + Query + ' ' + str(len(sequence) - left - right)
    return Insert, Delete, Mismatch, align, left, right, score, (len(Target) - Target.count('-')) / len(motif)

# tr = 'TAAAATAAACAAAATAAAATAAATTAAACAATTAAATTAAATAAAATAAATTAAA'
# motif = 'TAAAA'
# Insert, Delete, Mismatch, align, _, right, _, copy = segmented_global_align_algorithm(tr, motif, 0.15, 2, -5, -7, -3)
# print(align)
# print(copy)