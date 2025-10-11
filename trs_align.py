import math
from operator import itemgetter
import numpy as np
import regex
import scan_subread
import Segmented_Global_Alignment
from Bio.Align import PairwiseAligner
import utils


def can_merge_to_onetr(seq1, seq2, motif1, motif2, end, start, interrupted_seq):
    joint_motif = judge_motif_sim(motif1, motif2, seq1, seq2)
    if joint_motif != '' and 0 <= start - end - 1 <= 100:
        return joint_motif
    if len(motif1) < 7 and 0 <= start - end - 1 <= 100:
        pattern = f"({motif1}){{s<=1}}"
        matches = regex.finditer(pattern, str(motif2 + motif2))
        if any(matches):
            return can_standardize_motif(seq1, seq2, motif1)
    return ''


def trs_align_algorithm(n_candidate_trs, sub_read, p_indel, p_match, mat, mis, gap, ext, align_score, beta):
    if n_candidate_trs == []:
        return []
    qualified_trs = []
    n_candidate_trs.sort(key=itemgetter(1))
    after_merge_trs = []
    maybe_merge_list = [n_candidate_trs[0]]
    consensus_motif = ''
    for index, c_t in enumerate(n_candidate_trs[1:]):
        if maybe_merge_list[-1] == c_t:
            continue
        consensus_motif_test = can_merge_to_onetr(sub_read[maybe_merge_list[0][1]:maybe_merge_list[-1][2] + 1],
                                                  sub_read[c_t[1]:c_t[2] + 1], maybe_merge_list[-1][0], c_t[0],
                                                  maybe_merge_list[-1][2], c_t[1], sub_read[
                                                                                   maybe_merge_list[-1][2] + 1:max(
                                                                                       maybe_merge_list[-1][2] + 1,
                                                                                       c_t[1])])
        if consensus_motif_test != '':
            maybe_merge_list.append(c_t)
            consensus_motif = consensus_motif_test
            continue
        if index + 1 == len(n_candidate_trs[1:]):
            if len(maybe_merge_list) == 1:
                after_merge_trs.append(maybe_merge_list[0])
                maybe_merge_list = [c_t]
            else:
                final_merge_trs, early_qualified_trs = calculate_character_distance(maybe_merge_list, sub_read,
                                                                                    consensus_motif, p_indel, p_match,
                                                                                    mat, mis, gap, ext, align_score,
                                                                                    beta)
                after_merge_trs.extend(final_merge_trs)
                qualified_trs.extend(early_qualified_trs)
                maybe_merge_list = [c_t]
                consensus_motif = ''
            break
        i = 2
        while index + i < len(n_candidate_trs):
            # 右重复
            if maybe_merge_list[-1][2] < c_t[1]:
                # 重复部分均无法满足前合并
                if n_candidate_trs[index + i][1] > c_t[2]:
                    if len(maybe_merge_list) == 1:
                        after_merge_trs.append(maybe_merge_list[0])
                        maybe_merge_list = [c_t]
                    else:
                        final_merge_trs, early_qualified_trs = calculate_character_distance(maybe_merge_list, sub_read,
                                                                                            consensus_motif, p_indel,
                                                                                            p_match, mat, mis, gap, ext,
                                                                                            align_score, beta)
                        after_merge_trs.extend(final_merge_trs)
                        qualified_trs.extend(early_qualified_trs)
                        maybe_merge_list = [c_t]
                        consensus_motif = ''
                    break
                # 重复部分前合并成功
                consensus_motif_test = can_merge_to_onetr(sub_read[maybe_merge_list[0][1]:maybe_merge_list[-1][2] + 1],
                                                          sub_read[
                                                          n_candidate_trs[index + i][1]:n_candidate_trs[index + i][
                                                                                            2] + 1],
                                                          maybe_merge_list[-1][0], n_candidate_trs[index + i][0],
                                                          maybe_merge_list[-1][2], n_candidate_trs[index + i][1],
                                                          sub_read[
                                                          maybe_merge_list[-1][2] + 1:max(maybe_merge_list[-1][2] + 1,
                                                                                          n_candidate_trs[index + i][
                                                                                              1])])
                if consensus_motif_test != '':
                    maybe_merge_list.append(n_candidate_trs[index + i])
                    consensus_motif = consensus_motif_test
                    break
                # 重复部分前合并失败，但还处于重复之中
                i += 1
                continue
            # 左重复
            else:
                # 还处于重复之中，继续遍历
                if n_candidate_trs[index + i][1] < c_t[2]:
                    i += 1
                    continue
                # 刚好跳出重复，后合并成功
                consensus_motif_test = can_merge_to_onetr(sub_read[maybe_merge_list[0][1]:maybe_merge_list[-1][2] + 1],
                                                          sub_read[
                                                          n_candidate_trs[index + i][1]:n_candidate_trs[index + i][
                                                                                            2] + 1],
                                                          maybe_merge_list[-1][0], n_candidate_trs[index + i][0],
                                                          maybe_merge_list[-1][2], n_candidate_trs[index + i][1],
                                                          sub_read[
                                                          maybe_merge_list[-1][2] + 1:max(maybe_merge_list[-1][2] + 1,
                                                                                          n_candidate_trs[index + i][
                                                                                              1])])
                if consensus_motif_test != '':
                    maybe_merge_list.append(n_candidate_trs[index + i])
                    consensus_motif = consensus_motif_test
                    break
                # 刚好跳出重复，后合并失败
                if len(maybe_merge_list) == 1:
                    after_merge_trs.append(maybe_merge_list[0])
                    maybe_merge_list = [c_t]
                else:
                    final_merge_trs, early_qualified_trs = calculate_character_distance(maybe_merge_list, sub_read,
                                                                                        consensus_motif, p_indel,
                                                                                        p_match, mat, mis, gap, ext,
                                                                                        align_score, beta)
                    after_merge_trs.extend(final_merge_trs)
                    qualified_trs.extend(early_qualified_trs)
                    maybe_merge_list = [c_t]
                    consensus_motif = ''
                break

    if len(maybe_merge_list) == 1:
        after_merge_trs.append(maybe_merge_list[0])
    else:
        final_merge_trs, early_qualified_trs = calculate_character_distance(maybe_merge_list, sub_read, consensus_motif,
                                                                            p_indel, p_match, mat, mis, gap, ext,
                                                                            align_score, beta)
        after_merge_trs.extend(final_merge_trs)
        qualified_trs.extend(early_qualified_trs)

    for pos, a_m_t in enumerate(after_merge_trs):
        if pos > 0 and a_m_t[1] >= after_merge_trs[pos - 1][1] and a_m_t[2] <= after_merge_trs[pos - 1][2] and \
                after_merge_trs[pos - 1][0] in a_m_t[0] + a_m_t[0]:
            continue
        if a_m_t[2] - a_m_t[1] + 1 > 6000:
            insert, delete, mismatch, align, left, right, score, copy = Segmented_Global_Alignment.segmented_global_align_algorithm(
                sub_read[a_m_t[1]:a_m_t[2] + 1], a_m_t[0], p_indel, mat, mis, gap, ext, beta)
            if insert == None:
                continue
            length = a_m_t[2] - a_m_t[1] + 1 - right - left
            indel_rio = (insert + delete) / (length + delete)
            match_rio = (length - mismatch - insert) / (length - insert)
        else:
            l = round((a_m_t[2] - a_m_t[1] + 1) * (1 + p_indel + 0.001))
            seq = l // len(a_m_t[0]) * a_m_t[0] + a_m_t[0][:l % len(a_m_t[0])]

            # 创建 PairwiseAligner 对象
            aligner = PairwiseAligner()
            aligner.mode = 'local'  # 设置为局部比对模式
            # 设置比对参数
            aligner.match_score = mat
            aligner.mismatch_score = mis
            aligner.open_gap_score = gap
            aligner.extend_gap_score = ext
            # 执行局部比对
            try:
                alignments = aligner.align(seq, sub_read[a_m_t[1]:a_m_t[2] + 1])
                if alignments:
                    delete = alignments[0].length - alignments[0].coordinates[1][-1] + alignments[0].coordinates[1][0]
                    insert = alignments[0].counts()[0] - delete
                    mismatch = alignments[0].counts()[2]
                    length = alignments[0].length
                    score = alignments[0].score
                    align = alignments[0].format()
                    copy = (alignments[0].coordinates[0][-1]-alignments[0].coordinates[0][0]) / len(a_m_t[0])
                    left = alignments[0].coordinates[1][0]
                    right = a_m_t[2] - a_m_t[1] + 1 - alignments[0].coordinates[1][-1]
                else:
                    continue
            except OverflowError as e:
                continue

            indel_rio = (insert + delete) / (length + delete)
            match_rio = (length - mismatch - insert) / (length - insert)

        if length >= 25 and indel_rio <= p_indel and match_rio >= p_match and score >= align_score:
            qualified_trs.append(
                (a_m_t[0], a_m_t[1] + left, a_m_t[2] - right, indel_rio, match_rio, score, align, copy))

            if 100 <= left:
                left_candidate = scan_subread.filter_chaotic_repeats(
                    [(a_m_t[1], a_m_t[1] + left, a_m_t[1], a_m_t[1] + left, 1)], sub_read, len(a_m_t[0]))
                if left_candidate != []:
                    after_merge_trs.append(left_candidate[0])

            if right >= 100:
                right_candidate = scan_subread.filter_chaotic_repeats(
                    [(a_m_t[2] - right + 1, a_m_t[2], a_m_t[2] - right + 1, a_m_t[2], 1)], sub_read, len(a_m_t[0]))
                if right_candidate != []:
                    after_merge_trs.append(right_candidate[0])
        elif length >= 25 and try_change_motif(score, indel_rio, match_rio, align_score, p_indel, p_match):
            motif_p = utils.tri_gram_model(sub_read[a_m_t[1]:a_m_t[2] + 1], len(a_m_t[0]))
            if motif_p == None:
                continue
            for m_p in motif_p:
                if a_m_t[0] in m_p + m_p:
                    continue
                if a_m_t[2] - a_m_t[1] + 1 > 6000:
                    insert_p, delete_p, mismatch_p, align_p, left_p, right_p, score_p, copy_p = Segmented_Global_Alignment.segmented_global_align_algorithm(
                        sub_read[a_m_t[1]:a_m_t[2] + 1], m_p, p_indel, mat, mis, gap, ext, beta)
                    if insert == None:
                        continue
                    length_p = a_m_t[2] - a_m_t[1] + 1 - right_p - left_p
                    indel_rio_p = (insert_p + delete_p) / (length_p + delete_p)
                    match_rio_p = (length_p - mismatch_p - insert_p) / (length_p - insert_p)
                else:
                    l = round((a_m_t[2] - a_m_t[1] + 1) * (1 + p_indel + 0.001))
                    seq_p = l // len(m_p) * m_p + m_p[:l % len(m_p)]

                    # 创建 PairwiseAligner 对象
                    aligner_p = PairwiseAligner()
                    aligner_p.mode = 'local'  # 设置为局部比对模式
                    # 设置比对参数
                    aligner_p.match_score = mat
                    aligner_p.mismatch_score = mis
                    aligner_p.open_gap_score = gap
                    aligner_p.extend_gap_score = ext
                    # 执行局部比对
                    try:
                        alignments_p = aligner_p.align(seq_p, sub_read[a_m_t[1]:a_m_t[2] + 1])
                        if alignments_p:
                            delete_p = alignments_p[0].length - alignments_p[0].coordinates[1][-1] + \
                                       alignments_p[0].coordinates[1][0]
                            insert_p = alignments_p[0].counts()[0] - delete_p
                            mismatch_p = alignments_p[0].counts()[2]
                            length_p = alignments_p[0].length
                            score_p = alignments_p[0].score
                            align_p = alignments_p[0].format()
                            copy_p = alignments_p[0].coordinates[0][-1] / len(m_p)
                            left_p = alignments_p[0].coordinates[1][0]
                            right_p = a_m_t[2] - a_m_t[1] + 1 - alignments_p[0].coordinates[1][-1]
                        else:
                            continue
                    except OverflowError as e:
                        continue

                    indel_rio_p = (insert_p + delete_p) / (length_p + delete_p)
                    match_rio_p = (length_p - mismatch_p - insert_p) / (length_p - insert_p)

                if length_p >= 25 and indel_rio_p <= p_indel and match_rio_p >= p_match and score_p >= align_score:
                    qualified_trs.append((m_p, a_m_t[1] + left_p, a_m_t[2] - right_p, indel_rio_p, match_rio_p, score_p,
                                          align_p, copy_p))
    return qualified_trs


def can_shorten_distance(seq, motif):
    max_windows = [16, 16, 16, 16, 15, 15, 14, 14]
    motif_mark_indexes = utils.get_motif_marks(seq, motif)
    visit = [0] * len(seq)
    for m_p in motif_mark_indexes:
        visit[m_p:m_p + len(motif)] = [1] * len(motif)
    if len(motif) in [1, 4, 5, 6, 7]:
        window_ones_count = sum(visit[:20])
        max_ones_count = window_ones_count

        for i in range(1, len(visit) - 19):
            # 更新窗口中1的数量
            window_ones_count += visit[i + 19] - visit[i - 1]
            if window_ones_count > max_ones_count:
                max_ones_count = window_ones_count
    elif len(motif) in [2, 3]:
        pattern = f"({motif}){{s<=1}}"
        matches = regex.finditer(pattern, str(seq))
        positions = [match.start() for match in matches]
        for m_p in positions:
            for i in range(len(motif)):
                visit[m_p + i] = visit[m_p + i] + (1 - visit[m_p + i]) * (len(motif) - 1) / len(motif)
        windows = sum(visit[:20])
        window_ones_count = math.floor(windows)
        max_ones_count = window_ones_count

        for i in range(1, len(visit) - 19):
            # 更新窗口中1的数量
            windows += visit[i + 19] - visit[i - 1]
            window_ones_count = math.floor(windows)
            if window_ones_count > max_ones_count:
                max_ones_count = window_ones_count
    else:
        pattern = f"({motif}){{s<=3}}"
        matches = regex.finditer(pattern, str(seq))
        positions = [match.start() for match in matches]
        for m_p in positions:
            for i in range(len(motif)):
                visit[m_p + i] = visit[m_p + i] + (1 - visit[m_p + i]) * (len(motif) - 3) / len(motif)
        windows = sum(visit[:20])
        window_ones_count = math.floor(windows)
        max_ones_count = window_ones_count

        for i in range(1, len(visit) - 19):
            # 更新窗口中1的数量
            windows += visit[i + 19] - visit[i - 1]
            window_ones_count = math.floor(windows)
            if window_ones_count > max_ones_count:
                max_ones_count = window_ones_count

    if max_ones_count < max_windows[len(motif) - 1]:
        return False
    return True


def can_standardize_motif(seq1, seq2, motif1):
    motifs = utils.tri_gram_model(seq1 + seq2, len(motif1))
    if motifs == None:
        return ''
    start_probe = [24, 24, 24, 21, 24, 24, 24, 22]
    end_probe = [24, 24, 24, 22, 23, 23, 23, 22]
    to_align_seq = seq2 if motifs[0] in motif1 + motif1 else seq1
    if can_shorten_distance(to_align_seq, motifs[0]) == False:
        return ''
    probe = 25 // len(motifs[0]) * motifs[0] + motifs[0][:25 % len(motifs[0])]
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    # 设置比对参数
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    start_index = -1
    end_index = -1
    for i in range(len(to_align_seq) // 5 - 3):
        probe_align = aligner.align(probe, to_align_seq[i * 5:i * 5 + 20])
        if probe_align:
            pass
        else:
            continue
        if probe_align[0].score >= start_probe[len(motifs[0]) - 1]:
            start_index = max(i * 5 - 5, 0)
            break
    if start_index == -1:
        return ''
    for j in range((len(to_align_seq) - start_index - 1) // 5 - 3):
        probe_align = aligner.align(probe, to_align_seq[len(to_align_seq) - j * 5 - 20:len(to_align_seq) - j * 5])
        if probe_align:
            pass
        else:
            continue
        if probe_align[0].score >= end_probe[len(motifs[0]) - 1]:
            end_index = min(len(to_align_seq) - j * 5 + 4, len(to_align_seq) - 1)
            break
    if end_index == -1:
        return ''
    if end_index - start_index < 24:
        return ''
    return motifs[0]


def judge_motif_sim(m1, m2, s1, s2):
    if len(m1) < 7 and m1 in m2 + m2:
        return m1
    else:
        pattern = f"({m1}){{s<=1}}"
    matches = regex.finditer(pattern, str(m2 + m2))
    if any(matches):
        motifs = utils.tri_gram_model(s1 + s2, len(m1))
        if motifs == None:
            return ''
        return motifs[0]
    else:
        return ''


def try_change_motif(score, indel, match, align_score, p_indel, p_match):
    if align_score - score > 5 and indel <= p_indel and match >= p_match:
        return False
    if score >= align_score and indel - p_indel > 0.025 and match >= p_match:
        return False
    if score >= align_score and indel <= p_indel and p_match - match > 0.025:
        return False
    return True


def calculate_character_distance(maybe_merge_list_partial, sub_read, consensus_motif, p_indel, p_match, mat, mis, gap,
                                 ext, align_score, beta):
    final_merge_trs = []
    final_merge_list = []
    early_qualified_trs = []
    last = 0
    for inde, m_m_p in enumerate(maybe_merge_list_partial):
        if m_m_p[2] - m_m_p[1] + 1 > 6000:
            insert, delete, mismatch, align, left, right, score, copy = Segmented_Global_Alignment.segmented_global_align_algorithm(
                sub_read[m_m_p[1]:m_m_p[2] + 1], m_m_p[0], p_indel, mat, mis, gap, ext, beta)
            if insert == None:
                continue
            length = m_m_p[2] - m_m_p[1] + 1 - right - left
        else:
            seq = sub_read[m_m_p[1]:m_m_p[2] + 1]
            l = round(len(seq) * (1 + p_indel + 0.005))
            tr = l // len(consensus_motif) * consensus_motif + consensus_motif[:l % len(consensus_motif)]
            # 创建 PairwiseAligner 对象
            aligner = PairwiseAligner()
            aligner.mode = 'local'  # 设置为局部比对模式
            # 设置比对参数
            aligner.match_score = mat
            aligner.mismatch_score = mis
            aligner.open_gap_score = gap
            aligner.extend_gap_score = ext
            # 执行局部比对
            try:
                alignments = aligner.align(tr, seq)
                if alignments:
                    delete = alignments[0].length - alignments[0].coordinates[1][-1] + alignments[0].coordinates[1][0]
                    insert = alignments[0].counts()[0] - delete
                    mismatch = alignments[0].counts()[2]
                    length = alignments[0].length
                    score = alignments[0].score
                    align = alignments[0].format()
                    copy = (alignments[0].coordinates[0][-1]-alignments[0].coordinates[0][0]) / len(consensus_motif)
                    left = alignments[0].coordinates[1][0]
                    right = len(seq) - alignments[0].coordinates[1][-1]
                else:
                    continue
            except OverflowError as e:
                continue

        indel_rio = (insert + delete) / (length + delete)
        match_rio = (length - mismatch - insert) / (length - insert)
        if inde != 0 and left + last <= 70:
            if final_merge_list == []:
                final_merge_list.append(maybe_merge_list_partial[inde - 1])
                if early_qualified_trs and maybe_merge_list_partial[inde - 1][1] <= early_qualified_trs[-1][1] and \
                        maybe_merge_list_partial[inde - 1][2] >= early_qualified_trs[-1][2]:
                    del early_qualified_trs[-1]
            final_merge_list.append(m_m_p)
        elif inde == 0:
            if length >= 25 and indel_rio <= p_indel and match_rio >= p_match and score >= align_score:
                early_qualified_trs.append(
                    (m_m_p[0], m_m_p[1] + left, m_m_p[2] - right, indel_rio, match_rio, score, align, copy))
            elif length >= 25 and try_change_motif(score, indel_rio, match_rio, align_score, p_indel, p_match):
                motif_p = utils.tri_gram_model(sub_read[m_m_p[1]:m_m_p[2] + 1], len(m_m_p[0]))
                if motif_p == None:
                    continue
                for m_p in motif_p:
                    if m_m_p[0] in m_p + m_p:
                        continue
                    if m_m_p[2] - m_m_p[1] + 1 > 6000:
                        insert_p, delete_p, mismatch_p, align_p, left_p, right_p, score_p, copy_p = Segmented_Global_Alignment.segmented_global_align_algorithm(
                            sub_read[m_m_p[1]:m_m_p[2] + 1], m_p, p_indel, mat, mis, gap, ext, beta)
                        if insert == None:
                            continue
                        length_p = m_m_p[2] - m_m_p[1] + 1 - right_p - left_p
                        indel_rio_p = (insert_p + delete_p) / (length_p + delete_p)
                        match_rio_p = (length_p - mismatch_p - insert_p) / (length_p - insert_p)
                    else:
                        l_p = round((m_m_p[2] - m_m_p[1] + 1) * (1 + p_indel + 0.001))
                        seq_p = l_p // len(m_p) * m_p + m_p[:l_p % len(m_p)]

                        # 创建 PairwiseAligner 对象
                        aligner_p = PairwiseAligner()
                        aligner_p.mode = 'local'  # 设置为局部比对模式
                        # 设置比对参数
                        aligner_p.match_score = mat
                        aligner_p.mismatch_score = mis
                        aligner_p.open_gap_score = gap
                        aligner_p.extend_gap_score = ext
                        # 执行局部比对
                        try:
                            alignments_p = aligner_p.align(seq_p, sub_read[m_m_p[1]:m_m_p[2] + 1])
                            if alignments_p:
                                delete_p = alignments_p[0].length - alignments_p[0].coordinates[1][-1] + \
                                           alignments_p[0].coordinates[1][0]
                                insert_p = alignments_p[0].counts()[0] - delete_p
                                mismatch_p = alignments_p[0].counts()[2]
                                length_p = alignments_p[0].length
                                score_p = alignments_p[0].score
                                align_p = alignments_p[0].format()
                                copy_p = alignments_p[0].coordinates[0][-1] / len(m_p)
                                left_p = alignments_p[0].coordinates[1][0]
                                right_p = m_m_p[2] - m_m_p[1] + 1 - alignments_p[0].coordinates[1][-1]
                            else:
                                continue
                        except OverflowError as e:
                            continue

                        indel_rio_p = (insert_p + delete_p) / (length_p + delete_p)
                        match_rio_p = (length_p - mismatch_p - insert_p) / (length_p - insert_p)

                    if length_p >= 25 and indel_rio_p <= p_indel and match_rio_p >= p_match and score_p >= align_score:
                        early_qualified_trs.append((
                            m_p, m_m_p[1] + left_p, m_m_p[2] - right_p, indel_rio_p, match_rio_p,
                            score_p, align_p, copy_p))
        elif final_merge_list == []:
            if length >= 25 and indel_rio <= p_indel and match_rio >= p_match and score >= align_score:
                early_qualified_trs.append(
                    (m_m_p[0], m_m_p[1] + left, m_m_p[2] - right, indel_rio, match_rio, score, align, copy))
            elif length >= 25 and try_change_motif(score, indel_rio, match_rio, align_score, p_indel, p_match):
                motif_p = utils.tri_gram_model(sub_read[m_m_p[1]:m_m_p[2] + 1], len(m_m_p[0]))
                if motif_p == None:
                    continue
                for m_p in motif_p:
                    if m_m_p[0] in m_p + m_p:
                        continue
                    if m_m_p[2] - m_m_p[1] + 1 > 6000:
                        insert_p, delete_p, mismatch_p, align_p, left_p, right_p, score_p, copy_p = Segmented_Global_Alignment.segmented_global_align_algorithm(
                            sub_read[m_m_p[1]:m_m_p[2] + 1], m_p, p_indel, mat, mis, gap, ext, beta)
                        if insert == None:
                            continue
                        length_p = m_m_p[2] - m_m_p[1] + 1 - right_p - left_p
                        indel_rio_p = (insert_p + delete_p) / (length_p + delete_p)
                        match_rio_p = (length_p - mismatch_p - insert_p) / (length_p - insert_p)
                    else:
                        l_p = round((m_m_p[2] - m_m_p[1] + 1) * (1 + p_indel + 0.001))
                        seq_p = l_p // len(m_p) * m_p + m_p[:l_p % len(m_p)]

                        # 创建 PairwiseAligner 对象
                        aligner_p = PairwiseAligner()
                        aligner_p.mode = 'local'  # 设置为局部比对模式
                        # 设置比对参数
                        aligner_p.match_score = mat
                        aligner_p.mismatch_score = mis
                        aligner_p.open_gap_score = gap
                        aligner_p.extend_gap_score = ext
                        # 执行局部比对
                        try:
                            alignments_p = aligner_p.align(seq_p, sub_read[m_m_p[1]:m_m_p[2] + 1])
                            if alignments_p:
                                delete_p = alignments_p[0].length - alignments_p[0].coordinates[1][-1] + \
                                           alignments_p[0].coordinates[1][0]
                                insert_p = alignments_p[0].counts()[0] - delete_p
                                mismatch_p = alignments_p[0].counts()[2]
                                length_p = alignments_p[0].length
                                score_p = alignments_p[0].score
                                align_p = alignments_p[0].format()
                                copy_p = alignments_p[0].coordinates[0][-1] / len(m_p)
                                left_p = alignments_p[0].coordinates[1][0]
                                right_p = m_m_p[2] - m_m_p[1] + 1 - alignments_p[0].coordinates[1][-1]
                            else:
                                continue
                        except OverflowError as e:
                            continue

                        indel_rio_p = (insert_p + delete_p) / (length_p + delete_p)
                        match_rio_p = (length_p - mismatch_p - insert_p) / (length_p - insert_p)

                    if length_p >= 25 and indel_rio_p <= p_indel and match_rio_p >= p_match and score_p >= align_score:
                        early_qualified_trs.append(
                            (m_p, m_m_p[1] + left_p, m_m_p[2] - right_p, indel_rio_p, match_rio_p, score_p, align_p,
                             copy_p))
        else:
            final_merge_trs.append((consensus_motif, final_merge_list[0][1], final_merge_list[-1][2]))
            final_merge_list = []
            if length >= 25 and indel_rio <= p_indel and match_rio >= p_match and score >= align_score:
                early_qualified_trs.append(
                    (m_m_p[0], m_m_p[1] + left, m_m_p[2] - right, indel_rio, match_rio, score, align, copy))
            elif length >= 25 and try_change_motif(score, indel_rio, match_rio, align_score, p_indel, p_match):
                motif_p = utils.tri_gram_model(sub_read[m_m_p[1]:m_m_p[2] + 1], len(m_m_p[0]))
                if motif_p == None:
                    continue
                for m_p in motif_p:
                    if m_m_p[0] in m_p + m_p:
                        continue
                    if m_m_p[2] - m_m_p[1] + 1 > 6000:
                        insert_p, delete_p, mismatch_p, align_p, left_p, right_p, score_p, copy_p = Segmented_Global_Alignment.segmented_global_align_algorithm(
                            sub_read[m_m_p[1]:m_m_p[2] + 1], m_p, p_indel, mat, mis, gap, ext, beta)
                        if insert == None:
                            continue
                        length_p = m_m_p[2] - m_m_p[1] + 1 - right_p - left_p
                        indel_rio_p = (insert_p + delete_p) / (length_p + delete_p)
                        match_rio_p = (length_p - mismatch_p - insert_p) / (length_p - insert_p)
                    else:
                        l_p = round((m_m_p[2] - m_m_p[1] + 1) * (1 + p_indel + 0.001))
                        seq_p = l_p // len(m_p) * m_p + m_p[:l_p % len(m_p)]

                        # 创建 PairwiseAligner 对象
                        aligner_p = PairwiseAligner()
                        aligner_p.mode = 'local'  # 设置为局部比对模式
                        # 设置比对参数
                        aligner_p.match_score = mat
                        aligner_p.mismatch_score = mis
                        aligner_p.open_gap_score = gap
                        aligner_p.extend_gap_score = ext
                        # 执行局部比对
                        try:
                            alignments_p = aligner_p.align(seq_p, sub_read[m_m_p[1]:m_m_p[2] + 1])
                            if alignments_p:
                                delete_p = alignments_p[0].length - alignments_p[0].coordinates[1][-1] + \
                                           alignments_p[0].coordinates[1][0]
                                insert_p = alignments_p[0].counts()[0] - delete_p
                                mismatch_p = alignments_p[0].counts()[2]
                                length_p = alignments_p[0].length
                                score_p = alignments_p[0].score
                                align_p = alignments_p[0].format()
                                copy_p = alignments_p[0].coordinates[0][-1] / len(m_p)
                                left_p = alignments_p[0].coordinates[1][0]
                                right_p = m_m_p[2] - m_m_p[1] + 1 - alignments_p[0].coordinates[1][-1]
                            else:
                                continue
                        except OverflowError as e:
                            continue

                        indel_rio_p = (insert_p + delete_p) / (length_p + delete_p)
                        match_rio_p = (length_p - mismatch_p - insert_p) / (length_p - insert_p)

                    if length_p >= 25 and indel_rio_p <= p_indel and match_rio_p >= p_match and score_p >= align_score:
                        early_qualified_trs.append(
                            (m_p, m_m_p[1] + left_p, m_m_p[2] - right_p, indel_rio_p, match_rio_p, score_p, align_p,
                             copy_p))
        last = right

    if final_merge_list != []:
        final_merge_trs.append((consensus_motif, final_merge_list[0][1], final_merge_list[-1][2]))
    return final_merge_trs, early_qualified_trs


def simple_local_comparison_algorithm(seq, motif, p_indel, mat, mis, gap):
    l = round(len(seq) * (1 + p_indel + 0.005))
    tr = l // len(motif) * motif + motif[:l % len(motif)]
    score_matrix = np.zeros((len(seq) + 1, len(tr) + 1))
    traceback_matrix = np.zeros((len(seq) + 1, len(tr) + 1), dtype=object)

    max_score = 0
    max_pos = None

    for i in range(1, len(seq) + 1):
        for j in range(1, len(tr) + 1):
            match = score_matrix[i - 1, j - 1] + (mat if seq[i - 1] == tr[j - 1] else mis)
            delete = score_matrix[i, j - 1] + gap
            insert = score_matrix[i - 1, j] + gap
            score_matrix[i, j] = max(0, match, delete, insert)

            # 记录得分最高的位置
            if score_matrix[i, j] > max_score:
                max_score = score_matrix[i, j]
                max_pos = (i, j)

            # 偏替换,插入和删除偏删除
            if score_matrix[i, j] == match:
                traceback_matrix[i, j] = (i - 1, j - 1)
            elif score_matrix[i, j] == delete:
                traceback_matrix[i, j] = (i, j - 1)
            elif score_matrix[i, j] == insert:
                traceback_matrix[i, j] = (i - 1, j)

    i, j = max_pos
    align_mark = []
    target = []
    query = []
    delete = 0
    insert = 0
    mismatch = 0
    original_i = i
    original_j = j
    while score_matrix[i, j] != 0:
        ii, jj = traceback_matrix[i, j]
        if ii + 1 == i and jj + 1 == j:
            if seq[ii] == tr[jj]:
                align_mark.append('|')
            else:
                align_mark.append('.')
                mismatch += 1
            target.append(tr[jj])
            query.append(seq[ii])
        elif ii + 1 == i and jj == j:
            align_mark.append('-')
            query.append(seq[ii])
            target.append('-')
            insert += 1
        elif ii == i and jj + 1 == j:
            align_mark.append('-')
            query.append('-')
            target.append(tr[jj])
            delete += 1
        i = ii
        j = jj
    align_mark = ''.join(list(reversed(align_mark)))
    target = ''.join(list(reversed(target)))
    query = ''.join(list(reversed(query)))
    align = []
    target_start = j
    query_start = i
    align_start = 0
    for k in range(math.ceil(len(align_mark) / 60)):
        target_local = target[align_start:min(len(align_mark), align_start + 60)]
        query_local = query[align_start:min(len(align_mark), align_start + 60)]
        align_local = align_mark[align_start:min(len(align_mark), align_start + 60)]
        target_local = 'target' + (13 - len(str(target_start))) * ' ' + str(target_start) + ' ' + target_local
        query_local = 'query' + (14 - len(str(query_start))) * ' ' + str(query_start) + ' ' + query_local
        align_local = (19 - len(str(align_start))) * ' ' + str(align_start) + ' ' + align_local
        if k + 1 < math.ceil(len(align_mark) / 60):
            target_local += '\n'
            query_local += '\n\n'
            align_local += '\n'
            target_start += (60 - target_local.count('-'))
            query_start += (60 - query_local.count('-'))
            align_start += 60
        else:
            target_local += (' ' + str(original_j) + '\n')
            align_local += (' ' + str(len(align_mark)) + '\n')
            query_local += (' ' + str(original_i) + '\n' + '\n')
        align.append(target_local)
        align.append(align_local)
        align.append(query_local)
    align = ''.join(align)
    left = i
    right = len(seq) - original_i
    return insert, delete, mismatch, align, left, right, max_score
