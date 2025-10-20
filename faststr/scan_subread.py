import math
import regex
from sklearn.cluster import DBSCAN
import numpy as np
from Bio.Align import PairwiseAligner

from . import utils

def cluster_anchors(anchors, anchor_mark_list, indel, cover):
    if anchors == []:
        return []
    fuzzy_repeating_region = []
    difinition_rough_judgment = [19, 19, 19, 20, 18, 18, 17, 17]
    start_rought_judgment = [17, 18, 18, 18, 16, 16, 15, 15]
    end_rought_judgment = [17, 18, 18, 17, 15, 15, 14, 14]
    start_expand = [5, 5, 5, 5, 5, 5, 10, 10]
    end_expand = [5, 5, 5, 5, 5, 5, 10, 10]
    fuzzy_tr_structure_feature = [12, 12, 12, 12, 12, 12, 14, 12]
    fuzzy_tr_structure_threshold = [0, 8, 0, 0, 7.5, 0, 12.5, 0]
    fuzzy_tr_anchor_density = [0, 0.68, 0, 0, 0.64, 0, 0.68, 0]
    fuzzy_tr_anchor_boundary_density = [0, 14, 0, 0, 12, 0, 13, 0]
    anchor_dis = anchors[0][1] - 1
    dbscan_model = DBSCAN(eps=14, min_samples=14, n_jobs=1).fit(anchors)
    # 获取聚类标签
    labels = dbscan_model.labels_
    # 获取长度达到25,覆盖度达到(1 - indel) * cover的簇
    unique_labels = np.unique(labels)
    anchors = np.array(anchors)
    for label in unique_labels:
        if label == -1:  # 忽略噪声点
            continue
        cluster_indices = np.where(labels == label)[0]
        cluster_points = anchors[cluster_indices]
        cluster_first = min(item[0] for item in cluster_points)
        cluster_last = max(item[0] for item in cluster_points)
        l = cluster_last - cluster_first + 1
        anchor_density = len(cluster_points) / l
        if l >= 25 and anchor_density >= min(0.81, (1 - indel) * cover):
            start_index = -1
            end_index = -1
            logo = 0
            current_density = sum(anchor_mark_list[cluster_first:cluster_last + 1][:25])
            if logo == 0:
                for i in range(1, len(anchor_mark_list[cluster_first:cluster_last + 1]) - 25 + 1):
                    # 更新窗口内的区密度
                    current_density = current_density - anchor_mark_list[cluster_first:cluster_last + 1][i - 1] + \
                                      anchor_mark_list[cluster_first:cluster_last + 1][i + 25 - 1]
                    if start_index == -1:
                        if current_density >= start_rought_judgment[anchor_dis]:
                            start_index = max(i - start_expand[anchor_dis], 0)
                    if current_density >= difinition_rough_judgment[anchor_dis]:
                        logo = 1
                        break
            if logo == 0:
                continue
            current_density = sum(anchor_mark_list[cluster_first:cluster_last + 1][-25:])
            if current_density >= end_rought_judgment[anchor_dis]:
                end_index = len(anchor_mark_list[cluster_first:cluster_last + 1]) - 1
            else:
                for i in range(len(anchor_mark_list[cluster_first:cluster_last + 1]) - 2, start_index + 23, -1):
                    current_density = current_density - anchor_mark_list[cluster_first:cluster_last + 1][i + 1] + \
                                      anchor_mark_list[cluster_first:cluster_last + 1][i - 25 + 1]
                    if current_density >= end_rought_judgment[anchor_dis]:
                        end_index = min(i + end_expand[anchor_dis],
                                        len(anchor_mark_list[cluster_first:cluster_last + 1]) - 1)
                        break
            fuzzy_repeating_region.append(
                (cluster_first, cluster_last, cluster_first + start_index, cluster_first + end_index, 1))
        elif l >= 25:
            current_density = sum(anchor_mark_list[cluster_first:cluster_last + 1][:25])
            start_index = -1
            end_index = -1
            fuzzy_tr_start_index = -1
            fuzzy_tr_end_index = -1
            logo = 0
            if current_density >= start_rought_judgment[anchor_dis]:
                start_index = 0
            if current_density >= fuzzy_tr_anchor_boundary_density[anchor_dis]:
                fuzzy_tr_start_index = 0
            if current_density >= difinition_rough_judgment[anchor_dis]:
                logo = 1
            if current_density >= end_rought_judgment[anchor_dis]:
                end_index = min(29, len(anchor_mark_list[cluster_first:cluster_last + 1]) - 1)
            if current_density >= fuzzy_tr_anchor_boundary_density[anchor_dis]:
                fuzzy_tr_end_index = len(anchor_mark_list[cluster_first:cluster_last + 1]) - 1
            Q_support_domain_mod5 = [current_density]
            for i in range(1, len(anchor_mark_list[cluster_first:cluster_last + 1]) - 25 + 1):
                current_density = current_density - anchor_mark_list[cluster_first:cluster_last + 1][i - 1] + \
                                  anchor_mark_list[cluster_first:cluster_last + 1][i + 25 - 1]
                if start_index == -1:
                    if current_density >= start_rought_judgment[anchor_dis]:
                        start_index = max(i - start_expand[anchor_dis], 0)
                if fuzzy_tr_start_index == -1:
                    if current_density >= fuzzy_tr_anchor_boundary_density[anchor_dis]:
                        fuzzy_tr_start_index = i
                if logo == 0 and current_density >= difinition_rough_judgment[anchor_dis]:
                    logo = 1
                if current_density >= end_rought_judgment[anchor_dis]:
                    end_index = min(i + 24 + end_expand[anchor_dis],
                                    len(anchor_mark_list[cluster_first:cluster_last + 1]) - 1)
                if current_density >= fuzzy_tr_anchor_boundary_density[anchor_dis]:
                    fuzzy_tr_end_index = i + 24
                if i % 5 == 0:
                    Q_support_domain_mod5.append(current_density)
            Q_variance = np.mean((np.array(Q_support_domain_mod5) - fuzzy_tr_structure_feature[anchor_dis]) ** 2)
            if Q_variance >= fuzzy_tr_structure_threshold[anchor_dis]:
                if logo == 1:
                    fuzzy_repeating_region.append(
                        (cluster_first, cluster_last, cluster_first + start_index, cluster_first + end_index, 1))
            else:
                '''
                干扰区间的过滤
                '''
                if anchor_density >= fuzzy_tr_anchor_density[anchor_dis]:
                    fuzzy_repeating_region.append(
                        (cluster_first, cluster_last, cluster_first + fuzzy_tr_start_index,
                         cluster_first + fuzzy_tr_end_index, 0))
            # if logo == 0:
            #     for i in range(1, len(anchor_mark_list[cluster_first:cluster_last + 1]) - 25 + 1):
            #         # 更新窗口内的区密度
            #         current_density = current_density - anchor_mark_list[cluster_first:cluster_last + 1][i - 1] + \
            #                           anchor_mark_list[cluster_first:cluster_last + 1][i + 25 - 1]
            #         if start_index == -1:
            #             if current_density >= start_rought_judgment[anchor_dis]:
            #                 start_index = max(i - start_expand[anchor_dis], 0)
            #         if current_density >= difinition_rough_judgment[anchor_dis]:
            #             logo = 1
            #             break
            # if logo == 0:
            #     continue
            # current_density = sum(anchor_mark_list[cluster_first:cluster_last + 1][-25:])
            # if current_density >= end_rought_judgment[anchor_dis]:
            #     end_index = len(anchor_mark_list[cluster_first:cluster_last + 1]) - 1
            # else:
            #     for i in range(len(anchor_mark_list[cluster_first:cluster_last + 1]) - 2, start_index + 23, -1):
            #         current_density = current_density - anchor_mark_list[cluster_first:cluster_last + 1][i + 1] + \
            #                           anchor_mark_list[cluster_first:cluster_last + 1][i - 25 + 1]
            #         if current_density >= end_rought_judgment[anchor_dis]:
            #             end_index = min(i + end_expand[anchor_dis],
            #                             len(anchor_mark_list[cluster_first:cluster_last + 1]) - 1)
            #             break

    return fuzzy_repeating_region


def filter_chaotic_repeats(fuzzy_repeats, sub_read, n):
    candidate_trs = []
    max_windows = [18, 18, 17.5, 17, 18, 18, 17, 17]
    start_probe = [24, 23, 23, 20, 22, 23, 23, 22]
    end_probe = [24, 23, 22, 20, 22, 23, 22, 21]
    consistent_density = [0, 1.8, 0, 0, 5, 0, 4.8, 0]
    perfect_count = [0, 25, 0, 0, 15, 0, 10, 0]
    to_check_large_probe_trs = []
    for f_r in fuzzy_repeats:
        to_scan_sequence = sub_read[f_r[2]:f_r[3] + 1]
        if n < 5:
            motifs = utils.tri_gram_model(to_scan_sequence, n)
        else:
            motifs = utils.high_motif_detection(to_scan_sequence, n)
        if motifs == None:
            continue
        if f_r[4] == 1:
            for mot in motifs:
                # 过滤掉明显不符合串联重复结构的序列————小探针粗过滤
                motif_mark_indexes = utils.get_motif_marks(to_scan_sequence, mot)
                visit = [0] * len(to_scan_sequence)
                for m_p in motif_mark_indexes:
                    visit[m_p:m_p + n] = [1] * n
                if len(mot) in [1, 4, 6]:
                    window_ones_count = sum(visit[:20])
                    max_ones_count = window_ones_count

                    for i in range(1, len(visit) - 19):
                        # 更新窗口中1的数量
                        window_ones_count += visit[i + 19] - visit[i - 1]
                        if window_ones_count > max_ones_count:
                            max_ones_count = window_ones_count
                elif len(mot) in [2, 3]:
                    pattern = f"({mot}){{s<=1}}"
                    matches = regex.finditer(pattern, str(to_scan_sequence))
                    positions = [match.start() for match in matches]
                    for m_p in positions:
                        for i in range(len(mot)):
                            visit[m_p + i] = visit[m_p + i] + (1 - visit[m_p + i]) * (len(mot) - 1) / len(mot)
                    windows = sum(visit[:20])
                    window_ones_count = math.floor(windows)
                    max_ones_count = window_ones_count

                    for i in range(1, len(visit) - 19):
                        # 更新窗口中1的数量
                        windows += visit[i + 19] - visit[i - 1]
                        window_ones_count = math.floor(windows)
                        if window_ones_count > max_ones_count:
                            max_ones_count = window_ones_count
                elif len(mot) == 5:
                    pattern = f"({mot}){{s<=2}}"
                    matches = regex.finditer(pattern, str(to_scan_sequence))
                    positions = [match.start() for match in matches]
                    for m_p in positions:
                        for i in range(len(mot)):
                            visit[m_p + i] = visit[m_p + i] + (1 - visit[m_p + i]) * (len(mot) - 2) / len(mot)
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
                    pattern = f"({mot}){{s<=3}}"
                    matches = regex.finditer(pattern, str(to_scan_sequence))
                    positions = [match.start() for match in matches]
                    for m_p in positions:
                        for i in range(len(mot)):
                            visit[m_p + i] = visit[m_p + i] + (1 - visit[m_p + i]) * (len(mot) - 3) / len(mot)
                    windows = sum(visit[:20])
                    window_ones_count = math.floor(windows)
                    max_ones_count = window_ones_count

                    for i in range(1, len(visit) - 19):
                        # 更新窗口中1的数量
                        windows += visit[i + 19] - visit[i - 1]
                        window_ones_count = math.floor(windows)
                        if window_ones_count > max_ones_count:
                            max_ones_count = window_ones_count

                if max_ones_count // 0.5 * 0.5 < max_windows[len(mot) - 1]:
                    continue
                else:
                    to_check_large_probe_trs.append((f_r[0], f_r[1], mot))
        else:
            for mot in motifs:
                # 标记映射
                if n == 2:
                    pattern = r"|".join(['.' + mot[0] + mot[1], mot[0] + mot[1] + '.'])
                    matches = regex.finditer(pattern, to_scan_sequence)
                    motif_mark_indexes = [match.start() for match in matches]
                elif n == 5:
                    four_mers = set()
                    for i in range(len(mot)):
                        four_mers.add((mot + mot)[i:i + 4])
                    _, motif_mark_indexes = utils.find_most_frequent_substring_given(to_scan_sequence, four_mers, 4)
                elif n == 7:
                    matches = regex.finditer(f"({mot}){{s<=2}}", to_scan_sequence)
                    motif_mark_indexes = [match.start() for match in matches]

                # 获取标记密度
                motif_mark_gap = []
                for ind, m in enumerate(motif_mark_indexes[1:]):
                    motif_mark_gap.append(max(0, m - motif_mark_indexes[ind] - n))

                # 获取TR区间
                fluctuations_sum = 0
                fuzzy_tr_start = -1
                should_start = 0
                for g_index, g in enumerate(motif_mark_gap[:-perfect_count[len(mot) - 1]]):
                    if g_index < should_start:
                        continue
                    if g >= 1 and fuzzy_tr_start == -1:
                        continue
                    if g < 1 and fuzzy_tr_start == -1:
                        # Catch = [max(catch_g - 1, 0) for catch_g in G[g_index:g_index + 15]]
                        fluctuations_sum = sum(motif_mark_gap[g_index:g_index + perfect_count[len(mot) - 1]])
                        if fluctuations_sum / perfect_count[len(mot) - 1] <= consistent_density[len(mot) - 1]:
                            fuzzy_tr_start = g_index
                    elif fuzzy_tr_start != -1:
                        fluctuations_sum = fluctuations_sum - max(motif_mark_gap[g_index - 1], 0) + max(
                            motif_mark_gap[g_index + perfect_count[len(mot) - 1] - 1], 0)
                        if fluctuations_sum / perfect_count[len(mot) - 1] <= consistent_density[len(mot) - 1]:
                            continue
                        else:
                            to_check_large_probe_trs.append((motif_mark_indexes[fuzzy_tr_start] + f_r[3],
                                                             motif_mark_indexes[
                                                                 g_index + perfect_count[len(mot) - 1] - 2] + f_r[3],
                                                             mot))
                            fuzzy_tr_start = -1
                            fluctuations_sum = 0
                            should_start = g_index + perfect_count[len(mot) - 1] - 1
                if fuzzy_tr_start != -1:
                    to_check_large_probe_trs.append(
                        (motif_mark_indexes[fuzzy_tr_start] + f_r[3], motif_mark_indexes[-1] + f_r[3], mot))

    for t_o_t in to_check_large_probe_trs:
        # 过滤掉明显不符合串联重复结构的序列————大探针细过滤，同时获取大致起始和结束位置
        probe = 25 // len(t_o_t[2]) * t_o_t[2] + t_o_t[2][:25 % len(t_o_t[2])]
        to_scan_sequence = sub_read[t_o_t[0]:t_o_t[1] + 1]
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        # 设置比对参数
        aligner.match_score = 2
        aligner.mismatch_score = -3
        aligner.open_gap_score = -5
        start_index = -1
        end_index = -1
        for i in range(len(to_scan_sequence) // 5 - 3):
            probe_align = aligner.align(probe, to_scan_sequence[i * 5:i * 5 + 20])
            if probe_align:
                pass
            else:
                continue
            if probe_align[0].score >= start_probe[len(t_o_t[2]) - 1]:
                start_index = max(i * 5 - 5, 0)
                break
        if start_index == -1:
            continue
        for j in range((len(to_scan_sequence) - start_index - 1) // 5 - 3):
            probe_align = aligner.align(probe, to_scan_sequence[
                                               len(to_scan_sequence) - j * 5 - 20:len(to_scan_sequence) - j * 5])
            if probe_align:
                pass
            else:
                continue
            if probe_align[0].score >= end_probe[len(t_o_t[2]) - 1]:
                end_index = min(len(to_scan_sequence) - j * 5 + 4, len(to_scan_sequence) - 1)
                break
        if end_index == -1:
            continue
        if end_index - start_index < 24:
            continue
        # 获取tr的大致边界
        candidate_trs.append((str(t_o_t[2]), t_o_t[0] + start_index, t_o_t[0] + end_index))

    return candidate_trs
