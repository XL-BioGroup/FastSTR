from . import utils
from . import scan_subread
from . import trs_align



def get_subread_trs(sub_read, m, r, g, e, p_indel, p_match, score, beta):
    # 获取anchors
    anchor_dict, anchor_mark_dict = utils.get_anchors(sub_read)
    # 获取模糊重复区间
    fuzzy_repeats_dict = {i: [] for i in range(1, 9)}
    num = 0
    for n, n_anchors in anchor_dict.items():
        fuzzy_repeats_dict[n] = scan_subread.cluster_anchors(list(n_anchors), anchor_mark_dict[n], p_indel, p_match)
    del anchor_dict
    # 获取模糊重复区间的motif及大致比对区间
    candidate_trs_dict = {i: [] for i in range(1, 9)}
    for n, n_fuzzy_repeats in fuzzy_repeats_dict.items():
        candidate_trs_dict[n] = scan_subread.filter_chaotic_repeats(n_fuzzy_repeats, sub_read, n)
    del fuzzy_repeats_dict
    # 合并中断tr序列,比对备选序列，得到合格序列,得到一条sub_read上的最终tr集
    qualified_trs_list = []
    # num = 0
    for n, n_candidate_trs in candidate_trs_dict.items():
        qualified_trs_list.extend(trs_align.trs_align_algorithm(n_candidate_trs, sub_read, p_indel, p_match, m, r, g, e,
                                                                score, beta))
    del candidate_trs_dict
    return sorted(qualified_trs_list, key=lambda x: x[1])
