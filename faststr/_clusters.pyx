# cython: language_level=3
"""One-pass label statistics and compiled interval loops; thresholds unchanged.

Group labels stay in ascending order. DBSCAN still uses original anchor input
order. Early breaks, skipped first windows and every-fifth variance samples
are intentionally retained exactly. NumPy variance arithmetic is unchanged.
"""
import numpy as np
from .dbscan_1d import DBSCAN
from ._dbscan import labels_only

_INTP64 = np.dtype(np.intp).itemsize == 8

def _anchor_array(anchors):
    """Typed fill for built-in integer pairs, with the former NumPy fallback.

    Row order is unchanged. Representation guards are not detection gates.
    """
    cdef Py_ssize_t row, n = len(anchors)
    cdef Py_ssize_t[:, :] view
    cdef object item, x, period
    if type(anchors) is not list or n == 0 or not _INTP64:
        return np.asarray(anchors)
    array = np.empty((n, 2), dtype=np.intp)
    view = array
    for row in range(n):
        item = anchors[row]
        if type(item) is not tuple or len(item) != 2:
            return np.asarray(anchors)
        x, period = item
        if (type(x) is not int or type(period) is not int or
                x < -1099511627776 or x > 1099511627776 or
                period < 1 or period > 8):
            return np.asarray(anchors)
        view[row, 0] = x
        view[row, 1] = period
    return array

def cluster_anchors(anchors, anchor_mark_list, indel, cover):
    cdef Py_ssize_t row, label, groups, n
    cdef Py_ssize_t[:, :] points_view
    cdef Py_ssize_t[:] labels_view, counts_view, first_view, last_view
    cdef Py_ssize_t point, i, current_density
    cdef Py_ssize_t start_index, end_index, fuzzy_tr_start_index, fuzzy_tr_end_index
    cdef Py_ssize_t anchor_dis, logo
    cdef Py_ssize_t start_limit, end_limit, rough_limit, start_pad, end_pad, border_limit
    cdef list cluster_marks
    cdef Py_ssize_t[:] q_view
    cdef Py_ssize_t q_index
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
    # Cache unchanged integer thresholds; avoid Python boxing per position.
    start_limit = start_rought_judgment[anchor_dis]
    end_limit = end_rought_judgment[anchor_dis]
    rough_limit = difinition_rough_judgment[anchor_dis]
    start_pad = start_expand[anchor_dis]
    end_pad = end_expand[anchor_dis]
    border_limit = fuzzy_tr_anchor_boundary_density[anchor_dis]
    anchors = _anchor_array(anchors)
    labels = labels_only(anchors, eps=14, min_samples=14)
    # 获取长度达到25,覆盖度达到(1 - indel) * cover的簇
    # Convert only if needed for native indexing; the common integer array
    # already has np.intp dtype. DBSCAN validates the original input first.
    points_view = np.asarray(anchors, dtype=np.intp)
    labels_view = labels
    n = len(labels)
    groups = int(labels.max()) + 1
    group_stats = np.empty((3, groups), dtype=np.intp)
    counts, firsts, lasts = group_stats
    counts.fill(0)
    counts_view = counts
    first_view = firsts
    last_view = lasts
    for row in range(n):
        label = labels_view[row]
        if label < 0:
            continue
        point = points_view[row, 0]
        if counts_view[label] == 0:
            first_view[label] = point
            last_view[label] = point
        else:
            if point < first_view[label]:
                first_view[label] = point
            if point > last_view[label]:
                last_view[label] = point
        counts_view[label] += 1
    for label in range(groups):
        if counts_view[label] == 0:
            continue
        # Keep NumPy scalar coordinate arithmetic and density division as
        # in the original, and retain the identical NumPy variance formula.
        cluster_first = firsts[label]
        cluster_last = lasts[label]
        l = cluster_last - cluster_first + 1
        anchor_density = int(counts_view[label]) / l
        cluster_marks = anchor_mark_list[cluster_first:cluster_last + 1]
        if l >= 25 and anchor_density >= min(0.81, (1 - indel) * cover):
            start_index = -1
            end_index = -1
            logo = 0
            current_density = sum(cluster_marks[:25])
            if logo == 0:
                for i in range(1, len(cluster_marks) - 25 + 1):
                    # 更新窗口内的区密度
                    current_density = current_density - <Py_ssize_t>cluster_marks[i - 1] + \
                                      <Py_ssize_t>cluster_marks[i + 25 - 1]
                    if start_index == -1:
                        if current_density >= start_limit:
                            start_index = max(i - start_pad, 0)
                    if current_density >= rough_limit:
                        logo = 1
                        break
            if logo == 0:
                continue
            current_density = sum(cluster_marks[-25:])
            if current_density >= end_limit:
                end_index = len(cluster_marks) - 1
            else:
                for i in range(len(cluster_marks) - 2, start_index + 23, -1):
                    current_density = current_density - <Py_ssize_t>cluster_marks[i + 1] + \
                                      <Py_ssize_t>cluster_marks[i - 25 + 1]
                    if current_density >= end_limit:
                        end_index = min(i + end_pad,
                                        len(cluster_marks) - 1)
                        break
            fuzzy_repeating_region.append(
                (cluster_first, cluster_last, cluster_first + start_index, cluster_first + end_index, 1))
        elif l >= 25:
            current_density = sum(cluster_marks[:25])
            start_index = -1
            end_index = -1
            fuzzy_tr_start_index = -1
            fuzzy_tr_end_index = -1
            logo = 0
            if current_density >= start_limit:
                start_index = 0
            if current_density >= border_limit:
                fuzzy_tr_start_index = 0
            if current_density >= rough_limit:
                logo = 1
            if current_density >= end_limit:
                end_index = min(29, len(cluster_marks) - 1)
            if current_density >= border_limit:
                fuzzy_tr_end_index = len(cluster_marks) - 1
            # Same samples and dtype as np.array(the former integer list).
            # Keep NumPy subtraction, squaring, reduction and division intact.
            Q_support_domain_mod5 = np.empty(max(0, len(cluster_marks) - 25) // 5 + 1, dtype=np.intp)
            q_view = Q_support_domain_mod5
            q_view[0] = current_density
            q_index = 1
            for i in range(1, len(cluster_marks) - 25 + 1):
                current_density = current_density - <Py_ssize_t>cluster_marks[i - 1] + \
                                  <Py_ssize_t>cluster_marks[i + 25 - 1]
                if start_index == -1:
                    if current_density >= start_limit:
                        start_index = max(i - start_pad, 0)
                if fuzzy_tr_start_index == -1:
                    if current_density >= border_limit:
                        fuzzy_tr_start_index = i
                if logo == 0 and current_density >= rough_limit:
                    logo = 1
                if current_density >= end_limit:
                    end_index = min(i + 24 + end_pad,
                                    len(cluster_marks) - 1)
                if current_density >= border_limit:
                    fuzzy_tr_end_index = i + 24
                if i % 5 == 0:
                    q_view[q_index] = current_density
                    q_index += 1
            Q_variance = np.mean((Q_support_domain_mod5 - fuzzy_tr_structure_feature[anchor_dis]) ** 2)
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
            #     for i in range(1, len(cluster_marks) - 25 + 1):
            #         # 更新窗口内的区密度
            #         current_density = current_density - <Py_ssize_t>cluster_marks[i - 1] + \
            #                           <Py_ssize_t>cluster_marks[i + 25 - 1]
            #         if start_index == -1:
            #             if current_density >= start_limit:
            #                 start_index = max(i - start_pad, 0)
            #         if current_density >= rough_limit:
            #             logo = 1
            #             break
            # if logo == 0:
            #     continue
            # current_density = sum(cluster_marks[-25:])
            # if current_density >= end_limit:
            #     end_index = len(cluster_marks) - 1
            # else:
            #     for i in range(len(cluster_marks) - 2, start_index + 23, -1):
            #         current_density = current_density - <Py_ssize_t>cluster_marks[i + 1] + \
            #                           <Py_ssize_t>cluster_marks[i - 25 + 1]
            #         if current_density >= end_limit:
            #             end_index = min(i + end_pad,
            #                             len(cluster_marks) - 1)
            #             break

    return fuzzy_repeating_region
