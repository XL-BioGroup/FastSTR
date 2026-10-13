"""Exact original statistic scanners used when the native module is absent."""
import numpy as np
try:
    from ._motif_stats import scan_stats as _native_stats
except ImportError:
    _native_stats = None

def tri_stats_python(sequence):
    # 获取二元、三元计数和二元位置
    tri_count = {}
    bi_count = {'AA': 0, 'AG': 0, 'AC': 0, 'AT': 0, 'GA': 0, 'GG': 0, 'GC': 0, 'GT': 0, 'CA': 0, 'CG': 0, 'CC': 0,
                'CT': 0, 'TA': 0, 'TG': 0, 'TC': 0, 'TT': 0}
    bi_pos = {'AA': [], 'AG': [], 'AC': [], 'AT': [], 'GA': [], 'GG': [], 'GC': [], 'GT': [], 'CA': [], 'CG': [],
              'CC': [], 'CT': [], 'TA': [], 'TG': [], 'TC': [], 'TT': []}
    for i in range(len(sequence) - 2):
        if not set(sequence[i:i + 3]).issubset({'A', 'G', 'T', 'C'}):
            continue
        if sequence[i:i + 3] not in tri_count:
            tri_count[sequence[i:i + 3]] = 1
        else:
            tri_count[sequence[i:i + 3]] += 1
        bi_pos[sequence[i:i + 2]].append(i)
        bi_count[sequence[i:i + 2]] += 1
        if i == len(sequence) - 3:
            bi_pos[sequence[i + 1:i + 3]].append(i + 1)
            bi_count[sequence[i + 1:i + 3]] += 1

    return tri_count, bi_count, bi_pos, None

def high_stats_python(sequence):
    # 获取二元、三元计数和二元位置
    tri_count = {}
    bi_count = {'AA': 0, 'AG': 0, 'AC': 0, 'AT': 0, 'GA': 0, 'GG': 0, 'GC': 0, 'GT': 0, 'CA': 0, 'CG': 0, 'CC': 0,
                'CT': 0, 'TA': 0, 'TG': 0, 'TC': 0, 'TT': 0}
    bi_pos = {'AA': [], 'AG': [], 'AC': [], 'AT': [], 'GA': [], 'GG': [], 'GC': [], 'GT': [], 'CA': [], 'CG': [],
              'CC': [], 'CT': [], 'TA': [], 'TG': [], 'TC': [], 'TT': []}
    base_de_bruijn = np.zeros((4, 4))
    base_mark = {'A': 0, 'G': 1, 'C': 2, 'T': 3}
    base_annotation = {0: 'A', 1: 'G', 2: 'C', 3: 'T'}
    for i in range(len(sequence) - 2):
        if i > 0 and set(sequence[i - 1:i + 1]).issubset({'A', 'G', 'T', 'C'}):
            base_de_bruijn[base_mark[sequence[i - 1]]][base_mark[sequence[i]]] += 1
        if not set(sequence[i:i + 3]).issubset({'A', 'G', 'T', 'C'}):
            continue
        if sequence[i:i + 3] not in tri_count:
            tri_count[sequence[i:i + 3]] = 1
        else:
            tri_count[sequence[i:i + 3]] += 1
        bi_pos[sequence[i:i + 2]].append(i)
        bi_count[sequence[i:i + 2]] += 1
        if i == len(sequence) - 3:
            bi_pos[sequence[i + 1:i + 3]].append(i + 1)
            bi_count[sequence[i + 1:i + 3]] += 1
            base_de_bruijn[base_mark[sequence[i]]][base_mark[sequence[i + 1]]] += 1
            base_de_bruijn[base_mark[sequence[i + 1]]][base_mark[sequence[i + 2]]] += 1

    return tri_count, bi_count, bi_pos, base_de_bruijn

def scan_stats(sequence, graph):
    if _native_stats is not None:
        return _native_stats(sequence, graph)
    return high_stats_python(sequence) if graph else tri_stats_python(sequence)
