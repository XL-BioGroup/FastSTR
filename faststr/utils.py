import math
import numpy as np
import regex
from Bio import SeqIO


# 获取每条染色体从start_index到end_index上的碱基序列
def read_fasta(path, start_index=1, end_index=0):
    sequences = []
    with open(path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if end_index == 0:
                record.seq = record.seq.upper()[start_index - 1:]
            else:
                record.seq = record.seq.upper()[start_index - 1:end_index]

            sequences.append(record)
    return sequences


def make_sub_reads(gene_sequence, read_length, overlap_length):
    sub_reads = []
    if len(gene_sequence) <= read_length:
        return [gene_sequence.seq]
    # 滑动窗口提取sub_read
    for start in range(0, (len(gene_sequence) // (read_length - overlap_length) + 1) * (read_length - overlap_length),
                       read_length - overlap_length):
        if start >= len(gene_sequence):
            break
        end = start + read_length
        subsequence = gene_sequence.seq[start:min(end, len(gene_sequence))]
        sub_reads.append(subsequence)

    return sub_reads


def get_anchors(sub_read):
    anchor_dict = {i: set() for i in range(1, 9)}
    anchor_mark_dict = {i: [0] * len(sub_read) for i in range(1, 9)}
    same_bases = [[], [], [], []]
    for p, c in enumerate(sub_read):
        if c == 'A':
            same_bases[0].append(p)
        elif c == 'G':
            same_bases[1].append(p)
        elif c == 'C':
            same_bases[2].append(p)
        elif c == 'T':
            same_bases[3].append(p)
    for s_b in same_bases:
        i = 0
        for j in range(len(s_b)):
            while s_b[j] - s_b[i] > 8:
                i += 1
            for k in range(i, j):
                anchor_dict[s_b[j] - s_b[k]].add((s_b[k], s_b[j] - s_b[k]))
                anchor_dict[s_b[j] - s_b[k]].add((s_b[j], s_b[j] - s_b[k]))
                anchor_mark_dict[s_b[j] - s_b[k]][s_b[k]] = 1
    return anchor_dict, anchor_mark_dict


def if_motif_is_tr(motif):
    if not motif:
        return False
    # 拼接字符串并移除首尾字符
    doubled_motif = (motif + motif)[1:-1]
    # 检查原字符串是否出现在移除首尾字符后的新字符串中
    return motif in doubled_motif


def tri_gram_model(sequence, n):
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

    # 根据种子获取备选motif
    sorted_seed = sorted(bi_count, key=lambda x: bi_count[x], reverse=True)
    if bi_count[sorted_seed[0]] < 2:
        return None
    # 如果motif长度为1
    if n == 1:
        return [sorted_seed[0][0]]
    # 如果motif长度为2
    if n == 2:
        i = 0
        while i < len(sorted_seed):
            if sorted_seed[i][0] == sorted_seed[i][1] and bi_count[sorted_seed[i]] >= len(sequence) / 2:
                return None
            if sorted_seed[i][0] != sorted_seed[i][1] and bi_count[sorted_seed[i]] > 0:
                return [sorted_seed[i]]
            i += 1
        return None
    # 如果motif长度大于2小于7
    motifs = set()
    # true_motif = []
    k = 0
    while True:
        max_seed = sorted_seed[k]
        bipos_set = set(bi_pos[max_seed])

        for pos in bi_pos[max_seed]:
            if pos + n in bipos_set:
                if if_motif_is_tr(sequence[pos:pos + n]):
                    continue
                if not set(sequence[pos:pos + n]).issubset({'A', 'G', 'T', 'C'}):
                    continue
                # if sequence[pos:pos + n] in m + m:
                #     true_motif.append(sequence[pos:pos + n])
                motifs.add(sequence[pos:pos + n])
        if motifs and k == 2:
            break
        k += 1
        if k > 15:
            return None

    # 计算各个motif的合理概率，得到最佳motif
    tri_gram = {}
    for tri, _ in tri_count.items():
        if bi_count[tri[:2]] == 0:
            continue
        tri_gram[tri] = tri_count[tri] / bi_count[tri[:2]]
    p_rationality = []
    for mot in motifs:
        p_r = 1
        for b in range(n - 2):
            p_r += math.log(tri_gram[mot[b:b + 3]])
        p_rationality.append(p_r)
    sorted_prationalty = sorted(p_rationality)
    max_prationalty = sorted_prationalty[-1]
    motifs = list(motifs)
    alternative_motif = [str(motifs[index]) for index, value in enumerate(p_rationality) if value == max_prationalty]
    final_motif = [alternative_motif[0]]
    if len(alternative_motif) == 1:
        return final_motif
    for a_m in alternative_motif[1:]:
        logo = 0
        for f_m in final_motif:
            if a_m in f_m + f_m:
                logo = 1
        if logo == 0:
            final_motif.append(a_m)
    return final_motif


def get_realscan_sequence(sequence):
    if len(sequence) < 40:
        return sequence
    elif len(sequence) < 800:
        noisy_bases = round(len(sequence) * 0.05)
        start_noisy_bases = noisy_bases // 2
        end_noisy_bases = noisy_bases - start_noisy_bases
        return sequence[start_noisy_bases:len(sequence) - end_noisy_bases]
    else:
        return sequence[20:len(sequence) - 20]


def get_motif_marks(sequence, motif):
    if len(motif) < 4:
        pattern = f"({motif}){{s<=0}}"
    elif len(motif) < 6:
        pattern = f"({motif}){{s<=1}}"
    else:
        pattern = f"({motif}){{s<=2}}"
    matches = regex.finditer(pattern, str(sequence))
    positions = [match.start() for match in matches]
    return positions


def  high_motif_detection(sequence, n):
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

    # 根据种子获取备选motif
    sorted_seed = sorted(bi_count, key=lambda x: bi_count[x], reverse=True)
    if bi_count[sorted_seed[0]] < 2:
        return None

    motifs = set()
    k = 0
    while True:
        max_seed = sorted_seed[k]
        bipos_set = set(bi_pos[max_seed])

        for pos in bi_pos[max_seed]:
            if pos + n in bipos_set:
                if if_motif_is_tr(sequence[pos:pos + n]):
                    continue
                if not set(sequence[pos:pos + n]).issubset({'A', 'G', 'T', 'C'}):
                    continue
                motifs.add(sequence[pos:pos + n])
        if motifs and k == 2:
            break
        k += 1
        if k > 15:
            return None

    # 计算各个motif的合理概率，得到最佳motif
    tri_gram = {}
    for tri, _ in tri_count.items():
        if bi_count[tri[:2]] == 0:
            continue
        tri_gram[tri] = tri_count[tri] / bi_count[tri[:2]]
    p_rationality = []
    for mot in motifs:
        p_r = 1
        for b in range(n - 2):
            p_r += math.log(tri_gram[mot[b:b + 3]])
        p_rationality.append(p_r)
    sorted_prationalty = sorted(p_rationality)
    max_prationalty = sorted_prationalty[-1]
    motifs_mar = []
    motifs = list(motifs)
    arfa = 0
    if n == 5:
        arfa = 0.3
    elif n == 6:
        arfa = 0.4
    elif n == 7:
        arfa = 0.8
    else:
        arfa = 1
    for inde, mot in enumerate(motifs):
        if max_prationalty - p_rationality[inde] <= arfa:
            motifs_mar.append(mot)

    # copy_vector_variance = []
    # # copy = round(len(sequence) / n)
    # copy = estimated_copy_number(sequence, n, base_de_bruijn)
    # for mot in motifs_mar:
    #     mot_copy_vector_variance = get_copy_vector(mot, n, copy, base_de_bruijn, base_mark)
    #     if mot_copy_vector_variance != -1:
    #         copy_vector_variance.append(mot_copy_vector_variance)
    #
    # if copy_vector_variance == []:
    #     return None
    # min_variance = min(copy_vector_variance)
    # alternative_motif = [str(motifs_mar[index]) for index, value in enumerate(copy_vector_variance) if
    #                      value == min_variance]
    #
    # final_motif = [alternative_motif[0]]
    # if len(alternative_motif) == 1:
    #     return final_motif
    # for a_m in alternative_motif[1:]:
    #     logo = 0
    #     for f_m in final_motif:
    #         if a_m in f_m + f_m:
    #             logo = 1
    #     if logo == 0:
    #         final_motif.append(a_m)

    copy = estimated_copy_number(sequence, n, base_de_bruijn)
    Frequency = []
    for mot in motifs_mar:
        Frequency.append(calculate_base_frequencye(mot, copy, base_de_bruijn, base_mark))

    max_frequency = max(Frequency)
    alternative_motif = []
    for index, value in enumerate(Frequency):
        if value == max_frequency:
            alternative_motif.append(str(motifs_mar[index]))

    final_motif = [alternative_motif[0]]
    if len(alternative_motif) == 1:
        return final_motif
    for a_m in alternative_motif[1:]:
        logo = 0
        for f_m in final_motif:
            if a_m in f_m + f_m:
                logo = 1
        if logo == 0:
            final_motif.append(a_m)
    if len(final_motif) > 3:
        copy = len(sequence) // n
        Frequency = []
        for mot in motifs_mar:
            Frequency.append(calculate_base_frequencye(mot, copy, base_de_bruijn, base_mark))

        max_frequency = max(Frequency)
        alternative_motif = []
        for index, value in enumerate(Frequency):
            if value == max_frequency:
                alternative_motif.append(str(motifs_mar[index]))

        final_motif = [alternative_motif[0]]
        if len(alternative_motif) == 1:
            return final_motif
        for a_m in alternative_motif[1:]:
            logo = 0
            for f_m in final_motif:
                if a_m in f_m + f_m:
                    logo = 1
            if logo == 0:
                final_motif.append(a_m)

    return final_motif


def estimated_copy_number(sequence, n, base_de_bruijn):
    copy = len(sequence) // n
    # 展开数组为一维
    flattened = base_de_bruijn.flatten()
    # 获取最大的前8个数
    largest_eight = np.sort(flattened)[-8:][::-1]
    large = largest_eight[0]
    for i in range(1, n):
        if i == n - 1 or largest_eight[i] < copy:
            small = math.ceil(sum(largest_eight[:n]) / n)
            break
        if large - largest_eight[i] >= copy:
            small = largest_eight[i]
            break

    common_divisor = large - small
    max_common_divisor = min(common_divisor, small)
    while max_common_divisor >= copy:
        if small != max_common_divisor:
            common_divisor = max(max_common_divisor, small) - min(small, max_common_divisor)
            small = min(small, max_common_divisor)
            max_common_divisor = min(common_divisor, small)
        else:
            if max_common_divisor % copy == 0:
                max_common_divisor = copy
                break
            max_common_divisor = math.ceil(max_common_divisor / round(max_common_divisor / copy))
            break
    return math.ceil(max_common_divisor * 0.5 + copy * 0.5)


def calculate_base_frequencye(motif, copy, base_de_bruijn, base_mark):
    base_de_bruijn = base_de_bruijn.copy()
    base_frequency = 0
    for i in range(1, len(motif)):
        if base_de_bruijn[base_mark[motif[i - 1]]][base_mark[motif[i]]] >= copy:
            base_frequency += copy
            base_de_bruijn[base_mark[motif[i - 1]]][base_mark[motif[i]]] -= copy
        else:
            if base_de_bruijn[base_mark[motif[i - 1]]][base_mark[motif[i]]] >= 0.6 * copy:
                base_frequency += base_de_bruijn[base_mark[motif[i - 1]]][base_mark[motif[i]]]
            else:
                base_frequency += 0
            if base_de_bruijn[base_mark[motif[i - 1]]][base_mark[motif[i]]] > 0:
                base_de_bruijn[base_mark[motif[i - 1]]][base_mark[motif[i]]] = 0

    return base_frequency


def get_copy_vector(motif, n, copy, base_de_bruijn, base_mark):
    motif_base_de_bruijn = np.zeros((4, 4))
    for i in range(1, n):
        motif_base_de_bruijn[base_mark[motif[i - 1]]][base_mark[motif[i]]] += 1

    copy_vector = []
    for r in range(4):
        for l in range(4):
            if motif_base_de_bruijn[r][l] > 0:
                copy_vector.append(base_de_bruijn[r][l] / motif_base_de_bruijn[r][l])

    return np.mean((np.array(copy_vector) - copy) ** 2)


def find_most_frequent_substring_given(long_seq, four_mer_marker, k):
    freq_map = {}  # 字符串 -> 频数
    positions_map = {}  # 字符串 -> 出现位置列表
    max_count = 0
    most_frequent_string = ""

    # Convert short_strings list to a set for O(1) lookup
    short_strings_set = set(four_mer_marker)

    # Traverse all substrings of length k using sliding window
    for i in range(len(long_seq) - k + 1):
        substring = long_seq[i:i + k]

        # Only count the substring if it is in the set of given short strings
        if substring in short_strings_set:
            # Update frequency map
            if substring in freq_map and i - positions_map[substring][-1] >= k:
                freq_map[substring] += 1
                positions_map[substring].append(i)
            elif substring not in freq_map:
                freq_map[substring] = 1
                positions_map[substring] = [i]

            # Update the most frequent substring if needed
            if freq_map[substring] > max_count:
                max_count = freq_map[substring]
                most_frequent_string = substring

    # Return the most frequent substring and its positions
    return most_frequent_string, positions_map[most_frequent_string]
