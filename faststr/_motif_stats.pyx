# cython: language_level=3
"""Integer-coded counts, preserving original key/position order and types."""
import numpy as np

# Immutable strings shared across calls; each result still owns its dicts
# and position lists. Preserve AGCT insertion order and key types.
_PAIR_KEYS = tuple(x + y for x in 'AGCT' for y in 'AGCT')

def scan_stats(sequence, bint graph):
    cdef str text = str(sequence)
    cdef Py_ssize_t length = len(text)
    cdef unsigned char[:] codes = bytearray(length)
    cdef Py_ssize_t bi[16]
    cdef Py_ssize_t tri[64]
    cdef Py_ssize_t first[64]
    cdef Py_ssize_t edges[16]
    cdef Py_ssize_t i, code, pair, pos
    cdef int a, b, c, previous
    cdef Py_UCS4 ch
    cdef list order = []
    cdef list positions = [[] for _ in range(16)]
    cdef tuple keys = _PAIR_KEYS
    for i in range(16):
        bi[i] = 0
        edges[i] = 0
    for i in range(64):
        tri[i] = 0
    for i in range(length):
        ch = text[i]
        if ch == 65: codes[i] = 0
        elif ch == 71: codes[i] = 1
        elif ch == 67: codes[i] = 2
        elif ch == 84: codes[i] = 3
        else: codes[i] = 4
    for i in range(length - 2):
        a, b, c = codes[i], codes[i + 1], codes[i + 2]
        if graph and i > 0:
            previous = codes[i - 1]
            if previous < 4 and a < 4:
                edges[previous * 4 + a] += 1
        if a == 4 or b == 4 or c == 4:
            continue
        code = a * 16 + b * 4 + c
        if tri[code] == 0:
            first[code] = i
            order.append(code)
        tri[code] += 1
        pair = a * 4 + b
        positions[pair].append(i)
        bi[pair] += 1
        # Preserve the original conditional final-pair handling exactly.
        if i == length - 3:
            pair = b * 4 + c
            positions[pair].append(i + 1)
            bi[pair] += 1
            if graph:
                edges[a * 4 + b] += 1
                edges[b * 4 + c] += 1
    tri_count = {}
    for code in order:
        pos = first[code]
        # Slice the original object: str stays str, Bio.Seq stays Bio.Seq.
        tri_count[sequence[pos:pos + 3]] = tri[code]
    bi_count = {keys[i]: bi[i] for i in range(16)}
    bi_pos = {keys[i]: positions[i] for i in range(16)}
    base_graph = np.array([edges[i] for i in range(16)], dtype=float).reshape(4, 4) if graph else None
    return tri_count, bi_count, bi_pos, base_graph
