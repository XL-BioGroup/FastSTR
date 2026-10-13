# cython: language_level=3
"""Order-preserving implementation of utils.get_anchors_python.

Keep base traversal A,G,C,T, then increasing right endpoint and increasing
left endpoint. Identical set insertion order preserves downstream DBSCAN
border ownership. Marks flag only the left endpoint of each equal-base pair.
No density thresholds or caller-visible data structures are changed.
"""
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def get_anchors(str sequence):
    cdef Py_ssize_t size = len(sequence)
    cdef Py_ssize_t p, i, j, k, left, right, distance
    cdef Py_UCS4 base
    cdef unsigned int bit
    cdef unsigned char[:] seen = bytearray(size)
    cdef object left_object, right_object
    cdef list buckets = [[], [], [], []]
    cdef list positions, marks
    cdef list sets = [None] + [set() for _ in range(8)]
    cdef list mark_lists = [None] + [[0] * size for _ in range(8)]
    cdef set anchors
    for p in range(size):
        base = sequence[p]
        if base == 65:  # A
            buckets[0].append(p)
        elif base == 71:  # G
            buckets[1].append(p)
        elif base == 67:  # C
            buckets[2].append(p)
        elif base == 84:  # T
            buckets[3].append(p)

    for positions in buckets:
        i = 0
        for j in range(len(positions)):
            right_object = positions[j]
            right = right_object
            while right - <Py_ssize_t>positions[i] > 8:
                i += 1
            for k in range(i, j):
                left_object = positions[k]
                left = left_object
                distance = right - left
                bit = 1 << (distance - 1)
                # Duplicate set.add calls cannot change iteration order.
                # Retain the first insertion in exactly the original order,
                # and reuse position int objects instead of boxing again.
                if not (seen[left] & bit):
                    anchors = sets[distance]
                    anchors.add((left_object, distance))
                    seen[left] = seen[left] | bit
                if not (seen[right] & bit):
                    anchors = sets[distance]
                    anchors.add((right_object, distance))
                    seen[right] = seen[right] | bit
                mark_lists[distance][left] = 1
    return ({p: sets[p] for p in range(1, 9)},
            {p: mark_lists[p] for p in range(1, 9)})
