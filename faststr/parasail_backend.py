"""Small local-alignment API adapter for the FastSTR call sites.

DP is parasail SIMD SW; Bio.Alignment is used only for path representation,
counts and text formatting. One optimum is returned, not all tied paths.
"""
from functools import lru_cache
import re
import numpy as np
from Bio.Align import Alignment

import parasail


@lru_cache(maxsize=32)
def matrix(match, mismatch):
    return parasail.matrix_create('ACGTRYSWKMBDHVNacgtryswkmbdhvn', match, mismatch, case_sensitive=True)


class LazyAlignment:
    def __init__(self, target, query, result, params, traced):
        self.target, self.query = target, query
        self.result, self.params, self.traced = result, params, traced
        self.score = max(0, result.score)
        self._alignment = None

    def materialize(self):
        if self._alignment is None:
            result = self.result if self.traced else run(self.target, self.query, *self.params, trace=True)
            assert result.score == self.score
            cigar = result.cigar
            x, y = cigar.beg_query, cigar.beg_ref
            coords = [[x], [y]]
            for n, op in re.findall(r'(\d+)([=XMID])', cigar.decode.decode()):
                n = int(n)
                # parasail query is our target (first sequence).
                if op in '=XMI': x += n
                if op in '=XMD': y += n
                # parasail can encode an unaligned prefix as leading I/D.
                # Local SW does not score it; exclude it from aligned columns.
                if len(coords[0]) == 1 and op in 'ID':
                    coords = [[x], [y]]
                    continue
                coords[0].append(x); coords[1].append(y)
            assert x == result.end_query + 1 and y == result.end_ref + 1
            self._alignment = Alignment([self.target, self.query], np.array(coords, dtype=np.int64))
            self._alignment.score = self.score
        return self._alignment

    @property
    def coordinates(self): return self.materialize().coordinates
    @property
    def length(self): return self.materialize().length
    def counts(self): return self.materialize().counts()
    def format(self): return self.materialize().format()


def run(target, query, opening, extension, scoring, trace):
    fn = parasail.sw_trace_striped_16 if trace else parasail.sw_striped_16
    result = fn(target, query, opening, extension, scoring)
    if result.saturated:
        fn = parasail.sw_trace_striped_32 if trace else parasail.sw_striped_32
        result = fn(target, query, opening, extension, scoring)
        if result.saturated: raise OverflowError('parasail 32-bit alignment saturated')
    return result


class PairwiseAligner:
    def __init__(self):
        self.mode = 'local'
        self.match_score, self.mismatch_score = 1, 0
        self.open_gap_score, self.extend_gap_score = 0, 0

    def align(self, target, query):
        if self.mode != 'local': raise ValueError('Only local mode is supported')
        scores = (self.match_score, self.mismatch_score, self.open_gap_score, self.extend_gap_score)
        if any(int(x) != x for x in scores): raise ValueError('Integer scores required')
        match, mismatch, opening, extension = map(int, scores)
        if match <= 0 or mismatch > 0 or opening > 0 or extension > 0:
            raise ValueError('Expected positive match, nonpositive mismatch and gaps')
        target, query = str(target), str(query)
        if not target or not query: raise ValueError('Empty sequence')
        alphabet = set('ACGTRYSWKMBDHVNacgtryswkmbdhvn')
        if not set(target + query) <= alphabet: raise ValueError('Unsupported nucleotide')
        params = (-opening, -extension, matrix(match, mismatch))
        # Existing probes only need scores and have extension=0.
        traced = extension != 0
        result = run(target, query, *params, trace=traced)
        if result.score <= 0: return []
        return [LazyAlignment(target, query, result, params, traced)]

    def score(self, target, query):
        """Same SW score and saturation fallback, without alignment wrappers."""
        if self.mode != 'local': raise ValueError('Only local mode is supported')
        params = _score_parameters(self.match_score, self.mismatch_score,
                                   self.open_gap_score, self.extend_gap_score)
        target, query = str(target), str(query)
        if not target or not query: raise ValueError('Empty sequence')
        _check_dna(target)
        _check_dna(query)
        return max(0, run(target, query, *params, trace=False).score)


@lru_cache(maxsize=32)
def _score_parameters(match, mismatch, opening, extension):
    scores = (match, mismatch, opening, extension)
    if any(int(x) != x for x in scores): raise ValueError('Integer scores required')
    match, mismatch, opening, extension = map(int, scores)
    if match <= 0 or mismatch > 0 or opening > 0 or extension > 0:
        raise ValueError('Expected positive match, nonpositive mismatch and gaps')
    return -opening, -extension, matrix(match, mismatch)


_INVALID_DNA = re.compile('[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]')


@lru_cache(maxsize=1024)
def _check_dna(sequence):
    if _INVALID_DNA.search(sequence): raise ValueError('Unsupported nucleotide')
