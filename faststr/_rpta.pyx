# cython: language_level=3
"""Exact minimum-diagonal sufficient statistic for the INTERNAL RPTA caller.

The public Python RPTA still returns every set. global_motif_multiple_align
only reads min(set): min(A|B)=min(min(A),min(B)) and min(A+1)=min(A)+1.
No scoring, transition, or tie preference changes. Integer guard lives in
the Python caller; keep signed 64-bit scores and native path directions.
"""
import numpy as np
from collections import deque
cimport cython

cdef class PathRow:
    cdef object owner
    cdef Py_ssize_t i
    def __getitem__(self, Py_ssize_t j):
        cdef PackedPath p = self.owner
        cdef unsigned char direction = p.directions[self.i,j]
        if direction == 0: return -1
        if direction == 1: return (self.i-1,j-1,0)
        if direction == 2: return (self.i-1,j-1,1)
        if direction == 3: return (self.i,j-1,-1)
        return (self.i-1,j,-1)

cdef class PackedPath:
    cdef unsigned char[:, :] directions
    def __getitem__(self, Py_ssize_t i):
        cdef PathRow r=PathRow()
        r.owner=self
        r.i=i
        return r

@cython.boundscheck(False)
@cython.wraparound(False)
def align_min(bytes standard, bytes sequence, long long mat, long long mis,
              long long gap, int logo=0):
    cdef Py_ssize_t n=len(sequence), m=len(standard), i,j
    cdef long long[:, :] scores=np.empty((2,m+1),dtype=np.int64)
    cdef long long[:, :] minima=np.zeros((2,m+1),dtype=np.int64)
    cdef unsigned char[:, :] path=np.zeros((n+1,m+1),dtype=np.uint8)
    cdef long long delete,insert,replace,best,minimum,left,up,diag
    cdef int current,previous,mask
    cdef bint match
    cdef PackedPath packed=PackedPath()
    cdef list all_scores=[],all_minima=[]
    for j in range(m+1):
        scores[0,j]=j*gap
        if j: path[0,j]=3
    if logo==1:
        all_scores.append([scores[0,j] for j in range(m+1)])
        all_minima.append([{0} for j in range(m+1)])
    for i in range(1,n+1):
        current=i%2
        previous=1-current
        scores[current,0]=i*gap
        path[i,0]=4
        for j in range(1,m+1):
            match=standard[j-1]==sequence[i-1]
            delete=scores[current,j-1]+gap
            insert=scores[previous,j]+gap
            replace=scores[previous,j-1]+(mat if match else mis)
            best=max(insert,delete,replace)
            scores[current,j]=best
            mask=0
            minimum=n+m+1
            left=minima[current,j-1]
            up=minima[previous,j]
            diag=minima[previous,j-1]+1
            if best==delete:
                minimum=left
                mask+=2
            if best==insert:
                minimum=min(minimum,up)
                mask+=4
            if best==replace:
                minimum=min(minimum,diag)
                mask+=1
            minima[current,j]=minimum
            # Exactly the original seven-way tie resolution.
            if mask==1:
                path[i,j]=2 if match else 1
            elif mask==2:
                path[i,j]=3
            elif mask==4:
                path[i,j]=4
            elif mask==3:
                path[i,j]=(2 if match else 1) if diag<=left else 3
            elif mask==5:
                path[i,j]=(2 if match else 1) if diag<=up else 4
            elif mask==6:
                path[i,j]=3 if left<=up else 4
            else:
                if minimum==diag:
                    path[i,j]=2 if match else 1
                elif minimum==left:
                    path[i,j]=3
                else:
                    path[i,j]=4
        if logo==1:
            all_scores.append([scores[current,j] for j in range(m+1)])
            all_minima.append([{minima[current,j]} for j in range(m+1)])
    packed.directions=path
    if logo==1: return all_scores,all_minima,packed
    current=n%2
    return ([scores[current,j] for j in range(m+1)],
            [{minima[current,j]} for j in range(m+1)],packed)

def format_packed(PackedPath align_path, seq, motif, Py_ssize_t stab_l):
    cdef Py_ssize_t i=len(seq),j=stab_l,length=0
    cdef Py_ssize_t seq_length=len(seq)-1,motif_length=len(motif)-1
    cdef unsigned char direction
    target,query,alignment=deque(),deque(),deque()
    while True:
        direction=align_path.directions[i,j]
        if direction==0: break
        length+=1
        if direction==1 or direction==2:
            query.appendleft(seq[seq_length]); seq_length-=1
            target.appendleft(motif[motif_length]); motif_length-=1
            if motif_length==-1: motif_length+=len(motif)
            alignment.appendleft('|' if direction==2 else '.')
            i-=1; j-=1
        else:
            if direction==4:
                query.appendleft(seq[seq_length]); seq_length-=1
                target.appendleft('-'); i-=1
            else:
                query.appendleft('-')
                target.appendleft(motif[motif_length]); motif_length-=1
                if motif_length==-1: motif_length+=len(motif)
                j-=1
            alignment.appendleft('-')
    return [''.join(target),''.join(alignment),''.join(query),stab_l,length,len(seq)-1-seq_length]
