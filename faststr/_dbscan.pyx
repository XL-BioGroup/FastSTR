# cython: language_level=3
"""Exact labels-only DBSCAN: stable position order, linear neighbor scans.

Component labels use earliest ORIGINAL core index, not coordinate order.
Ambiguous borders belong to the smaller original-traversal component label.
Dense counting sort is only a representation optimization, with stable argsort
fallback for sparse coordinates. Hyperparameters are passed through unchanged.
"""
import numpy as np
import numbers
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def labels_only(X, eps=14, min_samples=14):
    if not isinstance(eps,numbers.Real) or not np.isfinite(eps) or eps<=0:
        raise ValueError('eps must be finite and positive')
    if not isinstance(min_samples,numbers.Integral) or min_samples<1:
        raise ValueError('min_samples must be a positive integer')
    a=np.asarray(X)
    if a.ndim!=2 or a.shape[0]==0 or a.shape[1]!=2 or a.dtype.kind not in 'iu':
        raise ValueError('Expected nonempty integer N x 2 anchors')
    # Avoid unsigned overflow before conversion.
    if a.dtype.kind=='u' and np.any(a>2**40):
        raise ValueError('Coordinates outside supported range')
    cdef long long[:, :] points=np.asarray(a,dtype=np.int64)
    cdef Py_ssize_t n=len(a),i,j,left=0,right=0,k,count=0,component=-1
    cdef long long low=points[0,0],high=low,x,last=0,radius
    cdef long long coordinate_bound=1099511627776
    cdef long long minimum=min_samples
    cdef Py_ssize_t span,offset,first_index,current_label,prev_core=-1,next_core
    cdef Py_ssize_t[:] order,counts,cursor,core_component,firsts,label_map,labels,components_by_first
    for i in range(n):
        x=points[i,0]
        if points[i,1]!=points[0,1] or x>coordinate_bound or x<-coordinate_bound:
            raise ValueError('Unsupported anchor coordinates/period')
        low=min(low,x); high=max(high,x)
    if eps>2**40: raise ValueError('eps outside supported range')
    radius=int(eps)
    span=high-low+1
    if span<=max(1024,4*n):
        counts=np.zeros(span,dtype=np.intp)
        cursor=np.empty(span,dtype=np.intp)
        order=np.empty(n,dtype=np.intp)
        for i in range(n): counts[points[i,0]-low]+=1
        offset=0
        for i in range(span):
            cursor[i]=offset
            offset+=counts[i]
        for i in range(n):
            k=points[i,0]-low
            order[cursor[k]]=i
            cursor[k]+=1
    else:
        order=np.argsort(a[:,0],kind='stable').astype(np.intp,copy=False)
    core_component=np.full(n,-1,dtype=np.intp)
    firsts=np.full(n,n,dtype=np.intp)
    for i in range(n):
        x=points[order[i],0]
        while left<n and points[order[left],0]<x-radius: left+=1
        while right<n and points[order[right],0]<=x+radius: right+=1
        if right-left>=minimum:
            if component==-1 or x-last>radius: component+=1
            core_component[i]=component
            firsts[component]=min(firsts[component],order[i])
            last=x
    result=np.full(n,-1,dtype=np.intp)
    labels=result
    if component==-1: return result
    # Linear relabeling by first original core index, no extra argsort.
    components_by_first=np.full(n,-1,dtype=np.intp)
    label_map=np.empty(component+1,dtype=np.intp)
    for i in range(component+1): components_by_first[firsts[i]]=i
    count=0
    for i in range(n):
        if components_by_first[i]>=0:
            label_map[components_by_first[i]]=count
            count+=1
    next_core=0
    for i in range(n):
        while next_core<n and (next_core<i or core_component[next_core]<0): next_core+=1
        x=points[order[i],0]
        current_label=count
        if prev_core>=0 and x-points[order[prev_core],0]<=radius:
            current_label=label_map[core_component[prev_core]]
        if next_core<n and points[order[next_core],0]-x<=radius:
            current_label=min(current_label,label_map[core_component[next_core]])
        if current_label<count: labels[order[i]]=current_label
        if core_component[i]>=0: prev_core=i
    return result
