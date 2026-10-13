"""Exact DBSCAN for FastSTR's integer (position, constant-period) anchors.

No generic metric/weights support. Core components are connected in sorted
one-dimensional position order. Their labels are ordered by the earliest CORE
point in the original input, matching sklearn's traversal order. A border point
reachable from two components belongs to the earlier component, not necessarily
the nearest or leftmost component. Border points never connect components.
"""
import numbers
import numpy as np


class DBSCAN:
    def __init__(self, *, eps=0.5, min_samples=5, n_jobs=None):
        if not isinstance(eps, numbers.Real) or not np.isfinite(eps) or eps <= 0:
            raise ValueError('eps must be finite and positive')
        if not isinstance(min_samples, numbers.Integral) or min_samples < 1:
            raise ValueError('min_samples must be a positive integer')
        self.eps=float(eps)
        self.min_samples=int(min_samples)
        self.n_jobs=n_jobs  # This implementation itself is single-threaded.

    def fit(self, X, y=None, sample_weight=None):
        if sample_weight is not None:
            raise ValueError('Weighted DBSCAN is outside the FastSTR anchor contract')
        a=np.asarray(X)
        if a.ndim != 2 or a.shape[0] == 0 or a.shape[1] != 2:
            raise ValueError('Expected a nonempty N x 2 anchor array')
        if a.dtype.kind not in 'iu':
            raise ValueError('FastSTR anchors must have integer coordinates')
        if np.any(a[:,1] != a[0,1]):
            raise ValueError('All anchors in one fit must have the same period')
        # FastSTR windows are small; this bound avoids integer overflow and
        # loss of exact integer coordinates when applying a floating eps.
        if np.any(a[:,0] > 2**40) or np.any(a[:,0] < -(2**40)) or self.eps > 2**40:
            raise ValueError('Coordinates/eps exceed the supported exact range')
        n=len(a)
        # Integer coordinates imply integer distances; flooring eps avoids
        # rounding an almost-integer radius up when adding it to a coordinate.
        radius=int(self.eps)
        order=np.argsort(a[:,0],kind='stable')
        x=a[order,0].astype(np.int64,copy=False)
        lo=np.searchsorted(x,x-radius,side='left')
        hi=np.searchsorted(x,x+radius,side='right')
        core=(hi-lo)>=self.min_samples  # Includes the query point itself.
        self.n_features_in_=2
        self.core_sample_indices_=np.sort(order[core])
        self.components_=a[self.core_sample_indices_].copy()
        self.labels_=np.full(n,-1,dtype=np.intp)
        if not np.any(core):
            return self

        core_x=x[core]
        components=np.cumsum(np.r_[True,np.diff(core_x)>radius])-1
        component_count=int(components[-1])+1
        first_core=np.full(component_count,n,dtype=np.intp)
        np.minimum.at(first_core,components,order[core])
        label_map=np.empty(component_count,dtype=np.intp)
        label_map[np.argsort(first_core)]=np.arange(component_count)
        core_labels=label_map[components]

        # At most two disconnected core components can reach a 1-D point:
        # the components of its closest core on the left and on the right.
        right=np.searchsorted(core_x,x,side='left')
        left=right-1
        li=np.clip(left,0,len(core_x)-1)
        ri=np.clip(right,0,len(core_x)-1)
        left_ok=(left>=0)&((x-core_x[li])<=radius)
        right_ok=(right<len(core_x))&((core_x[ri]-x)<=radius)
        left_label=np.where(left_ok,core_labels[li],component_count)
        right_label=np.where(right_ok,core_labels[ri],component_count)
        labels=np.minimum(left_label,right_label)
        labels[labels==component_count]=-1
        self.labels_[order]=labels
        return self

    def fit_predict(self, X, y=None, sample_weight=None):
        return self.fit(X,y,sample_weight).labels_
