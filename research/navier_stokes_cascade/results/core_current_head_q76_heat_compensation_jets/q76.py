"""Positive angular and root envelopes; no evaluation of huge radii."""
from math import factorial

def kernel(k):
    return 4+sum(6**i*factorial(i+1) for i in range(1,k+1))

def target(k):
    K=kernel(k);L=2*K+K*K
    return max(25*L,50*L,10*K)

def root(k):
    d=[6400]
    for j in range(1,k+1):
        d.append(100*(target(j)+100*sum(d[l]*d[j-l] for l in range(1,j))))
    return d

def record():
    return {'K_through_4':[kernel(k) for k in range(5)],
            'D_through_3':root(3),'physical_band_gain':False}
