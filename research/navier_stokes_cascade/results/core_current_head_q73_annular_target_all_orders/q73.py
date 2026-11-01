"""Normalized coefficient recurrence for the varying annular quadratic map."""
from fractions import Fraction as Q

def root_bounds(k,Z,h,d0=Q(6,10**13),inv=2000,q=(100,)*5):
    if h<=0 or k<0: raise ValueError("positive radius and nonnegative order required")
    d=[d0]
    for n in range(1,k+1):
        total=Z*Q(1,h)**n
        for i in range(min(n,len(q)-1)+1):
            for j in range(n-i+1):
                l=n-i-j
                if j<n and l<n: total+=q[i]*d[j]*d[l]
        d.append(inv*total)
    return d

def record():
    return {'manufactured_order3':str(root_bounds(3,2,Q(1,8),Q(1,10),2,(1,2,3))[-1]),
            'requires_high_order_contraction':False,'physical_band_gate_proved':False}
