"""Finite positive coefficients of the real heat-composition majorant."""
from fractions import Fraction as F
from math import factorial

def heat_table(r,k):
    if min(r,k)<0: raise ValueError("negative order")
    J={(i,j):F(2*(1,2,1)[j],factorial(i)) for i in range(r+1)
       for j in range(min(k,2)+1) if i+j}
    out={(0,0):F(4)}; power={(0,0):F(1)}
    for a in range(1,r+k+1):
        nxt={}
        for (i,j),v in power.items():
            for (m,n),w in J.items():
                key=i+m,j+n
                if key[0]<=r and key[1]<=k:
                    nxt[key]=nxt.get(key,0)+v*w
        power=nxt
        for key,v in power.items():out[key]=out.get(key,0)+factorial(a+1)*v
    return out

def record():
    W=heat_table(3,2)
    return {'heat_control_y3_eta2':str(W[3,2]),'reference_cutoff_width':2,
            'heat_cutoff_width':'3/10','physical_propagator_established':False}
