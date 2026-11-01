"""Positive finite-series controls for the pulse energy root."""
from fractions import Fraction as Q

def square(v,n):
    return sum(v[i]*v[n-i] for i in range(n+1) if i<len(v) and n-i<len(v))

def energy_root(n,base,inverse,main,forcing,ends):
    # ends: (energy weight, absolute affine slope, angular intercept jets)
    a=[base]
    for j in range(1,n+1):
        cost=main*square(a,j)+forcing[j]
        for w,v,p in ends:
            z=[v*(a[i] if i<j else 0)+(p[i] if i<len(p) else 0) for i in range(j+1)]
            cost+=w*square(z,j)
        a.append(inverse*cost)
    return a

def record():
    return {'manufactured_coefficients':list(map(str,energy_root(4,Q(1),Q(2),Q(1,4),[0,1,2,3,4],[(Q(1,10),1,[0,1,2])]))),
            'heat_changes_cancel_only_in_total_energy':True,'physical_wave_proved':False}
