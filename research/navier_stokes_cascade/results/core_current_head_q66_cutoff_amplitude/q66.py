"""Finite normalized scalar majorants; multivariate proof in README."""
from math import factorial
import sympy as s

def cutoff_bound(j):
    if j<0: raise ValueError("negative order")
    return 1 if j==0 else 729*8**(j-1)*factorial(j-1)**3

def amplitude_majorant(A,B,afloor,bcap,delta,x,m):
    if afloor<=0 or delta<=0 or bcap<0 or m<0:
        raise ValueError("invalid floor or order")
    trunc=lambda p:s.series(p,x,0,m+1).removeO().expand()
    R=trunc(sum((A/afloor)**k for k in range(m+1)))
    V=trunc(A+((bcap+B)**2*R-bcap**2)/afloor)
    Z=trunc(sum(s.Rational(cutoff_bound(j),factorial(j))*(8*V/delta)**j for j in range(m+1)))
    W=trunc(sum((4*V/delta)**k for k in range(m+1)))
    return trunc(s.sqrt((2+delta/2)/afloor)*Z*W*R)

def record():
    return {'cutoff_raw_bounds':[cutoff_bound(j) for j in range(7)],
            'actual_state_jets_enclosed':False,'pressure_reset':False}
