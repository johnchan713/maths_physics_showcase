"""Exact counts and Pell controls; no numerical continuum enclosure."""
from fractions import Fraction as F


def split(m):
    if not isinstance(m, int) or m < 4:
        raise ValueError("the two shell estimates require m >= 4")
    return m//2+1,(m+1)//2-1


def pell(j):
    a,b=1,0
    for _ in range(j):
        a,b=a+2*b,a+b
    return a,b


def record():
    k=F(1,100000)
    target=3-5*k
    m=(target-2)*2/k
    assert m.denominator==1
    return {str(i):dict(radial_counts=list(split(i)), phi_extra=i+1,
                        tail_order=str(2+i*k/2)) for i in (int(m)-1,int(m))}
