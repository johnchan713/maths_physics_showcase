"""Indices of the raw norm estimate, without constructing giant constants."""
from fractions import Fraction as F


def plan(n, r, p=99995):
    if r not in (0,1) or n<0 or p<r+2:
        raise ValueError("requires at least one remaining divisor power")
    pairs=[(s,2*p-r-s) for s in (p+1,p-1)]
    return dict(shells=[dict(radial=s,aux=q,power=s-r-q) for s,q in pairs],
                amplitude=n+1+2*p-r,weight=n+p+2,slow=n+1,
                curl=n+5,outer=n,rate=str(2+p*F(1,100000)))


def record():
    return dict(diagonal=plan(0,0),shear=plan(0,1),curl_rate=str(3-F(2,100000)))
