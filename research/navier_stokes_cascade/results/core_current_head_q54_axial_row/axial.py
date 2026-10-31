"""Axial two-step rate ledger, keeping the unmultiplied tail separate."""
from fractions import Fraction as F
from math import ceil


def rates(k=F(1,100000)):
    if not 0 < k < F(1,10):
        raise ValueError("requires a small positive graph loss")
    z0 = 1-k
    z0_error = min(1,2-k)
    f1 = z0+1-2*k
    f1_error = min(z0_error+1-2*k,z0+1-k,2-2*k,2-k)
    z1_error = min(f1_error,2)  # direct bound for A1 f1 inside the inverse
    leading = f1+1-2*k
    remainder = min(z1_error+1-2*k,f1+1-k,3-2*k)
    p = ceil((leading-2)/k)
    return dict(leading=leading,remainder=remainder,velocity_error=z1_error,
                direct_tail=2,required=3-6*k,p=p,phi_extra=2*p+2)


def record():
    return {key:str(value) for key,value in rates().items()}
