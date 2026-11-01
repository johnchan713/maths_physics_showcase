"""Exact shell budgets for the final radial compact term."""
from fractions import Fraction as Q
K=Q(1,100000)
TARGET=3-5*K

def row(start,r):
    p=(TARGET-start)/K
    if p.denominator!=1 or p<3 or r not in (0,1):raise ValueError('Unsupported shell budget')
    p=int(p)
    return dict(p=p,extra_source=2*p-r,extra_slow=p+1,
                radial_auxiliary_pairs=[[s,2*p-r-s] for s in (p+1,p-1)])

def record():
    return dict(diagonal=row(Q(1),0),shear=row(Q(1),1),remainder=row(2-3*K,0),
                status='conditional-final-radial-covariance-split',
                physical_wave_replaced=False,pressure_truncated=False,
                actual_envelopes_certified=False,complete_residual_or_blowup=False)
