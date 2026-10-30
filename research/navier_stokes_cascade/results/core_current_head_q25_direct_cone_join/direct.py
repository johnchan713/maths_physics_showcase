"""Rational direct-cone bounds, relative to actual provider hypotheses."""
from fractions import Fraction as Q

def exact(x):
    if isinstance(x,(float,bool)) or not isinstance(x,(int,Q)):raise TypeError('Exact inputs required')
    return Q(x)

def bounds(gamma,delta,A,u,jmin,jmax,vstar):
    g,d,A,u,lo,hi,v=map(exact,(gamma,delta,A,u,jmin,jmax,vstar))
    if not 0<g<=1 or not 0<=d<=g/4 or min(A,u,lo,v)<=0 or hi<lo:raise ValueError('Direct cone hypotheses required')
    lower=2-4*d-2*d*d;upper=2+4*d+2*d*d;numerator=g-d*(2-g)
    return {'det_lower':lo*lo*A*u*lower,'numerator_lower':lo*A*u*numerator/(A+u),
        'beta':numerator*v/(hi*upper),'original_beta':g*v/(8*hi),
        'inverse_bound':4*(1+d)/(lower*lo*min(A,u)),
        'original_K0':6/(lo*min(A,u))}

def record():
    return {'status':'canonical-direct-cone-to-root-margins',
        'normalized_determinant_interval':['7/8','25/8'],
        'beta':'4 gamma vstar/(25 jmax)', 'K0':'6/(jmin min(A,u))',
        'infinity_basis_cost_retained':True,
        'certified_direct_cone_tail_implication':True,
        'historical_minimum_band_membership_proved':False,
        'unconditional_residual_or_blowup':False}
