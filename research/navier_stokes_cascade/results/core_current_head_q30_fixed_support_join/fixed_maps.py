"""Exact fixtures for fixed-profile gaps and finite slow response conditioning."""
from fractions import Fraction as Q

def gap(Xa,Xb,c0,c1):
    Xa,Xb,c0,c1=map(Q,(Xa,Xb,c0,c1))
    if not 0<Xa<c0<c1<Xb:raise ValueError('Fixed inner and outer plateau collars required')
    ra,rb=c0/Xa,Xb/c1
    return min(2*(ra-1)/(ra+1),2*(rb-1)/(rb+1))

def determinants(g,v):
    g,v=Q(g),Q(v)
    if not g>1 or not 0<v<1:raise ValueError('Original separated bumps and positive fixed lambda required')
    z=(g*g,v/(g*g),v)
    return {'angular':abs((z[1]-z[0])*(z[2]-z[0])*(z[2]-z[1])),
            'axial':g*(1-v)}

def hull(l,u,g):
    l,u,g=map(Q,(l,u,g))
    if not 0<l<u or g<=1:raise ValueError('Positive bump hull and geometric dilation required')
    return l,g*g*u

def record():
    return {'result':'fixed-choice-cutoff-and-slow-hull-join',
        'cutoff_gap':'min(log(c0/Xa),log(Xb/c1))>0; fixed physical X profile',
        'slow_hull':'[l0,g^2*u0] compactly inside x-image of I4',
        'slow_potential_hull_not_bumps_only':True,
        'axial_determinant':'g*(1-exp(-2*dpatch*lambda))>0',
        'new_profile_or_band_selected':False,'uniform_in_lambda_to_zero':False,
        'full_pressure_identified_with_leading':False,
        'actual_source_and_background_premises_closed':False,
        'unconditional_residual_or_blowup':False}
