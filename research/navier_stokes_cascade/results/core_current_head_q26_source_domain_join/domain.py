"""Exact finite dyadic fixtures for the unchanged local-domain implication."""
from fractions import Fraction as Q

def coverage(ell_star,qstar,q,ell):
    if any(type(n) is not int or n<0 for n in (ell_star,ell)):raise ValueError('Nonnegative band indices required')
    if any(isinstance(x,(float,bool)) or not isinstance(x,(int,Q)) for x in (qstar,q)):raise TypeError('Exact rational fixtures required')
    qstar,q=Q(qstar),Q(q)
    if not 0<qstar<=1 or q<=0:raise ValueError('Positive domain required')
    Qband=Q(1,2**ell);qv=min(qstar/64,Q(1,2**(ell_star+2)))
    if q>qv or not Q(1,4)<=q/Qband<=4:raise ValueError('Point and enlarged label must meet the certified domain')
    return {'accepted':ell>=ell_star,'background_upper':4*Qband,
            'background_inside':4*Qband<=qstar/4,'q_unchanged':q}

def record():
    return {'status':'fixed-q-local-source-domain-closure',
        'domain':'0<q<=min(qstar/64,2^(-ell_*-2))',
        'label_bound':'Q<=4q<=2^(-ell_*)','background_bound':'4Q<=16q<=qstar/4',
        'recursion_domain_enlargement_factor':1,
        'unrestricted_exterior_slow_box_source_bound':False,
        'historical_localization_identified':False,'source_norm_proved':False,
        'field_or_selector_changed':False,'unconditional_residual_or_blowup':False}
