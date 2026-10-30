"""Exact affine lift fixtures and conditional derivative budget."""
from fractions import Fraction as Q

def path(H,w,b,c,scale=1,offset=0,radius=1):
    H=tuple(map(Q,H));w,b,c,scale,offset,radius=map(Q,(w,b,c,scale,offset,radius))
    if len(H)!=2 or not 0<b<1 or c<=0 or scale<1:raise ValueError('invalid fixed chart')
    lam=(b*H[0]+H[1])/(1+b*b)
    v=(scale*lam-offset+radius)/c
    a=c/scale
    return (H[0]+a*(w-v)*b,H[1]+a*(w-v))

def jacobian(b):
    b=Q(b)
    if not 0<b<1:raise ValueError('0<b<1 required')
    return ((1/(1+b*b),-b/(1+b*b)),(-b/(1+b*b),b*b/(1+b*b)))

def jet_factor(aux_order,path_order,S):
    if any(type(n) is not int or n<0 for n in (aux_order,path_order)):
        raise ValueError('nonnegative integer orders required')
    if isinstance(S,(float,bool)) or Q(S)<1:raise ValueError('exact S>=1 required')
    return Q(2)**(aux_order+path_order)/Q(S)**path_order

def record():
    return {'affine_auxiliary_jacobian':'I-v_t lambda_t',
            'higher_auxiliary_map_derivatives':0,
            'mixed_jet_factor':'2^(k+r)*S^(-r)',
            'extra_source_derivative_orders':0,
            'Gaussian_evaluated_at_path_point':True,
            'liftwise_not_deck_averaged':True,
            'conditional_extended_to_path_transfer_proved':True,
            'actual_extended_source_bounds_proved':False,
            'unconditional_residual_or_blowup':False}
