"""Exact conditional physical transfer of the corrected complete representatives."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path
HERE=Path(__file__).resolve().parent

def load(name,relative):
    s=importlib.util.spec_from_file_location(name,HERE.parent/relative)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

diagonal=load('q16_diagonal','core_complete_potential_diagonal_join/diagonal.py')
KINDS={'complete_wave_velocity':'velocity','wave_potential':'potential',
       'mean_azimuthal_potential':'potential','direct_axisymmetric_swirl':'velocity',
       'complete_fixed_gauge_pressure_difference':'pressure','mean_radial_velocity':'velocity'}

def order(n,minimum=0):
    if type(n) is not int or n<minimum:raise ValueError('invalid integer order')
    return n

def rows(h,j=1,m=0):
    if isinstance(h,(float,bool)):raise ValueError('exact parameter required')
    h=Q(h)
    if not 0<h<Q(1,100):raise ValueError('0<h<1/100 required')
    order(j,1);order(m)
    k=Q(1,100000);delta=Q(1,2)-3*k;A=Q(1,2)+h
    normalization={'velocity':A,'pressure':2*A,'potential':h}
    derivative_cost=1+5*h/2
    target=h*delta*j-2*A-2*m
    result={}
    for name,alpha in diagonal.representative_orders(j).items():
        norm=normalization[KINDS[name]]
        exponent=h*alpha-norm-m*derivative_cost
        result[name]={'normalized_order':alpha,'normalization_loss':norm,
                      'physical_exponent':exponent,'target':target,
                      'margin':exponent-target,'rate_above_stage':alpha-delta*j}
    return result

def symbolic():
    return {'h':'exp(-8*(exp(64)+10))','kappa':'1/100000','delta':'1/2-3*kappa',
            'normalized_margins':{'complete_wave_velocity':'1','wave_potential':'3/2',
             'mean_azimuthal_potential':'3/2-kappa','direct_axisymmetric_swirl':'3/2-kappa',
             'complete_fixed_gauge_pressure_difference':'3/2-kappa','mean_radial_velocity':'5/2-kappa'},
            'physical_margin':'h*(alpha-delta*j)+(2*A-normalization)+m*(1-5*h/2)',
            'minimum_physical_margin':'h*(3/2-kappa)+m*(1-5*h/2)',
            'sqrt_weight_constant':'factorial(n)*2^n*w^(-2*n), n=ceil((p+s)/2)',
            'actual_normalized_coefficients_bound':False,'same_field_representative_theorem_proved':False}

def polynomial_log_bound(coefficients,K):
    """P(S+K)<=sum c_i(16+K)^i (1+|log q|)^(2d) for ell>=16."""
    cs=[Q(x) for x in coefficients];K=Q(K)
    if K<0 or not cs or any(x<0 for x in cs):raise ValueError('nonnegative finite polynomial required')
    return sum(c*(16+K)**i for i,c in enumerate(cs)),2*(len(cs)-1)
