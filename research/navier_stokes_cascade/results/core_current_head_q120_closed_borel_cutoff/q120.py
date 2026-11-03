"""Closed norms for one existing Borel cutoff; no replacement scales or field identity."""
from fractions import Fraction as Q
from math import factorial
from functools import lru_cache
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
current=load('q120_current','core_current_head_q119_closed_slow_chart/q119.py')
borel=load('q120_borel','core_borel_cutoff_certificate/budgets.py')
order=current.order
IDENTITY='implemented-squared-flat-Borel'

def cutoff_recipe(m):
    order(m)
    low=(1,64,4928)
    if m<=2:return dict(order=m,exact_upper=low[m],identity=IDENTITY,actual_field_identity_certified=False)
    return dict(order=m,scalar=82,geometric_base=4097,power=m,factorial_argument=m,
                factorial_power=3,identity=IDENTITY,actual_field_identity_certified=False)

def materialize(recipe,maximum=40):
    order(maximum);order(recipe['order'])
    if recipe['order']>maximum:raise ValueError('Keep the large cutoff factorials in the analytical recipe')
    if 'exact_upper' in recipe:return Q(recipe['exact_upper'])
    return Q(recipe['scalar'])*recipe['geometric_base']**recipe['power']*factorial(recipe['factorial_argument'])**recipe['factorial_power']

@lru_cache(None)
def implemented_bound(m):
    # A finite-order norm provider, not a check that an unspecified field uses chi.
    return materialize(cutoff_recipe(m),maximum=m)

def cutoff_provider(identity):
    if identity!=IDENTITY:raise ValueError('The norm theorem is for the same implemented Borel function only')
    return implemented_bound

def cutoff_binary_upper(m):
    """Certified log2 upper budget; evaluates no factorial or large power."""
    order(m)
    if m<=2:return (0,6,13)[m]
    return 7+13*m+3*m*borel.ceil_log2(Q(m))

def deferred_row(h,n,m,profile_binary_upper,binary_log_lower):
    order(m)
    if type(n) is not int or n<1 or m>n:raise ValueError('Original diagonal row required')
    order(profile_binary_upper);h,lo=map(borel.exact,(h,binary_log_lower))
    if not 0<h<Q(1,100) or lo<0:raise ValueError('Original positive h and supplied scale-log lower bound required')
    bits=1+cutoff_binary_upper(m+1)+m*borel.ceil_log2(1+2*n*h)+profile_binary_upper
    slack=n*h*lo-n-bits
    return dict(n=n,m=m,profile_binary_upper=profile_binary_upper,
                normalized_constant_binary_upper=bits,normalized_log_slack=slack,
                sufficient_new_row_gate=slack>=0,factorials_materialized=False,
                actual_profile_or_scale_inputs_certified=False,
                infinite_schedule_certified=False,scales_selected_or_changed=False)

def claim_bracket(m,claim):
    order(m);claim=borel.exact(claim)
    if claim<1:raise ValueError('Padded through-order cutoff bound required')
    lower=Q(1 if m==0 else 12);upper=implemented_bound(m)
    return dict(lower=lower,certified_upper=upper,
                impossible_by_midpoint=claim<lower,certified_by_upper=claim>=upper,
                undecided_by_these_bounds=lower<=claim<upper,
                actual_field_identity_certified=False)

def retained_row(h,n,m,profile,old_K,binary_log_lower):
    """Audit a supplied existing log2(a_n) lower bound. Never choose a scale."""
    order(m)
    if type(n) is not int or n<1 or m>n:raise ValueError('Original normalized diagonal row 0<=m<=n required')
    h,profile,old_K,lo=map(borel.exact,(h,profile,old_K,binary_log_lower))
    if not 0<h<Q(1,100) or profile<0 or old_K<1 or lo<0:raise ValueError('Original positive h, full tuple norm, padded old cutoff and supplied log required')
    new_K=implemented_bound(m+1)
    common=2*(1+2*n*h)**m*profile
    old_C=max(Q(1),common*old_K);new_C=max(Q(1),common*new_K)
    rho=max(Q(1),new_K/old_K)
    old_bits=borel.ceil_log2(old_C);new_bits=borel.ceil_log2(new_C);extra=borel.ceil_log2(rho)
    slack=n*h*lo-n-old_bits
    return dict(n=n,m=m,old_cutoff=old_K,new_cutoff=new_K,old_constant=old_C,new_constant=new_C,
                inflation=rho,extra_binary_slack=extra,old_binary_slack=slack,
                sufficient_constant_transfer_gate=slack>=extra,
                sufficient_direct_new_row_gate=n*h*lo>=n+new_bits,
                old_claim_impossible_by_midpoint=old_K<12,
                supplied_log2_lower=lo,actual_norm_or_scale_intervals_certified=False,
                infinite_schedule_certified=False,scales_selected_or_changed=False)

def record():
    return dict(status='closed-implemented-Borel-Euler-bound-and-retained-row-slack',
                gaussian='12^j*(j!)^3',reciprocal='81*2048^j*(j!)^3',
                step_raw='82*2048^j*(j!)^3',Euler_through='82*4097^m*(m!)^3',
                original_low_certificates=[1,64,4928],
                retained_scale_transfer='old_binary_slack >= ceil(log2(max(1,Knew/Kold)))',
                cutoff_binary_budget='7+13*m+3*m*ceil(log2(m)) for m>=3',
                actual_field_identity_or_infinite_schedule_certified=False,
                full_cone_pulse_residual_or_blowup=False)
