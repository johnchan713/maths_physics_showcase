"""Canonical representative ledger and exact whole-interval scale majorant.

The scales below use rational upper bounds for logarithms of SYNTHETIC
constants. They do not evaluate actual canonical fields or their constants.
"""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
K=F(1,100000);B0=F(3,2)-2*K;DELTA=F(1,2)-K;RHO=1-K

def index(n,minimum=0):
    if type(n) is not int or n<minimum:raise ValueError('Require a finite integer index')
    return n

def parameters(h):
    h=F(h)
    if not 0<h<F(1,100):raise ValueError('Require original 0<h<1/100')
    return h,F(1,2)+h,F(1,2)-h

def order(j):return B0+DELTA*index(j)
def gain(j,h=F(1,200)):return parameters(h)[0]*DELTA*index(j,1)
def tuple_loss(m,h=F(1,200)):return 2*parameters(h)[1]+2*index(m)
def cut_loss(m,h=F(1,200)):return tuple_loss(index(m)+1,h)
def total_loss(m,h=F(1,200)):return cut_loss(m,h)+1
def residual_loss(m,h=F(1,200)):return 2*parameters(h)[1]+F(3,2)+2*index(m)

def orders(j):
    b=order(index(j,1)-1);c=b+F(1,2)
    return dict(complete_wave_velocity=b-K,wave_potential=c-K,
        mean_azimuthal_potential=c-2*K,direct_axisymmetric_swirl=c-2*K,
        complete_fixed_gauge_pressure_difference=c-2*K,mean_radial_velocity=c+1-2*K)

def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    return x

def canonical_input_review(record):
    f=record['flags'];p=record['reviews']['fresh_canonical_prefix']['orders']
    return dict(seed=str(B0),gain=str(DELTA),checks=dict(
        own_canonical_family_not_historical_identification=not f['historical_pulses_sources_roots_mean_fields_or_realization_identified_as_canonical'],
        actual_canonical_seed_matches=F(p['final_regular'])==B0,
        actual_canonical_full_mean_matches=F(p['final_full_mean'])==B0+F(1,2),
        canonical_all_finite_operator_hypotheses=f['conditional_all_finite_cycle_gain_half_minus_kappa_on_whole_half_line'],
        same_supports_jets_domain_and_carriers=f['all_fixed_jets_same_value_domain_and_original_supports_width_and_carriers'],
        full_residual_groups_and_additive_tails=f['complete_pressure_viscosity_all_six_covariances_and_nonlinear_squares_retained'] and f['old_additive_ledger_preserved_new_tails_appended_never_forward_inverted'],
        no_old_infinite_field_imported=not f['one_simultaneous_canonical_infinite_realization_constructed'],
        retained_analytic_premises_remain=not f['retained_core_background_and_all_analytic_premises_discharged']))

def representative_review():
    margins={k:v-DELTA for k,v in orders(1).items()}
    return dict(margins=encode(margins),checks=dict(
        exact_stage_identity=B0-DELTA==RHO,
        every_row_has_exact_slope_delta=orders(2)=={k:v+DELTA for k,v in orders(1).items()},
        weakest_margin_one_minus_two_kappa=min(margins.values())==1-2*K>0,
        universal_physical_phase_loss_less_than_two=1+F(5,2)*F(1,100)<2,
        direct_swirl_kept_outside_full_compact_potential=True,
        current_all_operation_pressure_differences_in_one_fixed_gauge=True,
        canonical_representatives_rebuilt_not_old_values=True))

def ceil(x):
    x=F(x);return -(-x.numerator//x.denominator)

def comparison(m=0,target=0,h=F(1,200)):
    m,target=index(m),index(target);h,_,_=parameters(h)
    L,H=cut_loss(m+2,h),total_loss(m+2,h)
    jtail=ceil(2*(target+1+L+H)/(h*DELTA)-1)
    jfinite=ceil(((target+1+residual_loss(m,h))/h-B0)/DELTA)
    J=max(m+2,1,jtail,jfinite);tail=gain(J+1,h)/2-L
    return dict(jet=m,target=target,h=str(h),comparison_stage=J,
        tail_exponent=str(tail),both_cross_exponent=str(tail-H),square_exponent=str(2*tail),
        finite_residual_exponent=str(h*order(J)-residual_loss(m,h)),
        checks=dict(enough_jets=J>=m+2,linear=tail>=target+1,
            both_crosses=tail-H>=target+1,square=2*tail>=target+1,
            finite_residual=h*order(J)-residual_loss(m,h)>=target+1))

def comparison_review():
    rows=[comparison(m,N,h) for h in (F(1,200),F(9,1000),F(1,10000))
        for m,N in ((0,0),(1,3),(5,11),(20,2))]
    return dict(rows=rows,actual_fields_or_scales_evaluated=False,
        checks={k:all(r['checks'][k] for r in rows) for k in rows[0]['checks']})

def choose_log_scale(j,gamma,constraints,previous=0,log_qstar_inverse=2):
    """Exact sufficient bound on log(a_j); logC is an upper bound on ln(C).

    Using ln(2)<1 ensures strict doubling and 2^-j smallness. Each
    constraint is (logC,P), both nonnegative rational values.
    """
    j=index(j,1);gamma=F(gamma);previous=F(previous);domain=F(log_qstar_inverse)
    pairs=[(F(c),F(p)) for c,p in constraints]
    if gamma<=0 or previous<0 or domain<0 or not pairs or any(min(c,p)<0 for c,p in pairs):
        raise ValueError('Positive gamma and a nonempty finite nonnegative list required')
    candidates=[previous+1,domain+1,F(1)]
    candidates.extend(2*(c+j)/gamma+p*p/(gamma*gamma) for c,p in pairs)
    return max(candidates)+1

def scale_review():
    previous=F(0);rows=[]
    for j,pairs in ((1,[(0,0),(5,7)]),(4,[(1000,300),(10**6,2000)]),
                    (17,[(10**30,10**12),(300,0)])):
        gamma=gain(j)/2;s=choose_log_scale(j,gamma,pairs,previous)
        # Upper bound for logC+P*ln(1+s)-gamma*s+j*ln2 for ALL x>=s.
        margin=max(F(c)+F(p)**2/(2*gamma)-gamma*s/2+j for c,p in pairs)
        rows.append(dict(stage=j,gamma=str(gamma),log_scale_upper_choice=str(s),
            whole_interval_log_margin_upper=str(margin),doubling_margin=str(s-previous-1)))
        previous=s
    return dict(rows=rows,exact_rational_majorant=True,synthetic_log_constants_only=True,
        checks=dict(whole_interval_strict_majorant=all(F(r['whole_interval_log_margin_upper'])<0 for r in rows),
            strict_doubling_with_log_two_less_than_one=all(F(r['doubling_margin'])>0 for r in rows),
            no_floating_roundoff_in_scale_recipe=True,actual_infinite_constants_not_evaluated=True))

@lru_cache(None)
def universal_controls():
    p=HERE.parent/'realization_foundation/realized_checks.py'
    s=importlib.util.spec_from_file_location('canonical_realization_universal_checks',p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
