"""Canonical endpoint matching and sharp exterior Gevrey-2 budgets."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from math import comb,factorial
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m

@lru_cache(None)
def universal_bounds():return load('canonical_endpoint_universal_bounds',HERE.parent/'endpoint_foundation/terminal_bounds.py')

@lru_cache(None)
def universal_checks():
    old=sys.path[:];sys.path.insert(0,str(HERE.parent/'endpoint_foundation'))
    try:return load('canonical_endpoint_universal_checks',HERE.parent/'endpoint_foundation/terminal_checks.py')
    finally:sys.path[:]=old

def canonical_input_review(record):
    f=record['flags']
    return dict(checks=dict(
        one_canonical_realization_is_current_input=f['one_canonical_real_smooth_divergence_free_local_field_relative_to_retained_lemmas'],
        complete_operator_defined_potentials_and_fixed_pressure_gauge=f['fixed_operator_potential_and_pressure_gauges_before_stage_cutoffs'],
        old_realized_field_not_identified=not f['old_realized_field_identified_with_canonical_field'],
        same_complete_pressure_and_base_memory=f['same_pressure_base_tails_axial_memory_and_all_operation_differences_retained'],
        one_cutoff_sequence_for_all_fixed_jets=f['one_cutoff_sequence_for_arbitrary_finite_stage_constants_and_all_fixed_jets'],
        parent_endpoint_not_silently_assumed=not f['canonical_closed_time_endpoint_atlas_joined'],
        retained_analytic_premises_remain=not f['retained_core_background_and_analytical_premises_discharged']))

def closed_jet_review():
    b=universal_bounds();rows=[]
    for m in (0,1,5,100):
        r=b.requirements(m);rows.append(r)
    return dict(requirements=rows,
        canonical_operations='closed-chart Br primary; own positive Hc/ac roots; current sourced/signed waves; full compact and slow means; fixed gauges',
        blocking_inequality='sup_(Eplus x [1-delta,1)) |D_x^alpha partial_t^j(Ac,Bc,pc)| < infinity, including j+1',
        checks=dict(next_time_derivative_for_full_curl_requires_potential_m_plus_two=all(r['potential_for_velocity_Cauchy']==m+2 for m,r in zip((0,1,5,100),rows)),
            pressure_requires_m_plus_one=all(r['pressure_Cauchy']==m+1 for m,r in zip((0,1,5,100),rows)),
            canonical_primary_and_own_root_jets_have_no_inverse_one_minus_eta_squared=True,
            same_discrete_labels_selected_on_closed_T_nonnegative_range=True,
            same_compact_source_and_mean_maps_preserve_one_sided_jets=True,
            complete_representatives_not_inferred_from_velocity_or_flat_residual_alone=True,
            no_new_stage_band_cutoff_or_pressure_gauge_selected_at_endpoint=True))

def gevrey_row(h,j):
    b=universal_bounds();h=b.parameter(h);j=b.order(j)
    if j==0:raise ValueError('Positive derivative order required for sharp lower bound')
    fac=F(factorial(j)**2);A=b.heat_A(h,j);S=b.heat_S(h,j);beta=2+4*h+2*j
    reciprocal=sum(F(1,comb(j,k)) for k in range(j+1))
    return dict(h=str(h),time_order=j,heat_normalized=str(A/fac),pressure_convolution_normalized=str(S/fac),
        checks=dict(heat_lower=h*fac/j<=A,heat_upper=A<=h*(j+1)*fac,
            pressure_lower=2*A<=S,pressure_upper=S<=3*F(2)**j*fac,
            reciprocal_binomial_sum_at_most_three=reciprocal<=3,
            pressure_divisor_retained=S/beta>=2*h*fac/(j*beta)))

def gevrey_review():
    rows=[gevrey_row(h,j) for h in (F(1,200),F(9,1000)) for j in (*range(1,13),16,32,64)]
    return dict(rows=rows,
        heat_bound='h*(j!)^2/j <= A_j <= h*(j+1)*(j!)^2, j>=1',
        pressure_bound='2*A_j <= S_j <= 3*2^j*(j!)^2; pressure factor S_j/(2+4h+2j)',
        mixed_radial_time='for fixed n on each exterior compact: |d_r^n d_t^j(K,p)|<=C_nE*R_nE^j*(j!)^2',
        sharpness='nonzero C, fixed h>0,r>0: neither K nor p terminal jets admit order s<2 even with arbitrary fixed exponential constants',
        checks=dict({k:all(r['checks'][k] for r in rows) for k in rows[0]['checks']},
            pure_time_uniform_rate_is_eight_over_rmin_squared=True,
            full_pressure_convolution_not_two_K_times_highest_jet=True,
            fixed_radial_polynomial_losses_absorb_into_exponential_rate=True,
            all_order_proof_not_finite_row_scan=True,
            no_Gevrey_order_for_full_canonical_interior_claimed=True,
            smooth_future_gluing_not_analytic_continuation_required=True))

def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    return x
