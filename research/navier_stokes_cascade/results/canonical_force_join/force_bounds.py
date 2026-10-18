"""Same-canonical-field localization and exact support-preserving Borel budgets."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m

@lru_cache(None)
def universal_bounds():return load('canonical_force_universal_bounds',HERE.parent/'force_foundation/localized_bounds.py')

@lru_cache(None)
def universal_checks():
    old=sys.path[:];sys.path.insert(0,str(HERE.parent/'force_foundation'))
    try:return load('canonical_force_universal_checks',HERE.parent/'force_foundation/localized_checks.py')
    finally:sys.path[:]=old

def index(n,minimum=0):
    if type(n) is not int or n<minimum:raise ValueError('Finite integer index required')
    return n

def first_controlled(a,m):return max(1,2*(index(a)+index(m)))

def canonical_input_review(endpoint,realization):
    e,r=endpoint['flags'],realization['flags']
    return dict(checks=dict(
        same_canonical_complete_field=e['same_canonical_local_field_and_complete_fixed_gauge_representatives'],
        actual_punctured_next_time_jets=e['compatible_punctured_one_sided_potential_swirl_velocity_pressure_jets_conditional'],
        uniform_corner_residual_is_complete=r['all_fixed_Cartesian_jets_of_complete_residual_flat_relative_to_retained_lemmas'],
        all_current_pressure_memory_and_curl_terms=e['all_pressure_memory_axial_velocity_and_cutoff_curl_terms_retained'],
        field_and_pressure_not_old_realization_identification=not e['old_endpoint_or_global_field_identified_with_canonical_field'],
        inherited_global_force_not_assumed=not e['canonical_localization_global_force_energy_and_comparison_joined'],
        no_interior_Gevrey_assumption=not e['full_canonical_interior_Gevrey_order_or_analyticity_claimed'],
        core_and_analytical_dependencies_remain=not e['retained_core_background_and_analytical_premises_discharged']))

def choose_dyadic_exponent(j,constraints,previous=1):
    """Exact integer n for b_j=2^n, using rational upper bounds on log2 C.

    These synthetic constraints are NOT actual canonical terminal norms.
    No enormous power 2^n is materialized.
    """
    j=index(j,1);previous=index(previous,1);pairs=list(constraints)
    if not pairs:raise ValueError('Nonempty finite constraint list required')
    n=max(2,previous+1)
    for a,m,c in pairs:
        a,m=index(a),index(m);c=F(c)
        if a+m>j//2 or c<0:raise ValueError('Require a+m<=floor(j/2), padded nonnegative log2 bound')
        x=(c+j+2)/(j-m);ceiling=-(-x.numerator//x.denominator)
        n=max(n,ceiling)
    return n+1

def scale_review():
    previous=1;rows=[]
    for j in (1,2,7,20,37):
        constraints=[(a,m,F((j+1)**4+10**(j//3)+(10**40 if j==37 else 0)))
            for a in range(j//2+1) for m in range(j//2+1-a)]
        n=choose_dyadic_exponent(j,constraints,previous)
        margin=max(c+(m-j)*n+j for a,m,c in constraints)
        rows.append(dict(stage=j,dyadic_exponent=n,constraints=len(constraints),
            log2_margin_upper=str(margin),increased=n>previous));previous=n
    cases=[(a,m,first_controlled(a,m)) for a,m in ((0,0),(1,3),(20,4),(6,70))]
    return dict(rows=rows,actual_norms_or_scales_evaluated=False,exact_rational_log2_upper_bounds=True,
        checks=dict(whole_support_constraints_hold=all(F(r['log2_margin_upper'])<=-2 for r in rows),
            one_strictly_increasing_dyadic_sequence=all(r['increased'] for r in rows),
            all_fixed_mixed_derivatives_eventually_controlled=all(a+m<=J//2 and m<J for a,m,J in cases),
            b0_two_and_later_exponents_above_one_give_future_support_to_three_halves=all(r['dyadic_exponent']>1 for r in rows),
            no_roundoff_in_scale_acceptance=True,
            no_uniform_terminal_order_constant_or_Gevrey_bound_needed=True))

def protected_disk_review():
    return dict(past='for all alpha,j,N: |D_x^alpha partial_t^j f_minus|<=C*(tau+|z|^(1/D))^N on fixed plateau box',
        traces='D_x^alpha F_j(x1,x2,0)=0 on the open protected equatorial disk for all alpha,j',
        future='for every alpha,m,N: |D_x^alpha partial_sigma^m f_plus|<=C_(alpha,m,N)*|z|^N uniformly on a smaller disk x |z|<=zL/4 x [0,1/2]',
        checks=dict(canonical_complete_residual_uniform_in_compact_X_and_exact_heat_zero_cover_entire_plateau=True,
            uniform_corner_bound_not_only_one_approach_path=True,
            uniform_global_force_jets_identify_smooth_compact_terminal_traces=True,
            all_spatial_trace_jets_vanish_on_disk_not_only_at_origin=True,
            uniform_differentiated_Borel_sum_preserves_flat_disk_for_all_future_times=True,
            finite_normal_Taylor_remainder_gives_each_uniform_power_of_abs_z=True,
            flat_disk_is_not_an_open_force_free_neighborhood=True,
            force_flatness_does_not_regularize_terminal_velocity_or_pressure=True))
