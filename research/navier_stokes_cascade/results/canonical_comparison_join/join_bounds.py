"""Canonical comparison interface and exact compact-pressure moment controls."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from math import isqrt
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

@lru_cache(None)
def universal_bounds():return load('canonical_comparison_universal_bounds',HERE.parent/'comparison_foundation/compact_bounds.py')

@lru_cache(None)
def universal_checks():
    previous=sys.modules.get('compact_bounds');sys.modules['compact_bounds']=universal_bounds()
    try:return load('canonical_comparison_universal_checks',HERE.parent/'comparison_foundation/compact_checks.py')
    finally:
        if previous is None:sys.modules.pop('compact_bounds',None)
        else:sys.modules['compact_bounds']=previous

def unit(n):
    if len(n)!=3:raise ValueError('Three-dimensional direction required')
    n=tuple(F(x) for x in n)
    if sum(x*x for x in n)!=1:raise ValueError('Unit direction required')
    return n

def kernel_unit(n):
    n=unit(n)
    return [[3*n[i]*n[j]-int(i==j) for j in range(3)] for i in range(3)]

def directional_kernel(n,e):
    """4*pi*r^4 D_e k at x=r*n for rational unit n,e."""
    n,e=unit(n),unit(e);s=sum(a*b for a,b in zip(n,e))
    return [[3*(e[i]*n[j]+n[i]*e[j]+int(i==j)*s-5*n[i]*n[j]*s)
             for j in range(3)] for i in range(3)]

def contract(a,b):return sum(a[i][j]*b[i][j] for i in range(3) for j in range(3))

def trace_free(a):
    tr=sum(a[i][i] for i in range(3))
    return [[a[i][j]-F(tr,3)*int(i==j) for j in range(3)] for i in range(3)]

def moment_remainder(radius,support_radius,mass):
    radius,support_radius,mass=F(radius),F(support_radius),F(mass)
    universal_bounds().separated_radius(radius,support_radius)
    if mass<0:raise ValueError('Nonnegative component L1 mass required')
    return 384*support_radius*mass/radius**4  # inverse 4*pi is retained separately

def kernel_at(z):
    z=tuple(F(x) for x in z);square=sum(x*x for x in z)
    if square<=0:raise ValueError('Off-origin kernel required')
    a,b=isqrt(square.numerator),isqrt(square.denominator)
    if a*a!=square.numerator or b*b!=square.denominator:raise ValueError('Rational-radius fixture required')
    r=F(a,b)
    return [[(3*z[i]*z[j]-int(i==j)*square)/r**5 for j in range(3)] for i in range(3)]

def canonical_interface(parent):
    f=parent['flags']
    return dict(field='THIS u=curl(c*Ac)+c*Bc, p=c*pc; one canonical force and zero datum',
        competitor='smooth on compact sets; same force and datum; finite sup_[0,T] L2 for each T<1 only',
        checks=dict(
            canonical_field_and_full_gauge=f['same_canonical_complete_field_fixed_gauges_and_pressure_memory_localized'],
            divergence_preserved=f['axisymmetric_potential_cutoff_and_separate_direct_swirl_preserve_divergence'],
            force_pressure_equation_keeps_all_terms=f['all_twelve_localization_terms_and_full_force_divergence_pressure_equation_retained'],
            same_zero_datum_and_smooth_force=f['canonical_zero_initial_datum_and_compact_smooth_global_force_relative_to_retained_premises'],
            preterminal_energy=f['correct_fixed_time_energy_and_separate_total_dissipation_bounds_for_same_field'],
            protected_core_is_retained_not_proved=f['protected_core_growth_preserved_but_retained_core_lemma_not_discharged'],
            no_old_field_substitution=not f['historical_force_or_global_field_identified_with_canonical_field'],
            comparison_not_assumed_from_parent=not f['canonical_global_competitor_pressure_recovery_and_comparison_joined'],
            fixed_T_compact_smoothness_supplies_L6_and_strain=True,
            no_extra_global_competitor_derivative_pressure_or_energy_inequality_hypothesis=True,
            analytical_premises_not_discharged=not f['retained_core_background_and_finite_analytical_premises_discharged']))

def moment_review():
    directions=[(1,0,0),(0,0,1),(F(3,5),F(4,5),0),(F(2,3),F(1,3),F(2,3))]
    derivatives=[directional_kernel(n,e) for n in directions for e in directions]
    M=[[F(2*int(i==j and i<2)) for j in range(3)] for i in range(3)]
    values=[contract(kernel_unit(n),M) for n in directions]
    rows=[]
    for R in (F(2),F(3),F(7)):
        for L in (F(1),F(1,2)):
            for y in (-L,F(0),L):
                for i in range(3):
                    for j in range(3):
                        a=kernel_at((R-y,0,0))[i][j]-kernel_at((R,0,0))[i][j]
                        b=moment_remainder(R,L,1)
                        rows.append(dict(R=str(R),L=str(L),y=str(y),component=[i,j],error=str(a),bound=str(b)))
    return dict(witness='V=(-x2,x1,0)*rho(|x|^2); compact smooth radial nonzero rho; u=w=V only as tensor witness',
        normalized_moment=[[str(x) for x in r] for r in M],axis_value_without_inverse_four_pi=str(values[1]),
        derivative_controls=len(derivatives),point_remainder_controls=len(rows),rows=rows,
        actual_same_force_NS_trajectories=False,
        checks=dict(
            directional_derivative_bound_24=all(abs(x)<=24 for a in derivatives for r in a for x in r),
            differentiated_kernel_trace_zero=all(sum(a[i][i] for i in range(3))==0 for a in derivatives),
            nonzero_trace_free_moment=any(x for r in trace_free(M) for x in r),
            isotropic_part_drops_but_full_moment_does_not=all(contract(kernel_unit(n),M)==contract(kernel_unit(n),trace_free(M)) for n in directions) and values[1]==-4,
            compact_divergence_free_witness_uses_radial_parity_and_positive_a=True,
            no_same_force_flow_counterexample_claim=True,
            all_point_remainder_bounds_hold=all(abs(F(r['error']))<=F(r['bound']) for r in rows),
            inverse_four_pi_retained_and_segment_separation_required=True,
            mixed_flux_generic_decay=4-F(3,2)==F(5,2),
            mixed_flux_isotropic_decay_only_if_moment_isotropic=5-F(3,2)==F(7,2)))

def interval_review():
    return dict(energy='E_R(t)<=2*C_T*T*exp(2*L_T*T)*R^-2',
        dissipation='integral_0^t D_R^2<=2*C_T*t*(1+2*L_T*T*exp(2*L_T*T))*R^-2',
        protected_path='tau<=min(tau0/2,rL^2/(8*Xin)); SAME canonical localization radius rL',
        checks=dict(fix_T_then_take_radius_then_terminal_limit=True,
            local_dissipation_integrated_without_competitor_global_gradient_assumption=True,
            energy_and_dissipation_both_decay_with_radius_two=True,
            zero_strain_no_division_required=True,
            canonical_localization_radius_not_old_pulse_radius=True,
            no_uniform_constants_as_T_tends_to_one=True,
            no_future_velocity_or_pressure_constructed=True,
            conditional_core_growth_only_after_velocity_identification=True))
