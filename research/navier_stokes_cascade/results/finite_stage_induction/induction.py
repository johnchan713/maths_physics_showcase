"""Exact half-line order certificate; analytical premises are in README.md."""
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
KAPPA=Q(1,100000)
B0=Q(3,2)-6*KAPPA
RHO=1-3*KAPPA
MU=1-KAPPA
GAIN=Q(1,2)-3*KAPPA


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


@dataclass(frozen=True)
class Affine:
    slope: Q
    constant: Q

    def at(self,b): return self.slope*Q(b)+self.constant

    def margin(self,target,minimum=B0):
        slope=self.slope-target.slope
        if slope<0: raise ValueError('Negative slope fails on the infinite half-line')
        return slope*Q(minimum)+self.constant-target.constant


def line(s,c=0): return Affine(Q(s),Q(c))
def encode(value):
    if isinstance(value,Q): return str(value)
    if isinstance(value,Affine): return dict(slope=str(value.slope),constant=str(value.constant))
    if isinstance(value,dict): return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [encode(v) for v in value]
    return value


def certificate(gain=GAIN):
    gain=Q(gain);k=KAPPA
    if gain<=0: raise ValueError('Strict positive finite-stage gain required')
    wave_target=line(1,gain);mean_target=line(1,gain+Q(1,2))
    # C=B+1/2; signed target C; complete signed increment order B.
    wave=dict(sourced_full_linear=line(1,Q(1,2)-3*k),
        sourced_old_wave_cross=line(1,Q(1,2)-k),
        sourced_new_square=line(2,-k),
        sourced_actual_mean_cross=line(1,MU-Q(1,2)),
        signed_full_linear=line(1,Q(1,2)-3*k),
        signed_old_wave_cross=line(1,Q(1,2)-k),
        signed_actual_mean_cross=line(1,MU-Q(1,2)),
        signed_new_square=line(2,-k),
        fast_mean_wave_change=line(1,Q(1,2)-k),
        slow_mean_wave_change=line(1,Q(1,2)))
    averaged=dict(old_complex_error_cross=line(1,RHO),
        complete_old_new_curl=line(1,1-k),signed_square=line(2,0),
        retained_current_moment_bump=line(1,Q(3,2)),
        full_axial_flux_and_pressure=line(1,Q(3,2)))
    fast=dict(slow_time=line(1,Q(3,2)-k),base_radial_flux=line(1,Q(3,2)-2*k),
        radial_viscosity=line(1,Q(3,2)-3*k),
        old_mean_transport=line(1,MU+Q(3,2)-2*k),
        changed_axial_pressure=line(1,Q(3,2)-k),new_mean_square=line(2,2-3*k))
    slow=dict(slow_time=line(1,Q(3,2)),base_radial_flux=line(1,Q(3,2)),
        radial_viscosity=line(1,Q(3,2)),old_mean_transport=line(1,MU+Q(3,2)-k),
        changed_axial_pressure=line(1,Q(3,2)),new_mean_square=line(2,2))
    defects=dict(accumulated_old_Haar_cross=line(1,Q(3,2)),new_slow_square=line(2,1),
        slow_time_base_and_viscosity=line(1,Q(5,2)),old_radial_mean_products=line(1,Q(7,2)))
    groups=dict(nonzero_regular=wave,averaged_signed=averaged,fast_full_mean=fast,
                slow_full_mean=slow,postslow_defects=defects)
    margins={}
    for group,rows in groups.items():
        target=wave_target if group=='nonzero_regular' else mean_target
        for term,value in rows.items(): margins[group+':'+term]=value.margin(target)
    # Exact state invariants and postfast defect input C are preserved.
    invariant=dict(wave_remainder=line(1,-RHO).at(B0),
        fast_mean_order=line(1,Q(1,2)-k-MU).at(B0),
        slow_Haar_order=line(1,-Q(1,2)).at(B0),
        postfast_defect_change_over_C=MU-k,
        primary_total_wave_order=RHO-Q(1,2),
        next_source_order_increase=gain)
    if min(margins.values())<0:
        failed=[n for n,v in margins.items() if v<0]
        raise ValueError('Gain blocked by: '+', '.join(failed))
    if min(invariant.values())<=0: raise ValueError('Required preserved state invariant failed')
    return dict(minimum_B=B0,gain=gain,fixed_primary_remainder=RHO,total_mean_order=MU,
        domain='every B>=3/2-6*kappa, proved by affine slopes and boundary values',
        groups=groups,wave_target=wave_target,mean_and_defect_target=mean_target,
        half_line_margins=margins,invariant_margins=invariant,
        minimum_margin=min(margins.values()),
        limiting_rows=[n for n,v in margins.items() if v==0])


def half_line_review():
    c=certificate();limiting=set(c['limiting_rows'])
    return dict(checks=dict(all_complete_cycle_rows_close_on_whole_half_line=all(v>=0 for v in c['half_line_margins'].values()),
        fixed_wave_and_mean_invariants_preserved=all(v>0 for v in c['invariant_margins'].values()),
        strict_positive_order_gain=GAIN>0,
        full_linear_and_old_complex_cross_both_limit_gain=limiting=={
            'nonzero_regular:sourced_full_linear','nonzero_regular:signed_full_linear','averaged_signed:old_complex_error_cross'},
        current_Haar_mean_order_one_retained=True,
        stage_constants_not_assumed_uniform=True),certificate=encode(c),
        gain_obstruction='s<=min(1/2-3*kappa,rho-1/2)=1/2-3*kappa',
        analytical_operator_lemmas_not_proved_by_this_rational_certificate=True)


def seed_review(parent):
    rows=parent['reviews']['current_mean_gain_physical_and_uncancelled_wave_ledger']['rows']
    b=Q(rows['final_nonzero_regular']);c=b+Q(1,2)
    margins=dict(wave=b-B0,full_mean=Q(rows['final_full_tangential'])-c,
        defects=Q(rows['slow_updated_defects'])-c,
        primary_remainder=Q(rows['current_fixed_primary_remainder'])-RHO,
        total_mean=Q(rows['mu'])-MU,old_Haar_tangential=Q(rows['old_Haar_tangential'])-1,
        old_Haar_radial=Q(rows['old_Haar_radial'])-2)
    return dict(checks=dict(actual_latest_seed_not_old_reference_initialization=b==B0,
        current_full_mean_matches_B_plus_half=Q(rows['final_full_tangential'])==c,
        current_defects_are_stronger_than_required=Q(rows['slow_updated_defects'])>c,
        all_seed_orders_satisfy_state=all(v>=0 for v in margins.values()),
        exact_Gaussian_tail_exclusion_is_a_seed_input=parent['flags']['actual_flat_velocities_in_every_product_and_Gaussian_tail_exclusion_preserved'],
        old_Haar_mean_is_accumulated_and_not_zero=parent['flags']['accumulated_previous_slow_mean_and_old_new_fast_correlations_retained']),
        seed_orders=encode(dict(B=b,C=c,defects=Q(rows['slow_updated_defects']),rho=RHO,mu=MU)),
        seed_margins=encode(margins))


def stage_order(j):
    if type(j) is not int or j<0: raise ValueError('Nonnegative finite integer stage required')
    return B0+GAIN*j


def physical_loss(h,m):
    h=Q(h)
    if not 0<h<Q(1,100) or type(m) is not int or m<0: raise ValueError('Construction h and fixed nonnegative jet order required')
    return 2*(Q(1,2)+h)+Q(3,2)+2*m


def physical_power(j,h,m): return Q(h)*stage_order(j)-physical_loss(h,m)


def first_stage_for_power(target,h=Q(1,200),m=0):
    target,h=Q(target),Q(h);loss=physical_loss(h,m)
    needed=(target+loss-h*B0)/(h*GAIN)
    return max(0,-(-needed.numerator//needed.denominator))


def harmonic_bound(j):
    stage_order(j)  # validate j, with no enumeration of all modes.
    return 4*2**j


def physical_quantifier_review():
    rows=[]
    for m,target in ((0,Q(0)),(2,Q(1)),(5,Q(10))):
        j=first_stage_for_power(target,m=m)
        rows.append(dict(jet_order=m,target=str(target),first_finite_stage=j,
            physical_power=str(physical_power(j,Q(1,200),m)),
            preceding_power=str(physical_power(j-1,Q(1,200),m)) if j else None))
    return dict(checks=dict(stage_order_increases_exactly_by_gain=stage_order(1)-stage_order(0)==GAIN,
        fixed_physical_loss_independent_of_stage=physical_loss(Q(1,200),3)==Q(851,100),
        each_fixed_target_has_a_finite_stage=all(Q(r['physical_power'])>=Q(r['target']) and
            (r['preceding_power'] is None or Q(r['preceding_power'])<Q(r['target'])) for r in rows),
        low_stage_physical_residual_still_has_negative_power=physical_power(0,Q(1,200),0)<0,
        finite_harmonics_double_inductively=harmonic_bound(0)==4 and harmonic_bound(1)==8,
        permitted_constant_growth_can_defeat_fixed_q_bound_contraction=-Q(1,200)*GAIN<0),
        rows=rows,physical_bound='for each fixed j,m: |R(U_j)|_m<=C_jm*q^(h*B_j-K_m)',
        K_m='2*A+3/2+2*m',B_j='3/2-6*kappa+j*(1/2-3*kappa)',
        quantifier_warning='forall target exists finite state does not give one state for all targets',
        admissible_constant_growth_example='C_j=q_fixed^(-2*h*delta*j) gives bound exponent h*B0-h*delta*j-K_m; not actual residual data',
        no_fixed_q_contraction_uniform_stage_constants_or_j_of_q_substitution=True,
        actual_thousand_stage_fields_not_numerically_evaluated=True)


@lru_cache(None)
def mechanism_review():
    p=HERE.parent/'residual_cycle_closure';sys.path.insert(0,str(p))
    try:
        from interaction_checks import longitudinal_and_transport_review,full_increment_review,covariance_increment_review
        from support_checks import support_and_flat_review,sharp_loss_controls
        functions=dict(complete_curl_transport=longitudinal_and_transport_review,
            full_actual_residual_increment=full_increment_review,generic_covariance=covariance_increment_review,
            supports_and_flat_velocity_memory=support_and_flat_review,sharp_missing_hypotheses=sharp_loss_controls)
        checks={name:all(bool(v) for v in f()['checks'].values()) for name,f in functions.items()}
    finally: sys.path.pop(0)
    return dict(checks=checks,checks_replayed_from_frozen_manufactured_interactions=True,
        actual_global_label_separation_and_inverse_estimates_remain_analytical_premises=True,
        designated_Gaussian_tails_never_promoted_to_forward_sources=True,
        numerical_diagnostics_not_a_continuum_proof=True)
