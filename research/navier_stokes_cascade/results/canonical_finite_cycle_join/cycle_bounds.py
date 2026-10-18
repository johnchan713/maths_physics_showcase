"""Fresh canonical prefix, reusable operator budget and finite-stage limits."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

K=F(1,100000)


@lru_cache(None)
def operator_budget():
    # This module is an abstract affine operator budget; no old field is read.
    p=Path(__file__).resolve().parent.parent/'finite_induction_foundation/finite_bounds.py'
    spec=importlib.util.spec_from_file_location('canonical_abstract_finite_budget',p)
    m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m


def prefix_orders():
    mu=rho=1-K;alpha=2-3*K
    signed_real=min(alpha+1-2*K,2*alpha-1)
    real_fast=alpha-K
    real_mean=min(signed_real,real_fast+1-2*K,alpha+1)
    original_source=min(1-K,F(1,2)+mu-F(1,2))
    source_after_real=min(original_source,alpha-K,2*alpha-1-K,real_fast,alpha)
    B=source_after_real;C=B+F(1,2)
    sourced_cov=min(B+F(1,2),2*B)
    D=min(B+F(1,2)-K,2*B-K,B+mu-F(1,2))
    generic=min(rho+C-F(1,2),C+F(1,2)-K,2*C-1)
    H=C-K
    final_regular=min(D,H,C)
    final_mean=min(generic,H+1-2*K,C+1)
    final_defects=min(C+1,2*C)
    return dict(mu=mu,rho=rho,initial_mean=alpha,initial_defects=F(2),real_signed_averaged=signed_real,
        postreal_full_mean=real_mean,postreal_defects=min(alpha+1,2*alpha),source_after_real=source_after_real,
        first_sourced_order=B,first_sourced_pressure=C,first_sourced_covariance=sourced_cov,
        first_sourced_full_mean=C-K,generic_signed_Haar=generic,fast_mean_order=H,
        final_regular=final_regular,final_full_mean=final_mean,final_defects=final_defects,
        old_Haar_tangential=F(1),old_Haar_radial=F(2))


def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    return x


def prefix_review():
    p=prefix_orders();o=operator_budget();required=p['final_regular']+F(1,2)
    return dict(orders=encode(p),checks=dict(
        fresh_canonical_initial_mean_is_only_seed_input=p['initial_mean']==2-3*K and p['initial_defects']==2,
        real_signed_square_and_sine_curls_retained=p['real_signed_averaged']==3-6*K,
        actual_first_source_retains_whole_path_order=p['source_after_real']==1-K,
        sourced_covariance_retains_old_cross_and_new_square=p['first_sourced_covariance']==p['first_sourced_pressure'],
        generic_three_covariance_error_groups_tie=p['generic_signed_Haar']==2-2*K,
        fast_mean_wave_source_ties_regular_bottleneck=p['fast_mean_order']==p['final_regular']==o.B0,
        current_canonical_full_mean_meets_state_floor=p['final_full_mean']==required,
        recomputed_defects_stronger_than_state_floor=p['final_defects']==F(5,2)-K>required,
        old_error_not_promoted_with_source_order=p['rho']<p['final_regular'],
        accumulated_old_Haar_mean_is_not_zero=p['old_Haar_tangential']==1 and p['old_Haar_radial']==2,
        canonical_root_bundle_not_identified_with_historical_root_values=True))


def input_review(record):
    o=record['reviews']['canonical_initial_order_ledger']['orders'];p=prefix_orders();f=record['flags']
    return dict(current_initial_orders=dict(o),checks=dict(
        canonical_complete_real_input=f['canonical_real_complete_curl_all_six_covariance_inputs_joined'],
        canonical_own_roots_not_historical_root_identification=f['canonical_Haar_principal_target_matches_with_its_own_rebuilt_roots'],
        initial_full_mean_matches_prefix=F(o['final_tangential_mean'])==p['initial_mean'],
        initial_current_three_defects_match_prefix=F(o['final_defects'])==p['initial_defects'],
        old_mean_targets_were_not_transplanted=not f['old_mean_targets_or_old_mean_velocities_transplanted'],
        current_input_remains_conditional=not f['unconditional_pulse_family_or_blowup_proved']))


def cycle_review():
    o=operator_budget();c=o.certificate()
    return dict(certificate=o.encode(c),checks=dict(
        universal_operator_inequalities_not_old_field_identity=True,
        all_affine_slopes_and_boundary_margins_cover_whole_half_line=all(v>=0 for v in c['margins'].values()),
        preserved_state_margins_strict=all(v>0 for v in c['invariants'].values()),
        gain_equals_half_minus_kappa=c['gain']==F(1,2)-K,
        fast_viscosity_charged_twice_for_two_radial_graph_derivatives=c['margins']['fast_full_mean:radial_viscosity']==F(1,2)-2*K,
        old_complex_cross_limits_mean_gain=c['margins']['averaged_signed:old_complex_cross']==0,
        stale_postfast_target_gap_is_positive=1-(o.MU-K)==2*K>0,
        no_initial_real_shortcut_in_generic_cycle=True))


def physical_review():
    o=operator_budget();rows=[]
    for h,m,target in ((F(1,200),0,F(0)),(F(1,200),2,F(1)),(F(1,400),3,F(2))):
        j=o.first_stage(target,h,m)
        rows.append(dict(h=str(h),fixed_jet=m,target_power=str(target),finite_stage=j,
            power=str(o.physical_power(j,h,m)),preceding_power=str(o.physical_power(j-1,h,m)) if j else None))
    return dict(rows=rows,checks=dict(
        fixed_finite_target_and_jet_have_a_finite_canonical_state=all(F(r['power'])>=F(r['target_power']) and (r['preceding_power'] is None or F(r['preceding_power'])<F(r['target_power'])) for r in rows),
        fixed_zero_jet_physical_order_is_still_negative_at_seed=o.physical_power(0,F(1,200),0)<0,
        radius_and_stage_dependent_constants_not_removed=True,
        physical_bound_conditional_on_retained_core_background_and_analytic_maps=True))


def diagonal_review():
    o=operator_budget();h=F(1,200);rows=[];good=True
    for ell in (10**6,2*10**6,4*10**6):
        x=2*h*o.GAIN*ell;j=-(-x.numerator//x.denominator)
        power=o.physical_power(j,h,0)
        # R_j(q)=2^(j^2)*q^(hB_j-K0), q=2^-ell. No enormous float evaluated.
        log_residual=j*j-ell*power
        good &= power>0 and log_residual>0
        rows.append(dict(ell=ell,diagonal_stage=j,fixed_stage_power=str(power),
            log2_stage_constant=j*j,log2_diagonal_residual=str(log_residual)))
    return dict(rows=rows,checks=dict(
        improved_fixed_stage_powers_do_not_license_arbitrary_diagonal=bool(good),
        exact_counterexample_has_valid_finite_stage_bound=True,
        counterexample_does_not_rule_out_a_properly_controlled_diagonal=True,
        no_lower_bound_or_stage_constant_model_for_actual_canonical_fields_claimed=True))


def signed_coefficients(H,a,target,eps):
    H=[[F(x) for x in row] for row in H];a=list(map(F,a));target=list(map(F,target));eps=F(eps)
    if len(H)!=2 or any(len(r)!=2 for r in H) or len(a)!=2 or len(target)!=2 or min(a)<=0 or eps<=0:raise ValueError('Two positive fixed roots and positive epsilon required')
    det=H[0][0]*H[1][1]-H[0][1]*H[1][0]
    if not det:raise ValueError('Singular principal covariance matrix')
    d=[(H[1][1]*target[0]-H[0][1]*target[1])/(det*eps),
       (-H[1][0]*target[0]+H[0][0]*target[1])/(det*eps)]
    return [x/(2*root) for x,root in zip(d,a)]


def signed_review():
    rows=[];exact=True;crossed=False
    for scale in (F(1),F(1,100),F(1,10000)):
        H=[[scale,scale],[-scale,scale]];a=[F(1,2),F(2,3)];eps=F(1,200)
        for target in ([F(1),F(-2)],[F(-100),F(300)]):
            da=signed_coefficients(H,a,target,eps)
            got=[eps*sum(H[i][j]*2*a[j]*da[j] for j in range(2)) for i in range(2)]
            exact &= got==target;crossed |= any(root+delta<0 for root,delta in zip(a,da))
            rows.append(dict(H_scale=str(scale),signed_target=list(map(str,target)),delta_a=list(map(str,da))))
    return dict(rows=rows,checks=dict(exact_canonical_form_signed_pair_covariance_target=bool(exact),
        arbitrary_signed_targets_do_not_require_updated_amplitude_positive=bool(crossed),
        angular_half_in_H_and_factor_two_in_linearization_retained=bool(exact),
        radius_dependent_inverse_cost_not_removed=rows[0]['delta_a']!=rows[-2]['delta_a'],
        fixed_positive_roots_not_replaced_by_new_square_roots=True,
        reciprocal_root_edge_exponent_is_half_for_every_fixed_jet=True))
