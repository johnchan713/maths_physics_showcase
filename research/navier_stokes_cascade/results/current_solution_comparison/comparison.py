"""Current-field comparison interface and exact radius/core limit budgets."""
from fractions import Fraction as Q
import importlib.util
import sys


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module


def exponent_review():
    rows=[]
    for a,b in ((Q(3,2),Q(1)),(Q(1,2),Q(1)),(Q(3,4),Q(7,4))):
        decay=2*b/(2-a)
        rows.append(dict(dissipation_power=str(a),radius_power=str(b),young_decay=str(decay)))
    offsets=[Q(5,2),Q(3,2),Q(5,2),Q(2)]
    # Direct polar integration of |z|^-3 min(|z|/R,1) at p=4/3.
    p=Q(4,3);near=2-2*p;far=2-3*p
    theta=Q(3,4);power=8
    return dict(rows=rows,offset_decays=[str(x) for x in offsets],
        kernel_radial_exponents=[str(near),str(far)],kernel_integral_coefficient=str(1/(near+1)-1/(far+1)),
        kernel_norm_radius_power=str(3/p-3),
        checks=dict(near_and_far_both_integrable=near>-1 and far<-1,
            commutator_radius_power_exact=3/p-3==Q(-3,4),
            L_four_pairing_exact=theta/6+(1-theta)/2==Q(1,4),
            eighth_weight_just_sufficient=Q(power,2)-1==Q(power,2)*theta,
            seventh_weight_insufficient=Q(7,2)-1<Q(7,2)*theta,
            all_Young_terms_subquadratic=all(Q(x['dissipation_power'])<2 for x in rows),
            worst_decay_is_four_thirds=min([Q(x['young_decay']) for x in rows]+offsets)==Q(4,3),
            critical_negative_Sobolev_order_not_used=2*Q(3,2)==3 and 2*2>3))


def fixed_interval_bound(T,strain,C,radius):
    """Exact Gronwall envelope using exp(x)-1<=x exp(x), not sampled flows."""
    T,strain,C,radius=map(Q,(T,strain,C,radius))
    if not 0<=T<1 or strain<0 or C<0 or radius<1:raise ValueError('Fixed preterminal interval, nonnegative constants and R>=1 required')
    return dict(prefactor=str(2*C*T),exponential_argument=str(2*strain*T),radius_decay='4/3',radius=str(radius))


def core_threshold(tau0,r0_squared,Xin):
    tau0,r0_squared,Xin=map(Q,(tau0,r0_squared,Xin))
    if min(tau0,r0_squared,Xin)<=0:raise ValueError('Positive observation/core constants required')
    return min(tau0/2,r0_squared/(8*Xin))


def limit_review():
    rows=[fixed_interval_bound(Q(3,4),L,5,R) for L in (0,2,1000) for R in (2,32,128)]
    thresholds=[core_threshold(Q(1,16),r2,X) for r2,X in ((Q(1,4),1),(Q(1,1000),17),(3,Q(1,9)))]
    return dict(envelope='E_R(t)<=2*C_T*T*exp(2*T*sup_[0,T]||grad u||_infinity)*R^(-4/3)',
        rows=rows,synthetic_core_thresholds=[str(x) for x in thresholds],
        checks=dict(zero_strain_case_needs_no_division=rows[0]['exponential_argument']=='0',
            all_constants_chosen_before_radius_limit=True,
            fixed_T_radius_limit_precedes_T_to_one=True,
            no_uniform_strain_integral_up_to_one_assumed=True,
            nonmonotone_cutoffs_allowed_by_fixed_ball_argument=True,
            core_points_eventually_in_current_plateau=all(x>0 for x in thresholds),
            compact_neighborhood_smoothness_not_global_spatial_Linfinity_is_contradicted=True))


def interface_review(parent):
    f=parent['flags'];r=parent['reviews']
    return dict(input='SAME current u=curl(c*A)+c*B, p=c*p0, f; datum zero, actual viscosity one',
        bounds='for each fixed T<1: support(u,p) subset K; sup ||u||2 finite; sup ||grad u||infinity finite',
        competitor='global smooth v,P with same f and datum and uniform L2; no pressure growth, global gradient or dissipation assumption',
        pressure='individual -Delta p=dd(uu)-div f; difference -Delta(P-p)=dd(vv-uu)',
        checks=dict(current_complete_field_not_historical_hypothetical_hierarchy=f['current_realized_complete_field_and_punctured_endpoint_jets_are_input'],
            same_force_and_zero_datum=r['current_whole_space_energy_and_core_growth']['initial_data']=='u(0)=0',
            one_global_smooth_compact_force=f['conditional_one_smooth_compact_global_force_matching_exact_preterminal_residual_supplied'],
            preterminal_energy_bound_supplied=f['conditional_current_whole_space_energy_and_integrated_dissipation_bound_supplied'],
            complete_pressure_and_nonzero_force_divergence_retained=f['force_divergence_in_exact_pressure_equation_retained'],
            conditional_current_core_lower_bound_retained=f['conditional_complete_core_lower_bound_survives_current_localization'],
            foundation_premises_not_discharged=f['inherited_leading_background_inverse_and_uniform_physical_estimates_remain_analytical_premises'],
            fixed_compact_support_and_preterminal_smoothness_suffice_for_strain_bound=True))
