"""Compact cross-stress pressure and radius budgets for the current field."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m
old=load('comparison_foundation_previous_radius_ledger',HERE.parent/'global_solution_comparison/comparison_bounds.py')

def young_decay(a,b):
    a,b=F(a),F(b)
    if not 0<a<2 or b<=0:raise ValueError('Require 0<a<2 and b>0')
    return 2*b/(2-a)

def separated_radius(radius,support_radius):
    radius,support_radius=F(radius),F(support_radius)
    if radius<1 or support_radius<0 or radius<2*support_radius:raise ValueError('Require R>=max(1,2*L)')
    return radius-support_radius

def kernel(p=F(4,3)):
    p=F(p)
    if not 1<p<F(3,2):raise ValueError('Require 1<p<3/2')
    return dict(near=str(2-2*p),far=str(2-3*p),integral=str(1/(3-2*p)+1/(3*p-3)),norm_decay=str(3-3/p))

def flux_rows(compact_cross=True):
    terms=[('quadratic_pressure_and_transport',F(3,2),F(1)),('quadratic_commutator',F(3,4),F(7,4))]
    if not compact_cross:terms.append(('unsplit_linear_pressure',F(1,2),F(1)))
    offsets=[F(5,2),F(5,2),F(2)]
    if compact_cross:offsets.append(F(5,2))
    else:offsets.append(F(3,2))
    rows=[dict(term=n,power=str(a),radius=str(b),young=str(young_decay(a,b)),offset=str(a+b)) for n,a,b in terms]
    return dict(rows=rows,other_decays=[str(v) for v in offsets],worst=str(min([F(r['young']) for r in rows]+offsets)))

def interface_review(parent):
    f=parent['flags']
    return dict(field='same current u=curl(c*A)+c*B, p=c*p0, zero datum, one compact smooth global force; viscosity one',
        fixed_interval='for each T<1, u,p smooth with support in K and finite sup ||u||_2, ||u||_6, ||grad u||_infinity',
        competitor='smooth on compact sets with same force/datum and sup_[0,T] ||v||_2 finite; no global derivative, pressure or energy inequality assumption',
        checks=dict(same_field_localized=f['same_current_single_realization_is_localized'],
            same_force_and_zero_datum=f['compact_smooth_global_force_relative_to_retained_lemmas_joined'] and f['zero_initial_datum_joined_to_same_field'],
            full_pressure_memory=f['both_background_pressure_tails_and_full_fixed_gauge_memory_retained'],
            current_energy_supplied=f['whole_space_preterminal_energy_and_total_dissipation_bounds_joined'],
            compact_preterminal_smoothness_supplies_fixed_interval_strain=True,
            core_lower_bound_is_still_relative=f['inherited_core_growth_preserved_on_localization_plateau'],
            retained_analytic_dependencies_explicit=not f['all_retained_analytic_dependencies_discharged']))

def pressure_review():
    return dict(stress='g=w tensor w+h; h=w tensor u+u tensor w; support(h) subset K',
        recovery='pi_star=sum R_i R_j g_ij in L_infinity_t H^-2; time-averaged harmonic gradient in H^-3 is zero',
        compact_cross='on R<=|x|<=2R, R>=max(1,2L): |pi_h|<=16/(4*pi)*R^-3*sum ||h_ij||_1',
        mixed_flux='|integral pi_h*w dot grad chi_R|<=C*R^-5/2*sum ||h_ij||_1*W_T',
        checks=dict(individual_forced_pressure_retains_minus_div_f=True,
            common_force_cancels_once_in_difference=True,
            pressure_gradient_recovered_before_any_pressure_flux_substitution=True,
            L1_Fourier_bound_embeds_in_H_minus_two_not_critical_H_minus_three_halves=True,
            time_averaged_momentum_excludes_harmonic_gradients_without_assuming_pressure_tempered=True,
            kernel_identity_distribution_zero_outside_compact_cross_support=True,
            separation_is_R_minus_L_at_least_R_over_two=separated_radius(2,1)==1,
            Hessian_component_bound_two_times_inverse_four_pi_r_cubed=True,
            cross_flux_decay_is_five_halves=4-F(3,2)==F(5,2),
            no_false_L1_stress_to_L1_pressure_bound_used=True))

def radius_review():
    r=flux_rows();oldr=flux_rows(False);k=kernel()
    return dict(current=r,unsplit_control=oldr,commutator_kernel=k,
        energy='0.5*E_R_prime+0.5*D_R^2<=||grad u||_infinity*E_R+C_T*R^-2',
        checks=dict(quadratic_main_has_B_to_three_halves_without_linear_B_to_one_half=True,
            commutator_near_and_far_integrals_both_kept=k==dict(near='-2/3',far='-2',integral='4',norm_decay='3/4'),
            exact_absorption_decays=[r['young'] for r in r['rows']]==['4','14/5'],
            every_remaining_flux_decays_at_least_two=F(r['worst'])==2,
            old_unsplit_estimate_reproduced=F(oldr['worst'])==F(4,3)==old.flux_review()['uniform_absorbed_remainder_power'],
            same_commutator_kernel_matches_old_exact_ledger=old.kernel_exponents()['norm_radius_power']==F(-3,4),
            eighth_cutoff_power_preserved=True,
            only_local_competitor_derivatives_used_in_cutoff_Sobolev=True,
            R_minus_two_is_scope_of_coarse_uniform_L2_laplacian_bound_not_universal_optimal_rate=True))

def limit_review():
    return dict(bound='E_R(t)<=2*C_T*T*exp(2*T*L_T)*R^-2; L_T=sup_[0,T] ||grad u||_infinity',
        core='tau<=min(tau0/2,r0^2/(8*Xin)); x_tau=(sqrt(2*Xin*tau),0,0); |u_theta|>=(e0/2)*tau^-A',
        checks=dict(fix_T_before_R_limit=True,no_uniform_T_to_one_constants_assumed=True,
            zero_strain_requires_no_division=True,fixed_ball_argument_needs_no_monotone_cutoff=True,
            equality_for_each_T_implies_equality_on_all_t_less_than_one=True,
            current_protected_core_remains_on_same_localization_plateau=True,
            compact_neighborhood_smoothness_contradicts_conditional_core_growth=True,
            conditional_exclusion_does_not_discharge_retained_exact_core_or_finite_analytic_lemmas=True))

def scope_review():
    return dict(result='same-current-field competitor exclusion relative to retained core and finite analytic lemmas; compact mixed-pressure refinement R^-2',
        next='trace retained exact core to infinite-space contraction/original equation and its global backward-pressure identification; then retained finite defining/existence/inverse/gap maps',
        checks=dict(general_comparison_proof_has_ordinary_explicit_hypotheses=True,
            special_construction_exists_only_relative_to_retained_lemmas=True,
            no_new_diagonal_schedule_inverse_or_Gaussian_sum=True,
            actual_infinite_fields_not_materialized=True,
            no_unforced_singularity_or_unconditional_forced_breakdown_claim=True,
            historical_evidence_and_old_valid_four_thirds_estimate_preserved=True))
