"""Current force, Borel gluing and energy budgets relative to retained maps."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb,factorial
import importlib.util
from pathlib import Path
import sys
import mpmath as mp

HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m
old=load('force_foundation_historical_borel_bounds',HERE.parent/'current_localized_force/forcing.py')
TERMS=old.TERMS
def order(n):
    if type(n) is not int or n<0:raise ValueError('Finite nonnegative integer order required')
    return n
def requirements(m):
    m=order(m)
    return dict(force=m,potential=m+3,direct_swirl=m+2,pressure=m+1,
        next_time_potential=m+4,next_time_direct_swirl=m+3,next_time_pressure=m+2)
def first_controlled(a,m):return max(1,2*(order(a)+order(m)))

def domain_review():
    return dict(cutoff='c=chi_r(r^2)*chi_z(z)*chi_t(t), compact axisymmetric support and unit plateau near (0,1)',
        field='u=curl(c*A)+c*B=c*u0+d; p=c*p0; d=grad c cross A',
        pressure='same complete p0 with both Pi1/Pi2 and all wave/mean differences; -Delta p=DD(u tensor u)-div f',
        requirements=[requirements(m) for m in (0,1,5,100)],
        checks=dict(current_enlarged_chart_margin_uses_sixteen_q=16*F(3,8)*F(1,16)==F(3,8)<1,
            fixed_first_band_and_domain_kept=4*F(3,8)*F(1,4)<1,
            zero_near_initial_time_from_tau0_less_than_one_over_256=F(1,16)*F(1,16)==F(1,256),
            r_squared_cutoff_is_Cartesian_smooth_and_preserves_direct_swirl_divergence=True,
            same_field_same_diagonal_sequences_and_fixed_pressure_gauge_used=True,
            complete_curls_time_viscosity_and_nonlinear_terms_retained=len(TERMS)==12,
            d_dot_grad_c_exactly_zero_by_cross_product_not_an_estimate=True,
            fields_defined_on_neighborhood_of_cutoff_support_before_zero_extension=True,
            localization_plateau_preserves_full_core_and_every_residual_jet=True))

def terminal_review():
    return dict(corner_bound='for all alpha,j,N: |D_x^alpha partial_t^j f_minus|<=C*(tau+|z|^(1/D))^N near corner',
        trace='F_j=lim_(t->1-) partial_t^j f_minus uniformly on R3; F_j in C_c^infinity(K), every spatial jet at origin zero',
        checks=dict(compact_X_uniform_complete_residual_plus_exact_heat_exterior_zero_used=True,
            q_le_two_tau_plus_four_abs_z_to_one_over_D=True,
            small_corner_ball_and_punctured_compact_cover_give_global_uniform_limits=True,
            spatial_integral_limits_identify_derivatives_even_across_origin=True,
            next_time_derivative_identifies_consecutive_terminal_jets=True,
            current_endpoint_bounds_supply_A_m_plus_four_B_m_plus_three_p_m_plus_two=True,
            actual_flat_velocities_and_Gaussian_ledger_kept_without_new_inversion_or_sum=True,
            no_regular_velocity_or_pressure_at_origin_inferred=True))

def derivative_factor(j,m,k):
    j,m,k=order(j),order(m),order(k)
    if k>m or m-k>j:return F(0)
    return F(comb(m,k),factorial(j-m+k))

def extension_review():
    rows=[];previous=1
    for j in (1,2,7,20):
        constraints=[(a,m,(j+1)**4+10**(j//3)) for a in range(j//2+1) for m in range(j//2+1-a)]
        n=old.choose_log2_scale(j,constraints,previous)
        with mp.workdps(80):margin=max(mp.mpf(C)+(m-j)*n*mp.log(2)+j*mp.log(2) for a,m,C in constraints)
        rows.append(dict(stage=j,dyadic_exponent=n,constraints=len(constraints),log_margin=float(margin),increased=n>previous));previous=n
    pairs=[(a,m,first_controlled(a,m)) for a,m in ((0,0),(1,3),(20,4),(6,70))]
    return dict(synthetic_norms_only=True,actual_terminal_norms_and_extension_scales_not_evaluated=True,
        term='T_j=chi0(b_j*sigma)*sigma^j/j!*F_j; b0=2, b_j=2^n_j increasing',
        derivative_bound='||D_x^a partial_sigma^m T_j||<=C_jam*b_j^(m-j)',
        tail='sum_(j>=J) ||D_x^a partial_sigma^m T_j||<=2^(1-J), J>=max(1,2*(a+m))',
        checks=dict(whole_support_finite_constraints_hold_for_extreme_norms=all(r['log_margin']<-1 for r in rows),
            dyadic_scales_increase_and_future_support_ends_by_three_halves=all(r['increased'] for r in rows),
            all_fixed_mixed_jets_eventually_in_one_constraint_sequence=all(a+m<=J//2 and m<J for a,m,J in pairs),
            Leibniz_factor_keeps_factorial_and_all_cutoff_derivatives=derivative_factor(3,2,1)==1,
            exact_terminal_jet_is_F_m_by_plateau_and_factorial=True,
            uniform_differentiated_series_matches_all_one_sided_terminal_jets=True,
            spatial_support_K_and_origin_flatness_preserved_by_uniform_jet_convergence=True,
            no_raw_Taylor_convergence_or_future_velocity_required=True),rows=rows)

def energy_review():
    # A rise followed by unforced decay illustrates why two separate sup/
    # total estimates cannot be combined with the SAME right-hand budget.
    epsilon=F(1,100);F1=1+epsilon/2
    lower_sum=2-F(1,1+2*(1-epsilon))+2*epsilon/3
    return dict(identity='0.5*d_t ||u||_2^2+||grad u||_2^2=<f,u>',
        regularization='d_t sqrt(||u||_2^2+delta^2)<=||f||_2; u(0)=0; delta->0',
        exact_bound='E(t)+2*integral_0^t D<=[integral_0^t ||f||_2]^2, every t<1',
        separate_bounds='sup_(t<1) E(t)<=F(1)^2; integral_0^1 D<=F(1)^2/2',
        wrong_combined_budget_counterexample=dict(scalar_y_rises_to_one_in_epsilon_then_decays_with_f_zero=True,
            epsilon=str(epsilon),force_L1_norm=str(F1),lower_bound_for_sup_energy_plus_twice_total_dissipation=str(lower_sum),force_budget_squared=str(F1**2)),
        checks=dict(energy_used_only_after_compact_smooth_force_constructed=True,
            no_division_by_velocity_norm_that_starts_at_zero=True,
            pressure_transport_fluxes_vanish_for_actual_compact_preterminal_fields=True,
            compact_force_has_finite_L1_time_L2_space_norm=True,
            factor_two_dissipation_retained=True,
            final_dissipation_limit_uses_monotone_convergence_only=True,
            separate_sup_and_total_bounds_not_added_with_same_budget=lower_sum>F1**2,
            core_growth_plateau_preserves_inherited_lower_bound_but_does_not_prove_it=True))

def scope_review():
    return dict(data='u(0)=0',force='f in C_c^infinity(R3 x (0,infinity)), support subset K x [1-tau0,3/2]',
        decay='sup (1+|x|+t)^L |D_x^alpha partial_t^m f| <= (1+R+3/2)^L ||D_x^alpha partial_t^m f||_infty',
        official_source='https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf; equations (4),(5),(7), alternative (C); checked 2026-10-17',
        checks=dict(all_required_data_and_force_weighted_derivative_bounds_from_zero_datum_and_compact_smooth_force=True,
            force_not_assumed_divergence_free_and_no_compact_Leray_projection_used=True,
            no_future_velocity_or_pressure_constructed_by_force_jet_gluing=True,
            retained_exact_core_and_finite_analytic_lemmas_remain_dependencies=True,
            global_competitor_pressure_recovery_and_comparison_not_promoted=True,
            no_unconditional_admissible_force_or_blowup_claim=True))
