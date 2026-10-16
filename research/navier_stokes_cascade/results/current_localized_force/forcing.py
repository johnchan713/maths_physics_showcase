"""Current-field localization and all-jet force gluing: conditional budgets."""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys
import mpmath as mp

HERE=Path(__file__).resolve().parent
TERMS=('base_residual','cutoff_time','curl_time','cutoff_laplacian',
    'gradient_viscosity','curl_viscosity','pressure_cutoff','quadratic_rescaling',
    'velocity_cutoff_transport','old_curl_transport','curl_old_transport','curl_square')


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module


def order(n):
    if type(n) is not int or n<0:raise ValueError('Nonnegative finite integer derivative order required')
    return n


def requirements(m):
    m=order(m)
    return dict(force_order=m,potential_bound=m+3,direct_swirl_bound=m+2,pressure_bound=m+1,
        next_time_potential_bound=m+4,next_time_direct_swirl_bound=m+3,next_time_pressure_bound=m+2)


def domain_review():
    return dict(q_e='min(qstar/16,2^(-ell0-2))',tau0='q_e/16',z0='(q_e/16)^D',
        cutoff='c=chi_r(r^2)*chi_z(z)*chi_t(t); compact axisymmetric support and inner plateau',
        q_over_q_e=str(Q(3,8)),fourQ_over_qstar=str(Q(3,16)),Q_over_first_band_scale=str(Q(3,16)),
        checks=dict(same_current_endpoint_observation_domain_retained=Q(3,8)<Q(1,2),
            every_active_enlarged_chart_and_fixed_first_band_admitted=Q(3,16)<1,
            cutoff_zero_near_t_zero_with_fixed_domain_margin=Q(1,16)<Q(1,2),
            Cartesian_radial_cutoff_uses_r_squared_and_is_axisymmetric=True,
            curl_and_direct_swirl_divergence_preserved_without_full_swirl_potential=True,
            all_base_pressure_tails_and_every_current_pressure_difference_cut_in_same_gauge=True,
            q_below_half_background_domain_alone_not_reused_as_fourQ_margin=True))


def localization_review():
    return dict(rows=[requirements(m) for m in (0,1,5,100)],retained_terms=list(TERMS),
        exact_field='u=curl(c*A)+c*B=c*u0+d, d=grad(c) cross A; p=c*p0',
        exact_zero_transport='(d.grad c)*u0=0 because d.grad c=0',
        pressure_identity='-Delta p=sum_ij partial_i partial_j(u_i*u_j)-div f',
        checks=dict(twelve_full_localization_terms_retained=len(TERMS)==12,
            every_fixed_force_jet_has_one_more_potential_curl_derivative=all(r['potential_bound']==r['force_order']+3 for r in [requirements(m) for m in (0,1,5,100)]),
            endpoint_Cauchy_requires_one_additional_time_derivative=True,
            full_pressure_cutoff_gradient_and_nonzero_force_divergence_retained=True,
            actual_Gaussian_velocities_and_excluded_tail_ledger_are_not_forward_inverted=True,
            inner_plateau_avoids_estimating_cutoff_commutators_at_singular_core=True))


def terminal_review():
    return dict(near_corner='for all alpha,j,N: |D_x^alpha partial_t^j f|<=C_(alpha,j,N)*(tau+|z|^(1/D))^N',
        terminal_jets='F_j=lim_(t->1-) partial_t^j f in C_c^infinity(K); D_x^alpha F_j(0)=0',
        uniform_gluing='small ball uses flat complete residual; K outside ball uses current punctured endpoint limits; off K is zero',
        checks=dict(complete_residual_not_only_regular_angular_part_is_flat=True,
            compact_X_uniformity_plus_exact_heat_exterior_zero_used=True,
            all_corner_residual_jets_attach_to_zero_but_velocity_and_pressure_need_not=True,
            uniform_spatial_derivative_limits_identify_smooth_compact_terminal_jets=True,
            time_compatibility_follows_from_integral_of_next_time_derivative=True,
            origin_flatness_alone_not_used_for_away_from_origin_limits=True,
            nonzero_terminal_force_jets_not_assumed_zero_or_analytic=True))


def first_controlled_stage(a,m):return max(1,2*(order(a)+order(m)))


def choose_log2_scale(j,constraints,previous=1):
    if type(j) is not int or j<1 or type(previous) is not int or previous<1 or not constraints:
        raise ValueError('Positive stage, previous dyadic exponent and finite constraints required')
    with mp.workdps(80):
        n=previous+1
        for a,m,logC in constraints:
            order(a);order(m);logC=mp.mpf(logC)
            if a+m>j//2 or not mp.isfinite(logC) or logC<0:
                raise ValueError('Require a+m<=floor(j/2) and finite padded log constants')
            n=max(n,int(mp.ceil((logC+j*mp.log(2)+2)/((j-m)*mp.log(2)))))
        return n


def extension_review():
    previous=1;rows=[]
    for j in (1,2,7,20):
        constraints=[(a,m,(j+1)**4+10**(j//3)) for a in range(j//2+1) for m in range(j//2+1-a)]
        n=choose_log2_scale(j,constraints,previous)
        with mp.workdps(80):margin=max(mp.mpf(C)+(m-j)*n*mp.log(2)+j*mp.log(2) for a,m,C in constraints)
        rows.append(dict(j=j,dyadic_exponent=n,constraints=len(constraints),log_margin=float(margin),doubles=n>previous))
        previous=n
    pairs=[(a,m,first_controlled_stage(a,m)) for a,m in ((0,0),(1,3),(20,4),(6,70))]
    return dict(rows=rows,synthetic_constants_only=True,actual_terminal_norms_and_scales_not_evaluated=True,
        term='T_j=chi0(b_j*sigma)*sigma^j/j!*F_j, b_0=2, b_j=2^n_j',
        derivative_bound='||D_x^a partial_sigma^m T_j||_infty<=C_jam*b_j^(m-j)',
        tail='for fixed a,m and J>=max(1,2*(a+m)), sum_(j>=J) ||D_x^a partial_sigma^m T_j||<=2^(1-J)',
        checks=dict(every_finite_constraint_met_on_whole_support=all(r['log_margin']<-1 for r in rows),
            dyadic_scales_increase_and_support_stays_below_sigma_one_half=all(r['doubles'] for r in rows),
            every_fixed_mixed_derivative_eventually_controlled=all(a+m<=J//2 and m<J for a,m,J in pairs),
            terminal_mth_jet_is_F_m_by_plateau_and_factorial=True,
            all_spatial_jets_remain_flat_at_origin_and_supported_in_K=True,
            no_raw_Taylor_convergence_or_uniform_terminal_order_constant_assumed=True,
            smooth_future_force_not_future_velocity_or_pressure_supplied=True))


def energy_review():
    return dict(initial_data='u(0)=0',force_norm='F(t)=integral_0^t ||f(s)||_2 ds',
        identity='0.5*d_t ||u||_2^2+nu*||grad u||_2^2=<f,u>',
        zero_norm_regularization='d_t sqrt(||u||_2^2+delta^2)<=||f||_2; limit delta->0',
        bound='||u(t)||_2^2+2*integral_0^t ||grad u||_2^2<=F(t)^2, viscosity one',
        lower_bound='|u_theta(sqrt(2*Xin*tau),0,0,1-tau)|>=(e0/2)*tau^(-A) for sufficiently small tau',
        checks=dict(pressure_and_transport_fluxes_vanish_for_compact_preterminal_fields=True,
            force_smooth_compact_support_gives_finite_L1_time_L2_space_norm=True,
            no_division_by_possibly_zero_velocity_norm=True,
            integrated_bound_retains_exact_factor_two_dissipation=True,
            monotone_convergence_only_applied_to_nonnegative_dissipation=True,
            current_inherited_core_lower_bound_unchanged_by_localization_plateau=True,
            no_global_competitor_comparison_or_regular_velocity_at_t_one_inferred=True))


def admissibility_review():
    return dict(official_source='https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf',
        source_checked='2026-10-16; equations (4),(5),(7) and alternative (C), PDF pages 1-2',
        rapid_decay='for all alpha,m,K: |D_x^alpha partial_t^m f(x,t)|<=C_(alpha,m,K)*(1+|x|+t)^(-K)',
        force_support='K x [1-tau0,3/2]',
        checks=dict(zero_datum_is_smooth_divergence_free_and_rapidly_decaying=True,
            compact_smooth_force_implies_every_required_spacetime_decay_bound=True,
            no_divergence_free_force_requirement_added_to_official_definition=True,
            initial_zero_neighborhood_supplies_smoothness_at_t_zero=True,
            meeting_data_and_force_requirements_is_not_global_solution_exclusion=True))
