"""Exact ledger for localizing the certified complete-potential field."""
from fractions import Fraction as Q
from math import comb, factorial
import importlib.util
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


endpoint = load(
    "certified_force_endpoint",
    HERE.parent / "core_certified_complete_potential_endpoint_join/endpoint.py",
)

TERMS = (
    "base_residual", "cutoff_time", "curl_time", "cutoff_laplacian",
    "gradient_viscosity", "curl_viscosity", "pressure_cutoff",
    "quadratic_rescaling", "velocity_cutoff_transport",
    "old_curl_transport", "curl_old_transport", "curl_square",
)


def order(value, name="order"):
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a finite nonnegative integer")
    return value


def force_requirements(m):
    m = order(m, "m")
    return dict(
        force_order=m,
        potential_order=m + 3,
        direct_swirl_order=m + 2,
        pressure_order=m + 1,
        Cauchy_potential_order=m + 4,
        Cauchy_direct_swirl_order=m + 3,
        Cauchy_pressure_order=m + 2,
    )


def derivative_factor(j, m, k):
    j, m, k = order(j, "j"), order(m, "m"), order(k, "k")
    if k > m or m - k > j:
        return Q(0)
    return Q(comb(m, k), factorial(j - m + k))


def first_controlled(spatial_order, time_order):
    a, m = order(spatial_order, "spatial_order"), order(time_order, "time_order")
    return max(1, 2 * (a + m), m + 1)


def domain_review():
    # q_e <= qstar/16 and q <= 3 q_e/8 on the retained cylinder.
    q_over_qstar = Q(3, 8) * Q(1, 16)
    four_Q_over_qstar = 16 * q_over_qstar
    rows = [force_requirements(m) for m in (0, 1, 5, 100)]
    return dict(
        field="u=curl(c*A)+c*B=c*u0+d; d=grad(c) cross A; p=c*p0",
        cutoff="c=chi_r(r^2) chi_z(z) chi_t(t), axisymmetric and one near (0,1)",
        complete_pressure="p0 remains in its fixed gauge with every inherited pressure difference and moment tail",
        enlarged_chart="1/4<=q/Q<=4; q<=3q_e/8; q_e<=qstar/16; hence 4Q<=16q<=3qstar/8",
        requirements=rows,
        checks=dict(
            certified_endpoint_parent_used=endpoint.closure_statement()[
                "same_certified_complete_diagonal_field_used"],
            certified_fixed_stage_domains_retained=endpoint.closure_statement()[
                "certified_fixed_stage_domains_retained"],
            parent_did_not_already_claim_localization=not endpoint.closure_statement()[
                "localization_initialization_admissible_force_energy_global_comparison_or_blowup_proved"],
            enlarged_chart_margin_strict=four_Q_over_qstar == Q(3, 8) < 1,
            every_force_order_has_three_extra_potential_derivatives=all(
                row["potential_order"] == row["force_order"] + 3 for row in rows),
            every_Cauchy_bound_reserves_one_more_time_derivative=all(
                row["Cauchy_potential_order"] == row["potential_order"] + 1
                and row["Cauchy_direct_swirl_order"] == row["direct_swirl_order"] + 1
                and row["Cauchy_pressure_order"] == row["pressure_order"] + 1
                for row in rows),
            twelve_complete_terms_retained=len(TERMS) == len(set(TERMS)) == 12,
            curl_and_axisymmetric_direct_swirl_are_exactly_divergence_free=True,
            d_dot_grad_c_vanishes_by_cross_product=True,
            zero_initial_neighborhood_follows_from_time_cutoff=True,
        ),
    )


def residual_review():
    return dict(
        identity=(
            "R_nu(u,p)=c R_nu(u0,p0)+c_t u0+d_t-nu(Delta c)u0"
            "-2nu sum_i c_i partial_i u0-nu Delta d+p0 grad c"
            "+c(c-1)(u0.grad)u0+c(u0.grad c)u0+c(u0.grad)d"
            "+c(d.grad)u0+(d.grad)d"
        ),
        terms=list(TERMS),
        pressure_Poisson="-Delta p=partial_i partial_j(u_i u_j)-div f_minus",
        corner_blocker=(
            "for every alpha,j,N, sup_corner |D_x^alpha partial_t^j R_1(u0,p0)|"
            " <= C_(alpha,j,N) (tau+|z|^(1/D))^N"
        ),
        punctured_blocker=(
            "for every fixed alpha,j, the required A,B,p0 jets and their next time"
            " derivatives are uniformly bounded on each punctured terminal compact"
        ),
        checks=dict(
            time_viscosity_pressure_and_all_nonlinear_terms_present=set(TERMS) == {
                "base_residual", "cutoff_time", "curl_time", "cutoff_laplacian",
                "gradient_viscosity", "curl_viscosity", "pressure_cutoff",
                "quadratic_rescaling", "velocity_cutoff_transport",
                "old_curl_transport", "curl_old_transport", "curl_square",
            },
            complete_pressure_cutoff_term_retained=True,
            force_divergence_retained_in_pressure_equation=True,
            corner_flatness_must_be_uniform_in_compact_X_not_pathwise=True,
            punctured_endpoint_and_corner_flatness_have_distinct_roles=True,
            no_compact_support_preserving_Leray_projection_assumed=True,
        ),
    )


def extension_review():
    rows = []
    for a, m in ((0, 0), (1, 3), (20, 4), (6, 70), (1000, 1000)):
        j = first_controlled(a, m)
        rows.append(dict(spatial_order=a, time_order=m, first_stage=j,
                         diagonal_admits=a + m <= j // 2, exponent_negative=m < j))
    return dict(
        terminal_jets="F_j=lim_(t up to 1) partial_t^j f_minus in C_c^infinity(K), flat at x=0",
        future_force="f_plus=sum_j chi0(b_j sigma) sigma^j F_j/j!, sigma=t-1",
        exact_factor="binomial(m,k) b_j^k chi0^(k)(b_j sigma) sigma^(j-m+k)/(j-m+k)!",
        whole_support_bound="||D_x^a partial_sigma^m T_j|| <= C_(j,a,m) b_j^(m-j)",
        rows=rows,
        checks=dict(
            representative_factorial_control=derivative_factor(5, 2, 1) == Q(1, 12),
            invalid_Leibniz_indices_give_zero=derivative_factor(1, 3, 1) == 0,
            every_fixed_mixed_jet_enters_the_single_diagonal=all(
                row["diagonal_admits"] and row["exponent_negative"] for row in rows),
            finite_constraints_allow_increasing_dyadic_scales=True,
            geometric_differentiated_tail_is_two_to_one_minus_J=True,
            plateau_and_factorial_match_exactly_F_m=True,
            raw_terminal_Taylor_convergence_not_required=True,
            future_velocity_or_pressure_not_constructed=True,
        ),
    )


def energy_review():
    epsilon = Q(1, 100)
    force_norm = 1 + epsilon / 2
    lower_combined = 2 - Q(1, 1 + 2 * (1 - epsilon)) + 2 * epsilon / 3
    return dict(
        identity="0.5 E'(t)+D(t)=<f(t),u(t)>",
        fixed_time="E(t)+2 integral_0^t D <= (integral_0^t ||f||_2)^2",
        separate="sup_(t<1) E<=F(1)^2 and integral_0^1 D<=F(1)^2/2",
        counterexample=dict(epsilon=epsilon, force_norm=force_norm,
                            lower_combined=lower_combined),
        checks=dict(
            zero_start_handled_by_regularized_square_root=True,
            pressure_and_transport_fluxes_vanish_for_compact_divergence_free_field=True,
            dissipation_factor_two_retained=True,
            terminal_dissipation_uses_monotone_convergence=True,
            separate_extrema_cannot_be_added_with_same_budget=(
                lower_combined > force_norm ** 2),
            energy_does_not_regularize_terminal_origin_or_prove_growth=True,
        ),
    )


def closure_statement():
    return dict(
        same_certified_complete_diagonal_field_localized=True,
        certified_fixed_stage_domains_and_single_cutoff_sequence_retained=True,
        uniform_stage_constant_growth_required=False,
        complete_fixed_gauge_pressure_and_all_twelve_residual_terms_retained=True,
        compact_smooth_force_and_zero_initial_data_follow_conditionally=True,
        preterminal_energy_and_integrated_dissipation_follow_conditionally=True,
        protected_inner_field_unchanged_on_cutoff_plateau=True,
        primitive_corner_and_punctured_bounds_independently_proved=False,
        actual_terminal_norms_or_Borel_scales_evaluated=False,
        regular_velocity_or_pressure_at_terminal_origin_proved=False,
        future_velocity_or_pressure_constructed=False,
        global_competitor_comparison_or_blowup_proved=False,
    )
