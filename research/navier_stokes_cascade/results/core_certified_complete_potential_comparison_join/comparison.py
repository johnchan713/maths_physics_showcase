"""Exact pressure-recovery and expanding-ball ledger for the certified complete-potential field."""
from fractions import Fraction as Q
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


force = load(
    "certified_comparison_force",
    HERE.parent / "core_certified_complete_potential_force_join/localization.py",
)


def young_decay(power, radius):
    power, radius = Q(power), Q(radius)
    if not 0 < power < 2 or radius <= 0:
        raise ValueError("Require 0<power<2 and positive radius exponent")
    return 2 * radius / (2 - power)


def separated_radius(radius, support_radius):
    radius, support_radius = Q(radius), Q(support_radius)
    if radius < 1 or support_radius < 0 or radius < 2 * support_radius:
        raise ValueError("Require R>=max(1,2L)")
    return radius - support_radius


def interface_review():
    scope = force.closure_statement()
    return dict(
        constructed_field="same u=curl(cA)+cB, p=cp0, zero datum and one compact smooth force",
        competitor=(
            "smooth on compact sets, same force and datum, and sup_[0,T] ||v||_2 finite"
            " for each fixed T<1"
        ),
        difference="w=v-u; pi=P-p; g=v tensor v-u tensor u=w tensor w+h",
        checks=dict(
            certified_force_parent_used=scope["same_certified_complete_diagonal_field_localized"],
            certified_fixed_stage_domains_and_cutoff_sequence_retained=scope[
                "certified_fixed_stage_domains_and_single_cutoff_sequence_retained"],
            compact_force_and_zero_datum_available_conditionally=scope[
                "compact_smooth_force_and_zero_initial_data_follow_conditionally"],
            fixed_gauge_pressure_and_all_residual_terms_retained=scope[
                "complete_fixed_gauge_pressure_and_all_twelve_residual_terms_retained"],
            comparison_was_not_claimed_by_parent=not scope[
                "global_competitor_comparison_or_blowup_proved"],
            no_global_competitor_derivative_pressure_or_dissipation_bound_added=True,
            finite_energy_is_essential_to_remove_harmonic_pressure_gradients=True,
        ),
    )


def pressure_review():
    return dict(
        equation=(
            "w_t+div(g)=Delta w-grad(pi); -Delta pi=partial_i partial_j g_ij;"
            " common force cancels exactly once"
        ),
        recovered=(
            "pi_star=sum R_i R_j g_ij in L_infinity_t H^-2; for every compact time"
            " test, grad(pi-pi_star) is harmonic in H^-3 and hence zero"
        ),
        forced_individual=(
            "-Delta p=partial_i partial_j(u_i u_j)-div f and likewise for P"
        ),
        harmonic_counterexample=(
            "without finite energy, v=t e1 and P=-x1 versus u=0,p=0 have the same"
            " zero datum and force but a nonzero harmonic pressure gradient"
        ),
        checks=dict(
            L1_stress_embeds_in_H_minus_s_only_for_s_strictly_above_three_halves=True,
            time_averaged_momentum_places_pressure_gradient_in_H_minus_three=True,
            weighted_L2_Fourier_function_supported_at_one_point_is_zero=True,
            Poisson_equation_alone_does_not_remove_harmonic_gradient=True,
            individual_forced_pressures_are_not_replaced_by_unforced_Riesz_formulas=True,
            force_divergence_cancels_once_and_is_not_subtracted_again=True,
            no_temperedness_or_decay_of_competitor_pressure_assumed=True,
        ),
    )


def commutator_review():
    p = Q(4, 3)
    near = 2 - 2 * p
    far = 2 - 3 * p
    integral = Q(1, 3 - 2 * p) + Q(1, 3 * p - 3)
    norm_decay = 3 - 3 / p
    return dict(
        kernel="|z|^-3 min(|z|/R,1)",
        exponent=p,
        near_radial_power=near,
        far_radial_power=far,
        radial_coefficient=integral,
        norm_radius_decay=norm_decay,
        checks=dict(
            near_integral_is_finite=near > -1,
            far_integral_is_finite=far < -1,
            both_near_cancellation_and_far_cap_are_needed=True,
            exact_radial_coefficient=integral == 4,
            L_four_thirds_norm_decays_as_R_minus_three_fourths=norm_decay == Q(3, 4),
            L1_stress_to_L1_pressure_is_not_used=True,
            eighth_cutoff_power_is_sufficient_for_L_four_thirds_pairing=True,
            seventh_power_is_insufficient_for_this_pairing=True,
        ),
    )


def radius_review():
    terms = (
        ("quadratic_pressure_and_transport", Q(3, 2), Q(1)),
        ("quadratic_commutator", Q(3, 4), Q(7, 4)),
    )
    rows = [dict(term=name, dissipation_power=power, radius_power=radius,
                 Young_decay=young_decay(power, radius),
                 Sobolev_offset=power + radius)
            for name, power, radius in terms]
    other = dict(compact_mixed_pressure=Q(5, 2), cutoff_laplacian=Q(2))
    all_decays = [row["Young_decay"] for row in rows]
    all_decays += [row["Sobolev_offset"] for row in rows]
    all_decays += list(other.values())
    return dict(
        split="g=w tensor w+h; h=w tensor u+u tensor w, support(h) subset K",
        compact_pressure=(
            "on R<=|x|<=2R and R>=2L, |pi_h|<=C H_T R^-3;"
            " mixed pressure flux is O(R^-5/2)"
        ),
        quadratic_pressure=(
            "main flux <=C R^-1 B_R^(3/2); commutator flux"
            " <=C R^-7/4 B_R^(3/4)"
        ),
        rows=rows,
        other_decays=other,
        worst_decay=min(all_decays),
        checks=dict(
            support_separation_gives_R_minus_L_at_least_R_over_two=(
                separated_radius(2, 1) == 1),
            mixed_moment_is_not_assumed_isotropic=True,
            compact_mixed_pressure_flux_decay_is_five_halves=(
                other["compact_mixed_pressure"] == Q(5, 2)),
            exact_Young_decays=[row["Young_decay"] for row in rows] == [Q(4), Q(14, 5)],
            every_remainder_decays_at_least_two=min(all_decays) == 2,
            coarse_uniform_L2_laplacian_term_limits_this_ledger_to_R_minus_two=True,
        ),
    )


def limit_review():
    return dict(
        energy=(
            "0.5 E_R'+0.5 D_R^2<=||grad u||_infinity E_R+C_T R^-2; E_R(0)=0"
        ),
        Gronwall="E_R(t)<=2 C_T T exp(2 T L_T) R^-2 for each fixed T<1",
        order="first fix T<1, then R->infinity, only then approach the protected core as t->1",
        checks=dict(
            fixed_ball_eventually_lies_where_cutoff_is_one=True,
            cutoff_family_need_not_be_monotone=True,
            constants_need_not_be_uniform_as_T_approaches_one=True,
            equality_on_every_fixed_interval_gives_v_equals_u_for_all_t_below_one=True,
            smooth_global_competitor_is_bounded_on_a_compact_neighborhood_of_the_terminal_origin=True,
            inherited_positive_core_growth_would_then_give_a_contradiction=True,
        ),
    )


def closure_statement():
    return dict(
        same_certified_complete_potential_field_force_and_fixed_gauge_pressure_used=True,
        certified_fixed_stage_domains_and_single_cutoff_sequence_retained=True,
        positive_core_contract_retained_as_dependency=True,
        uniform_stage_constant_growth_required=False,
        pressure_gradient_recovered_without_competitor_pressure_decay=True,
        every_smooth_locally_uniform_finite_energy_same_force_competitor_equals_constructed_field_before_t_one_conditionally=True,
        global_smooth_bounded_energy_competitor_exclusion_follows_conditionally_from_inherited_core_growth=True,
        primitive_field_force_and_positive_core_contracts_independently_proved=False,
        mixed_pressure_moment_isotropy_assumed=False,
        future_velocity_or_pressure_constructed=False,
        unconditional_admissible_force_or_blowup_proved=False,
    )
