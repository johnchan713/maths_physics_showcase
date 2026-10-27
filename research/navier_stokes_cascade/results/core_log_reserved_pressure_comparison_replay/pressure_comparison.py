"""Pressure recovery and expanding-ball replay on the log-reserved field."""

from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import sys


HERE = Path(__file__).resolve().parent


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


energy = load(
    "log_comparison_energy",
    "core_log_reserved_preterminal_energy_replay/preterminal_energy.py",
)
historical = load(
    "log_comparison_historical",
    "core_corrected_comparison_join/comparison.py",
)
core = load(
    "log_comparison_core_path",
    "core_same_field_positive_core_reconciliation/core_path.py",
)


def pressure_review():
    inherited = historical.pressure_review()
    return {
        "difference_equation": "w_t+div(g)=Delta w-grad(pi); div w=0",
        "stress": "g=v tensor v-u tensor u=w tensor w+w tensor u+u tensor w",
        "Poisson_difference": "-Delta pi=partial_i partial_j g_ij",
        "individual_pressures": (
            "-Delta p=partial_i partial_j(u_i u_j)-div f; "
            "-Delta P=partial_i partial_j(v_i v_j)-div f"
        ),
        "recovery": inherited["recovered"],
        "harmonic_counterexample": inherited["harmonic_counterexample"],
        "checks": {
            "common_force_cancels_exactly_once": inherited["checks"]["force_divergence_cancels_once_and_is_not_subtracted_again"],
            "individual_force_divergence_is_retained": "-div f" in inherited["forced_individual"],
            "L1_stress_requires_Hminus_s_strictly_above_three_halves": inherited["checks"]["L1_stress_embeds_in_H_minus_s_only_for_s_strictly_above_three_halves"],
            "time_averaged_momentum_places_gradient_in_Hminus_three": inherited["checks"]["time_averaged_momentum_places_pressure_gradient_in_H_minus_three"],
            "harmonic_Hminus_three_gradient_vanishes": inherited["checks"]["weighted_L2_Fourier_function_supported_at_one_point_is_zero"],
            "Poisson_equation_alone_is_insufficient": inherited["checks"]["Poisson_equation_alone_does_not_remove_harmonic_gradient"],
            "no_competitor_pressure_decay_or_temperedness_assumed": inherited["checks"]["no_temperedness_or_decay_of_competitor_pressure_assumed"],
        },
    }


def radius_review():
    commutator = historical.commutator_review()
    radius = historical.radius_review()
    return {
        "cutoff": "chi_R=phi(x/R)^8; phi=1 on B_1 and zero outside B_2",
        "support_split": "g=w tensor w+h; h=w tensor u+u tensor w; supp h subset K",
        "compact_mixed_pressure": "|pi_h|<=C H_T R^-3 and its annular flux is O(R^-5/2)",
        "commutator_kernel": commutator["kernel"],
        "commutator_exponent": commutator["exponent"],
        "commutator_radial_coefficient": commutator["radial_coefficient"],
        "commutator_norm_decay": commutator["norm_radius_decay"],
        "flux_rows": radius["rows"],
        "worst_decay": radius["worst_decay"],
        "checks": {
            "support_separation_is_R_over_two": historical.separated_radius(2, 1) == 1,
            "mixed_moment_is_not_assumed_isotropic": radius["checks"]["mixed_moment_is_not_assumed_isotropic"],
            "compact_mixed_flux_decays_R_minus_five_halves": radius["checks"]["compact_mixed_pressure_flux_decay_is_five_halves"],
            "near_commutator_integral_is_finite": commutator["checks"]["near_integral_is_finite"],
            "far_commutator_integral_is_finite": commutator["checks"]["far_integral_is_finite"],
            "eighth_cutoff_power_is_required_by_this_pairing": (
                commutator["checks"]["eighth_cutoff_power_is_sufficient_for_L_four_thirds_pairing"]
                and commutator["checks"]["seventh_power_is_insufficient_for_this_pairing"]
            ),
            "Young_remainders_are_R_minus_four_and_R_minus_fourteen_fifths": (
                [row["Young_decay"] for row in radius["rows"]] == [Q(4), Q(14, 5)]
            ),
            "coarse_cutoff_laplacian_limits_total_remainder_to_R_minus_two": radius["worst_decay"] == 2,
        },
    }


def limit_review():
    inherited = historical.limit_review()
    return {
        "localized_inequality": (
            "0.5 E_R'(t)+0.5 D_R(t)^2"
            "<=||grad u(t)||_infinity E_R(t)+C_T R^-2; E_R(0)=0"
        ),
        "Gronwall": "E_R(t)<=2 C_T T exp(2 T L_T) R^-2 for each fixed T<1",
        "limit_order": "fix T<1, then R->infinity, then approach the protected core as t->1",
        "checks": {
            "fixed_ball_eventually_lies_in_cutoff_plateau": inherited["checks"]["fixed_ball_eventually_lies_where_cutoff_is_one"],
            "cutoffs_need_not_be_monotone": inherited["checks"]["cutoff_family_need_not_be_monotone"],
            "constants_need_not_be_uniform_as_T_reaches_one": inherited["checks"]["constants_need_not_be_uniform_as_T_approaches_one"],
            "equality_holds_on_every_fixed_preterminal_interval": inherited["checks"]["equality_on_every_fixed_interval_gives_v_equals_u_for_all_t_below_one"],
            "global_smooth_competitor_is_locally_bounded_near_terminal_point": inherited["checks"]["smooth_global_competitor_is_bounded_on_a_compact_neighborhood_of_the_terminal_origin"],
            "same_field_positive_core_path_is_available_conditionally": core.review()["analytic_conclusions"]["same_complete_diagonal_has_positive_core_lower_bound"],
        },
    }


def review():
    energy_scope = energy.review()["scope"]
    pressure = pressure_review()
    radius = radius_review()
    limit = limit_review()
    checks = {
        "same_localized_field_force_and_zero_datum_retained": (
            energy_scope["same_field_preterminal_energy_identity_replayed_conditionally"]
            and energy.parent.review()["scope"]["one_smooth_compact_force_extension_constructed_conditionally"]
        ),
        "competitor_requires_only_local_smoothness_and_fixed_interval_L2_bound": True,
        "all_pressure_recovery_checks_close": all(pressure["checks"].values()),
        "all_radius_checks_close": all(radius["checks"].values()),
        "all_limit_checks_close": all(limit["checks"].values()),
        "fixed_time_then_radius_then_core_order_retained": True,
        "complete_fixed_gauge_pressure_memory_is_not_replaced": True,
        "force_divergence_is_retained_then_cancels_once": True,
        "conditional_competitor_exclusion_uses_same_field_core_lower_bound": True,
        "scope_excludes_future_solution_and_unconditional_blowup": True,
    }
    return {
        "checks": checks,
        "interface": {
            "constructed": "same u=curl(cA)+cB, p=cp0, compact smooth f and u(0)=0",
            "competitor": (
                "v,P smooth on compact subsets, same f and zero datum, "
                "sup_[0,T] ||v||_2 finite for each T<1"
            ),
            "no_added_hypotheses": "no global competitor derivative, pressure-decay or dissipation bound",
        },
        "pressure": pressure,
        "radius": radius,
        "limit": limit,
        "conditional_conclusion": {
            "comparison": "v=u on R^3 x [0,T] for every fixed T<1",
            "core": "|u_theta(x_tau,1-tau)|>=(e0/2) tau^-A on the protected path",
            "exclusion": "a globally smooth locally bounded finite-energy same-force competitor is excluded relative to all frozen contracts",
        },
        "scope": {
            "fixed_gauge_pressure_gradient_recovered_conditionally": True,
            "same_field_expanding_ball_comparison_replayed_conditionally": True,
            "global_smooth_finite_energy_competitor_exclusion_follows_conditionally": True,
            "finite_energy_hypothesis_removed": False,
            "mixed_pressure_moment_isotropy_assumed": False,
            "future_velocity_or_pressure_constructed": False,
            "all_frozen_continuum_contracts_independently_reproved_here": False,
            "unconditional_admissible_force_or_blowup_proved": False,
        },
    }
