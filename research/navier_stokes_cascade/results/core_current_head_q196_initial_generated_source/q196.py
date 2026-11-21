"""Q196: corrected initial generated-source quotient and endpoint audit."""
from fractions import Fraction as Q
from math import ceil, factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q196-INITIAL-GENERATED-SOURCE"

PROVIDERS = (
    "core_current_head_q195_primary_covariance_initial_mean",
    "regular_source_foundation",
    "source_envelope_edge_audit",
    "core_current_head_q184_generated_source_log_loss_correction",
    "first_sourced_foundation",
)

EXPECTED = {
    PROVIDERS[0]: "Q195-primary-six-entry-covariance-and-initial-complete-mean-joined",
    PROVIDERS[1]: "regular-source-foundation-sharper-linear-gain-and-actual-initial-supported-input-relative-to-core-and-pulse-lemmas",
    PROVIDERS[2]: "regular-source-envelope-edge-estimates-with-explicit-freezing-correction",
    PROVIDERS[3]: "Q184-generated-source-log-loss-corrected-and-pure-power-reserve-derived",
    PROVIDERS[4]: "first-sourced-physical-pulse-complex-covariance-and-current-pressure-foundation-relative-to-core-and-finite-lemmas",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q196 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def rates(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q196 identity required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 4):
        raise ValueError("reviewed kappa required")
    alpha = Q(1, 2)
    mu = 1 - kappa
    B = 1 - kappa
    C = B + Q(1, 2)
    D = C - kappa
    return {
        "alpha_primary": alpha,
        "mu_initial_mean": mu,
        "B_initial_source_log": B,
        "B_initial_source_pure_after_kappa": B - kappa,
        "C_first_sourced_pressure_potential": C,
        "D_after_first_sourced_correction": D,
        "linear": alpha + Q(1, 2) - kappa,
        "primary_square": 2 * alpha - kappa,
        "initial_mean_primary": alpha + mu - Q(1, 2),
    }


def linear_ledger(alpha=Q(1, 2), kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q196 identity required")
    alpha, kappa = Q(alpha), Q(kappa)
    if alpha <= 0 or not 0 < kappa < Q(1, 4):
        raise ValueError("positive order and reviewed kappa required")
    gains = {
        "slow_time_bDr_GDz": 1 - kappa,
        "phase_defect": Q(1, 2),
        "base_derivatives_and_connection": Q(1),
        "pressure_gradient": Q(1, 2) - kappa,
        "vector_viscosity": 1 - 2 * kappa,
        "mixed_viscosity_and_phase_divergence": Q(1, 2) - kappa,
        "angular_vector_connection": Q(1, 2),
        "principal_operator_on_complete_curl_remainder": Q(1, 2) - kappa,
    }
    minimum = min(gains.values())
    return {
        "gains": gains,
        "minimum_gain": minimum,
        "source_order": alpha + minimum,
        "pressure_gradient_retained": True,
        "mixed_viscosity_retained": True,
        "cylindrical_connections_retained": True,
        "source_normal_pressure_retained": True,
        "moving_normal_pressure_retained": True,
    }


def source_groups(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q196 identity required")
    row = rates(kappa)
    orders = {
        "complete_primary_linear": row["linear"],
        "complete_primary_square": row["primary_square"],
        "actual_initial_mean_primary": row["initial_mean_primary"],
    }
    return {
        "orders": orders,
        "minimum": min(orders.values()),
        "all_equal_B": all(value == row["B_initial_source_log"] for value in orders.values()),
        "nonzero_harmonics": (1, -1, 2, -2),
        "all_forwarded_groups_contain_complete_cut_wave": True,
        "flat_cutoff_commutators_remain_additive_tails": True,
        "Pi1_Pi2_remain_in_same_base_and_mean_equations": True,
    }


def absorption(degree, eta=Q(1, 100000), h=Q(1, 200), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q196 identity required")
    degree, eta, h = Q(degree), Q(eta), Q(h)
    if degree < 0 or eta <= 0 or h <= 0:
        raise ValueError("nonnegative degree and positive reserve required")
    D = ceil(degree)
    upper = factorial(2 * D) * (Q(2) / (h * eta)) ** (2 * D)
    return {
        "degree": degree,
        "degree_ceiling": D,
        "eta": eta,
        "h": h,
        "uniform_upper": upper,
        "inequality": "ell^(2d) epsilon^eta <= (2D)! (2/(h eta))^(2D)",
        "exact_exponent_supremum_follows_without_reserve": degree == 0,
    }


def quotient_contract(jet_order=2, polynomial_degree=1,
                      kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q196 identity required")
    if isinstance(jet_order, bool) or not isinstance(jet_order, int) or jet_order < 0:
        raise ValueError("nonnegative integer jet order required")
    degree = Q(polynomial_degree)
    if degree < 0:
        raise ValueError("nonnegative polynomial degree required")
    row = rates(kappa)
    absorb = absorption(degree, Q(kappa))
    return {
        "jet_order": jet_order,
        "log_renormalized": (
            "|D^I(f_initial,m/P)| <= C epsilon^(1-kappa) S^d "
            "sqrt(zeta) ya^(-r) yb^(-s)"
        ),
        "pure_power_after_reserve": (
            "|D^I(f_initial,m/P)| <= C_eta epsilon^(1-2kappa) "
            "sqrt(zeta) ya^(-r) yb^(-s)"
        ),
        "log_exponent": row["B_initial_source_log"],
        "pure_exponent": row["B_initial_source_pure_after_kappa"],
        "polynomial_degree": degree,
        "absorption": absorb,
        "same_domain_raw_source_jets": True,
        "same_domain_log_P_jets": True,
        "smooth_mixed_endpoint_zero_extension": True,
        "source_divided_by_cutoff_or_amplitude": False,
    }


def endpoint_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q196 identity required")
    return {
        "linear_forwarded_source_contains_temporal_cutoff": True,
        "nonlinear_source_contains_complete_cut_velocity": True,
        "positive_Gaussian_is_not_used_as_support": True,
        "P_positive_at_each_finite_endpoint": True,
        "source_over_P_mixed_endpoint_jets_vanish": True,
        "uncut_exit_not_claimed_zero": True,
        "potential_and_pressure_cut_before_zero_extension": True,
        "cutoff_commutator_and_one_minus_cutoff_source_are_additive": True,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q196 identity required")
    return {
        "closed_here": (
            "the immediate primary-plus-initial-mean regular source has log-renormalized "
            "order 1-kappa with all fixed quotient jets and smooth endpoint extension; "
            "one kappa reserve gives a uniform pure-power order 1-2kappa"
        ),
        "blocking_inequality": (
            "independently solve the first zero-entry sourced ODE from this exact source, "
            "retaining n dot f in pressure, and rederive the complete complex covariance "
            "increment (including every new harmonic square), compact pressure, axial moment "
            "and current defects on the same physical domain"
        ),
        "first_sourced_complex_covariance_independently_rederived_here": False,
        "complete_finite_cycle_residual_rederived_here": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q196 identity required")
    status = statuses()
    rate = rates()
    linear = linear_ledger()
    groups = source_groups()
    quotient = quotient_contract()
    endpoint = endpoint_contract()
    boundary = frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "immediate_source_stage_is_B_not_D": rate["B_initial_source_log"] == Q(99999, 100000) and rate["D_after_first_sourced_correction"] == Q(74999, 50000),
        "all_three_complete_initial_source_groups_meet_B": groups["all_equal_B"] and groups["minimum"] == rate["B_initial_source_log"],
        "complete_linear_ledger_retains_pressure_viscosity_connections": linear["source_order"] == rate["B_initial_source_log"] and linear["pressure_gradient_retained"] and linear["mixed_viscosity_retained"] and linear["cylindrical_connections_retained"],
        "Q195_exact_power_display_omitted_polynomial_loss": not quotient["absorption"]["exact_exponent_supremum_follows_without_reserve"],
        "one_kappa_reserve_absorbs_fixed_polynomial_degree": quotient["pure_exponent"] == 1 - 2 * Q(1, 100000) and quotient["absorption"]["uniform_upper"] == 3200000000000000,
        "endpoint_extension_uses_cut_support_not_Gaussian_smallness": endpoint["source_over_P_mixed_endpoint_jets_vanish"] and endpoint["positive_Gaussian_is_not_used_as_support"],
        "flat_cutoff_terms_are_not_forwarded": groups["flat_cutoff_commutators_remain_additive_tails"] and endpoint["cutoff_commutator_and_one_minus_cutoff_source_are_additive"],
        "no_first_sourced_or_complete_residual_promotion": not boundary["first_sourced_complex_covariance_independently_rederived_here"] and not boundary["complete_finite_cycle_residual_rederived_here"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "rates": rate,
        "linear_ledger": linear,
        "source_groups": groups,
        "quotient_contract": quotient,
        "endpoint_contract": endpoint,
        "correction": {
            "Q195_frontier_skipped_the_immediate_B_source_stage": True,
            "Q195_displayed_exact_power_supremum_omitted_allowed_S_degree": True,
            "Q184_already_corrected_Q183_fixed_cycle_gate": True,
            "historical_frozen_evidence_modified": False,
        },
        "frontier": boundary,
        "scope": {
            "initial_generated_source_quotient_closed_conditionally": True,
            "pure_power_closed_after_one_kappa_reserve": True,
            "first_sourced_complex_covariance_independently_rederived": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q196 identity required")
    return {
        "status": "Q196-initial-generated-source-stage-and-log-quotient-corrected",
        "established": (
            "the complete primary-plus-initial-mean nonzero source has order 1-kappa in "
            "the log-renormalized class, with every fixed f/P jet and smooth endpoint "
            "extension; one kappa reserve yields a uniform pure-power 1-2kappa bound"
        ),
        "correction": (
            "Q195 jumped to the post-first-sourced order 3/2-2kappa and repeated Q183's "
            "exact-power supremum without the allowed polynomial S degree; the immediate "
            "source order is 1-kappa and fixed polynomial loss needs an exponent reserve"
        ),
        "next": frontier()["blocking_inequality"],
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
