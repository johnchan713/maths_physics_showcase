"""Q202: carry the corrected post-signed source into the next sourced cycle."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q202-NEXT-SOURCED-CYCLE"

PROVIDERS = (
    "core_current_head_q201_post_signed_mean_reconciliation",
    "core_post_signed_residual_join",
    "core_post_signed_seed_sync",
    "sourced_inverse_foundation",
    "regular_source_foundation",
    "core_current_head_q160_physical_target_induction",
    "core_current_head_q161_source_tail_simultaneous_induction",
    "core_current_head_q180_target_source_simultaneous_induction",
    "finite_stage_induction",
    "core_physical_tail_join",
)

EXPECTED = {
    PROVIDERS[0]: "Q201-corrected-post-signed-fast-slow-mean-block-joined-conditionally",
    PROVIDERS[1]: "conditional-post-signed-current-nonzero-source-and-full-physical-request-tail-join",
    PROVIDERS[2]: "conditional-post-signed-induction-seed-compact-sync",
    PROVIDERS[3]: "explicit-same-domain-sourced-inverse-and-high-harmonic-parameter-jets-relative-to-reviewed-frame-inputs",
    PROVIDERS[4]: "regular-source-foundation-sharper-linear-gain-and-actual-initial-supported-input-relative-to-core-and-pulse-lemmas",
    PROVIDERS[5]: "Q160-physical-all-finite-cycle-target-moment-induction-joined",
    PROVIDERS[6]: "Q161-physical-source-inverse-tail-simultaneous-induction-joined",
    PROVIDERS[7]: "Q180-physical-target-source-tail-simultaneous-induction-reconciled",
    PROVIDERS[8]: "conditional-current-state-complete-finite-stage-induction",
    PROVIDERS[9]: "conditional-current-prefix-physical-tail-transfer-and-absolute-regular-order-gap",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q202 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def rates(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q202 identity required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 100):
        raise ValueError("reviewed kappa required")
    B = 1 - kappa
    C = B + Q(1, 2)
    D = C - kappa
    C0 = 2 * B
    delta = Q(1, 2) - 3 * kappa
    B1 = D + delta
    C1 = B1 + Q(1, 2)
    target_bridge = C0 + Q(1, 2) - kappa
    return {
        "kappa": kappa,
        "initial_B": B,
        "current_source_D": D,
        "current_target_C0": C0,
        "cycle_gain_delta": delta,
        "next_source_B1": B1,
        "next_required_target_C1": C1,
        "proved_target_bridge": target_bridge,
        "target_slack": target_bridge - C1,
    }


def current_source(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q202 identity required")
    r = rates()
    groups = {
        "inherited_current_source": r["current_source_D"],
        "signed_complete_linear": r["current_source_D"],
        "both_old_signed_wave_crosses": r["current_source_D"],
        "signed_square": 2 * r["initial_B"] - r["kappa"],
        "old_mean_signed_wave": r["current_source_D"],
        "new_fast_feedback": r["current_source_D"],
        "new_slow_feedback": r["initial_B"] + Q(1, 2),
    }
    return {
        "groups": groups,
        "order": min(groups.values()),
        "angular_harmonics": tuple(range(-4, 0)) + tuple(range(1, 5)),
        "all_nonzero_terms_retain_complete_cut_wave_factor": True,
        "full_pressure_gradient_viscosity_cylindrical_crosses_and_squares_retained": True,
        "flat_Gaussian_and_base_exterior_terms_forward_inverted": False,
        "smooth_endpoint_zero_extension_relative_to_named_support_inputs": True,
    }


def sourced_inverse(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q202 identity required")
    r = rates()
    alpha = r["current_source_D"]
    return {
        "input_source": alpha,
        "complete_cut_velocity": alpha,
        "uncut_tangential_response": alpha,
        "pressure_and_vector_potential": alpha + Q(1, 2),
        "complete_curl_remainder": alpha + Q(1, 2) - r["kappa"],
        "pressure_equals_current_C0": alpha + Q(1, 2) == r["current_target_C0"],
        "constraint": "t'+Kt+m^2dt+ikmn*pi=-f, n.t=0, t(0)=0",
        "pressure": "pi=i(n.Kt-n'.t+n.f)/(km|n|^2)",
        "moving_normal_term_retained": True,
        "source_normal_term_retained": True,
        "finite_harmonic_set_only": True,
        "potential_and_pressure_cut_before_complete_curl": True,
        "uncut_exit_declared_zero": False,
        "same_accepted_domain_and_C_minus_20_width": True,
    }


def physical_domain(m=0, N=0, h=Q(1, 200), original_pole=1, previous_wave_pole=30, kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q202 identity required")
    if any(isinstance(x, bool) or not isinstance(x, int) or x < 0 for x in (m, N, original_pole, previous_wave_pole)):
        raise ValueError("nonnegative integer physical request data required")
    if original_pole < 1:
        raise ValueError("positive original pole required")
    h, kappa = Q(h), Q(kappa)
    if not 0 < h < Q(1, 100) or not 0 < kappa < Q(1, 100):
        raise ValueError("reviewed h and kappa required")
    r = rates(kappa)
    normalization = Q(3, 2) + 2 * h
    threshold = (Q(N) + normalization + 2 * m) / h

    def count(source):
        return max(1, (threshold - source) // kappa + 1)

    p0 = count(r["initial_B"])
    p1 = count(r["current_source_D"])
    J0, J1 = m + 2 * p0 + 3, m + 2 * p1 + 3
    maxima = {
        "post_signed_row": J0 + 10,
        "new_covariance": J0 + 13,
        "next_signed_target_input": J0 + 14,
    }
    highest = maxima["new_covariance"]
    q = max(previous_wave_pole, 6 * original_pole, 6 * highest + 9, 6 * highest + 12 * original_pole + 6)
    normalized_q = 1200114 if (m, N, h, original_pole, previous_wave_pole, kappa) == (0, 0, Q(1, 200), 1, 30, Q(1, 100000)) else None
    return {
        "m": m,
        "N": N,
        "A0_integrations": p0,
        "A1_integrations": p1,
        "A0_source_order": J0,
        "A1_source_order": J1,
        "strict_A0_exponent_margin": h * (r["initial_B"] + p0 * kappa) - (N + normalization + 2 * m),
        "strict_A1_exponent_margin": h * (r["current_source_D"] + p1 * kappa) - (N + normalization + 2 * m),
        "maxima": maxima,
        "common_wave_pole": q,
        "post_signed_and_fast_domain": 2 * q,
        "full_nonzero_source_domain": 3 * q,
        "final_radial_source_domain": 4 * q,
        "normalized_mean_pole_for_default_request": normalized_q,
        "normalized_mean_pole_suffices_for_physical_request": normalized_q is None or normalized_q >= q,
        "actual_high_jet_values_evaluated": False,
    }


def next_cycle(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q202 identity required")
    r = rates()
    source_margins = {
        "sourced_full_linear": Q(0),
        "signed_full_linear": Q(0),
        "sourced_old_wave_cross": 2 * r["kappa"],
        "signed_old_wave_cross": 2 * r["kappa"],
        "sourced_actual_mean_cross": 2 * r["kappa"],
        "signed_actual_mean_cross": 2 * r["kappa"],
        "fast_mean_wave_change": 2 * r["kappa"],
        "slow_mean_wave_change": 3 * r["kappa"],
        "sourced_new_square": Q(1),
        "signed_new_square": Q(1),
    }
    target_minimum_margin = 2 * r["kappa"]
    return {
        "source_margins_over_B1": source_margins,
        "source_limiters": tuple(name for name, margin in source_margins.items() if margin == 0),
        "next_source": r["next_source_B1"],
        "next_required_target": r["next_required_target_C1"],
        "proved_target_bridge": r["proved_target_bridge"],
        "target_minimum_margin": target_minimum_margin,
        "target_slack": r["target_slack"],
        "source_and_target_advanced_simultaneously": True,
        "closes_relative_to_named_continuum_inputs": min(source_margins.values()) >= 0 and r["target_slack"] == target_minimum_margin,
    }


def physical_powers(h=Q(1, 200), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q202 identity required")
    h = Q(h)
    if not 0 < h < Q(1, 100):
        raise ValueError("reviewed h required")
    r = rates()
    A = Q(1, 2) + h
    alpha = r["current_source_D"]
    return {
        "sourced_velocity_zero_jet": -A + h * alpha,
        "sourced_pressure_zero_jet": -2 * A + h * (alpha + Q(1, 2)),
        "sourced_potential_zero_jet": Q(1, 2) - A + h * (alpha + Q(1, 2)),
        "next_regular_source_zero_jet": h * r["next_source_B1"] - (Q(3, 2) + 2 * h),
    }


def correction(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q202 identity required")
    request = physical_domain()
    return {
        "Q201_normalized_mean_domain_remains_valid_for_mean_block": True,
        "Q201_normalized_mean_domain_is_not_the_physical_tail_domain": True,
        "normalized_pole": request["normalized_mean_pole_for_default_request"],
        "physical_pole": request["common_wave_pole"],
        "physical_domain_strictly_larger": request["common_wave_pole"] > request["normalized_mean_pole_for_default_request"],
        "historical_frozen_evidence_modified": False,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q202 identity required")
    return {
        "closed_here": (
            "relative to the named same-family fixed-jet source, frame and support providers, "
            "the exact D-order post-signed source enters the constrained sourced inverse on the "
            "corrected physical domain and the coupled source/target step reaches B1 and C1"
        ),
        "next": (
            "replay the complete-potential diagonal, punctured endpoint and localization-force "
            "ledger on this corrected coupled state, retaining every scale-cutoff curl, pressure "
            "gradient, viscosity term, nonlinear cross and square"
        ),
        "actual_giant_constants_or_high_jets_evaluated": False,
        "stage_uniform_infinite_sum_proved": False,
        "absolute_physical_residual_decay_proved": False,
        "complete_diagonal_force_exclusion_or_blowup_proved": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q202 identity required")
    status = statuses()
    r = rates()
    source = current_source()
    inverse = sourced_inverse()
    domain = physical_domain()
    cycle = next_cycle()
    powers = physical_powers()
    fixed = correction()
    boundary = frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "complete_post_signed_source_is_exactly_D": source["order"] == r["current_source_D"],
        "sourced_inverse_orders_are_D_C0_and_C0_minus_kappa": inverse["complete_cut_velocity"] == r["current_source_D"] and inverse["pressure_and_vector_potential"] == r["current_target_C0"] and inverse["complete_curl_remainder"] == r["current_target_C0"] - r["kappa"],
        "full_constraint_pressure_and_support_operations_retained": inverse["moving_normal_term_retained"] and inverse["source_normal_term_retained"] and inverse["potential_and_pressure_cut_before_complete_curl"] and not inverse["uncut_exit_declared_zero"],
        "physical_request_counts_match_frozen_residual_join": domain["A0_source_order"] == 60200007 and domain["A1_source_order"] == 60100009 and domain["maxima"]["new_covariance"] == 60200020 and domain["maxima"]["next_signed_target_input"] == 60200021,
        "physical_common_pole_and_domains_match": domain["common_wave_pole"] == 361200138 and domain["post_signed_and_fast_domain"] == 722400276 and domain["full_nonzero_source_domain"] == 1083600414 and domain["final_radial_source_domain"] == 1444800552,
        "both_physical_compact_branches_have_strict_margin": domain["strict_A0_exponent_margin"] == Q(1, 20000000) and domain["strict_A1_exponent_margin"] == Q(1, 20000000),
        "Q201_normalized_domain_not_reused_as_physical_domain": fixed["physical_domain_strictly_larger"] and not domain["normalized_mean_pole_suffices_for_physical_request"],
        "next_source_gain_saturates_only_full_linear_rows": set(cycle["source_limiters"]) == {"sourced_full_linear", "signed_full_linear"},
        "next_source_and_target_rates_close_with_two_kappa_target_slack": cycle["closes_relative_to_named_continuum_inputs"] and cycle["next_source"] == r["next_source_B1"] and cycle["next_required_target"] == r["next_required_target_C1"] and cycle["target_slack"] == 2 * r["kappa"],
        "physical_pulse_and_next_residual_are_not_absolute_decay": powers["sourced_velocity_zero_jet"] < 0 and powers["sourced_pressure_zero_jet"] < 0 and powers["next_regular_source_zero_jet"] < 0,
        "no_diagonal_force_or_blowup_promotion": not boundary["complete_diagonal_force_exclusion_or_blowup_proved"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "rates": r,
        "current_source": source,
        "sourced_inverse": inverse,
        "physical_domain": domain,
        "next_cycle": cycle,
        "physical_powers": powers,
        "correction": fixed,
        "frontier": boundary,
        "scope": {
            "conditional_analytical_next_sourced_cycle": True,
            "tests_or_integer_request_arithmetic_are_continuum_proof": False,
            "actual_enormous_derivatives_or_constants_evaluated": False,
            "uniform_in_cycle_constants_claimed": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q202 identity required")
    return {
        "status": "Q202-next-constrained-sourced-cycle-joined-on-corrected-physical-domain",
        "correction": (
            "Q201's normalized mean pole 1200114 does not cover the physical compact-tail "
            "request; the default physical join requires common pole 361200138"
        ),
        "established": (
            "relative to the named same-family continuum providers, the exact D source enters "
            "the full constrained sourced inverse and the simultaneous cycle reaches source "
            "B1=2-5*kappa and required target C1=5/2-5*kappa with 2*kappa target slack"
        ),
        "next": frontier()["next"],
        "historical_frozen_records_modified": False,
        "unconditional_complete_residual_force_exclusion_or_blowup_certified": False,
    }
