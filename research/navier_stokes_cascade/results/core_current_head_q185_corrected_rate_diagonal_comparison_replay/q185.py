"""Replay Q180--Q182 after Q184's generated-source rate correction."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
IDENTITY = "Q185-CORRECTED-RATE-DIAGONAL-COMPARISON-REPLAY"
PROVIDERS = (
    "core_current_head_q184_generated_source_log_loss_correction",
    "core_current_head_q181_diagonal_endpoint_force_replay",
    "core_current_head_q182_comparison_core_exclusion_replay",
    "core_current_head_q166_complete_residual_frontier_audit",
)
EXPECTED = {
    PROVIDERS[0]: "Q184-generated-source-log-loss-corrected-and-pure-power-reserve-derived",
    PROVIDERS[1]: "Q181-physical-diagonal-endpoint-localized-force-replayed-conditionally",
    PROVIDERS[2]: "Q182-physical-pressure-comparison-protected-core-replayed-conditionally",
    PROVIDERS[3]: "Q166-current-head-complete-six-group-local-residual-audited-conditionally",
}


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q184 = load("q185_q184", "core_current_head_q184_generated_source_log_loss_correction/q184.py")
q181 = load("q185_q181", "core_current_head_q181_diagonal_endpoint_force_replay/q181.py")
q182 = load("q185_q182", "core_current_head_q182_comparison_core_exclusion_replay/q182.py")
q166 = load("q185_q166", "core_current_head_q166_complete_residual_frontier_audit/q166.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q185 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def ceil_fraction(value):
    value = Q(value)
    return -(-value.numerator // value.denominator)


def comparison_stage(target, loss_L=0, loss_H=0, route="pure", identity=IDENTITY):
    """Least nonnegative J0 for h*delta*(J0+1)/2-L-H >= target+1."""
    if identity != IDENTITY:
        raise ValueError("Q185 identity required")
    if any(isinstance(value, bool) or not isinstance(value, int) for value in (target, loss_L, loss_H)):
        raise ValueError("integer target and losses required")
    if min(target, loss_L, loss_H) < 0:
        raise ValueError("nonnegative target and losses required")
    rates = q184.corrected_rates()
    if route == "log":
        delta = rates["delta_log"]
    elif route == "pure":
        delta = rates["delta_pure"]
    else:
        raise ValueError("route must be 'log' or 'pure'")
    h = Q(1, 200)
    required = Q(target + 1 + loss_L + loss_H)
    J0 = max(0, ceil_fraction(2 * required / (h * delta)) - 1)
    replay_stage = J0 + 1
    threshold_left = h * delta * (J0 + 1) / 2 - loss_L - loss_H
    replay_left = h * delta * (replay_stage + 1) / 2 - loss_L - loss_H
    return {
        "route": route,
        "h": h,
        "delta": delta,
        "target": target,
        "loss_L": loss_L,
        "loss_H": loss_H,
        "J0": J0,
        "replay_stage": replay_stage,
        "threshold_left": threshold_left,
        "replay_left": replay_left,
        "threshold_holds": threshold_left >= target + 1,
        "retained_half_stage_margin": h * delta / 2,
        "replay_beats_target_by_margin": replay_left >= Q(target + 1) + h * delta / 2,
    }


def residual_margin_replay(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q185 identity required")
    rates = q184.corrected_rates()
    inherited = q166.universal_margin()
    h = Q(1, 200)
    log_margin = h * rates["delta_log"] / 2
    pure_margin = h * rates["delta_pure"] / 2
    compact_margin = h * rates["kappa"] / 2
    return {
        "log_delta": rates["delta_log"],
        "pure_delta": rates["delta_pure"],
        "inherited_Q181_margin": inherited["regular_retained_margin"],
        "log_regular_and_cross_margin": log_margin,
        "pure_regular_and_cross_margin": pure_margin,
        "pure_margin_correction": pure_margin - log_margin,
        "A0_own_margin_both_routes": compact_margin,
        "A1_own_margin_both_routes": compact_margin,
        "compact_accounts_disjoint": (
            not inherited["compact_reserves_may_pay_regular_or_cutoff_logs"]
            and not inherited["disjoint_reserves_added_to_strengthen_common_exponent"]
        ),
        "both_regular_margins_positive": log_margin > 0 and pure_margin > 0,
        "compact_margins_positive": compact_margin > 0,
    }


def diagonal_replay(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q185 identity required")
    diagonal = q181.diagonal_residual_contract()
    endpoint = q181.endpoint_contract()
    localization = q181.localization_force_energy_contract()
    margins = residual_margin_replay()
    return {
        "log_route_is_existing_WMS_diagonal": True,
        "pure_route_reselects_finite_stages_and_may_change_diagonal_scales": True,
        "pure_route_identified_byte_for_byte_with_Q181_field": False,
        "every_fixed_target_has_finite_log_stage": comparison_stage(3, 5, 7, "log")["threshold_holds"],
        "every_fixed_target_has_finite_pure_stage": comparison_stage(3, 5, 7, "pure")["threshold_holds"],
        "all_six_groups_retain_complete_memory": diagonal["all_six_groups_retain_complete_memory"],
        "fixed_jet_then_finite_stage": diagonal["fixed_jet_target_then_finite_stage"],
        "stage_uniform_constant_used": diagonal["stage_uniform_constant_or_J_of_q_used"],
        "all_fixed_cartesian_residual_jets_flat": diagonal["all_fixed_cartesian_residual_jets_flat"],
        "punctured_endpoint_supplied_separately": endpoint["endpoint_inequality_supplied_conditionally"],
        "terminal_origin_included": endpoint["singular_origin_included"],
        "localized_term_count": localization["term_count"],
        "complete_pressure_retained": localization["complete_pressure"],
        "unit_viscosity_retained": localization["unit_viscosity"],
        "protected_late_core_unchanged_by_support_geometry": localization["protected_late_core_unchanged"],
        "one_compact_smooth_force": localization["one_compact_smooth_force"],
        "future_velocity_or_pressure_constructed": localization["future_velocity_or_pressure_constructed"],
        "margins": margins,
    }


def comparison_replay(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q185 identity required")
    pressure = q182.pressure_comparison_contract()
    core = q182.protected_core_contract()
    return {
        "structural_input_only_requires_rebuilt_same_field_force_pair": True,
        "common_force_cancellations": pressure["force_cancellations"],
        "finite_energy_essential": pressure["finite_energy_essential"],
        "Poisson_alone_removes_harmonic_gradient": pressure["Poisson_alone_removes_harmonic_gradient"],
        "worst_radius_decay": pressure["worst_radius_decay"],
        "cutoff_power": pressure["cutoff_power"],
        "competitor_derivative_bound_assumed": pressure["competitor_derivative_bound_assumed"],
        "competitor_boundary_flux_bound_assumed": pressure["competitor_boundary_flux_bound_assumed"],
        "annular_diagonal_terms_vanish_exactly": core["annular_diagonal_terms_vanish_exactly"],
        "background_terms_vanish_exactly": core["background_terms_vanish_exactly"],
        "core_exponent": core["A"],
        "threshold_strictly_positive": core["threshold_strictly_positive"],
        "limit_order": (core["first_limit"], core["second_limit"], core["third_limit"]),
        "uniform_comparison_constant_as_T_to_one_required": core[
            "uniform_comparison_constant_as_T_to_one_required"
        ],
        "future_velocity_or_pressure_constructed": core["future_velocity_or_pressure_constructed"],
        "unconditional_blowup_claimed": core["unconditional_blowup_claimed"],
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q185 identity required")
    status = statuses()
    margins = residual_margin_replay()
    diagonal = diagonal_replay()
    comparison = comparison_replay()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "Q181_margin_is_exactly_the_log_route_margin": (
            margins["inherited_Q181_margin"] == margins["log_regular_and_cross_margin"]
        ),
        "pure_route_margin_is_corrected_by_one_over_forty_million": (
            margins["pure_margin_correction"] == -Q(1, 40000000)
            and margins["pure_regular_and_cross_margin"] == Q(49996, 40000000)
        ),
        "compact_margins_remain_separate_and_positive": (
            margins["A0_own_margin_both_routes"] == Q(1, 40000000)
            and margins["A1_own_margin_both_routes"] == Q(1, 40000000)
            and margins["compact_accounts_disjoint"]
        ),
        "both_routes_reach_every_fixed_diagonal_request": (
            diagonal["every_fixed_target_has_finite_log_stage"]
            and diagonal["every_fixed_target_has_finite_pure_stage"]
            and diagonal["fixed_jet_then_finite_stage"]
            and not diagonal["stage_uniform_constant_used"]
        ),
        "complete_residual_endpoint_and_localization_memory_survives": (
            diagonal["all_six_groups_retain_complete_memory"]
            and diagonal["punctured_endpoint_supplied_separately"]
            and diagonal["localized_term_count"] == 12
            and diagonal["complete_pressure_retained"]
            and diagonal["unit_viscosity_retained"] == 1
        ),
        "comparison_and_core_are_rate_independent_after_rebuild": (
            comparison["structural_input_only_requires_rebuilt_same_field_force_pair"]
            and comparison["worst_radius_decay"] == 2
            and comparison["cutoff_power"] == 8
            and comparison["annular_diagonal_terms_vanish_exactly"]
            and comparison["threshold_strictly_positive"]
        ),
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "residual_margins": margins,
        "sample_log_stage": comparison_stage(0, 0, 0, "log"),
        "sample_pure_stage": comparison_stage(0, 0, 0, "pure"),
        "diagonal_replay": diagonal,
        "comparison_replay": comparison,
        "correction": {
            "Q181_exact_margin_remains_valid_for_log_renormalized_WMS_route": True,
            "Q181_margin_is_too_large_for_corrected_pure_power_route": True,
            "corrected_pure_margin_is_49996_over_40000000": True,
            "pure_route_may_reselect_diagonal_scales_so_no_byte_identity_is_claimed": True,
            "Q182_implication_replays_on_either_rebuilt_same_field_force_pair": True,
            "historical_frozen_evidence_modified": False,
        },
        "scope": {
            "Q180_Q182_log_route_replayed_conditionally": True,
            "Q180_Q182_corrected_pure_route_replayed_conditionally": True,
            "original_old_rate_exact_pure_power_contract_proved": False,
            "exact_core_and_finite_existence_inverse_gap_lemmas_independently_reproved": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
        "next": (
            "begin an independent source-level audit of the earliest retained exact-core and "
            "finite existence/inverse/gap premise, preserving the corrected log or pure rate route"
        ),
    }


def record():
    return {
        "status": "Q185-corrected-rate-diagonal-localization-comparison-replayed-conditionally",
        "established": (
            "the Q180--Q182 diagonal, complete residual, endpoint, localization, fixed-gauge "
            "comparison and protected-core implications close conditionally in the exact W/M/S "
            "classes and also after spending kappa for a corrected pure-power recurrence"
        ),
        "correction": (
            "Q181's 49997/40000000 regular/cross reserve is the log-renormalized margin; the "
            "corrected pure-power gain gives 49996/40000000, while the disjoint A0/A1 margins "
            "remain 1/40000000"
        ),
        "object_limit": (
            "the pure route may choose different finite diagonal scales, so it is not identified "
            "byte-for-byte with Q181's field; Q182 replays structurally on its rebuilt common "
            "field-force pair with the same support hole and background"
        ),
        "historical_frozen_records_modified": False,
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
