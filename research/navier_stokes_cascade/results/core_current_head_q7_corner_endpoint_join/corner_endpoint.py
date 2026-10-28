"""Bind the Q6 residual to the selected corner and punctured endpoint ledger."""

from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
H = Q(1, 200)
KAPPA = Q(1, 100000)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, RESULTS / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


q6 = load(
    "q7_q6_complete_residual",
    "core_current_head_q6_complete_residual_join/complete_residual.py",
)
legacy = load(
    "q7_current_corner_endpoint",
    "core_current_log_reserved_corner_endpoint_join/transfer.py",
)
endpoint = load(
    "q7_current_endpoint_jets",
    "current_endpoint_jets/endpoint.py",
)


DIRECT_DEPENDENCIES = (
    "core_current_head_q6_complete_residual_join",
    "core_current_log_reserved_corner_endpoint_join",
    "current_endpoint_jets",
    "core_log_reserved_cartesian_corner_replay",
    "core_log_reserved_endpoint_compatibility_join",
    "core_corner_flatness_obstruction",
)

EXPECTED_STATUS = {
    "core_current_head_q6_complete_residual_join":
        "current-head-Q6-complete-six-group-local-residual-bound-conditionally",
    "core_current_log_reserved_corner_endpoint_join":
        "current-log-reserved-corner-and-punctured-endpoint-joined-conditionally",
    "current_endpoint_jets":
        "conditional-current-field-closed-time-atlas-and-punctured-terminal-jets",
    "core_log_reserved_cartesian_corner_replay":
        "log-reserved-six-group-residual-transferred-to-cartesian-corner-conditionally",
    "core_log_reserved_endpoint_compatibility_join":
        "log-reserved-corner-and-punctured-endpoint-inputs-joined-conditionally",
    "core_corner_flatness_obstruction":
        "closed-corner-transfer-proved-and-uniform-constant-obstruction-isolated",
}


def evidence(name):
    return json.loads((RESULTS / name / "evidence.json").read_text())


def order(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def selected_chart():
    chart = legacy.closed_chart(H)
    return {
        **chart,
        "h": str(H),
        "kappa": str(KAPPA),
        "same_Q6_parameter_choice": True,
        "alternate_h_diagnostic_promoted": False,
    }


@lru_cache(None)
def corner_target(jet=0, power=0, compact_logs=(0, 0),
                  regular_log=0, cutoff_log=0):
    jet, power = order(jet, "jet"), order(power, "power")
    row = q6.join(jet, power, compact_logs, regular_log, cutoff_log)
    exponent = Q(row["complete_minimum_exponent"])
    requested = Q(power + 1)
    chart = selected_chart()
    request = row["stage_request"]
    return {
        "jet": jet,
        "requested_cartesian_power": power,
        "comparison_stage": row["log_reserved_stage"],
        "stage_request": request,
        "closed_chart_exponent": str(exponent),
        "requested_strict_exponent": str(requested),
        "retained_exponent_margin": str(exponent - requested),
        "cartesian_factor": legacy.cartesian_constant(power + 1),
        "closed_chart": chart,
        "six_group_exponents": row["six_group_exponents"],
        "reserve_accounts": row["reserve_accounts"],
        "physical_memory": row["physical_memory"],
        "checks": {
            "Q6_controls_hold": all(row["checks"].values()),
            "source_and_corner_use_same_fixed_stage": request["stage"] == row["log_reserved_stage"],
            "extra_A0_A1_derivatives_reach_corner": all(
                request["raw_sources_with_extra_integration"][name]
                == request["raw_sources_before_extra_integration"][name] + 2
                for name in ("A0", "A1")
            ),
            "selected_terminal_chart_closes": chart["terminal_identity"] == 1,
            "closed_jacobian_floor_holds": chart["jacobian_floor"] >= Q(49, 50),
            "normalized_support_margin_positive": chart["normalized_support_margin"] > 0,
            "residual_beats_requested_cartesian_power": exponent >= requested,
            "all_six_Q6_groups_retained": tuple(row["six_group_exponents"]) == q6.RESIDUAL_GROUPS,
            "q_to_cartesian_factor_is_exact": legacy.cartesian_constant(power + 1) == 4 ** (power + 1),
            "no_scale_incompatible_floor_used": not chart["physical_support_floor_used"] and not chart["edge_weight_floor_used"],
            "only_selected_Q6_parameters_used": chart["same_Q6_parameter_choice"] and not chart["alternate_h_diagnostic_promoted"],
        },
    }


@lru_cache(None)
def endpoint_request(spatial_order=0, time_order=0, flat_power=0,
                     compact_logs=(0, 0), regular_log=0, cutoff_log=0):
    a = order(spatial_order, "spatial_order")
    j = order(time_order, "time_order")
    power = order(flat_power, "flat_power")
    total = a + j
    residual_jet = total + 1
    corner = corner_target(
        residual_jet, power, compact_logs, regular_log, cutoff_log
    )
    punctured = endpoint.requirements(total)
    requested = {
        "potential": total + 4,
        "direct_swirl": total + 3,
        "pressure": total + 2,
    }
    supplied = {
        "potential": punctured["future_force_Cauchy_potential_bound"],
        "direct_swirl": punctured["future_force_Cauchy_direct_swirl_bound"],
        "pressure": punctured["future_force_Cauchy_pressure_bound"],
    }
    endpoint_review = endpoint.endpoint_review()
    return {
        "spatial_order": a,
        "time_order": j,
        "total_mixed_order": total,
        "residual_Cauchy_jet": residual_jet,
        "comparison_stage": corner["comparison_stage"],
        "flat_power": power,
        "punctured_field_Cauchy_orders": supplied,
        "localization_Cauchy_orders": requested,
        "endpoint_blocking_inequality": endpoint_review["blocking_inequality"],
        "corner": corner,
        "checks": {
            "punctured_orders_match_localization_orders": supplied == requested,
            "potential_has_four_order_buffer": supplied["potential"] == total + 4,
            "direct_swirl_has_three_order_buffer": supplied["direct_swirl"] == total + 3,
            "pressure_has_two_order_buffer": supplied["pressure"] == total + 2,
            "residual_uses_one_additional_jet": residual_jet == total + 1,
            "Q6_corner_closes_requested_power": Q(corner["retained_exponent_margin"]) >= 0,
            "current_endpoint_controls_hold": all(endpoint_review["checks"].values()),
            "extra_compact_derivatives_reach_endpoint_request": corner["checks"]["extra_A0_A1_derivatives_reach_corner"],
        },
    }


@lru_cache(None)
def review():
    corners = [
        corner_target(0, 0),
        corner_target(1, 3, (7, 11), 13, 17),
        corner_target(5, 11, (37, 41), 43, 47),
    ]
    requests = [endpoint_request(a, j, power) for a, j, power in (
        (0, 0, 0), (1, 3, 5), (5, 4, 11), (3, 7, 2),
    )]
    statuses = {
        name: evidence(name).get("status") == status
        for name, status in EXPECTED_STATUS.items()
    }
    endpoint_review = endpoint.endpoint_review()
    checks = {
        "six_direct_provider_statuses_match": len(statuses) == 6 and all(statuses.values()),
        "Q6_is_direct_residual_parent": all(row["checks"]["Q6_controls_hold"] for row in corners),
        "selected_Q6_parameters_reach_closed_corner": all(row["checks"]["only_selected_Q6_parameters_used"] for row in corners),
        "extra_A0_A1_requests_are_not_dropped": all(row["checks"]["extra_A0_A1_derivatives_reach_corner"] for row in corners),
        "every_sample_residual_jet_is_cartesian_flat": all(all(row["checks"].values()) for row in corners),
        "punctured_and_localization_ledgers_match": all(row["checks"]["punctured_orders_match_localization_orders"] for row in requests),
        "every_sample_endpoint_request_closes": all(all(row["checks"].values()) for row in requests),
        "complete_Q6_memory_reaches_corner": tuple(corners[0]["physical_memory"]) == q6.PHYSICAL_MEMORY,
        "endpoint_uniform_physical_bound_remains_explicit_premise": endpoint_review["checks"]["inherited_all_order_physical_bounds_imply_needed_endpoint_inequality_conditionally"],
        "localization_force_energy_comparison_not_silently_replayed": True,
    }
    return {
        "checks": checks,
        "provider_status_checks": statuses,
        "sample_corner_targets": corners,
        "sample_endpoint_requests": requests,
        "blocking_inequality": endpoint_review["blocking_inequality"],
        "scope": {
            "exact_Q1_Q6_residual_is_cartesian_corner_flat_conditionally": True,
            "same_field_punctured_endpoint_ledger_is_joined_conditionally": True,
            "alternate_h_parameter_family_promoted": False,
            "endpoint_uniform_all_order_physical_bound_proved_here": False,
            "terminal_velocity_or_pressure_regular_at_origin_proved": False,
            "localization_force_energy_comparison_replayed_here": False,
            "hashes_or_tests_treated_as_continuum_proof": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
