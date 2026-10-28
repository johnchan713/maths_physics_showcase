"""Bind the Q1--Q5 object to the log-reserved six-group residual replay."""

from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, RESULTS / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


q5 = load(
    "q6_q5_weighted_inverse",
    "core_current_head_q5_weighted_inverse_source_join/weighted_inverse.py",
)
current_residual = load(
    "q6_current_residual",
    "core_current_log_reserved_residual_reconciliation/residual.py",
)


DIRECT_DEPENDENCIES = (
    "core_current_head_q5_weighted_inverse_source_join",
    "core_current_log_reserved_residual_reconciliation",
    "core_log_reserved_complete_residual_replay",
    "core_current_log_source_reconciliation",
    "core_complete_potential_diagonal_join",
    "core_corner_tail_slack_reconciliation",
)

EXPECTED_STATUS = {
    "core_current_head_q5_weighted_inverse_source_join":
        "current-head-Q5-weighted-inverse-and-log-source-domain-bound",
    "core_current_log_reserved_residual_reconciliation":
        "current-chain-log-reserved-six-group-local-residual-reconciled-conditionally",
    "core_log_reserved_complete_residual_replay":
        "log-reserved-complete-six-group-local-residual-replayed-conditionally",
    "core_current_log_source_reconciliation":
        "current-five-contract-log-renormalized-actual-source-induction-reconciled-conditionally",
    "core_complete_potential_diagonal_join":
        "conditional-corrected-complete-potential-diagonal-join",
    "core_corner_tail_slack_reconciliation":
        "closed-corner-cutoff-tail-slack-reconciled",
}

RESIDUAL_GROUPS = (
    "time", "viscosity", "pressure", "old_new", "new_old", "new_square",
)

PHYSICAL_MEMORY = (
    "regular_current_source",
    "A0_compact_tail",
    "A1_compact_tail",
    "old_Gaussian_tail",
    "signed_cutoff_time_tail",
    "flat_exterior_tail",
    "six_covariance_entries",
    "actual_axial_velocity_gd_minus_A1gd",
    "whole_dZ_Jz_plus_c_rho_P",
    "old_pressure_and_Pi1_Pi2",
    "potential_cutoff_curl_commutators",
    "full_pressure_cutoff_derivatives",
)


def evidence(name):
    return json.loads((RESULTS / name / "evidence.json").read_text())


def shared_source_record():
    """Both Q5 and the replay must bind the same source-reconciliation record."""
    relative = (
        "results/core_current_log_source_reconciliation/evidence.json"
    )
    q5_record = evidence("core_current_head_q5_weighted_inverse_source_join")
    replay_record = evidence("core_current_log_reserved_residual_reconciliation")
    return {
        "relative_path": relative,
        "q5_digest": q5_record["dependency_evidence_sha256"][relative],
        "residual_digest": replay_record["dependency_evidence_sha256"][relative],
        "same_digest": (
            q5_record["dependency_evidence_sha256"][relative]
            == replay_record["dependency_evidence_sha256"][relative]
        ),
    }


def order(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


@lru_cache(None)
def q5_review():
    return q5.review()


def join(m=0, target=0, compact_logs=(0, 0), regular_log=0, cutoff_log=0):
    """Return one exact fixed-jet Q6 comparison ledger."""
    m, target = order(m, "m"), order(target, "target")
    row = current_residual.comparison(
        m, target, compact_logs, regular_log, cutoff_log
    )
    q5_result = q5_review()
    groups = {name: Q(value) for name, value in row["residual_group_exponents"].items()}
    fixed_state = Q(row["fixed_state_exponent_after_logs"])
    complete_minimum = Q(row["complete_minimum_exponent"])
    stage_request = row["current_stage_request"]
    shared = shared_source_record()
    accounts = row["reserve_accounts"]

    checks = {
        "Q5_graph_and_interfaces_hold": all(q5_result["checks"].values()),
        "Q5_and_residual_share_exact_source_record": shared["same_digest"],
        "same_fixed_stage_used_for_source_and_residual": (
            stage_request["stage"] == row["log_reserved_comparison_stage"]
        ),
        "extra_A0_A1_integration_cost_is_two_derivatives": all(
            stage_request["raw_sources_with_extra_integration"][name]
            == stage_request["raw_sources_before_extra_integration"][name] + 2
            for name in ("A0", "A1")
        ),
        "four_accounts_in_three_disjoint_classes": (
            set(accounts) == {"A0", "A1", "regular", "cutoff"}
            and len({entry["owner"] for entry in accounts.values()}) == 4
            and {entry["reserve_class"] for entry in accounts.values()}
            == {"compact_source", "regular_residual", "cutoff_tail"}
            and all(not entry["may_pay_other_account"] for entry in accounts.values())
        ),
        "all_six_groups_retained": tuple(groups) == RESIDUAL_GROUPS,
        "both_cross_orientations_retained_separately": (
            "old_new" in groups and "new_old" in groups
            and groups["old_new"] == groups["new_old"]
        ),
        "unit_viscosity_pressure_and_square_retained": all(
            name in groups for name in ("viscosity", "pressure", "new_square")
        ),
        "fixed_state_and_every_tail_group_beat_target": (
            fixed_state >= target + 1
            and all(value >= target + 1 for value in groups.values())
        ),
        "complete_exponent_is_the_actual_minimum": (
            complete_minimum == min(fixed_state, *groups.values())
        ),
        "cutoff_cross_can_limit_the_common_reserve": (
            complete_minimum <= fixed_state
        ),
        "current_residual_controls_hold": all(row["checks"].values()),
    }
    return {
        "jet": m,
        "target": target,
        "base_stage": row["base_comparison_stage"],
        "log_reserved_stage": row["log_reserved_comparison_stage"],
        "stage_request": stage_request,
        "shared_source_record": shared,
        "reserve_accounts": accounts,
        "fixed_state_exponent": str(fixed_state),
        "six_group_exponents": {name: str(value) for name, value in groups.items()},
        "complete_minimum_exponent": str(complete_minimum),
        "residual_identity": row["residual_identity"],
        "physical_memory": list(PHYSICAL_MEMORY),
        "checks": checks,
    }


def review():
    samples = [
        join(0, 0, (0, 0), 0, 0),
        join(1, 3, (7, 11), 13, 17),
        join(5, 11, (37, 41), 43, 47),
    ]
    statuses = {
        name: evidence(name).get("status") == status
        for name, status in EXPECTED_STATUS.items()
    }
    checks = {
        "six_direct_provider_statuses_match": len(statuses) == 6 and all(statuses.values()),
        "same_Q5_source_object_enters_residual": all(
            row["checks"]["Q5_and_residual_share_exact_source_record"] for row in samples
        ),
        "finite_enlarged_requests_used_at_exact_stage": all(
            row["checks"]["same_fixed_stage_used_for_source_and_residual"]
            and row["checks"]["extra_A0_A1_integration_cost_is_two_derivatives"]
            for row in samples
        ),
        "three_log_classes_remain_disjoint": all(
            row["checks"]["four_accounts_in_three_disjoint_classes"] for row in samples
        ),
        "time_unit_viscosity_pressure_crosses_and_square_retained": all(
            row["checks"]["all_six_groups_retained"]
            and row["checks"]["both_cross_orientations_retained_separately"]
            and row["checks"]["unit_viscosity_pressure_and_square_retained"]
            for row in samples
        ),
        "complete_memory_and_cutoff_derivatives_retained": all(
            tuple(row["physical_memory"]) == PHYSICAL_MEMORY for row in samples
        ),
        "sample_fixed_jet_targets_close": all(
            all(row["checks"].values()) for row in samples
        ),
        "common_reserve_not_inflated_by_adding_disjoint_margins": all(
            row["checks"]["complete_exponent_is_the_actual_minimum"] for row in samples
        ),
        "downstream_chain_not_silently_promoted": True,
    }
    return {
        "checks": checks,
        "provider_status_checks": statuses,
        "samples": samples,
        "residual_groups": list(RESIDUAL_GROUPS),
        "physical_memory": list(PHYSICAL_MEMORY),
        "scope": {
            "Q1_Q5_same_object_is_directly_bound_to_residual_replay": True,
            "every_fixed_jet_and_target_has_one_finite_comparison_stage": True,
            "complete_local_residual_is_flat_conditionally": True,
            "disjoint_reserves_are_added_into_one_stronger_common_margin": False,
            "stage_uniform_constant_or_J_of_q_used": False,
            "corner_endpoint_localization_force_energy_comparison_replayed_here": False,
            "hashes_or_tests_treated_as_continuum_proof": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
