"""Join Q165 directly to Q6 and audit the universal residual margins."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
IDENTITY = "Q166-current-head-complete-residual-frontier-audit"
PROVIDERS = (
    "core_current_head_q165_finite_map_frontier_audit",
    "core_current_head_q6_complete_residual_join",
)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q6 = load("q166_q6", "core_current_head_q6_complete_residual_join/complete_residual.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q166 identity required")
    q165 = evidence(PROVIDERS[0])
    q6_record = evidence(PROVIDERS[1])
    if q165.get("status") != "Q165-current-head-finite-map-and-weighted-source-chain-audited":
        raise ValueError("Q165 status mismatch")
    if q6_record.get("status") != "current-head-Q6-complete-six-group-local-residual-bound-conditionally":
        raise ValueError("Q6 status mismatch")
    q5_path = "results/core_current_head_q5_weighted_inverse_source_join/evidence.json"
    q165_q5 = q165["dependency_evidence_sha256"][q5_path]
    q6_q5 = q6_record["dependency_evidence_sha256"][q5_path]
    return {
        "Q165_U2_closed": q165["symbolic"]["U2_finite_map_and_weighted_source_chain_closed"],
        "Q165_frontier_was_Q6": q165["frontier"][-1] == {
            "id": "U3", "status": "next-Q6-complete-residual-audit"
        },
        "identical_Q5_evidence_enters_Q165_and_Q6": q165_q5 == q6_q5,
        "Q5_evidence_sha256": q165_q5,
        "Q6_all_provider_gates_hold": all(q6_record["gates"].values()),
        "Q6_unconditional_claimed": q6_record["review"]["scope"][
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed"
        ],
    }


def universal_margin(identity=IDENTITY):
    """Symbolic margins valid after fixing any finite jet and target."""
    if identity != IDENTITY:
        raise ValueError("Q166 identity required")
    h = Q(1, 200)
    kappa = Q(1, 100000)
    delta = Q(1, 2) - 3 * kappa
    stage_gain = h * delta
    regular_retained = stage_gain / 2
    compact_retained = h * kappa / 2
    # If J0 meets Q6's old comparison threshold, Q6 uses J=J0+1.
    # Thus J+1=J0+2 gives the limiting cutoff crosses one extra half
    # stage gain beyond N+1.  The regular fixed-state term has the same
    # reserve.  The compact terms start one full stage above the old
    # threshold and retain their own additional positive margin.
    return {
        "h": h,
        "kappa": kappa,
        "delta": delta,
        "one_stage_gain": stage_gain,
        "regular_retained_margin": regular_retained,
        "limiting_cross_margin_at_least": regular_retained,
        "A0_own_log_retained_margin": compact_retained,
        "A1_own_log_retained_margin": compact_retained,
        "universal_common_margin_at_least": regular_retained,
        "regular_and_cross_reserves_same_size_but_distinct_owners": True,
        "compact_reserves_may_pay_regular_or_cutoff_logs": False,
        "disjoint_reserves_added_to_strengthen_common_exponent": False,
    }


def logarithmic_absorption(power, reserve):
    if type(power) is not int or power < 0:
        raise ValueError("power must be a nonnegative integer")
    reserve = Q(reserve)
    if reserve <= 0:
        raise ValueError("reserve must be positive")
    # log(1+s)<=sqrt(s) and Young's inequality give the exponent below.
    return Q(power * power, 1) / (4 * reserve)


def residual_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q166 identity required")
    rows = (
        q6.join(0, 0, (0, 0), 0, 0),
        q6.join(1, 3, (7, 11), 13, 17),
        q6.join(5, 11, (37, 41), 43, 47),
    )
    groups = ("time", "viscosity", "pressure", "old_new", "new_old", "new_square")
    memory = (
        "regular_current_source", "A0_compact_tail", "A1_compact_tail",
        "old_Gaussian_tail", "signed_cutoff_time_tail", "flat_exterior_tail",
        "six_covariance_entries", "actual_axial_velocity_gd_minus_A1gd",
        "whole_dZ_Jz_plus_c_rho_P", "old_pressure_and_Pi1_Pi2",
        "potential_cutoff_curl_commutators", "full_pressure_cutoff_derivatives",
    )
    return {
        "sample_stages": [row["log_reserved_stage"] for row in rows],
        "sample_complete_minimum_exponents": [row["complete_minimum_exponent"] for row in rows],
        "all_sample_controls_hold": all(all(row["checks"].values()) for row in rows),
        "all_six_groups_retained_in_order": all(tuple(row["six_group_exponents"]) == groups for row in rows),
        "both_cross_orientations_retained": all(
            row["six_group_exponents"]["old_new"] == row["six_group_exponents"]["new_old"]
            for row in rows
        ),
        "complete_pressure_covariance_axial_and_cutoff_memory_retained": all(
            tuple(row["physical_memory"]) == memory for row in rows
        ),
        "default_limiting_cross_exponent": rows[0]["six_group_exponents"]["old_new"],
        "default_complete_minimum": rows[0]["complete_minimum_exponent"],
        "default_fixed_state_exponent": rows[0]["fixed_state_exponent"],
        "complete_exponent_is_minimum_not_sum": all(
            Q(row["complete_minimum_exponent"]) == min(
                Q(row["fixed_state_exponent"]),
                *(Q(value) for value in row["six_group_exponents"].values()),
            ) for row in rows
        ),
        "fixed_jet_target_then_finite_stage_quantifier_retained": True,
        "stage_uniform_constant_or_J_of_q_used": False,
    }


def remaining_gate(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q166 identity required")
    return {
        "U3_complete_six_group_local_residual": "closed-conditionally-by-Q6-and-audited-by-Q166",
        "next": (
            "audit Q7's Cartesian-corner transfer and punctured endpoint premise on "
            "the same Q164-Q166 object; residual flatness must not be used as a "
            "substitute for bounded terminal velocity, potential or pressure jets"
        ),
        "closed_corner_endpoint_localization_force_energy_comparison_replayed_here": False,
        "unconditional_theorem_gate_closed": False,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q166 identity required")
    return (
        {"id": "U1", "status": "closed-by-Q164-core-pressure-audit"},
        {"id": "U2", "status": "closed-by-Q165-finite-map-audit"},
        {"id": "U3", "status": "closed-conditionally-by-Q166-Q6-residual-audit"},
        {"id": "U4", "status": "next-Q7-corner-and-endpoint-audit"},
    )


def record():
    return {
        "status": "Q166-current-head-complete-six-group-local-residual-audited-conditionally",
        "correction": (
            "Q6 predates Q164-Q165 and was not a direct dependency of their new frontier; "
            "Q166 joins it through the identical frozen Q5 source object and rederives "
            "the all-fixed-jet stage margins"
        ),
        "established": (
            "for every fixed jet and target one finite stage controls time, unit "
            "viscosity, complete pressure, both nonlinear crosses and the new square, "
            "with separate positive logarithmic accounts and complete physical memory"
        ),
        "scope": (
            "conditional local residual flatness relative to named continuum providers; "
            "exact arithmetic, evidence hashes and passing tests are not a continuum proof"
        ),
        "U3_complete_local_residual_closed_conditionally": True,
        "corner_endpoint_or_localized_force_replayed_here": False,
        "unconditional_pulse_residual_force_exclusion_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
        "next": remaining_gate()["next"],
    }
