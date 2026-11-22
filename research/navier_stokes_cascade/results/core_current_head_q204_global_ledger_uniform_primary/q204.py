"""Q204: join global ledgers and uniform primary pulses to Q203."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
IDENTITY = "Q204-GLOBAL-LEDGER-UNIFORM-PRIMARY"
PROVIDERS = (
    "core_current_head_q203_diagonal_endpoint_force_corrected_domain",
    "core_current_head_q148_coefficient_support_pressure",
    "core_current_head_q149_primary_pulse_join",
    "core_current_head_q152_physical_stress_feedback_join",
    "core_global_pressure_moment_recheck",
    "core_coefficient_support_recheck",
)
EXPECTED = {
    PROVIDERS[0]: "Q203-corrected-domain-diagonal-endpoint-localized-force-replayed-conditionally",
    PROVIDERS[1]: "Q147-final-profile-finite-coefficient-supports-and-reserved-pressure-obstruction-joined",
    PROVIDERS[2]: "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined",
    PROVIDERS[3]: "Q152-physical-principal-stress-gate-closed-and-feedback-D-joined",
    PROVIDERS[4]: "same-repaired-profile-global-pressure-moments-and-exterior-stress-identities-independently-rechecked",
    PROVIDERS[5]: "same-profile-finite-coefficient-supports-and-reserved-pressure-obstruction-rechecked",
}


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q203 = load("q204_q203", "core_current_head_q203_diagonal_endpoint_force_corrected_domain/q203.py")
q148 = load("q204_q148", "core_current_head_q148_coefficient_support_pressure/q148.py")
q149 = load("q204_q149", "core_current_head_q149_primary_pulse_join/q149.py")
q152 = load("q204_q152", "core_current_head_q152_physical_stress_feedback_join/q152.py")
global_checks = load("q204_global", "core_global_pressure_moment_recheck/proof_checks.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q204 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def same_object_chain(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q204 identity required")
    links = (
        ("core_current_head_q203_diagonal_endpoint_force_corrected_domain", "core_current_head_q202_next_sourced_cycle"),
        ("core_current_head_q202_next_sourced_cycle", "core_current_head_q180_target_source_simultaneous_induction"),
        ("core_current_head_q180_target_source_simultaneous_induction", "core_current_head_q161_source_tail_simultaneous_induction"),
        ("core_current_head_q161_source_tail_simultaneous_induction", "core_current_head_q153_physical_source_tail_diagonal_join"),
        ("core_current_head_q153_physical_source_tail_diagonal_join", "core_current_head_q152_physical_stress_feedback_join"),
        ("core_current_head_q152_physical_stress_feedback_join", "core_current_head_q149_primary_pulse_join"),
        ("core_current_head_q149_primary_pulse_join", "core_current_head_q148_coefficient_support_pressure"),
        ("core_current_head_q148_coefficient_support_pressure", "core_coefficient_support_recheck"),
        ("core_coefficient_support_recheck", "core_global_pressure_moment_recheck"),
    )
    out = {}
    for parent, child in links:
        child_path = f"results/{child}/evidence.json"
        out[f"{parent}_uses_{child}"] = (
            evidence(parent)["dependency_evidence_sha256"][child_path]
            == sha(PROJECT / child_path)
        )
    return out


def global_pressure_moments(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q204 identity required")
    r = global_checks.review()
    return {
        "renormalized_heat_moment": r["heat"]["R"],
        "raw_I_infinity_exists": r["heat"]["raw_I_infinity_exists"],
        "endpoint_to_global_formula": r["endpoint_to_global"]["difference_formula"],
        "common_future_tail_required": r["endpoint_to_global"]["common_future_tail"],
        "global_zero_defects": r["edits"]["global_zero_defects"],
        "functional_in_eta": r["edits"]["functional_in_eta"],
        "intermediate_pressure_memory_may_be_nonzero": r["edits"]["intermediate_pressure_memory_may_be_nonzero"],
        "same_pressure_gauge_gives_backward_tail": r["checks"]["same_pressure_gauge_gives_backward_tail"],
        "zero_five_defects_remove_all_exterior_memory": r["checks"]["zero_five_defects_remove_all_exterior_memory"],
        "pressure_eta_slope_defect_grows_linearly_in_X": r["checks"]["pressure_eta_slope_defect_grows_linearly_in_X"],
        "zero_U_removes_M_induced_radial_velocity": False,
        "all_exact_checks": all(r["checks"].values()),
    }


def finalized_profile():
    closure = {
        "stage": "same-profile-compact-repair-cone-and-global-pressure-join",
        "global_identities": q148.GLOBAL,
        "pressure": "Pi(X)=-integral_X^infty E^2/(2x) dx in the unchanged gauge",
        "same_exact_repair_and_inner_outer_value_cone_certified": True,
        "global_pressure_and_five_moment_identities_certified": True,
    }
    return q148.finalized_profile(closure, q148.IDENTITY)


def coefficient_supports(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q204 identity required")
    profile = finalized_profile()
    args = ((2, 3, 5, 7, 11), 8, 2, Q(3, 5), 7, Q(1, 1000))
    rows = {n: q148.finite_order(profile, n, *args, q148.IDENTITY) for n in (1, 2, 3)}
    pressure = q148.pressure_obstruction(10, 30, 3, Q(3, 7), Q(1, 5), q148.IDENTITY)
    return {
        "inverse_bound": rows[1]["inverse_bound"],
        "rows": rows[1]["rows"],
        "order_1_support": rows[1]["support"],
        "order_2_support": rows[2]["support"],
        "order_3_support": rows[3]["support"],
        "velocity_on_I4": rows[1]["velocity_on_I4"],
        "pressure_on_I4": rows[1]["pressure_on_I4"],
        "finite_order_only": rows[1]["finite_order_only"],
        "reserved_pressure_floor_sample": pressure["necessary_supremum_floor"],
        "pressure_floor_independent_of_K1_K2": pressure["independent_of"] == ("K1", "K2"),
        "pressure_floor_sharp_after_constant_cancellation": pressure["sharp_after_arbitrary_X_independent_cancellation"],
        "general_pressure_floor": "|chi(c1q)|q^(2h)mu^2(e^-6-e^-11)/(8rho3)",
    }


def schedule_and_primary(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q204 identity required")
    q148_symbolic = evidence(PROVIDERS[1])["symbolic"]
    joined = q149.close(q148_symbolic, q149.IDENTITY)
    schedule = q149.schedule_row(3, 2, Q(3, 200), 10, 1000, q149.IDENTITY)
    powers = q149.primary_powers(Q(1, 200), 0, 0, q149.IDENTITY)
    return {
        "Q148_schedule_obstruction_discharge": joined["schedule_correction"],
        "schedule_binary_margin_sample": schedule["binary_margin"],
        "schedule_strict_padded_rule_holds": schedule["strict_padded_rule_holds"],
        "scale_chosen_after_finite_Cnm": schedule["scale_is_chosen_after_finite_C_nm"],
        "preexisting_numeric_scale_certificate_required": schedule["preexisting_numeric_scale_certificate_required"],
        "uniform_primary_pulse_family_certified": joined["uniform_primary_pulse_family_certified"],
        "scale_domain": joined["scale_domain"],
        "velocity_bound": joined["velocity"],
        "pressure_bound": joined["pressure"],
        "potential_bound": joined["potential"],
        "supports": joined["supports"],
        "full_pressure_memory": joined["full_pressure_memory"],
        "canonical_and_physical_families_identified": joined["canonical_and_physical_primary_families_identified"],
        "complete_generated_hierarchy_certified": joined["complete_nonlinear_generated_pulse_hierarchy_certified"],
        "primary_zero_jet_powers": powers,
    }


def covariance_and_current_cycle(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q204 identity required")
    contracts = q152.provider_contracts()
    gate = q152.stress_gate()
    current = q203.q202.current_source()
    inverse = q203.q202.sourced_inverse()
    cycle = q203.q202.next_cycle()
    return {
        "same_physical_family": contracts["same_physical_family"],
        "principal_stress_gate_closed_exactly": gate["principal_gate_closed_exactly"],
        "complete_curl_margin": gate["complete_curl_margin"],
        "initial_and_dynamic_memory": (
            "all six covariance entries, radial and axial velocity, curl remainders, Pi1/Pi2, "
            "compact pressure and the whole derivative of c_rho P"
        ),
        "postfast_scalars": contracts["postfast_scalars"],
        "fresh_fast": contracts["fresh_fast"],
        "fresh_slow": contracts["fresh_slow"],
        "mean_feedback": contracts["mean_feedback"],
        "current_source_order": current["order"],
        "current_source_has_complete_operator_memory": current["full_pressure_gradient_viscosity_cylindrical_crosses_and_squares_retained"],
        "current_pressure_and_potential_order": inverse["pressure_and_vector_potential"],
        "moving_and_source_normal_pressure_retained": inverse["moving_normal_term_retained"] and inverse["source_normal_term_retained"],
        "next_source": cycle["next_source"],
        "next_target": cycle["next_required_target"],
        "target_slack": cycle["target_slack"],
        "cycle_closes_relative_to_named_inputs": cycle["closes_relative_to_named_continuum_inputs"],
    }


def diagonal_residual(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q204 identity required")
    diagonal = q203.complete_diagonal()
    local = q203.localization_force()
    return {
        "all_six_groups": diagonal["six_groups"],
        "complete_memory": diagonal["all_six_groups_retain_complete_memory"],
        "cutoff_curl_and_pressure_derivatives": diagonal["every_scale_cutoff_curl_and_pressure_derivative_retained"],
        "fixed_jet_then_finite_stage": diagonal["fixed_jet_then_finite_stage"],
        "stage_uniform_constant_used": diagonal["stage_uniform_constant_or_J_of_q_used"],
        "conditional_fixed_jet_flatness": diagonal["all_fixed_cartesian_residual_jets_flat_conditionally"],
        "localization_term_count": local["term_count"],
        "complete_localized_pressure": local["complete_pressure"],
        "unit_viscosity": local["unit_viscosity"],
    }


def precise_limits(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q204 identity required")
    primary = schedule_and_primary()["primary_zero_jet_powers"]
    sourced = q203.q202.physical_powers()
    return {
        "reserved_pressure_cannot_be_replaced_by_leading_pressure": True,
        "reserved_pressure_floor": "sup_I4|q^(2A)p_B-Pi0| >= |chi(c1q)|q^(2h)mu^2(e^-6-e^-11)/(8rho3)",
        "primary_velocity_zero_jet_power": primary["velocity"],
        "primary_pressure_zero_jet_power": primary["pressure"],
        "sourced_velocity_zero_jet_power": sourced["sourced_velocity_zero_jet"],
        "sourced_pressure_zero_jet_power": sourced["sourced_pressure_zero_jet"],
        "absolute_decay_certified": False,
        "uniform_primary_scale_bounds_certified": True,
        "uniform_generated_hierarchy_constants_certified": False,
        "Q202_constants_uniform_in_cycle": False,
        "tests_or_evidence_hashes_are_continuum_proof": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q204 identity required")
    status = statuses()
    chain = same_object_chain()
    global_rows = global_pressure_moments()
    coefficients = coefficient_supports()
    primary = schedule_and_primary()
    covariance = covariance_and_current_cycle()
    residual = diagonal_residual()
    limits = precise_limits()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "same_object_chain_reaches_global_ledgers": all(chain.values()),
        "five_global_pressure_moment_rows_close_in_one_gauge": (
            global_rows["global_zero_defects"] == ("M", "R", "J", "S", "Pi_axis+Cp")
            and global_rows["same_pressure_gauge_gives_backward_tail"]
            and global_rows["zero_five_defects_remove_all_exterior_memory"]
        ),
        "coefficient_supports_and_pressure_exceptions_are_exact": (
            coefficients["inverse_bound"] == 800
            and coefficients["order_1_support"] == {"velocity": "Xvel", "pressure": "Xv", "stress": "Xb"}
            and coefficients["order_2_support"] == {"velocity": "Xvel", "pressure": "Xv", "stress": "Xv"}
            and set(coefficients["order_3_support"].values()) == {"Xvel"}
            and coefficients["pressure_floor_independent_of_K1_K2"]
        ),
        "rule_defined_schedule_and_uniform_primary_close": (
            primary["schedule_strict_padded_rule_holds"]
            and primary["scale_chosen_after_finite_Cnm"]
            and not primary["preexisting_numeric_scale_certificate_required"]
            and primary["uniform_primary_pulse_family_certified"]
        ),
        "complete_covariance_pressure_and_axial_memory_reaches_current_cycle": (
            covariance["same_physical_family"]
            and covariance["principal_stress_gate_closed_exactly"]
            and covariance["current_source_has_complete_operator_memory"]
            and covariance["moving_and_source_normal_pressure_retained"]
            and covariance["cycle_closes_relative_to_named_inputs"]
        ),
        "complete_residual_ledger_reaches_diagonal_and_localization": (
            residual["complete_memory"]
            and residual["cutoff_curl_and_pressure_derivatives"]
            and residual["fixed_jet_then_finite_stage"]
            and not residual["stage_uniform_constant_used"]
            and residual["localization_term_count"] == 12
            and residual["complete_localized_pressure"]
            and residual["unit_viscosity"] == 1
        ),
        "remaining_limits_are_not_promoted": (
            limits["reserved_pressure_cannot_be_replaced_by_leading_pressure"]
            and limits["primary_velocity_zero_jet_power"] < 0
            and limits["sourced_velocity_zero_jet_power"] < 0
            and not limits["absolute_decay_certified"]
            and not limits["uniform_generated_hierarchy_constants_certified"]
        ),
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "same_object_chain": chain,
        "global_pressure_moments": global_rows,
        "coefficient_reserved_supports": coefficients,
        "schedule_uniform_primary": primary,
        "complete_covariance_current_cycle": covariance,
        "diagonal_localized_residual": residual,
        "precise_limits": limits,
        "scope": {
            "global_ledgers_and_uniform_primary_joined_to_Q203_conditionally": True,
            "uniform_complete_generated_hierarchy_proved": False,
            "absolute_pulse_smallness_proved": False,
            "reserved_pressure_deleted": False,
            "terminal_origin_regularized": False,
            "future_velocity_or_pressure_constructed": False,
            "unconditional_complete_residual_force_exclusion_or_blowup": False,
        },
        "next": (
            "audit the sourced and signed pulse value, derivative and closed-support constants "
            "through every finite cycle on the selected physical family, deciding whether the "
            "stage-dependent bounds suffice directly or require a new hierarchy majorant before "
            "any unconditional complete-residual claim"
        ),
    }


def record():
    return {
        "status": "Q204-global-pressure-coefficient-uniform-primary-and-covariance-ledgers-joined-conditionally",
        "established": (
            "the same Q203 object carries the five unchanged-gauge global moment/pressure identities, "
            "finite coefficient and reserved-mean supports, the rule-defined Borel schedule, uniform "
            "primary bounds on every used band, complete covariance/axial/pressure memory and the "
            "six-group diagonal plus twelve-term localization ledger"
        ),
        "correction": (
            "Q148's displayed Borel inequality is discharged by the implemented scale rule after each "
            "finite seminorm; unevaluated giant constants prevent a numerical scale certificate but do "
            "not obstruct the logical background or uniform primary construction"
        ),
        "precise_limit": (
            "reserved pressure has a sharp positive two-point floor, the primary and first sourced "
            "zero-jet velocity/pressure powers do not certify absolute decay, and constants for the "
            "complete generated hierarchy remain stage-dependent"
        ),
        "historical_frozen_records_modified": False,
        "unconditional_complete_residual_force_exclusion_or_blowup_certified": False,
    }
