"""Reconcile Q174 with the pre-existing same-family canonical source route."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q175-CANONICAL-SOURCE-ROUTE-RECONCILIATION"
PROVIDERS = (
    "core_current_head_q174_q5_operator_source_audit",
    "source_envelope_edge_audit",
    "current_target_cycle_bridge",
    "canonical_finite_cycle_join",
    "current_residual_ledger_audit",
    "explicit_tail_complete_join",
    "foundation_bundle_route",
)
EXPECTED = {
    "core_current_head_q174_q5_operator_source_audit": "Q174-Q5-operator-source-audited-with-family-and-typed-domain-boundary",
    "source_envelope_edge_audit": "regular-source-envelope-edge-estimates-with-explicit-freezing-correction",
    "current_target_cycle_bridge": "conditional-current-target-cancellation-cycle-bridge",
    "canonical_finite_cycle_join": "conditional-canonical-seed-and-all-finite-cycle-operator-bridge",
    "current_residual_ledger_audit": "conditional-current-nonzero-residual-ledger-and-fixed-jet-tail-bound",
    "explicit_tail_complete_join": "conditional-explicit-tail-complete-canonical-residual-and-fresh-force-join",
    "foundation_bundle_route": "foundation-route-closed-with-independent-continuum-validation-frontier",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q175 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def polynomial_to_log(S_power, h, ell):
    """Convert S^p=ell^(2p) to a fixed logarithmic loss exactly."""
    if type(S_power) is not int or S_power < 0:
        raise ValueError("nonnegative integer S power required")
    h = Q(h)
    if h <= 0 or type(ell) is not int or ell < 1:
        raise ValueError("positive h and integer ell required")
    degree = 2 * S_power
    polynomial = Q(ell ** degree, 1)
    log_weight = (1 + h * ell) ** degree
    ratio = polynomial / log_weight
    return {
        "S_power": S_power,
        "log_degree": degree,
        "ratio": ratio,
        "uniform_bound": h ** (-degree),
        "bound_holds": ratio <= h ** (-degree),
    }


def typed_route(harmonic_grammar, current_targets, canonical_rebuild,
                retained_analytic_inputs, independent_validation):
    for value in (harmonic_grammar, current_targets, canonical_rebuild,
                  retained_analytic_inputs, independent_validation):
        if type(value) is not bool:
            raise ValueError("boolean route flags required")
    conditional = harmonic_grammar and current_targets and canonical_rebuild and retained_analytic_inputs
    return {
        "conditional_typed_source_route_closed": conditional,
        "independently_validated_typed_source_route": conditional and independent_validation,
        "physical_to_canonical_transfer_used": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q175 identity required")
    status = statuses()
    q174 = evidence(PROVIDERS[0])
    envelope = evidence(PROVIDERS[1])
    target = evidence(PROVIDERS[2])
    cycle = evidence(PROVIDERS[3])
    residual = evidence(PROVIDERS[4])
    complete = evidence(PROVIDERS[5])
    foundation = evidence(PROVIDERS[6])

    ef = envelope["flags"]
    tf = target["flags"]
    cf = cycle["flags"]
    rf = residual["flags"]
    xf = complete["flags"]
    ff = foundation["flags"]
    route = typed_route(
        ef["generated_regular_source_grammar_has_full_path_and_entrance_inputs_relative_to_current_target_estimates"],
        tf["all_half_line_target_recurrence_relative_to_reviewed_current_map_inputs"],
        cf["fresh_canonical_initial_state_rebuilt_through_signed_sourced_and_mean_prefix"] and cf["universal_operator_estimates_reapplied_to_current_canonical_inputs"],
        ff["named_foundation_obligations_have_analytic_provider_routes"],
        ff["provider_continuum_theorems_independently_reproved_here"],
    )
    log_control = polynomial_to_log(7, Q(1, 200), 4096)
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "Q174_operator_and_typed_weight_correction_retained": q174["review"]["checks"]["compact_slow_maps_close_when_typed_inputs_are_supplied"] and q174["review"]["correction"]["correct_domain_is_typed_and_log_renormalized"],
        "Q174_physical_family_transfer_rejection_retained": q174["review"]["checks"]["Q161_physical_bounds_cannot_be_substituted_into_canonical_chain"],
        "source_grammar_supplies_harmonic_half_weight_inputs_conditionally": ef["all_fixed_mixed_normalized_source_jets_have_Bell_majorants"] and ef["signed_reciprocal_root_half_weight_cancellation_quantified"],
        "current_target_bridge_supplies_full_weight_targets_and_moments_conditionally": tf["current_target_estimate_derived_by_exact_state_difference_not_raw_field_bound"] and tf["generic_signed_target_cancellation_keeps_all_six_covariances_and_changed_pressure"],
        "canonical_cycle_is_rebuilt_on_its_own_family": cf["one_fixed_canonical_H_and_positive_primary_roots_for_all_signed_maps"] and not cf["historical_pulses_sources_roots_mean_fields_or_realization_identified_as_canonical"],
        "canonical_conditional_typed_route_closes_without_family_transfer": route["conditional_typed_source_route_closed"] and not route["physical_to_canonical_transfer_used"],
        "polynomial_S_cost_is_a_fixed_logarithmic_degree": log_control["log_degree"] == 14 and log_control["bound_holds"],
        "complete_pressure_memory_is_retained": cf["same_full_base_Pi1_Pi2_axial_velocity_and_whole_moving_pressure_moment"] and xf["all_six_covariances_axial_velocity_complete_pressure_and_Pi1_Pi2_retained"],
        "complete_finite_residual_is_already_conditionally_joined": rf["conditional_complete_finite_residual_estimate_relative_to_listed_inputs"] and xf["complete_finite_residual_keeps_viscosity_cutoffs_both_crosses_squares_and_tails"],
        "foundation_route_is_acyclic_and_same_object": ff["provider_route_is_acyclic"] and ff["same_selected_maps_roots_pressure_gauge_and_width_retained"],
        "independent_continuum_validation_remains_open": not route["independently_validated_typed_source_route"] and not ff["formal_or_external_independent_certificate"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "typed_route": route,
        "log_control": {k: str(v) if isinstance(v, Q) else v for k, v in log_control.items()},
        "correction": {
            "Q174_same_family_source_frontier_was_stale": True,
            "Q174_family_nonidentification_obstruction_remains_valid_for_transfer": True,
            "canonical_route_closes_by_fresh_same_family_rebuild_not_transfer": True,
            "Q174_instruction_to_replay_Q6_was_already_satisfied_conditionally": True,
            "historical_Q174_evidence_modified": False,
        },
        "conditional_conclusion": {
            "harmonic_bound": "for each fixed stage and jet, the canonical f/P source has epsilon power, zeta^(1/2), finite edge poles and a finite logarithmic degree",
            "mean_target_bound": "for each fixed stage and jet, the recomputed canonical Xi target has epsilon power, full zeta, exact tangential moments and complete pressure memory",
            "residual": "the same canonical family has a conditional complete finite residual bound with viscosity, cutoffs, both crosses, squares and all current interactions",
        },
        "scope": {
            "conditional_canonical_typed_source_membership_closed_relative_to_named_analytic_providers": True,
            "provider_theorems_independently_reproved_here": False,
            "uniform_in_stage_contraction_or_uncontrolled_diagonal_proved": False,
            "historical_first_band_or_physical_family_identified": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_proved": False,
        },
        "next": "independently validate the compact modulation/exact-repair theorem with its continuation and matching-root inputs on the selected profile, then recheck its annular/exterior and regional-cone joins; do not repeat the already routed canonical source or Q6 residual checkpoint",
    }


def record():
    return {
        "status": "Q175-canonical-typed-source-route-reconciled-and-independent-frontier-restored",
        "established": "the pre-existing same-family canonical source grammar, current-target recurrence and finite-cycle bridge conditionally close Q174's two typed source gates, and the complete residual was already conditionally joined",
        "correction": "Q174 correctly rejected transfer from the inherited physical family but incorrectly left the canonical same-family typed source route open; that route is rebuilt independently on canonical Hc, ac and current means",
        "obstruction": "the remaining frontier is independent continuum validation of the routed foundation theorems, beginning with compact modulation/exact repair and its continuation, matching-root and cone-join inputs, not another family-identification or residual replay",
        "historical_frozen_records_modified": False,
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
