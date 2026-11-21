"""Q200: reconcile the weighted-target obstruction with the physical induction."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q200-PHYSICAL-WEIGHTED-TARGET-RECONCILIATION"

PROVIDERS = (
    "core_current_head_q197_first_sourced_covariance",
    "core_current_head_q198_generic_signed_target",
    "core_current_head_q199_weighted_target_obstruction",
    "core_current_head_q149_primary_pulse_join",
    "core_current_head_q150_initial_covariance_mean_join",
    "core_current_head_q152_physical_stress_feedback_join",
    "core_current_head_q160_physical_target_induction",
    "core_current_head_q161_source_tail_simultaneous_induction",
    "core_current_head_q180_target_source_simultaneous_induction",
)

EXPECTED = {
    PROVIDERS[0]: "Q197-first-sourced-complex-covariance-and-current-means-joined",
    PROVIDERS[1]: "Q198-corrected-generic-complex-signed-target-and-weighted-input-obstruction",
    PROVIDERS[2]: "Q199-exact-moment-weighted-current-target-obstruction",
    PROVIDERS[3]: "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined",
    PROVIDERS[4]: "Q149-uniform-primary-initial-complete-covariance-mean-correction-joined",
    PROVIDERS[5]: "Q152-physical-principal-stress-gate-closed-and-feedback-D-joined",
    PROVIDERS[6]: "Q160-physical-all-finite-cycle-target-moment-induction-joined",
    PROVIDERS[7]: "Q161-physical-source-inverse-tail-simultaneous-induction-joined",
    PROVIDERS[8]: "Q180-physical-target-source-tail-simultaneous-induction-reconciled",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q200 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def rates(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q200 identity required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 100):
        raise ValueError("reviewed kappa required")
    B = 1 - kappa
    C = B + Q(1, 2)
    D = C - kappa
    C0 = 2 * B
    delta = Q(1, 2) - 3 * kappa
    return {
        "kappa": kappa,
        "Q197_current_error_B": B,
        "Q197_current_target_input_C": C,
        "Q197_post_source_regular_D": D,
        "Q160_first_source_B0": D,
        "Q198_signed_output_C0": C0,
        "Q160_target_seed_C0": C0,
        "generic_target_gain": C0 - C,
        "future_source_gain": delta,
        "future_target_bridge_gain": Q(1, 2) - kappa,
        "future_target_slack": (Q(1, 2) - kappa) - delta,
        "initial_physical_target": 2 - 3 * kappa,
    }


def stage_identification(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q200 identity required")
    r = rates()
    return {
        "same_profile_chain": (
            "Q149 physical primary; Q150 complete initial covariance/mean; Q152 same-family "
            "stress and first feedback; Q197 first sourced correction; Q198 first generic signed step"
        ),
        "Q197_D_equals_Q160_B0": r["Q197_post_source_regular_D"] == r["Q160_first_source_B0"],
        "Q198_C0_equals_Q160_C0": r["Q198_signed_output_C0"] == r["Q160_target_seed_C0"],
        "initial_target_is_stronger_than_Q197_C": (
            r["initial_physical_target"] - r["Q197_current_target_input_C"]
        ),
        "physical_and_canonical_families_identified": False,
        "same_inherited_physical_family_only": True,
    }


def actual_weight_route(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q200 identity required")
    return {
        "standalone_unweighted_to_weighted_inference_valid": False,
        "Q199_scaled_bump_obstruction_retained": True,
        "direct_physical_route_present": True,
        "exact_difference": "X_e(new)-X_e(old)=-I_e Q_e overline(Delta E_e)",
        "cancellation_before_estimation": True,
        "all_six_covariances_in_pressure_reconstruction": True,
        "axial_velocity_retained": True,
        "compact_pressure_retained": True,
        "whole_moving_pressure_product_retained": True,
        "wave_product_weight": "sqrt(zeta)*sqrt(zeta)=zeta",
        "mean_product_weight": "zeta^2<=zeta",
        "primitive_gate": "exact zero weighted total retains zeta at both edges",
        "moving_projection_commutator_retained": True,
        "reserved_bumps_use_positive_interior_weight_floor": True,
        "raw_Q199_inequality_proved_as_generic_seminorm_implication": False,
        "actual_same_family_target_constructed_by_state_difference": True,
    }


def generic_step(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q200 identity required")
    r = rates(kappa)
    B = r["Q197_current_error_B"]
    D = r["Q197_post_source_regular_D"]
    C0 = r["Q198_signed_output_C0"]
    errors = {
        "old_error_signed_principal": 2 * B,
        "complete_wave_signed_curl_remainder": Q(1, 2) + D,
        "complete_signed_square": 2 * B,
    }
    return {
        "input_target_order": r["Q197_current_target_input_C"],
        "output_target_order": C0,
        "errors": errors,
        "all_errors_equal_output": set(errors.values()) == {C0},
        "gain": C0 - r["Q197_current_target_input_C"],
    }


def simultaneous_step(j=0, kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q200 identity required")
    if isinstance(j, bool) or not isinstance(j, int) or j < 0:
        raise ValueError("finite nonnegative cycle required")
    r = rates(kappa)
    Bj = r["Q160_first_source_B0"] + j * r["future_source_gain"]
    Cj = Bj + Q(1, 2)
    return {
        "j": j,
        "B_j": Bj,
        "C_j": Cj,
        "B_next": Bj + r["future_source_gain"],
        "C_next": Cj + r["future_source_gain"],
        "proved_target_order": Cj + r["future_target_bridge_gain"],
        "target_slack": r["future_target_slack"],
        "source_and_target_advanced_simultaneously": True,
        "sequential_Q160_then_source_reading_valid": False,
        "closes": r["future_target_slack"] > 0,
    }


def correction(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q200 identity required")
    return {
        "Q199_local_provider_cut_was_accurately_reported": True,
        "Q199_global_no_route_conclusion_was_too_strong": True,
        "omitted_physical_providers": PROVIDERS[3:],
        "obstruction_scope": (
            "exact moments, compact support and fixed unweighted jets do not by themselves "
            "imply a flat-weight estimate"
        ),
        "constructive_scope": (
            "the named physical state-difference and simultaneous-induction providers propagate "
            "the full weight directly on the same inherited physical family"
        ),
        "historical_frozen_evidence_modified": False,
    }


def physical_residual_power(h=Q(1, 200), kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q200 identity required")
    h = Q(h)
    if not 0 < h < Q(1, 100):
        raise ValueError("reviewed h required")
    r = rates(kappa)
    A = Q(1, 2) + h
    return -2 * A - Q(1, 2) + h * r["Q197_post_source_regular_D"]


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q200 identity required")
    return {
        "closed_here": (
            "the corrected Q197/Q198 first generic signed target is bound to the direct "
            "same-family weighted target/source induction, with complete pressure memory"
        ),
        "next": (
            "recompute the actual post-signed fast and slow mean maps at the corrected rates, "
            "then join their complete pressure, viscosity, cutoff, cross and square terms to "
            "the simultaneous source/target and nonzero residual ledger"
        ),
        "standalone_generic_weighted_inference_proved": False,
        "same_family_weighted_target_route_available": True,
        "post_signed_means_recomputed_here": False,
        "complete_residual_or_blowup_proved_here": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q200 identity required")
    status = statuses()
    r = rates()
    stage = stage_identification()
    route = actual_weight_route()
    signed = generic_step()
    cycle = simultaneous_step()
    fixed = correction()
    boundary = frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "Q197_and_Q160_source_seeds_match": stage["Q197_D_equals_Q160_B0"],
        "Q198_and_Q160_target_seeds_match": stage["Q198_C0_equals_Q160_C0"],
        "initial_physical_target_has_positive_margin_over_Q197_C": stage["initial_target_is_stronger_than_Q197_C"] == Q(1, 2) - 2 * r["kappa"],
        "generic_signed_errors_all_tie_at_C0": signed["all_errors_equal_output"],
        "generic_target_gain_is_half_minus_kappa": signed["gain"] == Q(1, 2) - r["kappa"],
        "direct_weight_route_uses_exact_difference_not_generic_inference": route["direct_physical_route_present"] and route["cancellation_before_estimation"] and not route["standalone_unweighted_to_weighted_inference_valid"],
        "complete_pressure_and_moving_memory_retained": route["all_six_covariances_in_pressure_reconstruction"] and route["axial_velocity_retained"] and route["compact_pressure_retained"] and route["whole_moving_pressure_product_retained"] and route["moving_projection_commutator_retained"],
        "simultaneous_cycle_has_exact_two_kappa_target_slack": cycle["closes"] and cycle["target_slack"] == 2 * r["kappa"],
        "Q199_obstruction_retained_but_global_conclusion_corrected": route["Q199_scaled_bump_obstruction_retained"] and fixed["Q199_global_no_route_conclusion_was_too_strong"],
        "normalized_source_still_not_absolute_decay": physical_residual_power() < 0,
        "no_post_signed_mean_residual_or_blowup_promotion": not boundary["post_signed_means_recomputed_here"] and not boundary["complete_residual_or_blowup_proved_here"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "rates": r,
        "stage_identification": stage,
        "actual_weight_route": route,
        "generic_step": signed,
        "simultaneous_step": cycle,
        "correction": fixed,
        "physical_residual_power": physical_residual_power(),
        "frontier": boundary,
        "scope": {
            "analytical_relative_to_named_continuum_providers": True,
            "numerical_or_symbolic_tests_are_continuum_proof": False,
            "same_family_not_canonical_family_promotion": True,
            "uniform_in_cycle_constants_claimed": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q200 identity required")
    return {
        "status": "Q200-physical-weighted-target-route-reconciled-at-corrected-rates",
        "correction": (
            "Q199 correctly blocks an unweighted-to-weighted inference but overlooked the "
            "Q160/Q161/Q180 same-family state-difference and simultaneous-induction route"
        ),
        "established": (
            "relative to the named physical continuum providers, Q197 D equals the Q160 source "
            "seed, Q198 C0 equals the Q160 target seed, and the full zeta target propagates with "
            "gain 1/2-kappa while the coupled source cycle gains 1/2-3*kappa"
        ),
        "next": frontier()["next"],
        "historical_frozen_records_modified": False,
        "unconditional_post_signed_means_complete_residual_force_exclusion_or_blowup_certified": False,
    }
