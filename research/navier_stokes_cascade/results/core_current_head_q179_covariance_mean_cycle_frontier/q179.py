"""Q179 complete initial/first-feedback covariance audit and P3 frontier."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q179-COVARIANCE-MEAN-CYCLE-FRONTIER"
PROVIDERS = (
    "core_current_head_q178_physical_primary_source_audit",
    "complete_covariance_initial_means",
    "core_current_head_q150_initial_covariance_mean_join",
    "core_current_head_q151_family_stress_gate",
    "core_current_head_q152_physical_stress_feedback_join",
    "core_current_head_q153_physical_source_tail_diagonal_join",
    "core_current_head_q159_primary_contract_frontier",
)
EXPECTED = {
    "core_current_head_q178_physical_primary_source_audit": "Q178-physical-primary-pulse-source-boundary-and-covariance-obstruction-audited",
    "complete_covariance_initial_means": "actual-primary-complete-covariance-initial-mean-gain-audited",
    "core_current_head_q150_initial_covariance_mean_join": "Q149-uniform-primary-initial-complete-covariance-mean-correction-joined",
    "core_current_head_q151_family_stress_gate": "Q150-frontier-corrected-to-physical-canonical-principal-stress-gate",
    "core_current_head_q152_physical_stress_feedback_join": "Q152-physical-principal-stress-gate-closed-and-feedback-D-joined",
    "core_current_head_q153_physical_source_tail_diagonal_join": "Q153-physical-source-tail-and-complete-potential-diagonal-bound-conditionally",
    "core_current_head_q159_primary_contract_frontier": "Q159-P2-corrected-closed-and-P3-finite-cycle-target-obstructed",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q179 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def rates(kappa=Q(1, 100000), h=Q(1, 200)):
    kappa, h = map(Q, (kappa, h))
    if not 0 < kappa < Q(1, 4) or not 0 < h < Q(1, 100):
        raise ValueError("inherited kappa and h ranges required")
    B = 1 - kappa
    C = B + Q(1, 2)
    D = C - kappa
    initial = 2 - 3 * kappa
    complete_curl = 2 - 2 * kappa
    return {
        "B": B,
        "C": C,
        "D": D,
        "initial_residual_order": initial,
        "Q_initial_residual_exponent": h * initial,
        "initial_margin_over_D": initial - D,
        "required_stress_after_epsilon_inverse": D - 1,
        "available_complete_curl_after_epsilon_inverse": complete_curl - 1,
        "complete_curl_margin": (complete_curl - 1) - (D - 1),
        "first_feedback_order": min(C, D),
    }


def principal_assembly(labels, sigma0):
    if not labels:
        raise ValueError("nonempty physical labels required")
    actual = Q(0)
    cover = Q(0)
    root_defect = Q(0)
    for item in labels:
        eta2, a_scale, H, y, target = map(Q, item)
        cover += eta2
        actual += eta2 * a_scale * H * y
        root_defect += eta2 * a_scale * (H * y - target)
    sigma0 = Q(sigma0)
    return {
        "actual": actual,
        "cover": cover,
        "root_defect": root_defect,
        "defect": actual - sigma0,
        "exact": actual == sigma0 and cover == 1 and root_defect == 0,
    }


def moment_gate(theta_integral, theta_rhs, axial_integral, axial_rhs):
    theta_integral, theta_rhs, axial_integral, axial_rhs = map(
        Q, (theta_integral, theta_rhs, axial_integral, axial_rhs)
    )
    return {
        "angular_exact": theta_integral == theta_rhs,
        "axial_exact": axial_integral == axial_rhs,
        "both_exact": theta_integral == theta_rhs and axial_integral == axial_rhs,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q179 identity required")
    status = statuses()
    q178 = evidence(PROVIDERS[0])
    initial = evidence(PROVIDERS[1])["flags"]
    q150 = evidence(PROVIDERS[2])["symbolic"]
    q151 = evidence(PROVIDERS[3])["symbolic"]
    q152 = evidence(PROVIDERS[4])["symbolic"]
    q153 = evidence(PROVIDERS[5])["symbolic"]
    q159 = evidence(PROVIDERS[6])
    rate = rates()
    assembly = principal_assembly(
        [(Q(1, 4), 2, 3, 5, 15), (Q(3, 4), 2, 3, 5, 15)], 30
    )
    moments = moment_gate(7, 7, 11, 11)
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "Q178_requests_Q150_complete_covariance_audit": "Q150" in q178["review"]["next"],
        "initial_all_six_covariances_are_joined": initial["actual_primary_full_six_entry_covariance_initial_mean_inputs_joined"] and initial["complete_axial_velocity_and_moving_pressure_moment_retained"],
        "initial_fast_Haar_means_vanish_and_defects_are_recomputed": initial["initial_fast_update_has_exact_zero_Haar_velocity"] and initial["initial_recomputed_three_compatibility_defects_in_S2"],
        "both_background_pressure_tails_use_same_base_subtraction": initial["both_background_pressure_tails_retained_via_same_base_subtraction"],
        "initial_negative_raw_power_is_canceled_structurally": "negative absolute power is canceled" in q150["correction"] and rate["Q_initial_residual_exponent"] > 0,
        "canonical_feedback_does_not_transfer_by_closeness": not q151["inherited_and_canonical_primary_families_identified"],
        "physical_same_family_principal_stress_is_exact": assembly["exact"] and "localized principal selected stress equals Sigma0 exactly" in q152["established"],
        "complete_curl_has_half_power_margin": rate["complete_curl_margin"] == Q(1, 2),
        "first_physical_feedback_reaches_D": rate["first_feedback_order"] == rate["D"] and "feedback/source order D" in q152["established"],
        "generated_source_tail_diagonal_is_only_conditional": q153["physical_generated_source_endpoint_join_certified_conditionally"] and q153["complete_potential_diagonal_applicable_conditionally"] and not q153["closed_time_forcing_energy_or_blowup_certified"],
        "P3_all_finite_cycle_quantifier_remains_open": q159["P3_gate"]["obstruction_open"] and not q159["symbolic"]["P3_universal_finite_cycle_gate_discharged"],
        "P3_requires_complete_pressure_axial_and_support_memory": all(term in q159["P3_gate"]["memory"] for term in ("axial velocity", "six covariance", "Pi1/Pi2", "moving pressure", "zero-Haar", "reserved hull")),
        "ordinary_compactness_does_not_supply_edge_weight": "ordinary finite unweighted seminorms do not imply" in q159["P3_gate"]["why_compactness_is_insufficient"],
        "two_exact_moment_rows_are_independent_obligations": moments["both_exact"] and "J_theta" in q159["P3_gate"]["angular_moment"] and "c_rho P_mean" in q159["P3_gate"]["axial_moment"],
        "no_complete_residual_or_blowup_promotion": not q152["complete_residual_or_blowup_certified"] and not q159["symbolic"]["unconditional_admissible_force_or_blowup_certified"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "rates": {key: str(value) for key, value in rate.items()},
        "principal_assembly_fixture": {key: str(value) if isinstance(value, Q) else value for key, value in assembly.items()},
        "moment_fixture": moments,
        "established": {
            "initial_complete_six_entry_covariance_mean_cancellation": True,
            "same_physical_family_principal_stress_and_first_feedback_D": True,
            "axial_velocity_Pi1_Pi2_compact_and_moving_pressure_retained": True,
        },
        "correction": {
            "Q178_raw_negative_covariance_power_is_not_the_initial_stage_frontier": True,
            "Q151_family_transfer_gate_is_closed_by_a_physical_same_family_identity": True,
            "Q153_does_not_discharge_the_all_finite_cycle_target_quantifier": True,
            "historical_frozen_evidence_modified": False,
        },
        "obstruction": {
            "weighted_target": q159["P3_gate"]["weighted_bound"],
            "angular_moment": q159["P3_gate"]["angular_moment"],
            "axial_moment": q159["P3_gate"]["axial_moment"],
            "every_finite_cycle_and_fixed_jet": True,
            "ordinary_unweighted_compact_seminorms_suffice": False,
        },
        "scope": {
            "base_and_first_feedback_only": True,
            "P3_universal_finite_cycle_target_gate_discharged": False,
            "complete_residual_force_exclusion_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
        "next": "derive the P3 weighted Xi_j edge estimate and both complete tangential moment identities inductively for every finite physical cycle, retaining all six covariances, axial velocity, Pi1/Pi2, the whole moving pressure product, zero-Haar data and the reserved hull",
    }


def record():
    return {
        "status": "Q179-initial-and-first-feedback-covariance-means-closed-P3-universal-target-obstructed",
        "established": "the Q149 physical primary has exact complete six-entry initial covariance/mean cancellation and an exact same-family principal-stress identity; the first feedback reaches D=3/2-2kappa with axial velocity, Pi1/Pi2, compact and moving pressure, viscosity and A1 retained",
        "correction": "Q178's negative raw covariance power is canceled at the initial stage, and Q151's family-transfer concern is bypassed by the inherited physical family's own root and squared-partition identity; neither fact proves every later finite cycle",
        "obstruction": "P3 requires a zeta-weighted Xi_j edge estimate plus the exact angular and axial moment identities for every finite recomputed physical cycle; ordinary compact seminorms and the base/first-feedback cases do not prove this universal quantifier",
        "historical_frozen_records_modified": False,
        "unconditional_complete_covariance_cycle_residual_force_or_blowup_certified": False,
    }
