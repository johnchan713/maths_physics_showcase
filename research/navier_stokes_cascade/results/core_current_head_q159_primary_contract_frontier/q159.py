"""Correct the Q158 P2 frontier and isolate the first genuine finite-cycle gate."""
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q159-same-physical-primary-contract-and-finite-cycle-frontier"
PROVIDERS = (
    "core_current_head_q158_physical_dependency_frontier",
    "core_current_head_q149_primary_pulse_join",
    "core_leading_background_borel_join",
    "physical_pulse_transfer",
    "core_current_head_q150_initial_covariance_mean_join",
    "core_current_head_q152_physical_stress_feedback_join",
    "source_envelope_edge_audit",
    "finite_stage_induction",
    "core_current_head_q153_physical_source_tail_diagonal_join",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def p2_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q159 identity required")
    return {
        "object": (
            "Q147 finalized profile, its rule-defined real Borel background, "
            "and Q149 inherited physical primary family"
        ),
        "coefficient_quantifier": (
            "for every fixed n and m<=n, C_(n,m)<infinity before b_n=log a_n is chosen"
        ),
        "scale_rule": (
            "b_n>1+max(b_(n-1)+log2,B+log2,(log C_(n,m)+n log2)/(nh), all physical rows)"
        ),
        "coefficient_bound": "C_(n,m) q^(nh)<2^(-n) on q<=a_n^(-1)",
        "pulse_quantifier": (
            "one derivative-independent later band ell0; for every fixed physical (a,b), "
            "finite constants uniform over every used ell>=ell0 and label"
        ),
        "velocity": "C_ab Q^(-(a+1)sx-b st) S^B sqrt(zeta) delta^(-B) P",
        "pressure": "C_ab Q^(-(a+2)sx-b st) S^B sqrt(zeta) delta^(-B) P",
        "potential": "C_ab Q^(-a sx-b st) S^B sqrt(zeta) delta^(-B) P",
        "support": (
            "fixed closed transverse/time supports; potential curls do not enlarge support; "
            "smooth zero extension and one-sided q=0 background jets"
        ),
        "complete_pressure_and_cutoff_curl_retained": True,
        "new_profile_or_scale_sequence_selected": False,
        "analytically_discharged": True,
    }


def p3_gate(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q159 identity required")
    return {
        "base_case": (
            "Q150 proves the initial six-entry covariance/mean correction on Q149's physical family"
        ),
        "first_feedback": (
            "Q152 proves the same-family principal stress identity and conditional D-order feedback"
        ),
        "universal_request": (
            "for every finite cycle j and fixed mixed jet I, the actual recomputed signed/full-"
            "covariance target Xi_j must obey one weighted edge bound and both exact moments"
        ),
        "weighted_bound": (
            "|D^I Xi_j|<=C_(j,I) epsilon^(alpha_j) zeta "
            "ya^(-p_(j,I)) yb^(-q_(j,I))"
        ),
        "angular_moment": "integral R^2 overline(E_theta,j)=epsilon partial_Z J_theta,j",
        "axial_moment": (
            "integral R overline(E_z,j)=epsilon partial_Z(J_z,j+c_rho P_mean,j)"
        ),
        "memory": (
            "current axial velocity, all six covariance entries, Pi1/Pi2, moving pressure, "
            "both velocity moments, zero-Haar data and the reserved hull"
        ),
        "why_compactness_is_insufficient": (
            "zeta vanishes at both flat edges and reciprocal signed roots lose half a weight; "
            "ordinary finite unweighted seminorms do not imply the displayed divisibility"
        ),
        "source_envelope_provider_marks_every_later_target_verified": False,
        "finite_stage_result_is_conditional_on_this_request": True,
        "obstruction_open": True,
    }


def corrected_frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q159 identity required")
    return (
        {"id": "P1", "status": "closed", "contract": "PREMOD value/root envelope"},
        {"id": "P2", "status": "closed-by-Q149-and-Q159-audit", "contract": "Borel background and physical primary pulses"},
        {"id": "P3", "status": "open-at-all-finite-cycle-target-quantifier", "contract": "complete covariance/mean targets"},
        {"id": "P4", "status": "downstream-conditional", "contract": "generated-source quotient and endpoint extension"},
        {"id": "P5", "status": "downstream-conditional", "contract": "weighted inverses and physical tails"},
        {"id": "P6", "status": "downstream-conditional", "contract": "complete diagonal, endpoint and localized force"},
        {"id": "P7", "status": "conditional-consequence", "contract": "comparison and core exclusion"},
    )


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q159 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "Q158-physical-dependency-frontier-audited-and-first-loose-provider-bound",
        PROVIDERS[1]: "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined",
        PROVIDERS[2]: "same-leading-profile-and-rule-defined-background-borel-contract-joined",
        PROVIDERS[3]: "actual-background-to-physical-primary-pulse-transfer-audited",
        PROVIDERS[4]: "Q149-uniform-primary-initial-complete-covariance-mean-correction-joined",
        PROVIDERS[5]: "Q152-physical-principal-stress-gate-closed-and-feedback-D-joined",
        PROVIDERS[6]: "regular-source-envelope-edge-estimates-with-explicit-freezing-correction",
        PROVIDERS[7]: "conditional-current-state-complete-finite-stage-induction",
        PROVIDERS[8]: "Q153-physical-source-tail-and-complete-potential-diagonal-bound-conditionally",
    }
    for name, status in expected.items():
        if data[name].get("status") != status:
            raise ValueError(f"Q159 provider status mismatch: {name}")
    borel = data[PROVIDERS[2]]["flags"]
    pulse = data[PROVIDERS[3]]["flags"]
    source = data[PROVIDERS[6]]["flags"]
    induction = data[PROVIDERS[7]]["flags"]
    if not all((
        borel["all_positive_background_coefficients_exist_at_every_fixed_order"],
        borel["rule_defined_borel_sum_identified_with_actual_background"],
        borel["full_pressure_and_curl_cutoff_commutator_retained"],
        pulse["actual_background_inputs_to_primary_pulse_theorem_discharged"],
        pulse["common_all_later_band_primary_pulse_values_and_fixed_jets_derived"],
        pulse["complete_primary_curl_velocity_and_wave_pressure_physically_transferred"],
        pulse["exact_closed_support_and_one_sided_boundary_conditions_retained"],
    )):
        raise ValueError("Q159 P2 provider contract missing")
    if source["every_actual_later_covariance_mean_target_verified_by_grammar_or_protocol"]:
        raise ValueError("Q159 P3 obstruction was silently promoted")
    if induction["independent_foundation_review_complete"]:
        raise ValueError("Q159 conditional induction was silently promoted")
    return {
        "all_statuses_match": True,
        "P2_same_object_provider_chain_complete": True,
        "P3_all_finite_cycle_target_gate_open": True,
        "Q153_downstream_join_remains_conditional": True,
    }


def record():
    return {
        "status": "Q159-P2-corrected-closed-and-P3-finite-cycle-target-obstructed",
        "correction": (
            "Q158 incorrectly listed P2 as a conditional import: Q149 already binds the "
            "rule-defined Borel theorem and physical-pulse transfer to the exact finalized profile"
        ),
        "established": (
            "every fixed coefficient seminorm precedes its Borel scale, and complete primary "
            "velocity, pressure, potential, fixed physical jets and supports hold uniformly on "
            "every used band of the same Q149 family"
        ),
        "obstruction": (
            "P3 still requires the weighted Xi_j edge estimate and both complete tangential "
            "moment identities for every finite recomputed physical cycle; Q150 and Q152 close "
            "the base and first feedback but do not discharge that universal quantifier"
        ),
        "scope": (
            "analytical provider audit and a precise remaining inequality; numerical tests and "
            "dependency closure are not continuum proofs of P3"
        ),
        "next": (
            "derive the P3 weighted target and moment identities by induction on the actual "
            "physical cycle expressions, retaining full axial velocity, six covariances, pressure "
            "memory and zero-Haar/reserved-hull supports before reusing Q153"
        ),
        "P2_analytically_discharged": True,
        "P3_universal_finite_cycle_gate_discharged": False,
        "unconditional_admissible_force_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
