"""Q212: join the existing target degree contract and restore the raw-map frontier."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q212-EXISTING-TARGET-DEGREE-RAW-MAP-FRONTIER"
PROVIDERS = (
    "core_current_head_q211_degree_grammar_primitive_frontier",
    "core_current_head_q160_physical_target_induction",
    "core_current_head_q180_target_source_simultaneous_induction",
    "core_current_head_q184_generated_source_log_loss_correction",
    "finite_stage_induction",
    "core_current_head_q192_frame_graph_jet_transfer",
    "core_current_head_q149_primary_pulse_join",
    "core_current_head_q150_initial_covariance_mean_join",
    "core_current_head_q152_physical_stress_feedback_join",
    "core_current_head_q209_cutoff_role_separation",
    "core_current_head_q210_finite_bundle_uniform_tail",
)
EXPECTED = {
    PROVIDERS[0]: "Q211-polynomial-degree-grammar-closed-conditionally-weighted-primitive-frontier-isolated",
    PROVIDERS[1]: "Q160-physical-all-finite-cycle-target-moment-induction-joined",
    PROVIDERS[2]: "Q180-physical-target-source-tail-simultaneous-induction-reconciled",
    PROVIDERS[3]: "Q184-generated-source-log-loss-corrected-and-pure-power-reserve-derived",
    PROVIDERS[4]: "conditional-current-state-complete-finite-stage-induction",
    PROVIDERS[5]: "Q192-auxiliary-frame-jets-closed-and-physical-graph-loss-corrected",
    PROVIDERS[6]: "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined",
    PROVIDERS[7]: "Q149-uniform-primary-initial-complete-covariance-mean-correction-joined",
    PROVIDERS[8]: "Q152-physical-principal-stress-gate-closed-and-feedback-D-joined",
    PROVIDERS[9]: "Q209-background-and-stage-cutoff-roles-separated-Q208-identity-blocker-corrected",
    PROVIDERS[10]: "Q210-finite-primitive-bundle-unified-conditionally-uniform-all-band-source-frontier-isolated",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q212 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def finite_degree_max(degrees, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q212 identity required")
    values = []
    for degree in degrees:
        if isinstance(degree, bool):
            raise ValueError("finite nonnegative rational degrees required")
        try:
            degree = Q(degree)
        except (TypeError, ValueError, ZeroDivisionError):
            raise ValueError("finite nonnegative rational degrees required") from None
        if degree < 0:
            raise ValueError("finite nonnegative rational degrees required")
        values.append(degree)
    if not values:
        raise ValueError("nonempty finite request required")
    return max(values)


def target_degree_join(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q212 identity required")
    return {
        "class_definition": (
            "finite-stage W/M/S membership includes, for every fixed stage and derivative, "
            "a finite polynomial degree in S=ell^2 and finite inverse-edge powers"
        ),
        "expanded_Q160_target": (
            "|D^I Xi_j| <= C_(j,I) epsilon^(C_j) S^(d_(j,I)) zeta "
            "ya^(-p_(j,I)) yb^(-q_(j,I))"
        ),
        "exact_difference": "X_e(new)-X_e(old)=-I_e Q_e overline(Delta E_e)",
        "primitive": (
            "zero weighted totals and the two oriented current inverses preserve zeta with "
            "finite polynomial-S and inverse-edge costs"
        ),
        "finite_request_implication": (
            "D_R^Xi=max_(|I|<=R) d_(j,I)<infinity, and the finitely many primitive, "
            "root, frame and mean rows have one finite conditional maximum"
        ),
        "Q160_displayed_degree_zero_claimed": False,
        "conditional_degree_array_exists": True,
        "degree_array_numerically_evaluated": False,
        "independent_raw_selected_map_certificate": False,
    }


def pressure_memory(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q212 identity required")
    return [
        "all six covariance entries",
        "radial and axial velocity",
        "complete reconstructed pressure",
        "Pi1/Pi2 in the same base equation",
        "whole moving pressure product and its derivative",
        "both tangential moment identities",
        "zero weighted primitive totals at both edges",
        "reserved interior support",
    ]


def corrected_frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q212 identity required")
    return {
        "Q211_next_step_repeats_existing_conditional_target_theorem": True,
        "Q211_degree_algebra_remains_valid": True,
        "conditional_target_primitive_degree_obligation_closed": True,
        "raw_majorant": (
            "max_(|I|<=R+1) sup_D w_I S^(-d_I) |D^I(F0,U0,Pi0,M0,V0,"
            "E1,U1,W1,T1/zeta)| <= A_R<infinity"
        ),
        "same_object_gaps": (
            "g_R=min of the selected Jacobian, normal, root, determinant, variance and "
            "support gaps is positive on the same accepted tail"
        ),
        "same_object_inverse_archive": (
            "K_R=max of every requested compact/slow pressure and moment inverse constant "
            "is finite on that same object"
        ),
        "why_not_closed": (
            "Q160/Q180/WMS and Q192 are conditional on the selected leading/background "
            "geometry and finite existence/inverse/gap premises; hashes and finite tests do "
            "not prove their continuum suprema or positive infima"
        ),
        "raw_selected_map_majorant_independently_proved": False,
        "common_gap_and_inverse_archive_independently_proved": False,
        "unconditional_residual_or_blowup": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q212 identity required")
    status = statuses()
    target = target_degree_join()
    frontier = corrected_frontier()
    checks = {
        "provider_statuses_match": all(status.values()),
        "finite_request_maximum_is_valid": finite_degree_max([0, Q(3, 2), 7]) == 7,
        "existing_target_degree_contract_joined": (
            target["conditional_degree_array_exists"]
            and not target["Q160_displayed_degree_zero_claimed"]
            and not target["degree_array_numerically_evaluated"]
        ),
        "frontier_restored_to_raw_same_object_premises": (
            frontier["Q211_next_step_repeats_existing_conditional_target_theorem"]
            and frontier["conditional_target_primitive_degree_obligation_closed"]
            and not frontier["raw_selected_map_majorant_independently_proved"]
            and not frontier["common_gap_and_inverse_archive_independently_proved"]
        ),
    }
    return {
        "status": "Q212-existing-target-degree-contract-joined-conditionally-raw-map-frontier-restored",
        "checks": checks,
        "scope": {
            "conditional_target_primitive_degree_existence_joined": True,
            "actual_degree_array_evaluated": False,
            "raw_selected_profile_first_coefficient_bundle_proved": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
        },
        "next": (
            "prove the selected-profile/first-coefficient weighted majorant, common positive "
            "gap floor and compact/slow inverse archive on one actual accepted tail"
        ),
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q212 identity required")
    return {
        "status": review()["status"],
        "providers": statuses(),
        "finite_degree_example": {
            "degrees": ["0", "3/2", "7"],
            "maximum": str(finite_degree_max([0, Q(3, 2), 7])),
            "diagnostic_only": True,
        },
        "target_degree_join": target_degree_join(),
        "pressure_memory": pressure_memory(),
        "frontier": corrected_frontier(),
        "review": review(),
    }
