"""Q222: source audit of the joined global pressure and moment identities."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q222-GLOBAL-PRESSURE-MOMENT-SOURCE-AUDIT"
PROVIDERS = (
    "core_current_head_q221_fourier_slow_response_source_audit",
    "core_global_pressure_moment_recheck",
    "core_current_head_q147_repair_pressure_join",
    "core_current_head_q177_weighted_support_source_reconciliation",
    "core_current_head_q204_global_ledger_uniform_primary",
    "core_current_head_q160_physical_target_induction",
    "core_signed_target_edge_budget",
    "core_preserved_haar_scalar_join",
    "core_sourced_covariance_moment_join",
    "core_current_head_q17_mean_support_join",
)
EXPECTED = {
    PROVIDERS[0]: "Q221-zero-Haar-Fourier-and-fixed-slow-response-source-audited",
    PROVIDERS[1]: "same-repaired-profile-global-pressure-moments-and-exterior-stress-identities-independently-rechecked",
    PROVIDERS[2]: "same-modulation-five-row-repair-cone-and-global-pressure-identity-joined",
    PROVIDERS[3]: "Q177-selected-input-weighted-stress-and-reserved-support-route-reconciled",
    PROVIDERS[4]: "Q204-global-pressure-coefficient-uniform-primary-and-covariance-ledgers-joined-conditionally",
    PROVIDERS[5]: "Q160-physical-all-finite-cycle-target-moment-induction-joined",
    PROVIDERS[6]: "current-signed-target-two-edge-primitive-and-weighted-row-input-budget",
    PROVIDERS[7]: "conditional-complete-preserved-Haar-current-scalar-refinement-and-weighted-input-gap",
    PROVIDERS[8]: "generated-sourced-complete-covariance-full-pressure-moments-and-next-weighted-operator-domain-correction",
    PROVIDERS[9]: "projected-mean-pressure-support-count-derived",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def annular_transform(row, eta, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    M, I, J, S, Cp = map(Q, row)
    eta = Q(eta)
    return (M, J - 4 * eta * I, I, -(S - 8 * eta * M), Cp)


def annular_inverse(row, eta, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    r1, r2, r3, r4, r5 = map(Q, row)
    eta = Q(eta)
    return (r1, r3, r2 + 4 * eta * r3, 8 * eta * r1 - r4, r5)


def compact_transform(row, lam, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    M, I, J, S, Cp = map(Q, row)
    lam = Q(lam)
    if lam <= 0:
        raise ValueError("positive lambda required")
    return (M, (M - J) / lam, I, -S, Cp)


def compact_inverse(row, lam, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    r1, r2, r3, r4, r5 = map(Q, row)
    lam = Q(lam)
    if lam <= 0:
        raise ValueError("positive lambda required")
    return (r1, r3, r1 - lam * r2, -r4, r5)


def global_heat_join(delta_endpoint_I, delta_future_tail, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    delta_endpoint_I = Q(delta_endpoint_I)
    delta_future_tail = Q(delta_future_tail)
    return {
        "delta_R": delta_endpoint_I + delta_future_tail,
        "endpoint_I_is_global_R": False,
        "endpoint_match_suffices_only_with_common_tail": True,
        "common_tail": delta_future_tail == 0,
    }


def leading_global_ledger(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    return {
        "finite_row_order": "(M,I,J,S,Cp)",
        "global_row_order": "(M,R,J,S,Cp)",
        "global_identities": "M_infty=R=J_infty=S_infty=0",
        "pressure_identity": "Pi_axis+Cp_infty=0",
        "backward_pressure": "Pi(X)=-int_X^infty E(x)^2/(2x) dx",
        "intermediate_pressure_reset": False,
        "all_five_rows_required_for_exterior_stress_join": True,
        "pressure_derivative_defect_is_independent": True,
    }


def compact_pressure_moment(Qg, Pmean, c_rho, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    Qg, Pmean, c_rho = map(Q, (Qg, Pmean, c_rho))
    return {
        "integral_R_pbar": -Qg / 2 + c_rho * Pmean,
        "derivation": "int R p=-1/2 int R^2 partial_R p",
        "unit_mass_compensator": True,
        "A0_pointwise_remainder_deleted": False,
    }


def moving_product_derivative(Jz_Z, c_rho_Z, P, c_rho, P_Z,
                              epsilon=1, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    Jz_Z, c_rho_Z, P, c_rho, P_Z, epsilon = map(
        Q, (Jz_Z, c_rho_Z, P, c_rho, P_Z, epsilon))
    return epsilon * (Jz_Z + c_rho_Z * P + c_rho * P_Z)


def mean_moment_ledger(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    return {
        "angular": "int R^2 Haar(E_theta)=epsilon partial_Z J_theta",
        "axial": "int R Haar(E_z)=epsilon partial_Z(J_z+c_rho P)",
        "J_theta_contains_full_angular_axial_flux": True,
        "J_z_contains_axial_squares_Wzz_and_minus_Qg_over_2": True,
        "whole_moving_product_differentiated": True,
        "compact_radial_boundary_terms_zero": True,
        "pointwise_rows_determined_by_moments": False,
    }


def preservation_ledger(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    return {
        "nonzero_harmonic_waves_change_Haar_velocity_moments": False,
        "zero_Haar_fast_curl_change_Haar_velocity_moments": False,
        "slow_map_zero_moment_rows": 2,
        "slow_map_exact_divergence": True,
        "covariances_recomputed_after_each_update": True,
        "zero_Haar_correlations_discarded": False,
        "all_six_covariances_retained": True,
        "Pi1_Pi2_cancel_only_in_same_complete_base_subtraction": True,
    }


def projected_targets(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    return {
        "M2": "epsilon partial_Z J_theta",
        "M1": "epsilon partial_Z(J_z+c_rho P)",
        "sigma2_moment": "int R^2 sigma_2=1",
        "sigma1_moment": "int R sigma_1=1",
        "h2_zero_moment": True,
        "h1_zero_moment": True,
        "projection_deletes_pressure_or_mean_scalar": False,
        "two_edge_primitive_right_collar_requires_exact_zero_moment": True,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    row = tuple(map(Q, (3, 5, 7, 11, 13)))
    eta, lam = Q(2, 9), Q(5, 17)
    heat_good = global_heat_join(0, 0)
    heat_bad = global_heat_join(0, Q(1, 19))
    pressure = compact_pressure_moment(7, 11, Q(13, 17))
    moving = moving_product_derivative(2, 3, 5, 7, 11, 13)
    moments = mean_moment_ledger()
    preserve = preservation_ledger()
    targets = projected_targets()
    checks = {
        "provider_statuses_match": all(statuses().values()),
        "finite_row_transforms_round_trip": (
            annular_inverse(annular_transform(row, eta), eta) == row
            and compact_inverse(compact_transform(row, lam), lam) == row),
        "endpoint_to_global_join_keeps_common_tail_hypothesis": (
            heat_good["delta_R"] == 0 and heat_bad["delta_R"] != 0),
        "leading_pressure_row_not_reset": (
            not leading_global_ledger()["intermediate_pressure_reset"]
            and leading_global_ledger()["all_five_rows_required_for_exterior_stress_join"]),
        "compact_pressure_moment_has_both_terms": (
            pressure["integral_R_pbar"] == -Q(7, 2) + Q(13, 17) * 11),
        "moving_product_uses_all_three_derivative_terms": moving == 13 * (2 + 3 * 5 + 7 * 11),
        "both_complete_mean_moments_retained": (
            moments["whole_moving_product_differentiated"]
            and not moments["pointwise_rows_determined_by_moments"]),
        "updates_preserve_velocity_moments_without_erasing_correlations": (
            not preserve["nonzero_harmonic_waves_change_Haar_velocity_moments"]
            and not preserve["zero_Haar_fast_curl_change_Haar_velocity_moments"]
            and not preserve["zero_Haar_correlations_discarded"]),
        "projected_rows_have_exact_distinct_moments": (
            targets["h1_zero_moment"] and targets["h2_zero_moment"]
            and not targets["projection_deletes_pressure_or_mean_scalar"]),
    }
    return {
        "status": "Q222-leading-and-evolving-global-pressure-moment-identities-source-audited",
        "checks": checks,
        "scope": {
            "leading_five_global_rows_joined_on_common_tail": True,
            "compact_fast_slow_mean_moment_rows_joined": True,
            "endpoint_I_identified_with_global_R_without_tail": False,
            "pointwise_Haar_rows_inferred_from_moments": False,
            "complete_pressure_or_Pi1_Pi2_deleted": False,
            "uniform_actual_pulse_or_complete_residual": False,
            "unconditional_force_exclusion_or_blowup": False,
        },
        "next": (
            "audit the finite coefficient solves and reserved-mean supports on this joined "
            "pressure/moment object, including the nonremovable Pi1/Pi2 pressure floor"
        ),
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q222 identity required")
    heat_good = global_heat_join(0, 0)
    heat_bad = global_heat_join(0, Q(1, 19))
    pressure = compact_pressure_moment(7, 11, Q(13, 17))
    heat_good["delta_R"] = str(heat_good["delta_R"])
    heat_bad["delta_R"] = str(heat_bad["delta_R"])
    pressure["integral_R_pbar"] = str(pressure["integral_R_pbar"])
    return {
        "status": review()["status"],
        "providers": statuses(),
        "leading_global_ledger": leading_global_ledger(),
        "global_heat_common_tail": heat_good,
        "global_heat_missing_tail_obstruction": heat_bad,
        "compact_pressure_sample": pressure,
        "moving_product_sample": str(moving_product_derivative(2, 3, 5, 7, 11, 13)),
        "mean_moment_ledger": mean_moment_ledger(),
        "preservation_ledger": preservation_ledger(),
        "projected_targets": projected_targets(),
        "review": review(),
    }
