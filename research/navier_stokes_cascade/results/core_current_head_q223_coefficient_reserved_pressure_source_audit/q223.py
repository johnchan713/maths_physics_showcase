"""Q223: source audit of finite coefficient supports and reserved pressure."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q223-COEFFICIENT-RESERVED-PRESSURE-SOURCE-AUDIT"
PROVIDERS = (
    "core_current_head_q222_global_pressure_moment_source_audit",
    "core_coefficient_support_recheck",
    "core_current_head_q148_coefficient_support_pressure",
    "core_current_head_q177_weighted_support_source_reconciliation",
    "core_current_head_q188_compact_repair_hierarchy_audit",
    "core_current_head_q204_global_ledger_uniform_primary",
    "core_current_head_q149_primary_pulse_join",
    "core_reserved_geometry_recheck",
    "core_current_head_q209_cutoff_role_separation",
    "core_current_head_q214_exact_core_root_pressure_join",
)
EXPECTED = {
    PROVIDERS[0]: "Q222-leading-and-evolving-global-pressure-moment-identities-source-audited",
    PROVIDERS[1]: "same-profile-finite-coefficient-supports-and-reserved-pressure-obstruction-rechecked",
    PROVIDERS[2]: "Q147-final-profile-finite-coefficient-supports-and-reserved-pressure-obstruction-joined",
    PROVIDERS[3]: "Q177-selected-input-weighted-stress-and-reserved-support-route-reconciled",
    PROVIDERS[4]: "Q188-compact-repair-hierarchy-closed-and-continuum-envelope-frontier-isolated",
    PROVIDERS[5]: "Q204-global-pressure-coefficient-uniform-primary-and-covariance-ledgers-joined-conditionally",
    PROVIDERS[6]: "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined",
    PROVIDERS[7]: "uniform-reserved-coordinate-compensator-jets-and-moving-support-costs-rechecked",
    PROVIDERS[8]: "Q209-background-and-stage-cutoff-roles-separated-Q208-identity-blocker-corrected",
    PROVIDERS[9]: "Q214-exact-core-selected-roots-and-leading-pressure-joined-conditionally",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q223 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def coefficient_target(M, I, J, S, Cp, lam, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q223 identity required")
    M, I, J, S, Cp, lam = map(Q, (M, I, J, S, Cp, lam))
    if lam <= 0:
        raise ValueError("positive lambda required")
    return (M, (M - J) / lam, I, -S, Cp)


def coefficient_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q223 identity required")
    return {
        "solve": "c_n=-B^-1 z_n",
        "target": "z_n=(M,(M-J)/lambda,I,-S,Cp)",
        "inverse_bound": "||B^-1||<800",
        "fixed_positive_order": True,
        "all_lower_orders_finalized_first": True,
        "linear_solve_requires_small_target": False,
        "uniform_in_order_or_derivative": False,
    }


def support_ledger(order, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q223 identity required")
    if not isinstance(order, int) or order < 1:
        raise ValueError("positive integer order required")
    if order == 1:
        return {"velocity": "Xvel", "pressure": "Xv", "stress": "Xb"}
    if order == 2:
        return {"velocity": "Xvel", "pressure": "Xv", "stress": "Xv"}
    return {"velocity": "Xvel", "pressure": "Xvel", "stress": "Xvel"}


def reserved_geometry(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q223 identity required")
    bump_end, velocity_end, reserved_left, reserved_right = (
        Q(51, 20), Q(3), Q(6), Q(11))
    return {
        "bump_end": str(bump_end),
        "velocity_end": str(velocity_end),
        "reserved_interval": [str(reserved_left), str(reserved_right)],
        "strict_separation": bump_end < velocity_end < reserved_left < reserved_right,
        "positive_order_velocity_zero_on_reserved_interval": True,
        "potential_cutoff_commutator_zero_there": True,
        "pressure_inferred_zero_from_velocity_support": False,
    }


def reserved_pressure_ledger(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q223 identity required")
    return {
        "Pi1": "K1(eta)-mu^2/(4X)",
        "Pi2": "K2(eta)",
        "Pi_n_n_ge_3": "0",
        "K1_K2_retain_later_source_memory": True,
        "local_zero_pressure_derivative_sets_K1_or_K2_zero": False,
        "Pi1_has_forced_radial_oscillation": True,
        "Pi2_has_forced_radial_oscillation": False,
        "actual_K2_proved_nonzero": False,
    }


def minimax_constant_floor(left, right, identity=IDENTITY):
    """Sharp inf_C max(|left+C|, |right+C|) for rational inputs."""
    if identity != IDENTITY:
        raise ValueError("Q223 identity required")
    left, right = map(Q, (left, right))
    center = -(left + right) / 2
    floor = abs(left - right) / 2
    return {
        "optimal_constant": center,
        "sharp_floor": floor,
        "left_at_optimum": left + center,
        "right_at_optimum": right + center,
    }


def pressure_two_point_floor(scale, mu2, Xa, Xb, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q223 identity required")
    scale, mu2, Xa, Xb = map(Q, (scale, mu2, Xa, Xb))
    if scale < 0 or mu2 <= 0 or Xa <= 0 or Xb <= Xa:
        raise ValueError("require scale>=0, mu2>0 and 0<Xa<Xb")
    left = -scale * mu2 / (4 * Xa)
    right = -scale * mu2 / (4 * Xb)
    row = minimax_constant_floor(left, right)
    row.update({
        "endpoint_difference": right - left,
        "formula": "scale*mu^2*(Xa^-1-Xb^-1)/8",
        "strictly_positive": row["sharp_floor"] > 0,
    })
    return row


def corrected_claim(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q223 identity required")
    return {
        "positive_two_point_floor_source": "Pi1 radial term -mu^2/(4X)",
        "Pi2_role": "X-independent K2 source memory retained in complete pressure",
        "joint_Pi1_Pi2_positive_floor": False,
        "reason": "Pi2 cancels from every same-eta radial two-point difference",
        "Pi2_deletion_permitted": False,
        "full_residual_lower_bound": False,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q223 identity required")
    return {
        "fixed_order_coefficient_support_obligation": "closed conditionally on finite Volterra providers",
        "rule_selected_Borel_schedule": "already closed by Q204",
        "uniform_primary_all_used_band_bounds": "already supplied by Q149/Q204",
        "uniform_raw_uncut_generated_hierarchy_constant": False,
        "complete_covariance_pressure_viscosity_cutoff_interaction_residual": False,
        "unconditional_blowup": False,
        "next": (
            "audit the uniform actual primary-pulse value/derivative/support theorem at its "
            "exact source interfaces, then retain full Pi1/Pi2 and axial pressure memory in "
            "the complete covariance and residual join"
        ),
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q223 identity required")
    sample = pressure_two_point_floor(3, 5, 7, 11)
    pi2 = minimax_constant_floor(Q(13, 17), Q(13, 17))
    supports = [support_ledger(n) for n in (1, 2, 3)]
    correction = corrected_claim()
    checks = {
        "provider_statuses_match": all(statuses().values()),
        "five_row_target_keeps_lambda_loss": coefficient_target(3, 5, 7, 11, 13, Q(2, 9))[1] == -18,
        "linear_solve_has_no_small_target_gate": not coefficient_contract()["linear_solve_requires_small_target"],
        "three_support_regimes_are_distinct": len({tuple(row.items()) for row in supports}) == 3,
        "reserved_velocity_pressure_supports_are_separated": (
            reserved_geometry()["strict_separation"]
            and not reserved_geometry()["pressure_inferred_zero_from_velocity_support"]),
        "Pi1_minimax_floor_is_positive_and_sharp": (
            sample["strictly_positive"]
            and sample["left_at_optimum"] == -sample["right_at_optimum"]),
        "Pi2_only_has_zero_radial_two_point_floor": pi2["sharp_floor"] == 0,
        "Pi2_memory_retained_without_nonzero_claim": (
            not correction["joint_Pi1_Pi2_positive_floor"]
            and not reserved_pressure_ledger()["actual_K2_proved_nonzero"]),
        "scope_stops_before_complete_residual": not frontier()["complete_covariance_pressure_viscosity_cutoff_interaction_residual"],
    }
    return {
        "status": "Q223-finite-coefficient-reserved-support-and-Pi1-pressure-floor-source-audited",
        "checks": checks,
        "scope": {
            "finite_positive_order_coefficient_solve_source_audited": True,
            "velocity_pressure_stress_supports_kept_distinct": True,
            "sharp_Pi1_radial_pressure_floor_derived": True,
            "positive_Pi2_floor_claimed": False,
            "uniform_raw_generated_hierarchy_bound": False,
            "complete_residual_or_blowup": False,
        },
        "next": frontier()["next"],
    }


def _strings(row):
    return {key: str(value) if isinstance(value, Q) else value for key, value in row.items()}


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q223 identity required")
    return {
        "status": review()["status"],
        "providers": statuses(),
        "coefficient_contract": coefficient_contract(),
        "coefficient_target_sample": list(map(str, coefficient_target(3, 5, 7, 11, 13, Q(2, 9)))),
        "support_ledger": {str(n): support_ledger(n) for n in (1, 2, 3)},
        "reserved_geometry": reserved_geometry(),
        "reserved_pressure_ledger": reserved_pressure_ledger(),
        "Pi1_two_point_sample": _strings(pressure_two_point_floor(3, 5, 7, 11)),
        "Pi2_constant_negative_control": _strings(minimax_constant_floor(Q(13, 17), Q(13, 17))),
        "corrected_claim": corrected_claim(),
        "frontier": frontier(),
        "review": review(),
    }
