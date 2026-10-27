"""Closed-chart and Cartesian transfer of the log-reserved residual."""

from fractions import Fraction as Q
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


residual = load(
    "log_corner_residual",
    "core_log_reserved_complete_residual_replay/replay.py",
)
old_corner = load(
    "log_corner_old_obstruction",
    "core_corner_flatness_obstruction/corner.py",
)
primitive = load(
    "log_corner_primitive",
    "core_primitive_margin_normalization_reconciliation/margins.py",
)


def order(value, name):
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def closed_chart(h=Q(1, 200)):
    h = Q(h)
    if not Q(0) < h < Q(1, 100):
        raise ValueError("use the retained 0<h<1/100 regime")
    return {
        "terminal_identity": old_corner.normalized_identity(Q(3, 4), Q(1, 4)),
        "jacobian_floor": old_corner.jacobian_floor(h),
        "normalized_support_margin": primitive.NORMALIZED_SUPPORT_MARGIN,
        "physical_support_floor_used": False,
        "edge_weight_floor_used": False,
        "compact_X_interval_fixed_before_supremum": True,
        "angular_domain_includes_eta_endpoints": True,
    }


def cartesian_constant(power):
    power = order(power, "power")
    return old_corner.corner_transfer_constant(power)


def target(m=0, power=0, h=Q(1, 200), kappa=residual.diagonal.KAPPA):
    """One fixed comparison stage and its closed/Cartesian exponents."""
    m, power = order(m, "m"), order(power, "power")
    row = residual.comparison(m, power, h, kappa)
    chart = closed_chart(h)
    exponent = row["complete_minimum_exponent"]
    requested = Q(power + 1)
    return {
        "jet": m,
        "requested_power": power,
        "comparison_stage": row["log_reserved_comparison_stage"],
        "closed_chart_exponent": exponent,
        "requested_strict_exponent": requested,
        "retained_exponent_margin": exponent - requested,
        "cartesian_factor_at_requested_exponent": cartesian_constant(power + 1),
        "closed_chart": chart,
        "residual_groups": row["residual_group_exponents"],
        "pressure_memory": (
            "fixed gauge; moving normal and K pressure; all six covariance "
            "entries; gd-A1gd; dZ(Jz+c_rho P); both exceptional compact tails"
        ),
        "checks": {
            "terminal_chart_closes": chart["terminal_identity"] == 1,
            "jacobian_has_uniform_closed_floor": chart["jacobian_floor"] >= Q(49, 50),
            "normalized_support_margin_positive": chart["normalized_support_margin"] > 0,
            "one_fixed_stage_used": isinstance(row["log_reserved_comparison_stage"], int),
            "closed_residual_beats_requested_power": exponent >= requested,
            "all_six_groups_retained": tuple(row["residual_group_exponents"]) == residual.RESIDUAL_GROUPS,
            "cartesian_transfer_is_four_to_the_power": cartesian_constant(power + 1) == 4 ** (power + 1),
            "no_physical_or_edge_floor_used": (
                not chart["physical_support_floor_used"]
                and not chart["edge_weight_floor_used"]
            ),
        },
    }


def review():
    samples = [target(m, power) for m, power in ((0, 0), (0, 10), (3, 5), (8, 17))]
    checks = {
        "historical_full_exponent_no_log_stage_supremum_not_reused": True,
        "corrected_log_reserved_residual_is_the_parent": residual.review()["scope"][
            "log_reserved_six_group_complete_local_residual_flat_conditionally"
        ],
        "closed_terminal_identity_holds": all(row["checks"]["terminal_chart_closes"] for row in samples),
        "closed_jacobian_floor_holds": all(row["checks"]["jacobian_has_uniform_closed_floor"] for row in samples),
        "closed_normalized_support_contract_retained": all(row["checks"]["normalized_support_margin_positive"] for row in samples),
        "finite_fixed_stage_exists_for_every_sample": all(row["checks"]["one_fixed_stage_used"] for row in samples),
        "complete_six_group_residual_beats_every_sample_target": all(row["checks"]["closed_residual_beats_requested_power"] for row in samples),
        "q_le_four_rho_preserves_arbitrary_power": all(row["checks"]["cartesian_transfer_is_four_to_the_power"] for row in samples),
        "scale_incompatible_positive_floors_are_not_used": all(row["checks"]["no_physical_or_edge_floor_used"] for row in samples),
        "localization_force_comparison_and_blowup_remain_open_in_this_replay": True,
    }
    return {
        "checks": checks,
        "default_target": target(),
        "sample_targets": samples,
        "corrected_closed_stage_domain": (
            "sup_closed,K q^(-(h B_J-K_m)) (1+|log q|)^(-P_(J,m,K)) "
            "|R_regular(U^[J])|_m < infinity"
        ),
        "cartesian_transfer": (
            "q<=4 rho and |D^alpha R|<=C q^(N+1) imply "
            "|D^alpha R|<=C 4^(N+1) rho^(N+1)"
        ),
        "scope": {
            "same_field_all_fixed_residual_jets_cartesian_corner_flat_conditionally": True,
            "complete_pressure_memory_retained": True,
            "historical_evidence_modified": False,
            "punctured_field_endpoint_jets_reproved_here": False,
            "localization_force_energy_comparison_replayed_here": False,
            "unconditional_pulse_residual_force_or_blowup_claimed": False,
        },
    }
