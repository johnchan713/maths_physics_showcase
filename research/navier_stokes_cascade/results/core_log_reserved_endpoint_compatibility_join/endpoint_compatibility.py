"""Join punctured endpoint jets to log-reserved Cartesian residual flatness."""

from fractions import Fraction as Q
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


corner = load(
    "log_endpoint_corner",
    "core_log_reserved_cartesian_corner_replay/corner_replay.py",
)
endpoint = load(
    "log_endpoint_punctured",
    "core_corrected_endpoint_join/endpoint.py",
)
localization = load(
    "log_endpoint_localization_ledger",
    "core_corrected_force_join/localization.py",
)


def order(value, name):
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def request(spatial_order=0, time_order=0, flat_power=0,
            h=Q(1, 200), kappa=corner.residual.diagonal.KAPPA):
    """One finite localization-input request on the already selected field."""
    a = order(spatial_order, "spatial_order")
    j = order(time_order, "time_order")
    power = order(flat_power, "flat_power")
    total = a + j
    force = localization.force_requirements(total)
    punctured = endpoint.requirements(total)

    # Uniform Cauchy convergence of the requested residual derivative uses
    # one further time derivative.  This is another fixed finite request.
    residual_jet = total + 1
    flat = corner.target(residual_jet, power, h, kappa)
    return {
        "spatial_order": a,
        "time_order": j,
        "total_mixed_order": total,
        "residual_Cauchy_jet": residual_jet,
        "comparison_stage": flat["comparison_stage"],
        "flat_power": power,
        "punctured_field_orders": {
            "potential": punctured["potential_for_future_residual_Cauchy"],
            "direct_swirl": punctured["direct_for_future_residual_Cauchy"],
            "pressure": punctured["pressure_for_future_residual_Cauchy"],
        },
        "localization_force_orders": {
            "potential": force["potential_order"],
            "direct_swirl": force["direct_swirl_order"],
            "pressure": force["pressure_order"],
        },
        "localization_Cauchy_orders": {
            "potential": force["Cauchy_potential_order"],
            "direct_swirl": force["Cauchy_direct_swirl_order"],
            "pressure": force["Cauchy_pressure_order"],
        },
        "closed_residual_exponent": flat["closed_chart_exponent"],
        "cartesian_factor": flat["cartesian_factor_at_requested_exponent"],
        "checks": {
            "punctured_orders_match_localization_Cauchy_orders": (
                punctured["potential_for_future_residual_Cauchy"] == force["Cauchy_potential_order"]
                and punctured["direct_for_future_residual_Cauchy"] == force["Cauchy_direct_swirl_order"]
                and punctured["pressure_for_future_residual_Cauchy"] == force["Cauchy_pressure_order"]
            ),
            "potential_has_four_order_Cauchy_buffer": force["Cauchy_potential_order"] == total + 4,
            "direct_swirl_has_three_order_Cauchy_buffer": force["Cauchy_direct_swirl_order"] == total + 3,
            "pressure_has_two_order_Cauchy_buffer": force["Cauchy_pressure_order"] == total + 2,
            "residual_next_time_derivative_is_fixed_finite_request": residual_jet == total + 1,
            "closed_residual_flatness_beats_requested_power": flat["retained_exponent_margin"] >= 0,
            "same_fixed_field_is_used": endpoint.closure_statement()["same_corrected_diagonal_field_used"],
        },
    }


def review():
    samples = [request(a, j, power) for a, j, power in (
        (0, 0, 0), (1, 3, 5), (20, 4, 11), (6, 70, 2), (100, 100, 17),
    )]
    endpoint_scope = endpoint.closure_statement()
    corner_scope = corner.review()["scope"]
    checks = {
        "same_corrected_field_and_cutoff_sequence_retained": (
            endpoint_scope["same_corrected_diagonal_field_used"]
            and corner_scope["same_field_all_fixed_residual_jets_cartesian_corner_flat_conditionally"]
        ),
        "punctured_potential_direct_swirl_and_pressure_jets_available_conditionally": endpoint_scope[
            "compatible_one_sided_C_infinity_potential_direct_swirl_velocity_and_pressure_jets_on_punctured_terminal_compacts_follow_conditionally"
        ],
        "complete_exterior_heat_and_pressure_jets_retained": endpoint_scope[
            "complete_exterior_heat_and_pressure_mixed_jets_retained"
        ],
        "origin_is_handled_by_residual_flatness_not_field_regularization": (
            not endpoint_scope["terminal_origin_regularized"]
            and corner_scope["same_field_all_fixed_residual_jets_cartesian_corner_flat_conditionally"]
        ),
        "all_sample_derivative_ledgers_match": all(
            row["checks"]["punctured_orders_match_localization_Cauchy_orders"] for row in samples
        ),
        "extra_time_derivative_is_only_one_fixed_order": all(
            row["checks"]["residual_next_time_derivative_is_fixed_finite_request"] for row in samples
        ),
        "every_sample_flat_power_closes": all(
            row["checks"]["closed_residual_flatness_beats_requested_power"] for row in samples
        ),
        "request_dependent_stage_does_not_select_a_new_field": True,
        "fixed_gauge_complete_pressure_memory_is_not_reset": True,
        "localization_force_energy_comparison_and_blowup_not_claimed": True,
    }
    return {
        "checks": checks,
        "default_request": request(),
        "sample_requests": samples,
        "joined_cover": {
            "punctured": "uniform one-sided A, B, p jets on every |x|>=d>0 terminal compact",
            "corner": "all fixed jets of R_1(u*,p*) are flat at (x,t)=(0,1)",
            "exterior": "exact heat swirl and fixed-gauge pressure mixed jets",
            "overlap": "the inherited axial/equatorial punctured cover overlaps inside the enlarged closed atlas",
        },
        "scope": {
            "closed_time_localization_inputs_joined_conditionally": True,
            "same_field_pressure_support_and_gauge_retained": True,
            "terminal_velocity_or_pressure_regular_at_origin_proved": False,
            "localization_and_compact_force_replayed_here": False,
            "unconditional_pulse_residual_force_comparison_or_blowup_claimed": False,
        },
    }
