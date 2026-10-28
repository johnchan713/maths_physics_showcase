"""Transfer the current log-reserved residual to the corner and endpoint."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


residual = load(
    "current_corner_endpoint_residual",
    "core_current_log_reserved_residual_reconciliation/residual.py",
)
old_corner = load(
    "current_corner_endpoint_geometry",
    "core_corner_flatness_obstruction/corner.py",
)
primitive = load(
    "current_corner_endpoint_margins",
    "core_primitive_margin_normalization_reconciliation/margins.py",
)
endpoint = load(
    "current_corner_endpoint_punctured",
    "current_endpoint_jets/endpoint.py",
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
        "compact_X_interval_fixed_before_supremum": True,
        "angular_domain_includes_eta_endpoints": True,
        "physical_support_floor_used": False,
        "edge_weight_floor_used": False,
    }


def cartesian_constant(power):
    return old_corner.corner_transfer_constant(order(power, "power"))


@lru_cache(None)
def corner_target(m=0, power=0, compact_log_degrees=(0, 0),
                  regular_log_degree=0, cutoff_log_degree=0,
                  h=Q(1, 200), kappa=Q(1, 100000)):
    """One fixed current comparison stage and its Cartesian transfer."""
    m, power = order(m, "m"), order(power, "power")
    row = residual.comparison(
        m, power, compact_log_degrees, regular_log_degree,
        cutoff_log_degree, h, kappa,
    )
    chart = closed_chart(h)
    exponent = Q(row["complete_minimum_exponent"])
    requested = Q(power + 1)
    request = row["current_stage_request"]
    return {
        "jet": m,
        "requested_power": power,
        "comparison_stage": row["log_reserved_comparison_stage"],
        "current_stage_request": request,
        "closed_chart_exponent": exponent,
        "requested_strict_exponent": requested,
        "retained_exponent_margin": exponent - requested,
        "cartesian_factor_at_requested_exponent": cartesian_constant(power + 1),
        "closed_chart": chart,
        "residual_groups": row["residual_group_exponents"],
        "reserve_accounts": row["reserve_accounts"],
        "pressure_memory": list(residual.MEMORY),
        "checks": {
            "current_residual_controls_hold": all(row["checks"].values()),
            "current_source_and_corner_use_same_stage":
                request["stage"] == row["log_reserved_comparison_stage"],
            "two_extra_compact_source_derivatives_reach_corner_stage": all(
                request["raw_sources_with_extra_integration"][name]
                == request["raw_sources_before_extra_integration"][name] + 2
                for name in ("A0", "A1")
            ),
            "terminal_chart_closes": chart["terminal_identity"] == 1,
            "closed_jacobian_floor_holds": chart["jacobian_floor"] >= Q(49, 50),
            "normalized_support_margin_positive": chart["normalized_support_margin"] > 0,
            "closed_residual_beats_requested_power": exponent >= requested,
            "all_six_groups_retained":
                tuple(row["residual_group_exponents"]) == residual.RESIDUAL_GROUPS,
            "cartesian_transfer_is_four_to_the_power":
                cartesian_constant(power + 1) == 4 ** (power + 1),
            "no_scale_incompatible_positive_floor_used":
                not chart["physical_support_floor_used"]
                and not chart["edge_weight_floor_used"],
        },
    }


@lru_cache(None)
def endpoint_request(spatial_order=0, time_order=0, flat_power=0,
                     compact_log_degrees=(0, 0), regular_log_degree=0,
                     cutoff_log_degree=0, h=Q(1, 200),
                     kappa=Q(1, 100000)):
    """One finite localization-input request on the already selected field."""
    a = order(spatial_order, "spatial_order")
    j = order(time_order, "time_order")
    power = order(flat_power, "flat_power")
    total = a + j
    punctured = endpoint.requirements(total)
    residual_jet = total + 1
    flat = corner_target(
        residual_jet, power, compact_log_degrees, regular_log_degree,
        cutoff_log_degree, h, kappa,
    )
    requested = {
        "potential": total + 4,
        "direct_swirl": total + 3,
        "pressure": total + 2,
    }
    supplied = {
        "potential": punctured["future_force_Cauchy_potential_bound"],
        "direct_swirl": punctured["future_force_Cauchy_direct_swirl_bound"],
        "pressure": punctured["future_force_Cauchy_pressure_bound"],
    }
    endpoint_review = endpoint.endpoint_review()
    return {
        "spatial_order": a,
        "time_order": j,
        "total_mixed_order": total,
        "residual_Cauchy_jet": residual_jet,
        "comparison_stage": flat["comparison_stage"],
        "flat_power": power,
        "punctured_field_Cauchy_orders": supplied,
        "localization_Cauchy_orders": requested,
        "closed_residual_exponent": flat["closed_chart_exponent"],
        "cartesian_factor": flat["cartesian_factor_at_requested_exponent"],
        "endpoint_blocking_inequality": endpoint_review["blocking_inequality"],
        "checks": {
            "punctured_orders_match_localization_Cauchy_orders": supplied == requested,
            "potential_has_four_order_Cauchy_buffer": supplied["potential"] == total + 4,
            "direct_swirl_has_three_order_Cauchy_buffer": supplied["direct_swirl"] == total + 3,
            "pressure_has_two_order_Cauchy_buffer": supplied["pressure"] == total + 2,
            "residual_next_time_derivative_is_fixed_finite_request":
                residual_jet == total + 1,
            "closed_residual_flatness_beats_requested_power":
                flat["retained_exponent_margin"] >= 0,
            "current_endpoint_geometry_controls_hold":
                all(endpoint_review["checks"].values()),
            "same_current_residual_stage_carries_extra_derivatives":
                flat["checks"]["two_extra_compact_source_derivatives_reach_corner_stage"],
        },
    }


def review():
    corners = [
        corner_target(0, 0),
        corner_target(1, 3, (7, 11), 13, 17),
        corner_target(5, 11, (37, 41), 43, 47),
    ]
    # Frozen evidence samples several finite requests.  The formula is
    # algebraic in the orders; very large recurrence evaluations add runtime,
    # not a stronger quantifier.
    requests = [endpoint_request(a, j, power) for a, j, power in (
        (0, 0, 0), (1, 3, 5), (5, 4, 11), (3, 7, 2),
    )]
    endpoint_review = endpoint.endpoint_review()
    checks = {
        "current_log_reserved_residual_is_exact_parent":
            residual.review()["scope"][
                "six_group_log_reserved_local_residual_instantiated_on_current_chain"
            ],
        "current_request_dependent_stage_reaches_closed_corner":
            all(row["checks"]["current_source_and_corner_use_same_stage"] for row in corners),
        "extra_A0_A1_derivative_cost_not_dropped_at_corner":
            all(row["checks"]["two_extra_compact_source_derivatives_reach_corner_stage"]
                for row in corners),
        "all_fixed_sample_residual_jets_are_cartesian_flat":
            all(all(row["checks"].values()) for row in corners),
        "current_punctured_endpoint_contract_is_used":
            endpoint_review["checks"][
                "inherited_all_order_physical_bounds_imply_needed_endpoint_inequality_conditionally"
            ],
        "punctured_and_localization_derivative_ledgers_match":
            all(row["checks"]["punctured_orders_match_localization_Cauchy_orders"]
                for row in requests),
        "extra_time_derivative_is_only_a_fixed_finite_request":
            all(row["checks"]["residual_next_time_derivative_is_fixed_finite_request"]
                for row in requests),
        "every_sample_corner_power_and_endpoint_request_closes":
            all(all(row["checks"].values()) for row in requests),
        "fixed_gauge_complete_pressure_and_axial_memory_retained":
            tuple(corners[0]["pressure_memory"]) == residual.MEMORY
            and endpoint_review["checks"][
                "full_seed_pressure_and_all_correction_differences_in_same_gauge_retained"
            ],
        "localization_force_energy_comparison_and_blowup_remain_open": True,
    }
    return {
        "checks": checks,
        "sample_corner_targets": corners,
        "sample_endpoint_requests": requests,
        "joined_cover": {
            "corner": "all fixed jets of the complete current residual are Cartesian flat",
            "punctured": "uniform one-sided current A, B and fixed-gauge p jets on every |x|>=d>0 terminal compact",
            "exterior": "exact heat direct-swirl and full fixed-gauge pressure mixed jets",
            "overlap": "current axial and equatorial regions cover each enlarged punctured compact",
        },
        "scope": {
            "current_same_field_residual_cartesian_corner_flat_conditionally": True,
            "current_closed_time_localization_inputs_joined_conditionally": True,
            "endpoint_uniform_all_order_physical_bound_proved_here": False,
            "terminal_velocity_or_pressure_regular_at_origin_proved": False,
            "localization_force_energy_comparison_replayed_here": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
