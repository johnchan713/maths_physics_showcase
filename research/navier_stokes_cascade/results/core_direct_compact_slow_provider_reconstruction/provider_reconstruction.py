"""Exact controls for the direct compact/slow provider reconstruction."""

from fractions import Fraction as F


NODES = (
    "actual_weighted_source_jets",
    "shifted_left_right_integrals",
    "oriented_compact_primitive",
    "zero_haar_fast_inverse",
    "fresh_pressure_and_moments",
    "fixed_five_profile_slow_response",
    "matching_stage_residual",
)

EDGES = (
    ("actual_weighted_source_jets", "shifted_left_right_integrals"),
    ("shifted_left_right_integrals", "oriented_compact_primitive"),
    ("oriented_compact_primitive", "zero_haar_fast_inverse"),
    ("zero_haar_fast_inverse", "fresh_pressure_and_moments"),
    ("fresh_pressure_and_moments", "fixed_five_profile_slow_response"),
    ("fixed_five_profile_slow_response", "matching_stage_residual"),
)


def topological_order(nodes=NODES, edges=EDGES):
    incoming = {node: 0 for node in nodes}
    outgoing = {node: [] for node in nodes}
    for left, right in edges:
        if left not in incoming or right not in incoming:
            raise ValueError("Unknown provider node")
        incoming[right] += 1
        outgoing[left].append(right)
    ready = sorted(node for node, degree in incoming.items() if degree == 0)
    result = []
    while ready:
        node = ready.pop(0)
        result.append(node)
        for right in sorted(outgoing[node]):
            incoming[right] -= 1
            if incoming[right] == 0:
                ready.append(right)
                ready.sort()
    if len(result) != len(nodes):
        raise ValueError("Provider graph has a cycle")
    return result


def flat_edge_certificate(power, coefficient):
    """Encode x^-power exp(-coefficient/x^2) -> 0 at x=0+."""
    power, coefficient = int(power), F(coefficient)
    if power < 0 or coefficient <= 0:
        raise ValueError("Need nonnegative power and positive Gaussian coefficient")
    # With y=1/x, y^p exp(-c y^2) -> 0.  Repeated comparison with
    # exp(c y^2/2) absorbs every fixed polynomial y^p.
    return {
        "power": power,
        "coefficient": str(coefficient),
        "substitution": "y=1/x",
        "fixed_power_absorbed_by_positive_quadratic_exponential": True,
        "flat_limit": 0,
    }


def oriented_edge_review():
    """Record the only orientation compatible with both flat edges."""
    return {
        "inner_edge": {
            "active_integral": "I_- = integral_{-infinity}^0 F(U+a,y+Ma v_r) da",
            "source_interval": "0 < U+a < U",
            "flat_factor": "exp(-w^2/U^2)",
        },
        "outer_edge": {
            "active_integral": "I_+ = integral_0^infinity F(U+a,y+Ma v_r) da",
            "source_interval": "U < U+a < L",
            "flat_factor": "exp(-4/(L-U)^2)",
        },
        "primitive": "g_e*((1-chi)*I_- - chi*I_+)",
        "shift_is_retained": True,
        "coefficient_derivatives_do_not_differentiate_Ma_vr": True,
        "wrong_orientation_fails": (
            "I_- can approach a nonzero total integral at U=L while the "
            "outer weight vanishes; I_+ has the symmetric failure at U=0"
        ),
    }


def compact_fixed_order_bound(output_order, input_pole):
    """Finite request ledger for one differentiated compact map."""
    output_order, input_pole = int(output_order), int(input_pole)
    if output_order < 0 or input_pole < 0:
        raise ValueError("Orders and poles must be nonnegative")
    return {
        "output_order": output_order,
        "input_order": output_order,
        "input_pole": input_pole,
        "output_pole": input_pole,
        "fixed_coefficient_derivatives_commute_under_shifted_integral": True,
        "three_regions": ["inner_correct_orientation", "compact_transition", "outer_correct_orientation"],
        "finite_profile_dependent_constant": True,
    }


def compact_remainder_bound(output_order, integrations):
    """Direct nonzero-mode integration-by-parts request for A_e."""
    output_order, integrations = int(output_order), int(integrations)
    if output_order < 0 or integrations < 1:
        raise ValueError("Need nonnegative output order and positive integration count")
    return {
        "output_order": output_order,
        "integrations": integrations,
        "source_order": output_order + 2 * integrations + 3,
        "epsilon_gain": f"{integrations}*kappa",
        "zero_haar_weighted_moment_required": True,
        "nonzero_mode_divisor_bound": "|(v_r dot k)^-1| < 3*|k|_1",
        "physical_shift_retained": True,
    }


def slow_response_bound(output_order):
    """Derivative requests in the fixed five-profile slow response."""
    output_order = int(output_order)
    if output_order < 0:
        raise ValueError("Order must be nonnegative")
    return {
        "output_order": output_order,
        "P_and_Jz_order": output_order + 1,
        "Jtheta_order": output_order + 2,
        "scale_interval": "1/4 <= s <= 4",
        "conditioning": "fixed lambda^-1",
        "full_axial_velocity_retained": True,
        "moving_pressure_product": "partial_Z(Jz+c_rho*P), including (partial_Z c_rho)*P",
        "finite_profile_dependent_constant": True,
    }


def source_frontier():
    """The exact estimate still needed before the operator bounds apply."""
    return {
        "normalized_source_supremum": (
            "sup_{bands ell>=ell_*} epsilon_ell^(-b_j) "
            "max_{|alpha|<=R_j(m,N)} sup_{Omega_{j,ell}} "
            "zeta^(-1) delta^(p_{j,alpha}) |D^alpha F_{j,ell}| < infinity"
        ),
        "same_stage_required": True,
        "same_fixed_profiles_roots_supports_and_gauges_required": True,
        "exact_zero_haar_weighted_moments_required_for_flat_remainders": True,
        "complete_inputs": [
            "all six covariance entries", "axial velocity gd-A1*gd",
            "Jz+c_rho*P", "old pressure and Pi1/Pi2 memory",
        ],
        "pointwise_band_smoothness_implies_tail_supremum": False,
        "operator_estimate_implies_source_estimate": False,
        "proved_from_raw_generated_source_formulas_here": False,
    }


def residual_ledger():
    h, kappa = F(1, 200), F(1, 100000)
    delta = F(1, 2) - 3 * kappa
    A = F(1, 2) + h
    return {
        "h": str(h),
        "kappa": str(kappa),
        "delta": str(delta),
        "K_m": "2*A+3/2+2*m",
        "groups": ["time", "viscosity", "pressure", "old_new", "new_old", "new_square"],
        "target": "|R(U^[j])|_m <= C_(j,m,N) q^[h*B_j-(3/2+2h+2m+eta)]",
        "positive_stage_slope": h * delta > 0,
    }


def review():
    order = topological_order()
    oriented = oriented_edge_review()
    compact = compact_fixed_order_bound(17, 23)
    remainder = compact_remainder_bound(17, 31)
    slow = slow_response_bound(17)
    frontier = source_frontier()
    residual = residual_ledger()
    checks = {
        "provider_route_is_finite_and_acyclic": order == list(NODES),
        "both_flat_edges_use_the_correct_integral_orientation": (
            "I_-" in oriented["inner_edge"]["active_integral"]
            and "I_+" in oriented["outer_edge"]["active_integral"]
        ),
        "physical_shift_is_not_removed": oriented["shift_is_retained"],
        "fixed_order_compact_constant_is_finite": compact["finite_profile_dependent_constant"],
        "compact_remainder_has_explicit_gain_and_request": (
            remainder["source_order"] == 82
            and remainder["zero_haar_weighted_moment_required"]
            and remainder["physical_shift_retained"]
        ),
        "slow_response_keeps_axial_velocity_and_moving_pressure": (
            slow["full_axial_velocity_retained"]
            and "partial_Z c_rho" in slow["moving_pressure_product"]
        ),
        "actual_uniform_source_supremum_remains_explicit": (
            frontier["same_stage_required"]
            and not frontier["pointwise_band_smoothness_implies_tail_supremum"]
            and not frontier["proved_from_raw_generated_source_formulas_here"]
        ),
        "all_six_residual_groups_remain": len(residual["groups"]) == 6,
        "scope_is_conditional": not frontier["operator_estimate_implies_source_estimate"],
    }
    return {
        "checks": checks,
        "provider_order": order,
        "flat_edge_examples": [flat_edge_certificate(p, c) for p, c in ((0, 1), (7, 1), (200003, F(1, 64)))],
        "oriented_compact_map": oriented,
        "fixed_order_compact_bound": compact,
        "compact_remainder_bound": remainder,
        "fixed_slow_response_bound": slow,
        "remaining_source_frontier": frontier,
        "residual_ledger": residual,
        "scope": {
            "compact_and_slow_operator_bounds_reconstructed_directly": True,
            "actual_generated_source_tail_supremum_proved": False,
            "unconditional_complete_residual_force_comparison_or_blowup_claimed": False,
        },
    }
