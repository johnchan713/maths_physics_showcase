"""Exact arithmetic controls for the scale-compatible corner/force replay."""

from fractions import Fraction as F


def comparison_stage(h, delta, k_loss, l_loss, target):
    """Least nonnegative J satisfying both closed-corner exponent gates."""
    h, delta = F(h), F(delta)
    if h <= 0 or delta <= 0 or target < 0 or k_loss < 0 or l_loss < 0:
        raise ValueError("positive slopes and nonnegative losses are required")
    j = 0
    while True:
        b_j = F(3, 2) - F(2, 100000) + j * delta
        fixed = h * b_j - k_loss
        tail = h * delta * (j + 1) / 2 - l_loss
        if fixed >= target + 1 and tail >= target + 1:
            return {"stage": j, "fixed_exponent": fixed, "tail_exponent": tail}
        j += 1


def logarithmic_half_gain_bound(gain, power):
    """Exponent in sup (1+s)^P exp(-gain*s/2) <= exp(P^2/gain)."""
    gain = F(gain)
    if gain <= 0 or power < 0:
        raise ValueError("gain must be positive and power nonnegative")
    return F(power * power, 1) / gain


def cartesian_transfer_constant(power):
    """q <= 4 s turns C q^power into (4^power C) s^power."""
    if power < 0:
        raise ValueError("power must be nonnegative")
    return 4 ** power


def support_distances(scales):
    """Physical distances from the normalized 3/8012 separation."""
    margin = F(3, 8012)
    values = [margin * F(scale) for scale in scales]
    if not values or any(scale <= 0 for scale in scales):
        raise ValueError("positive scales are required")
    return values


def force_jet_requirements(order):
    """Fixed representative orders required by an order-m force jet."""
    if order < 0:
        raise ValueError("order must be nonnegative")
    return {
        "potential": order + 3,
        "direct_swirl": order + 2,
        "pressure": order + 1,
        "next_time_derivative_reserved": True,
    }


def review():
    h = F(1, 200)
    delta = F(1, 2) - F(3, 100000)
    stage = comparison_stage(h, delta, 13, 17, 9)
    previous = stage["stage"] - 1
    previous_fixed = h * (F(3, 2) - F(2, 100000) + previous * delta) - 13
    previous_tail = h * delta * (previous + 1) / 2 - 17
    distances = support_distances([F(1, 2**n) for n in (4, 20, 100)])
    checks = {
        "comparison_stage_meets_both_gates": (
            stage["fixed_exponent"] >= 10 and stage["tail_exponent"] >= 10
        ),
        "comparison_stage_is_minimal": previous_fixed < 10 or previous_tail < 10,
        "half_gain_keeps_finite_log_constant": logarithmic_half_gain_bound(F(7, 100), 5) == F(2500, 7),
        "cartesian_transfer_keeps_requested_power": cartesian_transfer_constant(9) == 4**9,
        "physical_separation_has_no_uniform_floor": distances[-1] < F(1, 10**20),
        "force_jet_reserves_complete_orders": force_jet_requirements(6) == {
            "potential": 9,
            "direct_swirl": 8,
            "pressure": 7,
            "next_time_derivative_reserved": True,
        },
    }
    return {
        "parameters": {"h": str(h), "delta": str(delta)},
        "sample_stage_gate": {key: str(value) for key, value in stage.items()},
        "checks": checks,
        "analytic_conclusions": {
            "scale_compatible_certificate_suffices_for_closed_chart_constants": True,
            "selected_half_gain_absorbs_finite_logs_and_constants": True,
            "all_fixed_residual_jets_are_cartesian_corner_flat_relative_to_the_certified_provider_chain": True,
            "same_field_compact_force_and_preterminal_energy_follow_relative_to_punctured_endpoint_bounds": True,
        },
        "scope": {
            "positive_physical_support_floor_used": False,
            "positive_edge_weight_floor_used": False,
            "same_field_positive_core_lower_bound_independently_proved_here": False,
            "unconditional_blowup_or_competitor_exclusion_proved": False,
        },
    }
