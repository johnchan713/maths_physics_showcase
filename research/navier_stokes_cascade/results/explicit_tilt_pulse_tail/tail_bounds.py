"""Exact same-cone tilt selection and conservative dyadic tail arithmetic."""

from fractions import Fraction as F


def rat(value):
    if type(value) is not F:
        raise TypeError("exact Fraction input required")
    return value


def ceil_fraction(value):
    value = rat(value)
    return -((-value.numerator) // value.denominator)


def floor_fraction(value):
    value = rat(value)
    return value.numerator // value.denominator


def ceil_log2(value):
    """Least nonnegative m with value <= 2**m, using integers only."""
    value = rat(value)
    if value <= 0:
        raise ValueError("positive logarithm argument required")
    if value <= 1:
        return 0
    m = max(0, value.numerator.bit_length() - value.denominator.bit_length())
    while F(2**m) < value:
        m += 1
    while m and F(2 ** (m - 1)) >= value:
        m -= 1
    return m


def ceil_root(value, degree):
    """Least nonnegative integer n with n**degree >= value."""
    value = rat(value)
    if value < 0 or type(degree) is not int or degree <= 0:
        raise ValueError("nonnegative value and positive integer degree required")
    if value == 0:
        return 0
    lo, hi = 0, 1
    while F(hi**degree) < value:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if F(mid**degree) < value:
            lo = mid
        else:
            hi = mid
    return hi


def tilt_certificate(epsilon):
    """Select u^2=2/epsilon under Ract<=1-epsilon."""
    epsilon = rat(epsilon)
    if not 0 < epsilon <= 1:
        raise ValueError("need 0 < epsilon <= 1")
    u2 = 2 / epsilon
    kappa = u2 / (1 + u2)
    actual_floor = kappa - (1 - epsilon)
    d0 = epsilon / 2
    return {
        "epsilon": epsilon,
        "u_squared": u2,
        "kappa": kappa,
        "actual_deficit_floor": actual_floor,
        "chosen_d0": d0,
        "unused_slack": actual_floor - d0,
        "coarse_ratio_upper": (1 + 1 / u2) * (1 - epsilon),
    }


def tail_certificate(h, coefficient_cost, rb, omega, lpsi, cmesh, cpert,
                     dbox, dchart, epsilon):
    """A sufficient integer tail for the selected-cutoff mixed scale gate."""
    values = tuple(map(rat, (
        h, coefficient_cost, rb, omega, lpsi, cmesh, cpert,
        dbox, dchart, epsilon,
    )))
    h, coefficient_cost, rb, omega, lpsi, cmesh, cpert, dbox, dchart, epsilon = values
    if not 0 < h < 1 or min(coefficient_cost, rb, omega, lpsi, cpert,
                            dbox, dchart) <= 0 or cmesh < 0:
        raise ValueError("invalid positive gate data")
    tilt = tilt_certificate(epsilon)
    d0, kappa = tilt["chosen_d0"], tilt["kappa"]
    gamma = omega * omega * d0 / (9 * kappa)
    budget = min(
        rb / 8,
        omega * omega * d0 / (16 * lpsi),
        gamma / (128 * cpert),
    )
    m_required = ceil_log2(coefficient_cost / budget)
    ell_exponential = 2 + ceil_fraction(F(m_required, 1) / (2 * h))
    sixth_power_required = max(F(1), 16 * cmesh / gamma, dbox / dchart)
    ell_mesh = ceil_root(sixth_power_required, 6)
    ell = max(2, ell_exponential, ell_mesh)
    exponent = floor_fraction(2 * h * (ell - 2))
    scale_floor = F(2**exponent)
    error = coefficient_cost / scale_floor
    H = gamma / 8 - cmesh / F(ell**6)
    mixed_budget = scale_floor * min(
        rb / 8,
        omega * omega * d0 / (16 * lpsi),
        H / (8 * cpert),
    )
    return {
        "tilt": tilt,
        "gamma": gamma,
        "pointwise_budget": budget,
        "binary_exponent_required": m_required,
        "ell_exponential": ell_exponential,
        "sixth_power_required": sixth_power_required,
        "ell_mesh": ell_mesh,
        "ell_tail": ell,
        "scale_floor": scale_floor,
        "error_upper": error,
        "H_lower": H,
        "mixed_budget": mixed_budget,
        "checks": {
            "raw_deficit_is_valid": d0 <= tilt["actual_deficit_floor"],
            "strict_unused_raw_margin": tilt["unused_slack"] > 0,
            "coarse_same_tilt_ratio_is_strict": tilt["coarse_ratio_upper"] < 1,
            "mesh_half_budget": cmesh / F(ell**6) <= gamma / 16,
            "positive_H": H >= gamma / 16 > 0,
            "base_error_gate": error <= rb / 8,
            "deficit_error_gate": error <= omega * omega * d0 / (16 * lpsi),
            "current_value_error_gate": error <= gamma / (128 * cpert),
            "selected_cutoff_mixed_gate": coefficient_cost <= mixed_budget,
            "positive_chart_gate": F(ell**6) >= dbox / dchart,
        },
    }


def physical_exponents(h, spatial_order, time_order):
    """The inherited uniform primary-pulse physical Q losses."""
    h = rat(h)
    if type(spatial_order) is not int or type(time_order) is not int \
            or min(spatial_order, time_order) < 0:
        raise ValueError("nonnegative integer derivative orders required")
    sx = F(1, 2) + h / 2
    st = 1 + 3 * h / 2
    return {
        "sx": sx,
        "st": st,
        "velocity": (spatial_order + 1) * sx + time_order * st,
        "pressure": (spatial_order + 2) * sx + time_order * st,
        "potential": spatial_order * sx + time_order * st,
    }
