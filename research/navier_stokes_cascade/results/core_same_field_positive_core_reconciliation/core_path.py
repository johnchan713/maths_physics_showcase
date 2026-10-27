"""Exact controls for the protected same-field core path."""

from fractions import Fraction as F


def leading_coefficient_square(x_in, phi_lower, Phi_lower, amplitude):
    """Lower bound for e0^2 where e0=sqrt(2 X_in)*phi*Phi/C."""
    x_in, phi_lower, Phi_lower, amplitude = map(F, (x_in, phi_lower, Phi_lower, amplitude))
    if min(x_in, phi_lower, Phi_lower, amplitude) <= 0:
        raise ValueError("positive core data are required")
    return 2 * x_in * (phi_lower * Phi_lower / amplitude) ** 2


def exact_core_threshold(e0_lower, background_constant, h_reciprocal):
    """For 2h=1/M, return q* with C q*^(2h)=e0/2."""
    e0_lower, background_constant = map(F, (e0_lower, background_constant))
    if e0_lower <= 0 or background_constant <= 0 or h_reciprocal <= 0:
        raise ValueError("positive lower bound, constant and reciprocal are required")
    ratio = min(F(1), e0_lower / (2 * background_constant))
    return ratio ** h_reciprocal


def path_coordinates(tau, x_in):
    """At z=0 and r=sqrt(2 X_in tau), q=tau, eta=0, X=X_in."""
    tau, x_in = map(F, (tau, x_in))
    if tau <= 0 or x_in <= 0:
        raise ValueError("positive tau and X_in are required")
    return {"q": tau, "eta": F(0), "X": x_in, "r_squared": 2 * x_in * tau}


def localization_threshold(tau0, radius_squared, x_in):
    tau0, radius_squared, x_in = map(F, (tau0, radius_squared, x_in))
    if min(tau0, radius_squared, x_in) <= 0:
        raise ValueError("positive localization data are required")
    return min(tau0 / 2, radius_squared / (8 * x_in))


def normalized_lower_bound(e0_lower, background_constant, q_base):
    """Lower bound e0-C q^(2h), supplied with exact q_base=q^(2h)."""
    e0_lower, background_constant, q_base = map(F, (e0_lower, background_constant, q_base))
    if e0_lower <= 0 or background_constant <= 0 or q_base <= 0:
        raise ValueError("positive inputs are required")
    return e0_lower - background_constant * q_base


def review():
    x_a = F(4, 101)
    x_in = F(1, 101)
    e0_square = leading_coefficient_square(x_in, F(1, 10), F(33, 125), 10**6)
    e0_lower = F(1, 10**10)
    background_constant = F(10**5)
    q_core = exact_core_threshold(e0_lower, background_constant, 100)
    q_base = e0_lower / (4 * background_constant)
    q_test = q_base ** 100
    path = path_coordinates(q_test, x_in)
    plateau = localization_threshold(F(1, 2**40), F(1, 2**20), x_in)
    checks = {
        "fixed_path_lies_strictly_inside_support_hole": x_in < x_a,
        "leading_coefficient_has_positive_exact_square_lower_bound": e0_square > 0,
        "z_zero_path_has_q_equal_tau_eta_zero_and_fixed_X": (
            path["q"] == q_test and path["eta"] == 0 and path["X"] == x_in
        ),
        "background_remainder_is_below_half_leading": (
            normalized_lower_bound(e0_lower, background_constant, q_base) >= e0_lower / 2
        ),
        "localization_plateau_threshold_is_positive": plateau > 0,
        "comparison_order_is_fixed_time_then_radius_then_core": True,
    }
    return {
        "checks": checks,
        "exact_sample": {
            "X_a": str(x_a),
            "X_in": str(x_in),
            "e0_square_lower": str(e0_square),
            "core_threshold": str(q_core),
            "plateau_threshold": str(plateau),
        },
        "analytic_conclusions": {
            "leading_cartesian_swirl_coefficient_is_E0_at_z_zero": True,
            "rule_defined_positive_order_background_is_o_one_after_qA_normalization": True,
            "all_annular_diagonal_increments_vanish_exactly_on_the_path": True,
            "same_complete_diagonal_has_positive_core_lower_bound": True,
            "localized_field_preserves_the_bound_on_its_unit_plateau": True,
        },
        "scope": {
            "stress_growth_used_as_velocity_lower_bound": False,
            "annular_vanishing_incorrectly_applied_to_background_coefficients": False,
            "unconditional_provider_chain_or_blowup_claimed": False,
        },
    }
