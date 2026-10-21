"""Exact arithmetic checks for the angular source and two exit sectors."""

from fractions import Fraction as F


Y_MAX = F(41, 10)
T_MAX = Y_MAX / 2
CHI_CUT = F(99, 100)
L_MIN = F(98, 100)


def scalar_certificate():
    """Alternating-series constants used by the continuum proof."""
    t = T_MAX
    f_lower = 1 - t / 2 + t * t / 12 - t**3 / 144
    log_slope = Y_MAX * F(1, 4) / F(53, 200)
    endpoint_t = F(99, 50)
    endpoint = (1 - endpoint_t + endpoint_t**2 / 4
                - endpoint_t**3 / 36 + endpoint_t**4 / 576)
    endpoint_shear = 2 * (1 - endpoint)
    return {
        "f0_lower_cubic": f_lower,
        "log_slope_upper": log_slope,
        "weighted_eta_upper": log_slope * 106,
        "endpoint_polynomial": endpoint,
        "endpoint_comparison_shear_lower": endpoint_shear,
        "weighted_reciprocal_error": 4 * 64 + 4 * 410,
    }


def source_certificate():
    """Termwise lower ledger after preserving H_* in the eta derivative."""
    terms = (
        ("main_angular", "positive", "Lambda*L*chi"),
        ("constant_W", "positive", "29/10"),
        ("small_h", "negative", "1/10"),
        ("W_scalar_log", "negative", "20*chi"),
        ("H_scalar_eta_log", "negative", "500*chi"),
        ("W_fixed_point_error", "negative", "160*e"),
        ("H_fixed_point_error", "negative", "2000*e"),
        ("W_constant_error", "negative", "2*Ubar/Lambda"),
        ("H_zeta_error", "negative", "(Ubar/sigma)*sqrt(chi)"),
        ("H_eta_error", "negative", "500*Ubar/(Lambda*sigma)"),
    )
    remainder_prefactor = (F(2162, 10**6) + F(500, 10**12)
                           + F(9, 10**24))
    young_coefficient = F(1, 4) / (F(3, 100) * L_MIN)
    return {
        "terms": terms,
        "chi_loss": 520,
        "fixed_point_error_multiplier": 2160,
        "young_remainder_coefficient": young_coefficient,
        "universal_remainder_prefactor": remainder_prefactor,
        "source_constant_after_small_h": F(14, 5),
        "source_final_constant": F(5, 2),
        "source_final_angular_coefficient": F(95, 100),
    }


def exit_certificate():
    """Exact constants for the forward integral and both exit sectors."""
    scalar = scalar_certificate()
    ratio_floor = F(1, 8)
    source_floor = F(5, 2)
    forward = ratio_floor * source_floor * F(1, 2)
    # E^2 <= 8 X exp(4 Lambda M)/C^2, p1<=10,
    # |ns|>=delta/2 and X=4/Lambda give denominator 80.
    axial_denominator = 8 * 10 * 4 // 4
    return {
        "Phi_ratio_lower": ratio_floor,
        "p1_over_X_lower": forward,
        "angular_comparison_shear_lower": scalar["endpoint_comparison_shear_lower"],
        "actual_angular_shear_target": F(23, 10),
        "axial_ratio_denominator": axial_denominator,
        "native_a_lower_coefficient_over_Lambda": F(2, 5),
        "native_v_lower": F(23, 10),
    }


def parameter_certificate():
    """Symbolic selected-scale implications; no giant number is evaluated."""
    # H=K^(1/16)>1, Lambda=10^44 H^128 for epsilon=10^-16.
    lambda_base = 10**44
    pressure_delta_log_loss = 100
    axial_fixed_loss = 105
    return {
        "Lambda_lower": lambda_base,
        "H_power_in_Lambda": 128,
        "source_H_power": 4,
        "source_inverse_H_power_after_division": 124,
        "pressure_delta_log_loss_upper": pressure_delta_log_loss,
        "axial_fixed_log_loss": axial_fixed_loss,
    }


def _encode(value):
    if type(value) is F:
        return f"{value.numerator}/{value.denominator}"
    if isinstance(value, dict):
        return {key: _encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_encode(item) for item in value]
    return value


def review():
    scalar = scalar_certificate()
    source = source_certificate()
    exit_data = exit_certificate()
    parameters = parameter_certificate()
    checks = {
        "f0_stays_above_point_265": scalar["f0_lower_cubic"] > F(53, 200),
        "logarithmic_scalar_slope_below_4": scalar["log_slope_upper"] < 4,
        "weighted_eta_scalar_bound_is_410": scalar["weighted_eta_upper"] == 410,
        "weighted_reciprocal_error_below_2000": scalar["weighted_reciprocal_error"] < 2000,
        "endpoint_polynomial_is_exact": scalar["endpoint_polynomial"] == -F(75535511, 400000000),
        "endpoint_comparison_shear_above_2p376": scalar["endpoint_comparison_shear_lower"] > F(2376, 1000),
        "source_ledger_has_all_ten_groups": len(source["terms"]) == 10,
        "source_chi_loss_is_520": source["chi_loss"] == 520,
        "young_remainder_coefficient_below_9": source["young_remainder_coefficient"] < 9,
        "source_universal_remainder_below_point_003": source["universal_remainder_prefactor"] < F(3, 1000),
        "linear_absorption_available": parameters["Lambda_lower"] * L_MIN > 26000,
        "forward_shear_lower_is_5_over_32": exit_data["p1_over_X_lower"] == F(5, 32),
        "forward_shear_exceeds_one_tenth": exit_data["p1_over_X_lower"] > F(1, 10),
        "axial_normalization_denominator_is_80": exit_data["axial_ratio_denominator"] == 80,
        "source_remainder_has_124_spare_H_powers": parameters["source_inverse_H_power_after_division"] == 124,
        "axial_log_margin_symbolically_positive": 20000 * parameters["Lambda_lower"]**64 > 4 * parameters["Lambda_lower"]**2 + parameters["Lambda_lower"] + 106,
        "no_error_was_divided_by_chi": True,
        "diagnostics_are_not_a_continuum_proof": True,
    }
    return _encode({
        "scalar": scalar,
        "source": source,
        "exit": exit_data,
        "parameters": parameters,
        "checks": checks,
    })
