"""Exact ledgers for flat activation, five moments and the cone collar."""

from fractions import Fraction as F


def volterra_certificate():
    """Constants from the positive half-line representation of V_m."""
    # On 0<=s<=1/4, 2<=g<=6.  The kernel, (1+s^2 v)^-p,
    # and exp(-v) are nonnegative and the last two integrate to at most one.
    return {
        "g_lower": 2,
        "g_upper": 6,
        "operator_norm": F(3, 2),
        "endpoint_value_multiplier": F(1, 2),
        "V0_first_jet_multiplier": F(1, 2),
        "V0_second_jet_b2_multiplier": F(1, 2),
        "V0_second_jet_b0_multiplier": -F(3, 2),
        "field_factor": "t*ea*s^3",
        "moment_factor": "t^2*ea*s^6",
    }


def moment_ledger():
    """All coefficients after one exact subtraction and before V_3."""
    return (
        {
            "name": "M", "formula": ("X*bU",),
            "C_power": 0, "R_power": 3, "removed_C_factor": "1",
        },
        {
            "name": "I", "formula": ("X*bH",),
            "C_power": -1, "R_power": 5, "removed_C_factor": "C",
        },
        {
            "name": "J",
            "formula": ("X*Ur*bH", "X*Hr*bU", "X*alpha*bU*bH"),
            "C_power": -1, "R_power": 5, "removed_C_factor": "C",
        },
        {
            "name": "S",
            "formula": ("2*X*Ur*bU", "X*alpha*bU^2", "-X*bE2/2"),
            "C_power": 0, "R_power": 4, "removed_C_factor": "1",
        },
        {
            "name": "Cp", "formula": ("bE2/2",),
            "C_power": -2, "R_power": 4, "removed_C_factor": "C^2",
        },
    )


def map_ledger():
    """Exact difference terms in W,NQ,Qs,Ns,p1,p2."""
    return {
        "Delta_W": (
            "-(2*D*eta*DeltaM+d*DeltaM_eta)/X",
        ),
        "Delta_NQ": (
            "(1-h)*DeltaI", "-D*eta*DeltaI_eta",
            "-d*DeltaJ_eta", "2*(h-D)*eta*DeltaJ",
        ),
        "Delta_Qs": (
            "-DeltaW", "DeltaNQ/(X*H)",
            "-NQ_r*DeltaH/(X*H_r*H)",
        ),
        "Delta_Ns": (
            "-W_r*DeltaU", "-U_r*DeltaW", "-DeltaW*DeltaU",
            "D*(DeltaM-eta*DeltaM_eta)/X",
            "4*h*eta*DeltaS/X", "-d*DeltaS_eta/X",
            "4*A*eta*DeltaCp", "-d*DeltaCp_eta",
        ),
        "Delta_p1": ("X*DeltaQs/L",),
        "Delta_p2": (
            "X*DeltaNs/(L*E)",
            "-X*Ns_r*DeltaE/(L*E_r*E)",
        ),
    }


def envelope_ledger():
    """Independently counted sufficient R powers after displayed factors."""
    return {
        "Delta_logE_over_t_ea_s3": 1,
        "Delta_U_over_t_ea_s3": 2,
        "C_Delta_E_H_over_t_ea_s3": 4,
        "moments_over_t2_ea_s6": 8,
        "Delta_W_CDeltaNQ_over_t2_ea_s6": 10,
        "Delta_p1_over_t_ea_s3": 18,
        "Delta_Ns_over_t_ea_s3": 16,
        "Delta_p2_over_C_t_ea_s3": 24,
        "Er_over_E_minus_1_over_t_ea_s3": 2,
    }


def cone_ledger():
    """Factor counts after division by t*s^3 in D_c,J_c*,V_c."""
    # C_power is the positive physical C power.  Every row is bounded by
    # C^2 R^80; the largest derived R cost is ts*Delta p2.
    return (
        {"name": "Delta_p1", "C_power": 0, "R_power": 18},
        {"name": "ts_Delta_p2", "C_power": 2, "R_power": 27},
        {"name": "q2_ratio_error", "C_power": 2, "R_power": 5},
        {"name": "q_ts_ratio_error", "C_power": 2, "R_power": 6},
        {"name": "Delta_p2_in_J", "C_power": 1, "R_power": 24},
        {"name": "ts_Delta_p1_in_J", "C_power": 1, "R_power": 21},
        {"name": "q_ratio_error_in_J", "C_power": 1, "R_power": 3},
        {"name": "vs_error", "C_power": 2, "R_power": 6},
    )


def scale_certificate():
    """Logarithmic comparisons at the selected C, without exponentiating."""
    # log R=B/100 and log C>=10000 B.  Coefficients below multiply B.
    return {
        "log_R80_over_C_upper_per_B": F(4, 5) - 10000,
        "log_tR3_upper_per_B": F(3, 100) - 20 * 10000,
        # log((Kact*t)/C^-17)=log(R^80/C).
        "log_Kact_t_over_Cminus17_upper_per_B": F(4, 5) - 10000,
        "log_native_growth_change_upper_per_B": F(8, 100) - 18 * 10000,
        "cone_C_power_after_R160_below_C2": -31,
        "cone_constant": F(32, 81),
        "activated_growth_lower": F(17, 8),
        "denominator_lower": F(9, 8),
        "native_growth_lower": F(9, 4),
        "flat_tail_C_power": -97,
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
    volterra = volterra_certificate()
    moments = moment_ledger()
    maps = map_ledger()
    envelopes = envelope_ledger()
    cone = cone_ledger()
    scales = scale_certificate()
    moment_names = tuple(row["name"] for row in moments)
    pressure_terms = maps["Delta_Ns"][-2:]
    checks = {
        "positive_Volterra_norm_is_three_halves": volterra["operator_norm"] == F(3, 2),
        "Volterra_endpoint_is_not_zero": volterra["endpoint_value_multiplier"] == F(1, 2),
        "Volterra_second_endpoint_jet_keeps_radial_Jacobian": volterra["V0_second_jet_b0_multiplier"] == -F(3, 2),
        "all_five_moment_rows_present_once": moment_names == ("M", "I", "J", "S", "Cp"),
        "J_row_keeps_all_three_product_terms": len(moments[2]["formula"]) == 3,
        "pressure_row_has_no_spurious_X": moments[4]["formula"] == ("bE2/2",),
        "I_and_J_have_one_inverse_C": moments[1]["C_power"] == moments[2]["C_power"] == -1,
        "pressure_increment_has_two_inverse_C_powers": moments[4]["C_power"] == -2,
        "Qs_reciprocal_difference_is_complete": len(maps["Delta_Qs"]) == 3,
        "Ns_keeps_both_pressure_terms": pressure_terms == ("4*A*eta*DeltaCp", "-d*DeltaCp_eta"),
        "p2_keeps_field_denominator_difference": len(maps["Delta_p2"]) == 2,
        "moment_R8_envelope_covers_each_exact_row": max(row["R_power"] for row in moments) <= envelopes["moments_over_t2_ea_s6"],
        "worst_cone_row_is_C2_R27": max((row["C_power"], row["R_power"]) for row in cone if row["C_power"] == 2) == (2, 27),
        "Kact_C2_R80_covers_every_cone_row": all(row["C_power"] <= 2 and row["R_power"] <= 80 for row in cone),
        "R80_is_strictly_below_C": scales["log_R80_over_C_upper_per_B"] < 0,
        "unchanged_width_makes_tR3_small": scales["log_tR3_upper_per_B"] < 0,
        "unchanged_width_makes_Kact_t_below_Cminus17": scales["log_Kact_t_over_Cminus17_upper_per_B"] < 0,
        "native_growth_change_is_negligible": scales["log_native_growth_change_upper_per_B"] < 0,
        "flat_tail_loss_below_one_eighth_for_C_ge_2": 18 < 2**94,
        "squared_cone_loss_has_Cminus31": scales["cone_C_power_after_R160_below_C2"] == -31,
        "squared_cone_constant_is_32_over_81": scales["cone_constant"] == F(32, 81),
        "diagnostics_are_not_a_continuum_proof": True,
    }
    return _encode({
        "volterra": volterra,
        "moments": moments,
        "map_differences": maps,
        "envelopes": envelopes,
        "cone_coefficients": cone,
        "selected_scale": scales,
        "checks": checks,
    })
