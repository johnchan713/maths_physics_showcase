"""Exact controls for the independently rederived core fixed-point ledger."""

from fractions import Fraction as F
from math import comb, factorial


A0 = 256
D0_SYMBOL = "80/rho"
PRESSURE_PRIMITIVE = 80 * A0**3


def weight(n, k, rho=F(1)):
    if type(n) is not int or type(k) is not int or n < 0 or k < 0 or rho <= 0:
        raise ValueError("nonnegative integer indices and positive rho required")
    return F(factorial(k) * comb(n + k, k),
             20**n * (n + 1)**2 * (k + 1)**2) / rho**k


def radial_shift_ratio(n, k):
    return weight(n, k) / weight(n + 1, k)


def eta_radial_shift_ratio(n, k, rho=F(1)):
    return weight(n, k + 1, rho) / weight(n + 1, k, rho)


def mixed_allocation(i, j, nu):
    if min(i, j) < 0 or nu not in (1, 2):
        raise ValueError("invalid allocation indices")
    return F((i + 1) * j, (i + j + 1) * (i + j + nu))


def convolution_sum(index):
    if type(index) is not int or index < 0:
        raise ValueError("nonnegative integer required")
    return sum(F((index + 1)**2, (i + 1)**2 * (index - i + 1)**2)
               for i in range(index + 1))


def cauchy_multiplier(k):
    if type(k) is not int or k < 0:
        raise ValueError("nonnegative integer required")
    return F((k + 1)**2, 4**k)


def inverse_word_bound(k):
    if type(k) is not int or k < 0:
        raise ValueError("nonnegative integer required")
    return F(30720**k, factorial(k) * factorial(k + 1))


def ledger():
    """Every normalized source term after the common 1/L bound."""
    return (
        # R1 fixed-linear terms.
        {"name": "R1_Wstar_Phi", "group": "R1_linear", "constant": F(128)},
        {"name": "R1_h_seed_Phi", "group": "R1_linear", "constant": F(37, 50)},
        {"name": "R1_Wstar_DPhi", "group": "R1_linear", "constant": F(128)},
        {"name": "R1_Hstar_Phi_eta", "group": "R1_linear", "constant": F(128)},
        # R1 terms carrying Lambda^-1.  The zeta term records its M multiplier.
        {"name": "R1_B_average_u_Phi", "group": "R1_quadratic", "constant": F(4)},
        {"name": "R1_B_average_u_eta_Phi", "group": "R1_quadratic", "constant": F(10)},
        {"name": "R1_h_u_Phi", "group": "R1_quadratic", "constant": F(2, 25)},
        {"name": "R1_B_average_u_DPhi", "group": "R1_quadratic", "constant": F(4)},
        {"name": "R1_B_average_u_eta_DPhi", "group": "R1_quadratic", "constant": F(10)},
        {"name": "R1_d_u_Phi_eta", "group": "R1_quadratic", "constant": F(10)},
        {"name": "R1_d_u_zeta_Phi", "group": "R1_zeta", "constant": F(10)},
        # R2 fixed-linear and ordinary quadratic terms.
        {"name": "R2_seed_coefficient_u", "group": "R2_linear", "constant": F(186)},
        {"name": "R2_Wstar_Du", "group": "R2_linear", "constant": F(128)},
        {"name": "R2_Hstar_u_eta", "group": "R2_linear", "constant": F(128)},
        {"name": "R2_minus_2A_eta_u2", "group": "R2_quadratic", "constant": F(8)},
        {"name": "R2_B_average_u_Du", "group": "R2_quadratic", "constant": F(4)},
        {"name": "R2_B_average_u_eta_Du", "group": "R2_quadratic", "constant": F(10)},
        {"name": "R2_d_u_u_eta", "group": "R2_quadratic", "constant": F(10)},
        # All three pressure-memory terms.
        {"name": "R2_minus_4A_eta_p", "group": "R2_pressure", "constant": F(16)},
        {"name": "R2_d_p_eta", "group": "R2_pressure", "constant": F(10)},
        {"name": "R2_minus_2eta_Dp", "group": "R2_pressure", "constant": F(8)},
    )


def group_sums(rows=None):
    rows = ledger() if rows is None else rows
    out = {}
    for row in rows:
        out[row["group"]] = out.get(row["group"], F(0)) + row["constant"]
    return out


def majorant_coefficients():
    sums = group_sums()
    # The zeta multiplier is kept symbolic in the theorem.
    q1_constant = sums["R1_quadratic"] * A0**2
    q1_M = sums["R1_zeta"] * A0**2
    q2_ordinary = sums["R2_quadratic"] * A0**2
    q2_pressure = sums["R2_pressure"] * PRESSURE_PRIMITIVE * A0
    return {
        "R1_linear_A0": sums["R1_linear"] * A0,
        "R1_quadratic_A0_squared": q1_constant,
        "R1_zeta_A0_squared": q1_M,
        "R2_linear_A0": sums["R2_linear"] * A0,
        "R2_quadratic_A0_squared": q2_ordinary,
        "pressure_primitive_A0_cubed": PRESSURE_PRIMITIVE,
        "R2_pressure_A0_fourth": q2_pressure,
    }


def lipschitz_coefficients():
    major = majorant_coefficients()
    return {
        "R1_linear": major["R1_linear_A0"],
        "R1_quadratic_per_S": 2 * major["R1_quadratic_A0_squared"],
        "R1_zeta_per_M_S": 2 * major["R1_zeta_A0_squared"],
        "R2_linear": major["R2_linear_A0"],
        "R2_quadratic_per_S": 2 * major["R2_quadratic_A0_squared"],
        "R2_pressure_per_S": 2 * major["R2_pressure_A0_fourth"],
    }


def universal_contraction():
    # H=K^(1/16)>1 and epsilon=10^-16.  E*max(S,Fb,Lb)<=H/10^6,
    # Lambda=10^12 H^128 epsilon^-2.
    epsilon = F(1, 10**16)
    coefficient = epsilon * epsilon / 10**18
    return {
        "epsilon_squared_over_10_power_18": coefficient,
        "H_denominator_power": 127,
        "strict_universal_ceiling": F(1, 10**50),
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
    rows = ledger()
    sums = group_sums(rows)
    major = majorant_coefficients()
    lip = lipschitz_coefficients()
    universal = universal_contraction()
    expected_names = {
        "R1_Wstar_Phi", "R1_h_seed_Phi", "R1_Wstar_DPhi", "R1_Hstar_Phi_eta",
        "R1_B_average_u_Phi", "R1_B_average_u_eta_Phi", "R1_h_u_Phi",
        "R1_B_average_u_DPhi", "R1_B_average_u_eta_DPhi", "R1_d_u_Phi_eta",
        "R1_d_u_zeta_Phi", "R2_seed_coefficient_u", "R2_Wstar_Du",
        "R2_Hstar_u_eta", "R2_minus_2A_eta_u2", "R2_B_average_u_Du",
        "R2_B_average_u_eta_Du", "R2_d_u_u_eta", "R2_minus_4A_eta_p",
        "R2_d_p_eta", "R2_minus_2eta_Dp",
    }
    finite_convolutions = [convolution_sum(n) for n in range(65)]
    finite_allocations = [mixed_allocation(i, j, nu)
                          for i in range(20) for j in range(20) for nu in (1, 2)]
    cauchy = [cauchy_multiplier(k) for k in range(65)]
    inverse = [inverse_word_bound(k) for k in range(12)]
    checks = {
        "all_21_normalized_source_terms_present_once": len(rows) == 21 and {r["name"] for r in rows} == expected_names,
        "R1_linear_sum_is_19237_over_50": sums["R1_linear"] == F(19237, 50),
        "R1_quadratic_sum_is_952_over_25": sums["R1_quadratic"] == F(952, 25),
        "R1_zeta_multiplier_is_10": sums["R1_zeta"] == 10,
        "R2_linear_sum_is_442": sums["R2_linear"] == 442,
        "R2_quadratic_sum_is_32": sums["R2_quadratic"] == 32,
        "R2_pressure_sum_is_34": sums["R2_pressure"] == 34,
        "pressure_primitive_uses_three_algebra_products": major["pressure_primitive_A0_cubed"] == 80 * A0**3,
        "pressure_source_uses_fourth_algebra_power": major["R2_pressure_A0_fourth"] == 2720 * A0**4,
        "pressure_lipschitz_keeps_two_S_factor": lip["R2_pressure_per_S"] == 5440 * A0**4,
        "finite_convolution_controls_below_exact_A0": max(finite_convolutions)**2 < A0,
        "mixed_allocations_at_most_one": max(finite_allocations) <= 1,
        "cauchy_fixed_coefficient_multipliers_at_most_one": max(cauchy) <= 1,
        "inverse_word_recurrence_has_two_factorial_denominators": all(
            inverse[k + 1] / inverse[k] == F(30720, (k + 1) * (k + 2))
            for k in range(len(inverse) - 1)
        ),
        "universal_selected_formula_gives_10_minus_50": universal["epsilon_squared_over_10_power_18"] == F(1, 10**50),
        "remaining_H_power_is_127": universal["H_denominator_power"] == 127,
    }
    return _encode({
        "operator_constants": {
            "A0": A0,
            "D0": D0_SYMBOL,
            "maximum_sampled_convolution_product": max(finite_convolutions)**2,
            "maximum_sampled_mixed_allocation": max(finite_allocations),
            "maximum_sampled_cauchy_multiplier": max(cauchy),
            "finite_samples_are_not_all_index_proofs": True,
        },
        "source_terms": rows,
        "source_sums": sums,
        "majorant_coefficients": major,
        "lipschitz_coefficients": lip,
        "universal_contraction": universal,
        "checks": checks,
    })
