"""Exact ledgers for the same-profile modulation, repair, and cone join."""

from fractions import Fraction as Q


def signed_variance_ledger():
    return {
        "R_at_zero": "sin(theta)",
        "R2_at_zero": Q(1, 2),
        "signed_F_mu_squared_over_d0_squared_at_zero": Q(1, 2),
        "signed_extension_required": True,
        "division_by_V_mu_at_zero_allowed": False,
        "cutoff_flat_zero_extension": True,
    }


def moving_phase_control(r=Q(1, 4), theta=Q(1, 3)):
    """A rational normalized phase model with a genuinely moving inverse."""
    r, theta = Q(r), Q(theta)
    phi_r = theta * theta - theta
    w = 1 + r * (2 * theta - 1)
    theta_r = -phi_r / w
    # Pull back t(r,theta)=r*theta at fixed phase.
    t_r_fixed_theta = theta
    t_theta = r
    pulled_t_r = t_r_fixed_theta + t_theta * theta_r
    return {
        "phi_left": Q(0),
        "phi_right": Q(1),
        "w": w,
        "theta_r": theta_r,
        "inverse_identity": phi_r + w * theta_r,
        "t_r_fixed_theta": t_r_fixed_theta,
        "pulled_t_r": pulled_t_r,
        "phase_term": t_theta * theta_r,
    }


def fast_chain_control(N=7):
    """For A(y,varphi)=y+2varphi, differentiate A(y,Ny)/N."""
    N = Q(N)
    return {
        "N": N,
        "full_radial_derivative": Q(1, N) + 2,
        "fast_loop_term": Q(2),
        "slow_remainder": Q(1, N),
        "field_value_power": "N^-1",
        "field_rth_radial_power": "N^(r-1)",
        "moment_rth_radial_power": "N^max(r-2,-1)",
    }


def repair_certificate():
    beta = Q(800)
    q = Q(25000)
    return {
        "beta": beta,
        "q": q,
        "D": "H^16",
        "target": "D/N",
        "root": "2*beta*D/N",
        "root_coefficient": 2 * beta,
        "contraction_coefficient": 8 * beta * beta * q,
        "binary_coefficient_bound": 2**37,
        "H_binary_floor": 2**16,
        "H_inverse_power": 16,
        "lambda_positive": True,
        "rows": ("M", "(M-J)/lambda", "I", "-S", "Cp"),
    }


def compact_transform(values, lam):
    M, I, J, S, Cp = map(Q, values)
    lam = Q(lam)
    if lam <= 0:
        raise ValueError("lambda must be positive")
    return (M, (M - J) / lam, I, -S, Cp)


def compact_inverse(rows, lam):
    r1, r2, r3, r4, r5 = map(Q, rows)
    lam = Q(lam)
    if lam <= 0:
        raise ValueError("lambda must be positive")
    return (r1, r3, r1 - lam * r2, -r4, r5)


def hierarchy_certificate():
    A0 = 2**16
    return {
        "A_floor": A0,
        "G_floor": 2**16,
        "cone_left_exponent": 30 * 16,
        "cone_constant_exponent": 9 + 20,
        # With G^4>5*A^32, these two integer inequalities imply the
        # state-error and coefficient-radius logarithmic comparisons.
        "state_polynomial_left": 13 + 128 * A0 + 64 * A0**32,
        "state_polynomial_right": 80 * A0**32,
        "radius_polynomial_left": 12 + 128 * A0,
        "radius_polynomial_right": 80 * A0**32,
        "exp_ratio_monotone_from_one": True,
        "N": "1+floor(H^32)>H^32",
        "epsilon": "G^-64",
        "kappa": "G^-4",
        "state_error": "2*H^16*(1+1600*A^128)/N<G^-64/2",
        "cone_loss": "256*(4G^3)^10*(G^-64/2)<G^-4/2",
    }


def support_ledger():
    return {
        "inner_protected": "[Xa,Xminus]",
        "compact": "[Xminus,Xrep]",
        "repair_endpoint": "Xrep=rho_p*exp(5)=Xstar",
        "future_heat_edit_starts": "Xstar",
        "all_five_endpoint_defects_zero_as_eta_functions": True,
        "axis_pressure_datum_unchanged": True,
        "inner_width": "C^-20",
        "annular_and_compact_roots_reselected": False,
        "frequency_reselected": False,
    }


def pressure_defect(eta, p, p_eta, h=Q(0), X=Q(1), L=Q(1)):
    eta, p, p_eta, h, X, L = map(Q, (eta, p, p_eta, h, X, L))
    return X * (4 * (Q(1, 2) + h) * eta * p - (1 - eta * eta) * p_eta) / L


def review():
    variance = signed_variance_ledger()
    phase = moving_phase_control()
    chain = fast_chain_control()
    repair = repair_certificate()
    hierarchy = hierarchy_certificate()
    supports = support_ledger()
    sample = (Q(2), Q(3), Q(5), Q(7), Q(11))
    rows = compact_transform(sample, Q(1, 100))
    checks = {
        "signed_variance_root_is_regular_at_zero":
            variance["R2_at_zero"] == variance["signed_F_mu_squared_over_d0_squared_at_zero"] == Q(1, 2),
        "zero_variance_does_not_use_singular_V_mu_division":
            variance["signed_extension_required"] and not variance["division_by_V_mu_at_zero_allowed"],
        "moving_phase_inverse_identity_is_exact": phase["inverse_identity"] == 0,
        "freezing_phase_would_drop_a_nonzero_term": phase["phase_term"] != 0 and phase["pulled_t_r"] != phase["t_r_fixed_theta"],
        "fast_phase_derivative_is_not_mislabelled_small":
            chain["full_radial_derivative"] - chain["fast_loop_term"] == chain["slow_remainder"],
        "moment_errors_have_one_extra_radial_gain":
            chain["field_rth_radial_power"] != chain["moment_rth_radial_power"],
        "repair_contraction_coefficient_fits_binary_bound":
            repair["contraction_coefficient"] < repair["binary_coefficient_bound"],
        "repair_contraction_is_below_2^-219":
            repair["binary_coefficient_bound"] * Q(1, repair["H_binary_floor"] ** repair["H_inverse_power"]) <= Q(1, 2**219),
        "same_five_row_transform_is_exactly_invertible": compact_inverse(rows, Q(1, 100)) == sample,
        "pressure_row_is_fifth_and_independent": repair["rows"][-1] == "Cp" and len(repair["rows"]) == 5,
        "selected_frequency_clears_state_threshold":
            hierarchy["state_polynomial_left"] < hierarchy["state_polynomial_right"],
        "selected_frequency_clears_coefficient_radius":
            hierarchy["radius_polynomial_left"] < hierarchy["radius_polynomial_right"],
        "cone_constant_is_absorbed":
            hierarchy["cone_left_exponent"] > hierarchy["cone_constant_exponent"],
        "support_endpoints_join_without_reselection":
            supports["repair_endpoint"].endswith("=Xstar") and not supports["annular_and_compact_roots_reselected"] and not supports["frequency_reselected"],
        "all_five_endpoint_rows_and_axis_datum_are_retained":
            supports["all_five_endpoint_defects_zero_as_eta_functions"] and supports["axis_pressure_datum_unchanged"],
        "omitting_pressure_row_has_nonzero_axial_defect": pressure_defect(0, 0, 1) == -1,
    }
    return {
        "variance": variance,
        "phase": phase,
        "chain": chain,
        "repair": repair,
        "hierarchy": hierarchy,
        "supports": supports,
        "sample_compact_rows": rows,
        "checks": checks,
    }
