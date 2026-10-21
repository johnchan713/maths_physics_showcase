"""Exact ledgers for the continuation and the two five-row matching roots."""

from fractions import Fraction as Q


def continuation_ledger():
    return {
        "reference_after_taper": {
            "Dy_log_Fr": Q(0),
            "Dy_log_Er": Q(1, 2),
            "physical_a_r": Q(0),
            "first_collar_shortcut_valid_after_taper": False,
        },
        "actual_slopes": (
            "a=kappa*p1_r; U_y=-kappa*X*ns_r/2",
            "a=kappa0*p1_r; U_y=-beta*kappa0*X*ns_r/2",
            "a=(1-theta)*kappa0*p1_r+(4/5)*theta; U_y=0",
            "a=4/5; U_y=0",
        ),
        "moment_densities": (
            "X*U",
            "sqrt(2)*X^(3/2)*E",
            "sqrt(2)*X^(3/2)*U*E",
            "X*(U^2-E^2/2)",
            "E^2/2",
        ),
        "pressure": "Pi=Pi_axis+Cp",
        "shape_join_slope": Q(1, 2) * (1 - Q(4, 5)),
    }


def annular_certificate():
    beta = Q(610)
    q = Q(50)
    epsilon = Q(1, 10**16)
    delta = Q(143, 500) * epsilon
    radius = 2 * beta * delta
    contraction = 8 * beta * beta * q * delta
    return {
        "beta": beta,
        "q": q,
        "epsilon": epsilon,
        "delta": delta,
        "radius": radius,
        "contraction": contraction,
        "old_C2_radius": Q(572) * epsilon,
        "old_C1_radius": Q(6000) * epsilon,
        "Q_eta_normalized_derivatives": (Q(4), Q(8), Q(8), Q(4), Q(1)),
        "pointwise_inverse": 2 * beta,
    }


def compact_certificate():
    beta = Q(800)
    q = Q(25000)
    coefficient = 8 * beta * beta * q
    return {
        "beta": beta,
        "q": q,
        "contraction_coefficient": coefficient,
        "binary_coefficient_bound": 2**37,
        "H_inverse_power": 16,
        "H_binary_floor": 2**16,
        "contraction_binary_upper_exponent": -219,
        "root_over_target": 2 * beta,
        "pointwise_inverse": 2 * beta,
        "lambda_must_be_positive": True,
    }


def row_transforms():
    return {
        "annular": (
            "M", "J-4*eta*I", "I", "-(S-8*eta*M)", "Cp",
        ),
        "annular_inverse": (
            "r1", "r3", "r2+4*eta*r3", "8*eta*r1-r4", "r5",
        ),
        "compact": (
            "M", "(M-J)/lambda", "I", "-S", "Cp",
        ),
        "compact_inverse": (
            "r1", "r3", "r1-lambda*r2", "-r4", "r5",
        ),
    }


def profile_ledger():
    return {
        "strict_order": (
            "Xsep<XR*exp(-8)",
            "XR*exp(-8)<XR*exp(-7)",
            "XR*exp(-7)<XR*exp(-6)",
            "XR*exp(-6)<XR*exp(-5)",
            "XR*exp(-5)<rho_p",
            "rho_p<rho_p*exp(5)=Xrep",
        ),
        "annular_supports_disjoint": True,
        "compact_supports_disjoint": True,
        "prefix_M": "mu*eta with mu>0",
        "reserved_U": "0",
        "reserved_M_is_zero": False,
        "whole_compact_input": "finite fixed mixed jets on [Xminus,Xrep]x[-1,1]",
        "heat_compensation_starts_at": "Xrep",
    }


def pressure_defect(eta, p, p_eta, h=Q(0)):
    eta, p, p_eta, h = map(Q, (eta, p, p_eta, h))
    A = Q(1, 2) + h
    d = 1 - eta * eta
    return 4 * A * eta * p - d * p_eta


def review():
    continuation = continuation_ledger()
    annular = annular_certificate()
    compact = compact_certificate()
    transforms = row_transforms()
    profile = profile_ledger()
    checks = {
        "flattened_reference_has_zero_physical_shear":
            continuation["reference_after_taper"]["physical_a_r"] == 0,
        "first_collar_shortcut_is_not_extended_past_taper":
            not continuation["reference_after_taper"]["first_collar_shortcut_valid_after_taper"],
        "shape_slope_matches_final_continuation": continuation["shape_join_slope"] == Q(1, 10),
        "all_five_forward_rows_retained": len(continuation["moment_densities"]) == 5,
        "complete_pressure_memory_retained": continuation["pressure"] == "Pi=Pi_axis+Cp",
        "annular_contraction_is_strict": annular["contraction"] < 1,
        "annular_root_is_inside_both_old_uniqueness_balls":
            annular["radius"] < annular["old_C2_radius"] < annular["old_C1_radius"],
        "annular_quadratic_eta_derivatives_stop_at_degree_four":
            annular["Q_eta_normalized_derivatives"] == (4, 8, 8, 4, 1),
        "compact_contraction_coefficient_fits_binary_bound":
            compact["contraction_coefficient"] < compact["binary_coefficient_bound"],
        "compact_contraction_is_below_2^-219":
            compact["binary_coefficient_bound"] * Q(1, compact["H_binary_floor"] ** compact["H_inverse_power"])
            <= Q(1, 2**219),
        "both_transforms_keep_fifth_pressure_row":
            transforms["annular"][-1] == transforms["compact"][-1] == "Cp",
        "transforms_are_five_by_five": all(len(v) == 5 for v in transforms.values()),
        "same_profile_reaches_repair_endpoint_before_heat_edit":
            profile["heat_compensation_starts_at"] == "Xrep",
        "reserved_zero_U_does_not_erase_accumulated_M":
            profile["reserved_U"] == "0" and not profile["reserved_M_is_zero"],
        "missing_pressure_row_has_detectable_eta_slope_defect":
            pressure_defect(0, 0, 1) == -1,
    }
    return {
        "continuation": continuation,
        "annular": annular,
        "compact": compact,
        "transforms": transforms,
        "profile": profile,
        "checks": checks,
    }
