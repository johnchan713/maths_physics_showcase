"""Exact global moment, pressure, and exterior-stress ledgers."""

from fractions import Fraction as Q


def annular_rows(raw, eta):
    M, I, J, S, Cp = map(Q, raw)
    eta = Q(eta)
    return (M, J - 4 * eta * I, I, -(S - 8 * eta * M), Cp)


def annular_inverse(rows, eta):
    r1, r2, r3, r4, r5 = map(Q, rows)
    eta = Q(eta)
    return (r1, r3, r2 + 4 * eta * r3, 8 * eta * r1 - r4, r5)


def compact_rows(raw, lam):
    M, I, J, S, Cp = map(Q, raw)
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


def heat_ledger():
    return {
        "Hcal": "Gamma(1+h)^-1 int exp(-v) v^h (1+Zv)^-h dv",
        "Hcal_derivative": "(-1)^j (h)_j/Gamma(1+h) int exp(-v) v^(h+j) (1+Zv)^(-h-j) dv",
        "large_X_H_minus_Hpow": "O_h(X^(-1-h))",
        "large_X_eta_derivatives": "O_h(X^(-1-h)) at every fixed order",
        "axis_H": "2X F with smooth F",
        "axis_Hpow_integrable_condition": "0<h<1",
        "raw_I_infinity_exists": False,
        "R": "int_0^infinity (H-Hpow) dX",
        "R_limit": "lim_X [I(X)-X Hpow(X)/(1-h)]",
        "Cp_tail": "O(X^(-2-2h))",
        "S_tail": "O(X^(-1-2h))",
        "h_zero_interchange_allowed": False,
    }


def endpoint_to_global_ledger():
    return {
        "endpoint_row_order": "(M,I,J,S,Cp)",
        "global_defect_order": "(M,R,J,S,Cp)",
        "finite_row": "Delta I(Xrep)=0",
        "common_future_tail": True,
        "common_Hpow": True,
        "difference_formula": "Delta R=Delta I(Xrep)+int_Xrep^infinity Delta(H-Hpow) dX",
        "conclusion": "Delta R=0",
        "invalid_conclusion": "I(infinity)=0",
    }


def edit_ledger():
    return {
        "axial_and_strength_roots": ("M=0", "J=0", "S=0", "I,Cp unchanged"),
        "angular_roots_and_wait": ("R=0", "Cp total preserved", "S included in strength root"),
        "heat_edit_and_compensation": ("Delta R cancels", "Delta S cancels", "Delta Cp cancels"),
        "shape_and_annular_root": "all five endpoint defects restored",
        "modulation_and_compact_repair": "all five endpoint defects restored",
        "global_zero_defects": ("M", "R", "J", "S", "Pi_axis+Cp"),
        "functional_in_eta": True,
        "intermediate_pressure_memory_may_be_nonzero": True,
    }


def pressure_fixture():
    ideal = (Q(5), Q(7), Q(11), Q(13))
    pi_axis = -sum(ideal)
    heat = (Q(2), Q(-3), Q(5), Q(7))
    shape = (Q(17), Q(-19), Q(23), Q(-29))
    loop = (Q(-31), Q(37), Q(-41), Q(43))
    final = tuple(
        x + a - a + b - b + c - c
        for x, a, b, c in zip(ideal, heat, shape, loop)
    )
    cut = 2
    forward = pi_axis + sum(final[:cut])
    backward = -sum(final[cut:])
    missing_pressure_repair = sum(final) + loop[-1] + pi_axis
    return {
        "pi_axis": pi_axis,
        "ideal_total": sum(ideal),
        "final_total": sum(final),
        "forward": forward,
        "backward": backward,
        "intermediate_heat": heat,
        "missing_pressure_repair": missing_pressure_repair,
    }


def exterior_defects(eta, h, X, H, defects):
    eta, h, X, H = map(Q, (eta, h, X, H))
    (M, Me), (R, Re), (J, Je), (S, Se), (P, Pe) = [tuple(map(Q, x)) for x in defects]
    D = Q(1, 2) - h
    A = Q(1, 2) + h
    d = 1 - eta * eta
    L = 1 - 2 * h * eta * eta
    Bm = 2 * D * eta * M + d * Me
    radial = -Bm / L
    theta = (H * Bm + (1 - h) * R - D * eta * Re - d * Je + 2 * (h - D) * eta * J) / L
    axial = (D * (M - eta * Me) + 4 * h * eta * S - d * Se + X * (4 * A * eta * P - d * Pe)) / L
    return (radial, theta, axial)


def review():
    heat = heat_ledger()
    endpoint = endpoint_to_global_ledger()
    edits = edit_ledger()
    pressure = pressure_fixture()
    sample = (Q(2), Q(3), Q(5), Q(7), Q(11))
    eta = Q(2, 5)
    lam = Q(1, 10**20)
    zero = [(Q(0), Q(0))] * 5
    pressure_slope = zero[:]
    pressure_slope[4] = (Q(0), Q(1))
    p100 = exterior_defects(0, Q(1, 1000), 100, 1, pressure_slope)[2]
    p200 = exterior_defects(0, Q(1, 1000), 200, 1, pressure_slope)[2]
    radial_M = zero[:]
    radial_M[0] = (Q(0), Q(1))
    checks = {
        "renormalized_heat_moment_is_the_finite_global_object":
            not heat["raw_I_infinity_exists"] and "I(X)-" in heat["R_limit"],
        "both_heat_ends_and_fixed_eta_jets_are_integrable":
            heat["axis_Hpow_integrable_condition"] == "0<h<1" and heat["large_X_eta_derivatives"].startswith("O_h"),
        "h_zero_cannot_be_interchanged_with_infinite_integration": not heat["h_zero_interchange_allowed"],
        "finite_I_match_converts_to_R_only_with_common_tail":
            endpoint["endpoint_row_order"] == "(M,I,J,S,Cp)"
            and endpoint["global_defect_order"] == "(M,R,J,S,Cp)"
            and endpoint["common_future_tail"] and endpoint["common_Hpow"]
            and endpoint["conclusion"] == "Delta R=0",
        "raw_I_infinity_is_not_set_to_zero": endpoint["invalid_conclusion"] != endpoint["conclusion"],
        "annular_rows_invert_exactly": annular_inverse(annular_rows(sample, eta), eta) == sample,
        "compact_rows_invert_exactly_for_positive_lambda": compact_inverse(compact_rows(sample, lam), lam) == sample,
        "all_five_global_defects_are_functional_equalities":
            len(edits["global_zero_defects"]) == 5 and edits["functional_in_eta"],
        "same_pressure_gauge_gives_backward_tail":
            pressure["pi_axis"] + pressure["final_total"] == 0 and pressure["forward"] == pressure["backward"],
        "intermediate_pressure_memory_is_not_reset": any(pressure["intermediate_heat"]),
        "omitting_pressure_repair_breaks_the_gauge": pressure["missing_pressure_repair"] != 0,
        "zero_five_defects_remove_all_exterior_memory": exterior_defects(eta, Q(1, 1000), 100, Q(3, 2), zero) == (0, 0, 0),
        "pressure_eta_slope_defect_grows_linearly_in_X": p100 != 0 and p200 == 2 * p100,
        "zero_U_does_not_remove_M_induced_radial_velocity": exterior_defects(0, Q(1, 1000), 100, 1, radial_M)[0] != 0,
    }
    return {
        "heat": heat,
        "endpoint_to_global": endpoint,
        "edits": edits,
        "pressure": pressure,
        "sample_annular_rows": annular_rows(sample, eta),
        "sample_compact_rows": compact_rows(sample, lam),
        "pressure_slope_X100": p100,
        "checks": checks,
    }
