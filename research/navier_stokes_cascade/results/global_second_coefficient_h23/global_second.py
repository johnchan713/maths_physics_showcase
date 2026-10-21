#!/usr/bin/env python3
"""Exact assembly of the full-support second repair and global H2,3 bound."""

from fractions import Fraction as F
from math import comb


def rat(value):
    if isinstance(value, (bool, float)):
        raise TypeError("Exact rational inputs required")
    return F(value)


def serial(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(item) for item in value]
    return value


def regular_v_bound(k, h, angular_order, u_next):
    """Equation (13) of the ordered full-support moment audit."""
    h, u_next = map(rat, (h, u_next))
    if type(k) is not int or k < 0 or angular_order not in (2, 3):
        raise ValueError("Coefficient order and angular order 2 or 3 required")
    if not 0 < h <= F(1, 100) or u_next < 0:
        raise ValueError("Retained h range and nonnegative U bound required")
    ell = 1 - 2 * h
    r2 = 1 / ell + 6 * h / ell**2 + 16 * h**2 / ell**3
    r3 = r2 + 16 * h**2 / ell**3 + 64 * h**3 / ell**4
    reciprocal = r2 if angular_order == 2 else r3
    D = F(1, 2) - h
    return reciprocal * (4 + 4 * (D + 2 * h * k) + 4 * (angular_order + 1)) * u_next


def shifted_second_source(h, nu, x_velocity, leading_u5, first_u3):
    """Full-support Q10,Q11, including both shifted axial-viscosity moments."""
    h, nu, x_velocity, leading_u5, first_u3 = map(
        rat, (h, nu, x_velocity, leading_u5, first_u3)
    )
    if nu < 0 or x_velocity <= 0 or min(leading_u5, first_u3) < 0:
        raise ValueError("Admissible support and nonnegative coefficient bounds required")
    A, D = F(1, 2) + h, F(1, 2) - h
    lam = 2 * h
    ell = 1 - 2 * h
    r1 = 1 / ell + 4 * h / ell**2
    r2 = 1 / ell + 6 * h / ell**2 + 16 * h**2 / ell**3
    u = (leading_u5, first_u3)
    v = (
        regular_v_bound(0, h, 2, leading_u5),
        regular_v_bound(1, h, 2, first_u3),
    )
    previous_v3 = regular_v_bound(0, h, 3, leading_u5)
    mixed = u[0] * v[1] + u[1] * v[0]
    quadratic = 2 * v[0] * v[1]
    shifted, totals = [], []
    for weight in (0, 1):
        volume = x_velocity ** (weight + 1) / F(weight + 1)
        beta = weight
        axial = (
            nu
            * r1
            * r2
            * (6 + 4 * abs(beta - D))
            * (12 + 4 * abs(beta))
            * volume
            * previous_v3
        )
        if weight == 0:
            total = (
                r1
                * volume
                * ((lam + 4 * D) * v[1] + (6 + 4 * abs(lam - A)) * mixed)
                + volume * quadratic / 2
                + 2 * nu * v[1]
                + axial
            )
        else:
            total = (
                r1
                * volume
                * ((1 + lam + 4 * D) * v[1] + (6 + 4 * abs(lam + D)) * mixed)
                + volume * quadratic / 2
                + axial
            )
        shifted.append(axial)
        totals.append(total)
    return {
        "u_Ceta2": u,
        "v_Ceta2": v,
        "previous_v_Ceta3": previous_v3,
        "ordered_mixed_B": mixed,
        "ordered_quadratic_C": quadratic,
        "shifted_axial_viscosity": tuple(shifted),
        "Q10": totals[0],
        "Q11": totals[1],
        "R1": r1,
        "R2": r2,
    }


def first_field_values(local_first, repair_c13, radial_sqrt_lower, bump_value=320):
    """Keep F1's inverse-radius repair cost separate from the U1 bound."""
    local_first, repair_c13, radial_sqrt_lower, bump_value = map(
        rat, (local_first, repair_c13, radial_sqrt_lower, bump_value)
    )
    if min(local_first, repair_c13, bump_value) < 0 or radial_sqrt_lower <= 0:
        raise ValueError("Nonnegative envelopes and positive radial square-root lower bound required")
    return {
        "U1": local_first + bump_value * repair_c13,
        "F1": local_first + bump_value * repair_c13 * max(F(1), 1 / radial_sqrt_lower),
        "bump_value": bump_value,
    }


def coefficient_cutoff_jets(delta, x_keep, normalized=(1, 729, 2916, 58320)):
    """Physical X jets through three, retaining logarithmic Jacobian terms."""
    delta, x_keep = map(rat, (delta, x_keep))
    normalized = tuple(map(rat, normalized))
    if delta <= 0 or x_keep <= 0 or len(normalized) != 4 or normalized[0] != 1:
        raise ValueError("Positive cutoff scales and four normalized jets with L0=1 required")
    if min(normalized) < 0:
        raise ValueError("Nonnegative cutoff jet bounds required")
    l0, l1, l2, l3 = normalized
    return (
        l0,
        l1 / (delta * x_keep),
        (l2 / delta**2 + l1 / delta) / x_keep**2,
        (l3 / delta**3 + 3 * l2 / delta**2 + 2 * l1 / delta) / x_keep**3,
    )


def product_jets(cutoff, inner):
    cutoff, inner = tuple(map(rat, cutoff)), tuple(map(rat, inner))
    if len(cutoff) != 4 or len(inner) != 4 or min((*cutoff, *inner)) < 0:
        raise ValueError("Four nonnegative cutoff and inner jets required")
    return tuple(
        sum(F(comb(order, j)) * cutoff[j] * inner[order - j] for j in range(order + 1))
        for order in range(4)
    )


def finalized_first_jets(inner_first, cutoff, repair_c13, p_lower, radial_sqrt_lower):
    """Finalize the first F,U radial rows used by the second collar envelope."""
    repair_c13 = rat(repair_c13)
    if repair_c13 < 0:
        raise ValueError("Nonnegative first repair coefficient required")
    cut = product_jets(cutoff, inner_first)
    bumps = repair_bump_jets(p_lower, radial_sqrt_lower)
    U = tuple(cut[j] + repair_c13 * bumps["U"][j] for j in range(4))
    F = tuple(cut[j] + repair_c13 * bumps["E_to_F"][j] for j in range(4))
    return {"cut": cut, "U": U, "F": F, "mixed_Ceta3": max((*U, *F))}


def second_inner_jets(core_factors, collar_affine, t_collar, m2_core, m2_collar):
    """Evaluate the nested core/collar J2,r formula before the second cutoff."""
    core_factors = tuple(map(rat, core_factors))
    collar_affine = tuple(tuple(map(rat, row)) for row in collar_affine)
    t_collar, m2_core, m2_collar = map(rat, (t_collar, m2_core, m2_collar))
    if (
        len(core_factors) != 4
        or len(collar_affine) != 4
        or any(len(row) != 2 for row in collar_affine)
        or min((*core_factors, *(item for row in collar_affine for item in row))) < 0
        or t_collar < 1
        or min(m2_core, m2_collar) < 0
    ):
        raise ValueError("Four nonnegative core/collar rows, T>=1 and Picard bounds required")
    powers = (0, 0, 2, 4)
    return tuple(
        core_factors[r] * m2_core
        + t_collar**powers[r] * (collar_affine[r][0] * m2_collar + collar_affine[r][1])
        for r in range(4)
    )


def second_discrepancies(x_cut, x_velocity, inner_value, f0, u0, f1, u1, q10, q11):
    """Five ordinary second-order discrepancies with fixed lower squares."""
    values = tuple(map(rat, (x_cut, x_velocity, inner_value, f0, u0, f1, u1, q10, q11)))
    x_cut, x_velocity, inner_value, f0, u0, f1, u1, q10, q11 = values
    if min(x_cut, x_velocity) <= 0 or min(values[2:]) < 0:
        raise ValueError("Positive supports and nonnegative bounds required")
    return (
        x_cut * inner_value,
        x_cut**2 * inner_value,
        2 * f0 * x_cut * inner_value + x_velocity * f1**2 + q10 / 2,
        x_cut**2 * (u0 + f0) * inner_value + x_velocity**2 * u1 * f1,
        2 * u0 * x_cut * inner_value
        + f0 * x_cut**2 * inner_value
        + x_velocity * u1**2
        + x_velocity**2 * f1**2 / 2
        + q11 / 2,
    )


def inverse_800(rho, lam, geometric_k, e_star, discrepancies):
    rho, lam, geometric_k, e_star = map(rat, (rho, lam, geometric_k, e_star))
    moments = tuple(map(rat, discrepancies))
    if min(rho, lam, geometric_k, e_star) <= 0 or len(moments) != 5 or min(moments) < 0:
        raise ValueError("Positive geometry and five nonnegative discrepancies required")
    m1, m2, m3, m4, m5 = moments
    w1 = 4 / e_star
    rows = (
        m1 / rho,
        (m1 / rho + geometric_k * w1 * m4) / lam,
        geometric_k * m2,
        w1 * m5 / rho,
        w1 * m3,
    )
    return {"W1": w1, "rows": rows, "coefficient_Ceta1": 800 * max(rows)}


def repair_bump_jets(p_lower, radial_sqrt_lower, bounds=(320, 121600, 220000000, 18870643200000)):
    """Actual width-1/10 U and relative-E-to-F radial jets through order three."""
    p_lower, radial_sqrt_lower = map(rat, (p_lower, radial_sqrt_lower))
    b0, b1, b2, b3 = map(rat, bounds)
    if min(p_lower, radial_sqrt_lower) <= 0 or min(b0, b1, b2, b3) < 0:
        raise ValueError("Positive support lower bounds and nonnegative bump jets required")
    U = (
        b0,
        b1 / p_lower,
        (b2 + b1) / p_lower**2,
        (b3 + 3 * b2 + 2 * b1) / p_lower**3,
    )
    E_to_F = (
        b0 / radial_sqrt_lower,
        (b1 + b0 / 2) / (p_lower * radial_sqrt_lower),
        (b2 + 2 * b1 + 3 * b0 / 4) / (p_lower**2 * radial_sqrt_lower),
        (b3 + F(9, 2) * b2 + F(23, 4) * b1 + F(15, 8) * b0)
        / (p_lower**3 * radial_sqrt_lower),
    )
    return {"logarithmic": (b0, b1, b2, b3), "U": U, "E_to_F": E_to_F}


def f_to_e_jets(f_jets, radial_sqrt_lower, radial_sqrt_upper):
    """Convert finalized F jets to E=sqrt(2X)F on the positive support."""
    f_jets = tuple(map(rat, f_jets))
    radial_sqrt_lower, radial_sqrt_upper = map(
        rat, (radial_sqrt_lower, radial_sqrt_upper)
    )
    if len(f_jets) != 4 or min(f_jets) < 0:
        raise ValueError("Four nonnegative F jets required")
    if not 0 < radial_sqrt_lower <= radial_sqrt_upper:
        raise ValueError("Ordered positive square-root radius bounds required")
    radius = (
        radial_sqrt_upper,
        1 / radial_sqrt_lower,
        1 / radial_sqrt_lower**3,
        3 / radial_sqrt_lower**5,
    )
    E = tuple(
        sum(F(comb(order, j)) * radius[j] * f_jets[order - j] for j in range(order + 1))
        for order in range(4)
    )
    return {"radius_jets": radius, "E": E}


def finalized_h23(
    h, nu, x_velocity, leading_u5, local_first, repair_c13,
    radial_sqrt_lower, radial_sqrt_upper, x_cut, inner_second_jets, cutoff_jets,
    leading_f1, leading_u1, rho, lam, geometric_k, e_star, p_lower,
):
    first = first_field_values(local_first, repair_c13, radial_sqrt_lower)
    source = shifted_second_source(h, nu, x_velocity, leading_u5, first["U1"])
    cut = product_jets(cutoff_jets, inner_second_jets)
    discrepancies = second_discrepancies(
        x_cut, x_velocity, cut[0], leading_f1, leading_u1,
        first["F1"], first["U1"], source["Q10"], source["Q11"],
    )
    repair = inverse_800(rho, lam, geometric_k, e_star, discrepancies)
    bumps = repair_bump_jets(p_lower, radial_sqrt_lower)
    final_u = tuple(
        cut[j] + repair["coefficient_Ceta1"] * bumps["U"][j]
        for j in range(4)
    )
    final_f = tuple(
        cut[j] + repair["coefficient_Ceta1"] * bumps["E_to_F"][j]
        for j in range(4)
    )
    conversion = f_to_e_jets(final_f, radial_sqrt_lower, radial_sqrt_upper)
    return {
        "first_fields": first,
        "source": source,
        "cut_second_jets": cut,
        "discrepancies": discrepancies,
        "repair": repair,
        "repair_bumps": bumps,
        "final_U2_jets": final_u,
        "final_F2_jets": final_f,
        "F_to_E": conversion,
        "H23": max((*final_u, *conversion["E"])),
    }


def selected_formula():
    return {
        "leading": "K0,5=H^262 bounds every required leading Ceta2,Ceta3,Ceta4 row",
        "first": "U1<=K1,3*=J13+320C13; F1<=J13+320C13 max(1,rp^-1)",
        "source": "Q10,Q11 are the full [0,Xv] k=1 formulas with both shifted axial-viscosity terms",
        "inner": "J2,r from the nested second core/collar; P2,r=sum binom(r,j) k_j J2,r-j",
        "collar_envelope": "Tcol=max(1,H0_mixedJ3, finalized_first_FU_mixedCeta3); J2,r=Ac,r M2core+Tcol^p_r(Acol,r M2col+Bcol,r)",
        "repair": "C21=800 max(mu1/rho3,(mu1/rho3+Kg(4/e*)mu4)/lambda,Kg mu2,(4/e*)mu5/rho3,(4/e*)mu3)",
        "final": "finalize U2,F2 jets; convert E2=sqrt(2X)F2 by the retained Rlo,Rhi product rule; H2,3*=max_r<=3(U2_r,E2_r)",
        "cutoff": "k1,k2,k3 retain all Delta^-j logarithmic-Jacobian terms with L=(1,729,2916,58320)",
    }


def review():
    cutoff = coefficient_cutoff_jets(F(1, 8), 2)
    first_jets = finalized_first_jets((2, 3, 5, 7), cutoff, 17, 31, F(1, 2))
    t_collar = max(F(1), 101, first_jets["mixed_Ceta3"])
    inner_second = second_inner_jets(
        (2, 3, 5, 7), ((11, 13), (17, 19), (23, 29), (31, 37)),
        t_collar, 41, 43,
    )
    out = finalized_h23(
        F(1, 200), 1, 7, 11, 13, 17, F(1, 2), 3, 19,
        inner_second, cutoff, 23, 29, 31, F(1, 1000), 37, 41, 31,
    )
    inviscid = shifted_second_source(F(1, 200), 0, 7, 11, out["first_fields"]["U1"])
    return serial({
        "selected_formula": selected_formula(),
        "exact_fixture": {
            "cutoff": cutoff,
            "finalized_first_jets": first_jets,
            "Tcollar": t_collar,
            "inner_second_jets": inner_second,
            "assembled": out,
        },
        "checks": {
            "first_F_and_U_are_not_conflated": out["first_fields"]["F1"] > out["first_fields"]["U1"],
            "both_shifted_axial_viscosity_moments_retained": (
                out["source"]["Q10"] > inviscid["Q10"]
                and out["source"]["Q11"] > inviscid["Q11"]
            ),
            "both_ordered_cross_terms_retained": out["source"]["ordered_mixed_B"] > 0,
            "lower_first_squares_enter_repair": (
                out["discrepancies"][2] > out["source"]["Q10"] / 2
                and out["discrepancies"][4] > out["source"]["Q11"] / 2
            ),
            "cutoff_third_jet_keeps_lower_jacobian_terms": cutoff[3] == F(58320 * 8**3 + 3 * 2916 * 8**2 + 2 * 729 * 8, 2**3),
            "inverse_uses_full_W1": out["repair"]["W1"] == F(4, 41),
            "actual_third_bump_jet_is_used": out["repair_bumps"]["logarithmic"][3] == 18870643200000,
            "second_collar_uses_finalized_first_radial_jets": t_collar == first_jets["mixed_Ceta3"],
            "nested_second_core_and_collar_are_both_present": inner_second[0] == 2 * 41 + 11 * 43 + 13,
            "disjoint_repair_has_no_factor_five": out["final_U2_jets"][0] == out["cut_second_jets"][0] + out["repair"]["coefficient_Ceta1"] * 320,
            "F_is_converted_to_E_before_H23": out["F_to_E"]["E"][0] == 3 * out["final_F2_jets"][0],
            "H23_is_maximum_of_all_final_E_and_U_jets": out["H23"] == max((*out["final_U2_jets"], *out["F_to_E"]["E"])),
        },
    })


if __name__ == "__main__":
    import json
    print(json.dumps(review(), indent=2, sort_keys=True))
