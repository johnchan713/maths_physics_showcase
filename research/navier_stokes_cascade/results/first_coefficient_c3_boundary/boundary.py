#!/usr/bin/env python3
"""Exact order ledger and bounds at the global first-coefficient C_eta^3 boundary."""

from fractions import Fraction as F
from math import comb, factorial


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


def derivative_ladder(moment_order):
    """Minimal leading angular order for signed leading-source moments."""
    if type(moment_order) is not int or moment_order < 0:
        raise ValueError("moment_order must be a nonnegative integer")
    return {
        "target_signed_moment_order": moment_order,
        "a_b_moment_order": moment_order + 1,
        "regular_v_order": moment_order + 1,
        "leading_U_order": moment_order + 2,
        "reason": (
            "I0,I1 contain one eta derivative of a,b; v contains one eta "
            "derivative of the regular radial average of U"
        ),
    }


def reciprocal_cap(h):
    """Safe C_eta^m cap for 1/(1-2h eta^2), any m>=2, h<=1/100."""
    h = rat(h)
    if not 0 < h <= F(1, 100):
        raise ValueError("Retained 0<h<=1/100 required")
    # |2h eta^2|_m <= 8h in every factorial norm m>=2.
    return 1 / (1 - 8 * h)


def regular_v4_bound(h, leading_u5):
    """Bound the regular leading v in C_eta^4 from leading U in C_eta^5."""
    h, leading_u5 = map(rat, (h, leading_u5))
    if leading_u5 < 0:
        raise ValueError("Nonnegative leading U C_eta^5 bound required")
    D = F(1, 2) - h
    factor = reciprocal_cap(h) * (4 + 4 * D + 20)
    return {"factor": factor, "v4": factor * leading_u5}


def signed_source_c3_bounds(h, nu, Xv, leading_u5):
    """Full-support signed I0,I1 C_eta^3 bounds, not absolute-source bounds."""
    h, nu, Xv, leading_u5 = map(rat, (h, nu, Xv, leading_u5))
    if nu < 0 or Xv <= 0 or leading_u5 < 0:
        raise ValueError("Nonnegative viscosity/U bound and positive support required")
    A = F(1, 2) + h
    D = F(1, 2) - h
    R = reciprocal_cap(h)
    v4 = regular_v4_bound(h, leading_u5)["v4"]
    volume0 = Xv
    volume1 = Xv * Xv / 2
    q0 = (
        2 * nu * v4
        + R
        * volume0
        * (8 * D * v4 + (16 + 4 * A) * leading_u5 * v4)
        + volume0 * v4 * v4 / 2
    )
    q1 = (
        R
        * volume1
        * ((1 + 8 * D) * v4 + (16 + 4 * D) * leading_u5 * v4)
        + volume1 * v4 * v4 / 2
    )
    return {
        "reciprocal_cap": R,
        "v4": v4,
        "I0_Ceta3": q0,
        "I1_Ceta3": q1,
        "signed_full_support_only": True,
    }


def first_repair_moments_c3(xcut, seed, f0, u0, I0, I1):
    """Five ordinary first-repair discrepancy bounds in the C_eta^3 algebra."""
    xcut, seed, f0, u0, I0, I1 = map(
        rat, (xcut, seed, f0, u0, I0, I1)
    )
    if xcut <= 0 or min(seed, f0, u0, I0, I1) < 0:
        raise ValueError("Positive cutoff and nonnegative norm bounds required")
    return (
        xcut * seed,
        xcut * xcut * seed,
        2 * f0 * xcut * seed + I0 / 2,
        xcut * xcut * (u0 + f0) * seed,
        2 * u0 * xcut * seed + f0 * xcut * xcut * seed + I1 / 2,
    )


def linear_repair_c3_bound(rho, lam, geometric_K, e_star, moments):
    """Exact inverse-800 first coefficient repair bound at angular order three."""
    rho, lam, geometric_K, e_star = map(
        rat, (rho, lam, geometric_K, e_star)
    )
    moments = tuple(map(rat, moments))
    if (
        rho <= 0
        or lam <= 0
        or geometric_K <= 0
        or e_star <= 0
        or len(moments) != 5
        or min(moments) < 0
    ):
        raise ValueError("Positive geometry/amplitude and five nonnegative moments required")
    mu1, mu2, mu3, mu4, mu5 = moments
    W3 = 5 / e_star
    rows = (
        mu1 / rho,
        (mu1 / rho + geometric_K * W3 * mu4) / lam,
        geometric_K * mu2,
        W3 * mu5 / rho,
        W3 * mu3,
    )
    return {"W3": W3, "normalized_rows": rows, "coefficient_Ceta3": 800 * max(rows)}


def local_c3_rows(c2_rows, delta):
    """Upgrade the established whole-inner strip join from C_eta^2 to C_eta^3."""
    delta = rat(delta)
    rows = tuple(tuple(map(rat, row)) for row in c2_rows)
    if not 0 < delta <= 1 or len(rows) != 4 or any(len(row) != 3 for row in rows):
        raise ValueError("Four established triple rows and a positive strip required")
    c2 = 1 + 2 / delta + 4 / delta**2
    c3 = c2 + 8 / delta**3
    return {
        "Ceta2_factor": c2,
        "Ceta3_factor": c3,
        "radial_Ceta3_bounds": tuple(
            tuple(c3 * value / c2 for value in row) for row in rows
        ),
    }


def cutoff_product(rows, cutoff_jets):
    """Eta-independent radial cutoff propagation through physical order three."""
    rows = tuple(tuple(map(rat, row)) for row in rows)
    cutoff_jets = tuple(map(rat, cutoff_jets))
    if (
        len(rows) != 4
        or any(len(row) != 3 or min(row) < 0 for row in rows)
        or len(cutoff_jets) != 4
        or min(cutoff_jets) < 0
    ):
        raise ValueError("Four complete nonnegative rows and cutoff jets required")
    return tuple(
        tuple(
            sum(
                F(comb(r, j)) * cutoff_jets[j] * rows[r - j][component]
                for j in range(r + 1)
            )
            for component in range(3)
        )
        for r in range(4)
    )


def c4_not_c5_countercontrol(n):
    """Bounds for n^-4 sin(n eta): uniform C4, fifth jet n/120 at eta=0."""
    if type(n) is not int or n < 1:
        raise ValueError("Positive integer frequency required")
    c4 = sum(F(n) ** (j - 4) / factorial(j) for j in range(5))
    return {
        "n": n,
        "Ceta4_upper": c4,
        "uniform_Ceta4_cap": F(65, 24),
        "normalized_fifth_at_zero": F(n, 120),
    }


def review():
    ladder = [derivative_ladder(order) for order in range(6)]
    source = signed_source_c3_bounds(F(1, 200), 1, 7, 3)
    moments = first_repair_moments_c3(11, 13, 17, 19, source["I0_Ceta3"], source["I1_Ceta3"])
    repair = linear_repair_c3_bound(2, F(1, 100), 23, 29, moments)
    c2_rows = ((1, 2, 3), (2, 3, 5), (3, 5, 7), (5, 7, 11))
    local = local_c3_rows(c2_rows, F(1, 4))
    cut = cutoff_product(local["radial_Ceta3_bounds"], (1, 729, 2916, 58320))
    counter = [c4_not_c5_countercontrol(n) for n in (1, 8, 64, 512)]
    return serial(
        {
            "derivative_ladder": ladder,
            "source_fixture": source,
            "repair_moments_fixture": moments,
            "repair_fixture": repair,
            "local_fixture": local,
            "cutoff_fixture": cut,
            "countercontrol": counter,
            "checks": {
                "signed_Ceta3_moments_require_leading_U_Ceta5": ladder[3]["leading_U_order"] == 5,
                "historical_Ceta1_count_requires_leading_U_Ceta3": ladder[1]["leading_U_order"] == 3,
                "reciprocal_neumann_cap_below_25_over_23": reciprocal_cap(F(1, 100)) == F(25, 23),
                "regular_v4_cost_below_29": regular_v4_bound(F(1, 100), 1)["factor"] < 29,
                "axis_viscosity_retained_in_I0": source["I0_Ceta3"] > signed_source_c3_bounds(F(1, 200), 0, 7, 3)["I0_Ceta3"],
                "weighted_radial_viscosity_has_no_axis_term": source["I1_Ceta3"] == signed_source_c3_bounds(F(1, 200), 0, 7, 3)["I1_Ceta3"],
                "linear_repair_uses_W3_five_over_e_star": repair["W3"] == F(5, 29),
                "local_third_angular_loss_strictly_exceeds_second": local["Ceta3_factor"] > local["Ceta2_factor"],
                "third_cutoff_jet_is_charged": cut[3][0] > cutoff_product(local["radial_Ceta3_bounds"], (1, 729, 2916, 0))[3][0],
                "C4_countercontrol_is_uniform": all(F(row["Ceta4_upper"]) <= F(65, 24) for row in serial(counter)),
                "C5_countercontrol_grows": counter[-1]["normalized_fifth_at_zero"] > counter[0]["normalized_fifth_at_zero"],
            },
        }
    )


if __name__ == "__main__":
    import json

    print(json.dumps(review(), indent=2, sort_keys=True))
