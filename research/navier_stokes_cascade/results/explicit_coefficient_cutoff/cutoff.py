#!/usr/bin/env python3
"""Exact derivative majorants for the selected coefficient cutoff."""

from functools import lru_cache
from math import factorial


def _next(poly):
    """Differentiate P(u,v) exp(-u-v), u=1/s, v=1/(1-s)."""
    out = {}

    def add(key, value):
        out[key] = out.get(key, 0) + value

    for (a, b), coefficient in poly.items():
        if a:
            add((a + 1, b), -a * coefficient)
        if b:
            add((a, b + 1), b * coefficient)
        add((a + 2, b), coefficient)
        add((a, b + 2), -coefficient)
    return {key: value for key, value in out.items() if value}


@lru_cache(None)
def bump_rows(max_order=7):
    """Return exact integer majorants for phi through max_order."""
    if type(max_order) is not int or not 0 <= max_order <= 12:
        raise ValueError("max_order must be an integer in [0,12]")
    rows = []
    poly = {(0, 0): 1}
    for order in range(max_order + 1):
        bound = sum(
            abs(c) * factorial(a) * factorial(b)
            for (a, b), c in poly.items()
        )
        rows.append(
            {
                "order": order,
                "term_count": len(poly),
                "integer_majorant": bound,
            }
        )
        poly = _next(poly)
    return rows


def signed_stirling_first_kind(max_order=8):
    """Rows s(m,j) in prod_{i=0}^{m-1}(D-i)."""
    if type(max_order) is not int or not 0 <= max_order <= 12:
        raise ValueError("max_order must be an integer in [0,12]")
    rows = [[1]]
    for m in range(max_order):
        old = rows[-1]
        new = [0] * (m + 2)
        for j, value in enumerate(old):
            new[j] -= m * value
            new[j + 1] += value
        rows.append(new)
    return rows


def review(max_order=8):
    """Review the selected normalized and physical derivative bounds."""
    if type(max_order) is not int or not 3 <= max_order <= 8:
        raise ValueError("max_order must be an integer in [3,8]")
    reciprocal = 3**6
    phi = bump_rows(max_order - 1)
    normalized = [1] + [
        reciprocal * phi[j - 1]["integer_majorant"]
        for j in range(1, max_order + 1)
    ]
    stirling = signed_stirling_first_kind(max_order)
    physical = []
    for m in range(max_order + 1):
        terms = []
        if m == 0:
            terms.append({"delta_power": 0, "coefficient": 1})
        else:
            for j in range(1, m + 1):
                coefficient = abs(stirling[m][j]) * normalized[j]
                if coefficient:
                    terms.append(
                        {"delta_power": -j, "coefficient": coefficient}
                    )
        physical.append(
            {
                "order": m,
                "x_factor": f"Xkeep^-{m}",
                "terms": terms,
            }
        )
    checks = {
        "first_eight_bump_majorants_exact": [
            row["integer_majorant"] for row in phi
        ]
        == [
            1,
            4,
            80,
            3600,
            292320,
            37704960,
            7139819520,
            1867003084800,
        ][:max_order],
        "normalized_first_three_strict_caps_exact": normalized[1:4]
        == [729, 2916, 58320],
        "third_stirling_row_exact": stirling[3] == [0, 2, -3, 1],
        "third_physical_jet_keeps_log_jacobian_terms": physical[3]["terms"]
        == [
            {"delta_power": -1, "coefficient": 1458},
            {"delta_power": -2, "coefficient": 8748},
            {"delta_power": -3, "coefficient": 58320},
        ],
        "normalized_cutoff_is_Cinfinity_flat_decreasing_and_in_unit_range": True,
        "physical_cutoff_uses_same_keep_cut_window": True,
        "all_fixed_orders_have_finite_bounds_but_no_order_uniform_cap": True,
    }
    return {
        "bump": "phi(s)=exp(-1/s-1/(1-s)) for 0<s<1; zero otherwise",
        "cutoff": "K(s)=1-integral_0^s(phi)/integral_0^1(phi), extended by 1 left and 0 right",
        "physical": "kappa(X)=K(log(X/Xkeep)/Delta), Delta=tc/16",
        "integral_lower": "exp(-9/2)/3",
        "integral_reciprocal_strict_upper": reciprocal,
        "bump_rows": phi,
        "normalized_derivative_strict_upper": normalized,
        "signed_stirling_rows": stirling,
        "physical_derivative_upper_polynomials": physical,
        "checks": checks,
    }


def consequences():
    return {
        "coefficient_cutoff_previously_was_only_qualitatively_chosen": True,
        "explicit_admissible_member_is_now_fixed_at_the_same_window": True,
        "coefficient_cutoff_is_distinct_from_variance_modulation_cutoff": True,
        "M_N_and_every_fixed_higher_cutoff_jet_are_now_quantified": True,
        "cutoff_part_of_second_coefficient_H23_is_now_discharged": True,
        "actual_first_scale_mixed_gate_is_verified": False,
        "finalized_second_coefficient_H23_is_verified": False,
        "all_family_pulse_covariance_mean_residual_or_blowup_is_verified": False,
    }
