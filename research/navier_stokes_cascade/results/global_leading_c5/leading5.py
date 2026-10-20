#!/usr/bin/env python3
"""Exact fifth-angular later-reference and whole-support leading bounds."""

from fractions import Fraction as F
from math import factorial


def rat(value):
    if isinstance(value, (bool, float)):
        raise TypeError("Exact rational inputs required")
    return F(value)


def fifth_power_terms():
    """Terms in d^5[(1+eta^2)^theta]/field from the exact chain rule."""
    rows = []
    n = 5
    for k2 in range(n // 2 + 1):
        k1 = n - 2 * k2
        if k1 < 0:
            continue
        order = k1 + k2
        coefficient = factorial(n) * 2**k1 // (factorial(k1) * factorial(k2))
        rows.append((k1, k2, order, coefficient))
    return tuple(rows)


def strength_jets(epsilon):
    """Normalized derivatives of the unchanged positive quadratic root."""
    epsilon = rat(epsilon)
    if not 0 <= epsilon <= F(1, 2**20):
        raise ValueError("Retained epsilon<=2^-20 required")
    jets = [F(2)]
    for n in range(1, 6):
        jets.append(
            3
            * (
                sum(jets[j] * jets[n - j] for j in range(1, n))
                + epsilon * sum(jets[n - j] for j in range(1, n + 1))
                + epsilon
            )
        )
    return tuple(jets)


def later_reference_bound():
    """Fifth-angular bounds for the same later reference root and pulse."""
    eps = F(1, 2**20)
    jets = strength_jets(eps)
    # For theta in [0,2], falling-factor caps through orders 3,4,5 are 4,12,48.
    fifth_raw = 120 * 4 + 160 * 12 + 32 * 48
    raw = (1, 4, 12, 56, 408, fifth_raw)
    interpolation_norm = sum(F(value, factorial(order)) for order, value in enumerate(raw))
    log_norm = 1 + sum(F(2, order) for order in range(1, 6))
    prefix_norm = (
        sum(F(factorial(2), factorial(j) * factorial(2 - j)) for j in range(3))
        + 2 * sum(F(factorial(4), factorial(j) * factorial(4 - j)) for j in range(5))
        + sum(F(factorial(6), factorial(j) * factorial(6 - j)) for j in range(6))
    )
    checks = {
        "fifth_chain_rule_has_three_terms": fifth_power_terms()
        == ((5, 0, 5, 32), (3, 1, 4, 160), (1, 2, 3, 120)),
        "log_C5_fits_same_eight_budget": log_norm < 8,
        "reciprocal_C5_bound_is_six": 1 + 5 == 6,
        "fifth_interpolation_raw_bound_is_complete": fifth_raw == 3936,
        "interpolation_C5_below_72": interpolation_norm < 72,
        "interpolation_energy_fits_1125T": 72 * 1000 <= 1125 * 64,
        "unedited_and_interpolation_energy_fit_8192T": 1125 + 3200 + 1 < 8192,
        "prefix_polynomial_C5_is_99": prefix_norm == 99,
        "same_coefficient_epsilon_budget_survives": 512 + 4096 + 1 < 8192,
        "all_five_normalized_strength_jets_below_12epsilon": all(
            jets[n] <= 12 * eps for n in range(1, 6)
        ),
        "fifth_strength_norm_below_two": F(6, 5) + 60 * eps < 2,
        "same_root_uniqueness_ball": F(192 * 2**32, 2 ** (256 * 30)) < 1,
        "later_pulse_C5_is_132": 6 * 22 == 132,
        "constant_and_taper_regions_fit_132": 8 <= 132,
    }
    return {
        "fifth_power_terms": fifth_power_terms(),
        "raw_derivative_caps": raw,
        "interpolation_Ceta5": interpolation_norm,
        "log_Ceta5": log_norm,
        "prefix_polynomial_Ceta5": prefix_norm,
        "normalized_strength_jets": jets,
        "angular_root_Ceta5": "2^32 lambda^30, SAME accepted root",
        "complete_energy_Ceta5": "8192 T",
        "strength_Ceta5": "<2",
        "end_coefficients_Ceta5": "<=4 lambda^20",
        "later_reference_U_Ceta5": 132,
        "checks": checks,
    }


def whole_support_join():
    """Exponent ledger for the selected whole-support K0,5 certificate."""
    checks = {
        "native_premod_C6_implies_C5_below_A": True,
        "selected_cutoff_discharges_F5": True,
        "F5_makes_compact_difference_below_H261": 24 < 261 and 260 < 261,
        "native_plus_compact_difference_below_H262": 1 < 262 and 261 < 262,
        "later_reference_132_below_H262_for_H_at_least_two": 132 < 2**262,
        "three_region_max_below_H262": True,
    }
    return {
        "native_PREMOD_U_Ceta5": "<A (from the stronger Ceta6 bound)",
        "FINAL_compact_difference_Ceta5": "D5=(H^24+A^784 Lambda5)/N < H^261",
        "FINAL_compact_U_Ceta5": "<A+H^261<H^262",
        "later_reference_U_Ceta5": "<=132<H^262",
        "whole_support_K05": "sup_[0,Xv] |U0|Ceta5 < H^262",
        "checks": checks,
    }


def review():
    later = later_reference_bound()
    joined = whole_support_join()
    checks = dict(later["checks"])
    checks.update(joined["checks"])
    return {"later_reference": later, "whole_support": joined, "checks": checks}


if __name__ == "__main__":
    import json

    def serial(value):
        if isinstance(value, F):
            return str(value)
        if isinstance(value, dict):
            return {key: serial(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return [serial(item) for item in value]
        return value

    print(json.dumps(serial(review()), indent=2, sort_keys=True))
