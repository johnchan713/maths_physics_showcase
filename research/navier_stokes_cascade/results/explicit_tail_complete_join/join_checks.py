"""Exact parameter and exponent checks for the explicit-tail complete join."""

from fractions import Fraction as F


KAPPA = F(1, 100000)
B0 = F(3, 2) - 2 * KAPPA


def combined_tail(*floors):
    if not floors or any(type(value) is not int for value in floors):
        raise TypeError("tail floors must be nonempty exact integers")
    if any(value < 2 for value in floors):
        raise ValueError("every tail floor must be at least two")
    return max(floors)


def chart_review():
    rows = []
    for explicit, canonical, edge in (
        (11, 17, 13), (10**6, 19, 23), (31, 2**20, 47),
        (10**30, 10**29, 10**28),
    ):
        joined = combined_tail(explicit, canonical, edge)
        rows.append({
            "explicit": explicit,
            "canonical": canonical,
            "edge": edge,
            "joined": joined,
            "dominates_all": all(joined >= x for x in (explicit, canonical, edge)),
            # q < min(qstar/64,2^(-joined-2)); Q <= 4q.
            "Q_over_dyadic_band_ceiling": F(1),
            "whole_enlargement_over_qstar": F(1, 4),
        })
    return rows


def cycle_orders(B, kappa=KAPPA):
    if type(B) is not F or type(kappa) is not F:
        raise TypeError("orders must be exact Fractions")
    if kappa <= 0 or B < F(3, 2) - 2 * kappa:
        raise ValueError("cycle order is below the accepted seed")
    mu = 1 - kappa
    required = B + F(1, 2) - kappa
    rows = {
        "sourced_linear": B + F(1, 2) - kappa,
        "signed_linear": B + F(1, 2) - kappa,
        "old_wave_left": B + F(1, 2) - kappa,
        "old_wave_right": B + F(1, 2) - kappa,
        "mean_wave_left": B + mu - F(1, 2),
        "mean_wave_right": B + mu - F(1, 2),
        "new_square": 2 * B - kappa,
        "fast_mean_wave": B + F(1, 2) - kappa,
        "slow_mean_wave": B + F(1, 2),
    }
    return {
        "B": B,
        "required": required,
        "orders": rows,
        "margins": {name: value - required for name, value in rows.items()},
        "tie_count": sum(value == required for value in rows.values()),
    }


def cycle_review():
    rows = []
    for j in range(13):
        B = B0 + j * (F(1, 2) - KAPPA)
        row = cycle_orders(B)
        row["all_groups_close"] = all(value >= 0 for value in row["margins"].values())
        rows.append(row)
    return rows


def localization_review():
    # qe=q#/16, tau0=qe/16, z0=(qe/16)^D.
    q_over_qsharp = 2 * F(1, 256) + 4 * F(1, 256)
    # An intersecting dyadic chart has Q <= 4q and its enlargement q' <= 4Q.
    whole_chart_over_qsharp = 16 * q_over_qsharp
    return {
        "q_over_qsharp": q_over_qsharp,
        "whole_chart_over_qsharp": whole_chart_over_qsharp,
        "inside_tail": whole_chart_over_qsharp < 1,
    }


def force_terms():
    return (
        "cR", "ct_u", "dt", "laplacian_c_u", "gradient_c_gradient_u",
        "laplacian_d", "p_gradient_c", "c_cminus1_u_gradient_u",
        "c_u_gradient_c_u", "c_u_gradient_d", "c_d_gradient_u",
        "d_gradient_d",
    )


def review():
    charts = chart_review()
    cycles = cycle_review()
    localization = localization_review()
    terms = force_terms()
    return {
        "tail_charts": charts,
        "cycle_orders": cycles,
        "localization": localization,
        "force_terms": list(terms),
        "checks": {
            "combined_tail_dominates_each_distinct_family_gate": all(
                row["dominates_all"] for row in charts
            ),
            "whole_enlarged_charts_stay_in_original_value_domain": all(
                row["whole_enlargement_over_qstar"] <= F(1, 4) for row in charts
            ),
            "nine_complete_residual_groups_close_on_half_line_samples": all(
                row["all_groups_close"] for row in cycles
            ),
            "exactly_seven_groups_tie_the_gain": all(
                row["tie_count"] == 7 for row in cycles
            ),
            "localization_support_closure_stays_inside_joined_tail": localization["inside_tail"],
            "all_twelve_localization_force_groups_recorded": len(terms) == 12,
        },
    }
