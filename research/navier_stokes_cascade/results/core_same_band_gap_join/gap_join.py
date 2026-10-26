#!/usr/bin/env python3
"""Exact algebra controls for the same-background accepted-tail gap join."""
from fractions import Fraction as F


NODES = (
    "leading_regional_gaps", "rule_defined_background", "finite_head",
    "weighted_stress", "accepted_tail", "current_value_cone",
    "full_phase", "reference_normal_floor", "same_band_contract",
)
EDGES = (
    ("leading_regional_gaps", "rule_defined_background"),
    ("rule_defined_background", "finite_head"),
    ("finite_head", "weighted_stress"),
    ("weighted_stress", "accepted_tail"),
    ("leading_regional_gaps", "current_value_cone"),
    ("accepted_tail", "current_value_cone"),
    ("rule_defined_background", "full_phase"),
    ("accepted_tail", "full_phase"),
    ("full_phase", "reference_normal_floor"),
    ("current_value_cone", "same_band_contract"),
    ("reference_normal_floor", "same_band_contract"),
)


def topological_order(nodes=NODES, edges=EDGES):
    incoming = {node: 0 for node in nodes}
    outgoing = {node: [] for node in nodes}
    if len(incoming) != len(nodes):
        raise ValueError("Duplicate node")
    for left, right in edges:
        if left not in incoming or right not in incoming:
            raise ValueError("Unknown node")
        incoming[right] += 1
        outgoing[left].append(right)
    ready = sorted(node for node, degree in incoming.items() if degree == 0)
    order = []
    while ready:
        node = ready.pop(0)
        order.append(node)
        for right in sorted(outgoing[node]):
            incoming[right] -= 1
            if incoming[right] == 0:
                ready.append(right)
                ready.sort()
    if len(order) != len(nodes):
        raise ValueError("Cycle")
    return order


def background_floors(f, fbar, t, d, gbar, omega, gamma):
    """Denominator-free floors after the certified value perturbation."""
    f, fbar, t, d, gbar, omega, gamma = map(F, (f, fbar, t, d, gbar, omega, gamma))
    if min(f, t, d, omega, gamma) <= 0 or min(fbar, gbar) < 0 or gamma > 1:
        raise ValueError("Invalid leading gap data")
    return {
        "F": f / 2,
        "minus_gtheta": t / 2,
        "D": d / 2,
        "lambda_squared": f * t * d / (4 * (gbar + 1) ** 2),
        "c_squared": d / (4 * (fbar + 1) * (gbar + 1)),
        "orientation": omega / 3,
        "cone_margin": omega ** 2 * gamma * (2 - gamma) / 16,
    }


def accepted_tail_margins(ell, h, log2_value_ratio, r0, length_required, epsilon, k_inverse):
    """Check supplied exact rows of the existing one-integer tail selector."""
    if type(ell) is not int or ell < 3:
        raise ValueError("Integer band at least three required")
    h, log2_value_ratio, r0, length_required, epsilon, k_inverse = map(
        F, (h, log2_value_ratio, r0, length_required, epsilon, k_inverse)
    )
    if min(h, r0, length_required, epsilon, k_inverse) <= 0 or h >= F(1, 100):
        raise ValueError("Invalid accepted-tail data")
    return {
        "value_log2": 2 * h * (ell - 2) - log2_value_ratio,
        "length": 2 * r0 * ell ** 2 - length_required,
        "slow": 1 - ell ** 4 * (epsilon + epsilon ** 2 + k_inverse),
    }


def reference_normal_floor(theta, beta):
    """Reference-relative normal floor and the retained coarse comparison."""
    theta, beta = F(theta), F(beta)
    if not 0 <= theta <= F(1, 32) or beta <= 0:
        raise ValueError("Invalid normal-floor data")
    reciprocal = 1 / ((1 - theta) ** 2 * beta)
    return {
        "amplitude_fraction": 1 - theta,
        "reciprocal": reciprocal,
        "explicit_31_over_32": 31 * F(1, 32),
        "coarse_reciprocal": 4 / beta,
        "sharper_than_coarse": reciprocal < 4 / beta,
    }


def chart_tail(qstar, ell_star):
    """The common q-tail keeps every enlarged active chart in the same domain."""
    qstar = F(qstar)
    if qstar <= 0 or qstar > 1 or type(ell_star) is not int or ell_star < 1:
        raise ValueError("Invalid chart-tail data")
    qv = min(qstar / 64, F(1, 2 ** (ell_star + 2)))
    return {
        "qv": qv,
        "largest_Q": 4 * qv,
        "largest_enlarged_q": 16 * qv,
        "inside_qstar": 16 * qv < qstar,
        "at_or_after_selected_band": 4 * qv <= F(1, 2 ** ell_star),
    }


def review():
    order = topological_order()
    cycle = False
    try:
        topological_order(edges=EDGES + (("same_band_contract", "leading_regional_gaps"),))
    except ValueError:
        cycle = True
    floors = background_floors(F(1, 5), 3, F(1, 7), F(1, 11), 4, F(1, 13), F(1, 17))
    ell = 400000
    epsilon = F(1, 2 ** 2000)
    margins = accepted_tail_margins(ell, F(1, 200), 100, F(1, 6057072), 1000,
                                    epsilon, F(1, 2 ** 1000))
    normal = reference_normal_floor(F(1, 32), F(1, 101))
    chart = chart_tail(F(1, 8), 20)
    checks = {
        "provider_route_is_acyclic": len(order) == len(NODES),
        "downstream_to_leading_cycle_is_rejected": cycle,
        "all_background_floors_are_positive": all(value > 0 for value in floors.values()),
        "one_integer_selector_rows_have_strict_margin": all(value > 0 for value in margins.values()),
        "reference_floor_is_at_least_31_over_32": normal["amplitude_fraction"] >= normal["explicit_31_over_32"],
        "reference_reciprocal_is_stronger_than_four_over_beta": normal["sharper_than_coarse"],
        "whole_enlarged_chart_stays_in_same_tail": chart["inside_qstar"] and chart["at_or_after_selected_band"],
    }
    return {
        "route": order,
        "controls": {
            "background_floors": {key: str(value) for key, value in floors.items()},
            "accepted_tail_margins": {key: str(value) for key, value in margins.items()},
            "normal_floor": {key: str(value) if isinstance(value, F) else value for key, value in normal.items()},
            "chart_tail": {key: str(value) if isinstance(value, F) else value for key, value in chart.items()},
        },
        "checks": checks,
    }
