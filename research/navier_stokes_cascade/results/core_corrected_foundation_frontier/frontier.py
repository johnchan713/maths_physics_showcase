#!/usr/bin/env python3
"""Exact dependency/frontier controls for the corrected construction."""
from fractions import Fraction as F


NODES = (
    "corrected_comparison", "pressure_identification", "value_domain",
    "fixed_point", "growth", "activation", "continuation_roots",
    "modulation_repair", "global_pressure", "coefficient_support",
    "primary_tail", "covariance_means", "dynamic_means", "slow_jet",
    "pressure_shape", "fixed_profile", "initial_high_jets",
    "finite_recurrence", "corrected_diagonal",
)

EDGES = (
    ("pressure_identification", "value_domain"),
    ("fixed_point", "growth"), ("growth", "activation"),
    ("activation", "continuation_roots"),
    ("continuation_roots", "modulation_repair"),
    ("modulation_repair", "global_pressure"),
    ("global_pressure", "coefficient_support"),
    ("coefficient_support", "primary_tail"),
    ("primary_tail", "covariance_means"),
    ("covariance_means", "dynamic_means"),
    ("dynamic_means", "slow_jet"),
    ("slow_jet", "pressure_shape"),
    ("pressure_shape", "fixed_profile"),
    ("fixed_profile", "initial_high_jets"),
    ("initial_high_jets", "finite_recurrence"),
    ("finite_recurrence", "corrected_diagonal"),
    ("corrected_diagonal", "corrected_comparison"),
)

PRIMITIVES = (
    "actual analytic leading/background and Borel identification",
    "same-band phase, normal, value and cone gaps",
    "canonical homogeneous ODE with full-path cutoffs and pressure",
    "fixed compensator/projector profiles and common supports",
    "exact compact/slow inverse maps on their declared weighted domains",
)


def topological_order(nodes=NODES, edges=EDGES):
    nodes = tuple(nodes)
    if len(set(nodes)) != len(nodes):
        raise ValueError("Duplicate node")
    incoming = {node: 0 for node in nodes}
    outgoing = {node: [] for node in nodes}
    for source, target in edges:
        if source not in incoming or target not in incoming:
            raise ValueError("Unknown endpoint")
        incoming[target] += 1
        outgoing[source].append(target)
    ready = sorted(node for node, degree in incoming.items() if degree == 0)
    order = []
    while ready:
        source = ready.pop(0)
        order.append(source)
        for target in sorted(outgoing[source]):
            incoming[target] -= 1
            if incoming[target] == 0:
                ready.append(target)
                ready.sort()
    if len(order) != len(nodes):
        raise ValueError("Dependency cycle")
    return order


def pressure_shape_bound(m, epsilon, H, K, J):
    """Coefficient in the fixed-profile pressure-shape estimate."""
    m = int(m)
    epsilon, H, K, J = map(F, (epsilon, H, K, J))
    if m < 0 or min(epsilon, H, K, J) < 0:
        raise ValueError("Nonnegative data required")
    return epsilon**5 * (m + 1) * (m + 2) * (m + 3) * (m + 4) * H * K * J / 2


def finite_conjunction_threshold(bounds):
    """One later integer tail for finitely many fixed lower bounds."""
    bounds = tuple(int(value) for value in bounds)
    if not bounds or min(bounds) < 0:
        raise ValueError("Finite nonnegative bounds required")
    return max(bounds)


def review():
    order = topological_order()
    cycle_rejected = False
    try:
        topological_order(edges=EDGES + (("corrected_comparison", "fixed_point"),))
    except ValueError:
        cycle_rejected = True
    b0 = pressure_shape_bound(0, F(1, 2), 2, 3, 5)
    b1 = pressure_shape_bound(1, F(1, 2), 2, 3, 5)
    return {
        "counts": {"nodes": len(NODES), "edges": len(EDGES), "primitive_contracts": len(PRIMITIVES)},
        "topological_order": order,
        "pressure_shape_controls": {
            "m0": str(b0), "m1": str(b1),
            "m1_over_m0": str(b1 / b0),
            "finite_tail_max": finite_conjunction_threshold((17, 4, 29, 11)),
        },
        "checks": {
            "provider_graph_is_acyclic": len(order) == len(NODES),
            "added_downstream_to_core_edge_is_rejected": cycle_rejected,
            "pressure_shape_bound_keeps_four_derivative_shifts": b1 / b0 == 5,
            "finite_value_gates_share_one_later_tail": finite_conjunction_threshold((17, 4, 29, 11)) == 29,
            "primitive_frontier_is_explicit_and_finite": len(PRIMITIVES) == 5,
        },
        "primitive_contracts": list(PRIMITIVES),
    }
