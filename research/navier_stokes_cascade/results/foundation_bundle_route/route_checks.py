"""Exact dependency and common-tail checks for the foundation route."""

from fractions import Fraction as F


NODES = (
    "core_contraction",
    "core_growth",
    "continuation",
    "reference_roots",
    "axial_root",
    "matching_roots",
    "modulation_repair",
    "annular_exterior",
    "regional_cone",
    "global_pressure",
    "weighted_stress",
    "background_coefficients",
    "global_coefficients",
    "explicit_primary_tail",
    "canonical_source_cycle",
    "complete_residual_force",
)


EDGES = (
    ("core_contraction", "core_growth"),
    ("core_growth", "continuation"),
    ("continuation", "matching_roots"),
    ("reference_roots", "axial_root"),
    ("axial_root", "matching_roots"),
    ("continuation", "modulation_repair"),
    ("matching_roots", "modulation_repair"),
    ("matching_roots", "annular_exterior"),
    ("modulation_repair", "annular_exterior"),
    ("core_growth", "regional_cone"),
    ("modulation_repair", "regional_cone"),
    ("annular_exterior", "regional_cone"),
    ("reference_roots", "global_pressure"),
    ("axial_root", "global_pressure"),
    ("matching_roots", "global_pressure"),
    ("regional_cone", "global_pressure"),
    ("global_pressure", "weighted_stress"),
    ("annular_exterior", "weighted_stress"),
    ("weighted_stress", "background_coefficients"),
    ("background_coefficients", "global_coefficients"),
    ("regional_cone", "explicit_primary_tail"),
    ("global_coefficients", "explicit_primary_tail"),
    ("explicit_primary_tail", "canonical_source_cycle"),
    ("weighted_stress", "canonical_source_cycle"),
    ("global_pressure", "canonical_source_cycle"),
    ("canonical_source_cycle", "complete_residual_force"),
)


def topological_order(nodes=NODES, edges=EDGES):
    """Return a stable topological order, rejecting missing nodes and cycles."""
    if len(set(nodes)) != len(nodes):
        raise ValueError("duplicate node")
    known = set(nodes)
    if any(a not in known or b not in known for a, b in edges):
        raise ValueError("edge has unknown endpoint")
    incoming = {node: set() for node in nodes}
    outgoing = {node: set() for node in nodes}
    for source, target in edges:
        incoming[target].add(source)
        outgoing[source].add(target)
    ready = [node for node in nodes if not incoming[node]]
    order = []
    while ready:
        node = ready.pop(0)
        order.append(node)
        for target in tuple(outgoing[node]):
            incoming[target].remove(node)
            if not incoming[target] and target not in order and target not in ready:
                ready.append(target)
        ready.sort(key=nodes.index)
    if len(order) != len(nodes):
        raise ValueError("dependency cycle")
    return tuple(order)


def combined_tail(*floors):
    """A finite conjunction of later-band gates is their maximum."""
    if not floors:
        raise ValueError("at least one floor required")
    if any(type(value) is not int or value < 0 for value in floors):
        raise TypeError("floors must be nonnegative integers")
    return max(floors)


def chart_domain(q_over_qstar=F(1, 64)):
    """Worst enlarged-chart ratios used by the certified later tail."""
    if type(q_over_qstar) is not F or q_over_qstar <= 0:
        raise TypeError("positive Fraction required")
    whole_band = 16 * q_over_qstar
    localized_q = F(3, 128)
    localized_chart = 16 * localized_q
    return {
        "unlocalized_whole_band": whole_band,
        "localized_q": localized_q,
        "localized_whole_chart": localized_chart,
    }


def review():
    order = topological_order()
    position = {node: index for index, node in enumerate(order)}
    domain = chart_domain()
    return {
        "node_count": len(NODES),
        "edge_count": len(EDGES),
        "topological_order": order,
        "all_edges_forward": all(position[a] < position[b] for a, b in EDGES),
        "tail_max_control": combined_tail(16, 31, 17, 29) == 31,
        "domain": domain,
        "domain_inside": all(value < 1 for value in (
            domain["unlocalized_whole_band"],
            domain["localized_q"],
            domain["localized_whole_chart"],
        )),
    }
