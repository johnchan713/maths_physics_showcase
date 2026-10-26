"""Exact controls for the compact/slow weighted-domain join."""
from fractions import Fraction as Q
from math import ceil, comb, factorial, floor


NODES = (
    "fixed_profile_projectors",
    "weighted_initial_rows",
    "native_haar_primitives",
    "shifted_compact_maps",
    "zero_haar_temporal_inverse",
    "fixed_slow_response",
    "first_real_refresh",
    "whole_support_extension",
    "compact_slow_domain_contract",
)
EDGES = (
    ("fixed_profile_projectors", "native_haar_primitives"),
    ("weighted_initial_rows", "native_haar_primitives"),
    ("weighted_initial_rows", "shifted_compact_maps"),
    ("native_haar_primitives", "first_real_refresh"),
    ("shifted_compact_maps", "first_real_refresh"),
    ("zero_haar_temporal_inverse", "first_real_refresh"),
    ("fixed_slow_response", "first_real_refresh"),
    ("first_real_refresh", "compact_slow_domain_contract"),
    ("whole_support_extension", "compact_slow_domain_contract"),
)


def topological_order(nodes=NODES, edges=EDGES):
    incoming = {node: 0 for node in nodes}
    outgoing = {node: [] for node in nodes}
    for left, right in edges:
        if left not in incoming or right not in incoming:
            raise ValueError("Unknown dependency node")
        incoming[right] += 1
        outgoing[left].append(right)
    ready = [node for node in nodes if incoming[node] == 0]
    result = []
    while ready:
        node = ready.pop(0)
        result.append(node)
        for right in outgoing[node]:
            incoming[right] -= 1
            if incoming[right] == 0:
                ready.append(right)
    if len(result) != len(nodes):
        raise ValueError("Dependency cycle")
    return result


def exact_positive(value):
    if not isinstance(value, (int, Q)) or value <= 0:
        raise ValueError("Positive exact datum required")
    return Q(value)


def symbolic_two_edge(pole, width_squared, shell_length, plateau_gap):
    """Return positive coefficient/exponent pairs; never materialize exp."""
    if type(pole) is not int or pole < 0:
        raise ValueError("Nonnegative integer pole required")
    a = exact_positive(width_squared)
    length = exact_positive(shell_length)
    gap = exact_positive(plateau_gap)
    if 2 * gap > length:
        raise ValueError("Plateau gap must fit both shell collars")
    r = min(Q(1), length / 4, gap)
    q = max(0, ceil(Q(pole - 3, 2)))

    def d(c):
        return sum(Q(comb(q, j) * factorial(j)) * (r * r / c) ** j for j in range(q + 1)) / (2 * c)

    # Safe positive symbolic branches of the retained two-edge theorem.
    edge = max(1, factorial(ceil(Q(pole, 2))) * (Q(2) / a) ** ceil(Q(pole, 2)))
    return (
        (d(a), Q(4) / (length - r) ** 2),
        (d(Q(4)), a / (length - r) ** 2),
        (length * edge, (a + 4) / (r * r)),
    )


def fast_request_ledger(output_order):
    if type(output_order) is not int or output_order < 0:
        raise ValueError("Nonnegative output order required")
    return {
        "temporal_tangential_source": output_order + 4,
        "temporal_axial_source": output_order + 4,
        "radial_potential_source": output_order + 5,
        "postfast_scalar": output_order + 8,
        "whole_axial_pressure_moment": output_order + 9,
    }


def slow_scale_powers(h=Q(1, 200)):
    h = exact_positive(h)
    a = Q(1, 2) + h
    return {"u_from_P": a, "u_from_Jz": a - 1, "gamma": a - Q(3, 2), "psi": a - 1}


def first_real_request_ledger(output_order=0, kappa=Q(1, 100000), nested_p=1, quotient_allowance=12, old_pole=30):
    if type(output_order) is not int or output_order < 0:
        raise ValueError("Nonnegative output order required")
    kappa = exact_positive(kappa)
    if not all(type(v) is int and v >= 0 for v in (nested_p, quotient_allowance, old_pole)):
        raise ValueError("Nonnegative integer ledger inputs required")
    b = 1 - kappa
    alpha = 2 - 3 * kappa
    source = alpha - kappa
    target = 3 - 6 * kappa
    p1 = floor((target - source) / kappa) + 1
    p0 = floor((target - b) / kappa) + 1
    a1 = output_order + 2 * p1 + 3
    a0 = output_order + 2 * p0 + 3
    prefast = a0 + 10
    covariance = a0 + 13
    target_input = a0 + 14
    common_pole = max(old_pole, 6 * nested_p, 6 * covariance + 9, 6 * covariance + quotient_allowance + 6)
    return {
        "B": b,
        "alpha": alpha,
        "D": source,
        "nu": target,
        "p1": p1,
        "p0": p0,
        "A1_raw_source": a1,
        "A0_raw_source": a0,
        "max_new_prefast_row": prefast,
        "max_new_covariance": covariance,
        "max_signed_target_input": target_input,
        "common_pole": common_pole,
    }


def strict_margin_review(kappa=Q(1, 100000)):
    row = first_real_request_ledger(kappa=kappa)
    source_margin = row["D"] - row["D"]
    viscosity_order = row["D"] + 1 - 2 * kappa
    viscosity_margin = viscosity_order - row["nu"]
    return {
        "source_margin": source_margin,
        "viscosity_margin": viscosity_margin,
        "strict_without_extra_gain": source_margin + viscosity_margin > 0,
        "required_inequality": "delta_source + delta_viscosity > 0",
    }


def pressure_product_review():
    # c=1+Z, P=2+3Z+Z^2: compare complete derivative with the stale-c shortcut.
    z = Q(2, 5)
    c, cp = 1 + z, 1
    p, pp = 2 + 3 * z + z * z, 3 + 2 * z
    return {"complete": cp * p + c * pp, "dropped_moving_compensator": c * pp}


def review():
    order = topological_order()
    kernel = symbolic_two_edge(12, Q(1, 100), Q(8), Q(1))
    ledger = first_real_request_ledger()
    pressure = pressure_product_review()
    powers = slow_scale_powers()
    checks = {
        "dependency_route_is_acyclic": set(order) == set(NODES),
        "same_choice_precedes_native_primitives": order.index("fixed_profile_projectors") < order.index("native_haar_primitives"),
        "all_symbolic_two_edge_branches_positive": all(c > 0 and e >= 0 for c, e in kernel),
        "symbolic_exponentials_are_not_materialized": all(isinstance(c, Q) and isinstance(e, Q) for c, e in kernel),
        "fast_map_retains_distinct_m_plus_4_and_m_plus_5_inputs": fast_request_ledger(7)["radial_potential_source"] == 12,
        "whole_pressure_moment_uses_m_plus_9": fast_request_ledger(7)["whole_axial_pressure_moment"] == 16,
        "slow_scale_powers_are_distinct": len(set(powers.values())) == 3 and powers["u_from_Jz"] == powers["psi"],
        "first_real_A1_request": ledger["A1_raw_source"] == 200001,
        "first_real_A0_request": ledger["A0_raw_source"] == 399995,
        "first_real_target_request": ledger["max_signed_target_input"] == 400009,
        "first_real_common_pole": ledger["common_pole"] == 2400066,
        "radial_viscosity_has_zero_margin": strict_margin_review()["viscosity_margin"] == 0,
        "generic_positive_route_has_no_strict_gain": not strict_margin_review()["strict_without_extra_gain"],
        "moving_compensator_derivative_is_nonzero": pressure["complete"] != pressure["dropped_moving_compensator"],
    }
    return {
        "checks": checks,
        "symbolic_kernel": [[str(c), str(e)] for c, e in kernel],
        "first_real_ledger": {k: str(v) for k, v in ledger.items()},
        "strict_margin": {k: str(v) for k, v in strict_margin_review().items()},
        "scope": "finite rational and symbolic controls; continuum mapping follows from the cited weighted source, fixed-profile, shifted-map and support estimates",
    }
