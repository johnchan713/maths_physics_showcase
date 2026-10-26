#!/usr/bin/env python3
"""Exact algebra controls for the canonical full-path ODE join."""
from fractions import Fraction as F
from math import factorial


NODES = (
    "same_background_tail", "full_phase_coefficients", "normal_floor",
    "tangent_propagator", "canonical_entrance", "homogeneous_primary",
    "sourced_inverse", "cut_potential_pressure", "full_path_contract",
)
EDGES = (
    ("same_background_tail", "full_phase_coefficients"),
    ("same_background_tail", "normal_floor"),
    ("full_phase_coefficients", "tangent_propagator"),
    ("normal_floor", "tangent_propagator"),
    ("tangent_propagator", "canonical_entrance"),
    ("canonical_entrance", "homogeneous_primary"),
    ("full_phase_coefficients", "sourced_inverse"),
    ("normal_floor", "sourced_inverse"),
    ("tangent_propagator", "sourced_inverse"),
    ("homogeneous_primary", "cut_potential_pressure"),
    ("sourced_inverse", "cut_potential_pressure"),
    ("cut_potential_pressure", "full_path_contract"),
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


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


def matvec(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def add(*vectors):
    return tuple(sum(values, F(0)) for values in zip(*vectors))


def scale(value, vector):
    return tuple(value * entry for entry in vector)


def projected_ode(n, nprime, matrix, forcing, tangent, k, mode, damping):
    """Evaluate the exact projected equation and pressure coefficient pi/i."""
    n = tuple(map(F, n)); nprime = tuple(map(F, nprime))
    matrix = tuple(tuple(map(F, row)) for row in matrix)
    forcing = tuple(map(F, forcing)); tangent = tuple(map(F, tangent))
    k, mode, damping = F(k), F(mode), F(damping)
    norm2 = dot(n, n)
    if norm2 <= 0 or k == 0 or mode == 0 or dot(n, tangent) != 0:
        raise ValueError("Nonzero carrier, normal and tangent datum required")
    kt = matvec(matrix, tangent)
    pressure_over_i = (dot(n, kt) - dot(nprime, tangent) + dot(n, forcing)) / (k * mode * norm2)
    projected_forcing = add(forcing, scale(-dot(n, forcing) / norm2, n))
    correction = scale((dot(n, kt) - dot(nprime, tangent)) / norm2, n)
    tangent_prime = add(scale(-1, kt), correction,
                        scale(-(mode * mode) * damping, tangent),
                        scale(-1, projected_forcing))
    # In the original equation i*k*m*n*pi = -k*m*n*(pi/i).
    original_residual = add(tangent_prime, kt,
                            scale((mode * mode) * damping, tangent),
                            scale(-k * mode * pressure_over_i, n), forcing)
    return {
        "pressure_over_i": pressure_over_i,
        "constraint_derivative": dot(nprime, tangent) + dot(n, tangent_prime),
        "original_residual": original_residual,
        "projected_forcing": projected_forcing,
    }


def ambient_lift(n, nprime, matrix, damping, profile_rate):
    """Check the normal row identity of the Gaussian-adapted ambient lift."""
    n = tuple(map(F, n)); nprime = tuple(map(F, nprime))
    matrix = tuple(tuple(map(F, row)) for row in matrix)
    damping, profile_rate = F(damping), F(profile_rate)
    norm2 = dot(n, n)
    h = scale(1 / norm2, n)
    n_k = tuple(sum(n[r] * matrix[r][c] for r in range(3)) for c in range(3))
    a0 = []
    for r in range(3):
        a0.append(tuple(-matrix[r][c] + h[r] * (n_k[c] - nprime[c]) for c in range(3)))
    a = tuple(tuple(a0[r][c] - (damping if r == c else 0) for c in range(3)) for r in range(3))
    ahat = tuple(tuple(a[r][c] + (damping + profile_rate) * h[r] * n[c]
                                for c in range(3)) for r in range(3))
    normal_row = tuple(nprime[c] + sum(n[r] * ahat[r][c] for r in range(3)) for c in range(3))
    return {
        "normal_row": normal_row,
        "expected_normal_row": scale(profile_rate, n),
        "tangent_agreement": matvec(ahat, (n[1], -n[0], F(0))) == matvec(a, (n[1], -n[0], F(0))),
    }


def product_jets(left, right):
    """Raw one-variable Leibniz jets."""
    if len(left) != len(right):
        raise ValueError("Equal jet lengths required")
    return tuple(sum(F(factorial(n), factorial(j) * factorial(n-j)) * F(left[j]) * F(right[n-j])
                     for j in range(n + 1)) for n in range(len(left)))


def cutoff_review():
    cutoff = (F(0), F(2), F(-3), F(5))
    potential = (F(7), F(11), F(13), F(17))
    pressure = (F(19), F(23), F(29), F(31))
    cut_potential = product_jets(cutoff, potential)
    cut_pressure = product_jets(cutoff, pressure)
    return {
        "cut_potential": cut_potential,
        "cut_pressure": cut_pressure,
        "potential_derivative_commutator": cut_potential[1] - cutoff[0] * potential[1],
        "pressure_derivative_commutator": cut_pressure[1] - cutoff[0] * pressure[1],
    }


def series_mul(left, right, order):
    return [sum(left[j] * right[n-j] for j in range(n + 1)) for n in range(order + 1)]


def inverse_sqrt_series(a, b0, order):
    """Coefficients b of a*b^2=1 when the exact positive b0 is supplied."""
    a = list(map(F, a)) + [F(0)] * max(0, order + 1 - len(a))
    b = [F(b0)]
    if a[0] * b[0] * b[0] != 1:
        raise ValueError("Exact positive inverse square root required")
    for n in range(1, order + 1):
        trial = b + [F(0)]
        known = series_mul(a, series_mul(trial, trial, n), n)[n]
        b.append(-known / (2 * a[0] * b[0]))
    return b


def varying_growth_review():
    # Expand at v=L/2 with u=3/4,L=1: s=u(1+x), so 1+s^2 has square value 25/16 at x=0.
    u = F(3, 4)
    a = [1 + u*u, 2*u*u, u*u, F(0), F(0)]
    inv = inverse_sqrt_series(a, F(4, 5), 4)
    damping_factor = F(64, 125)  # (1+u^2)^(-3/2)
    p = [inv[j] - damping_factor * a[j] for j in range(5)]
    jets = [p[j] * factorial(j) for j in range(5)]
    frozen = [jets[0]] + [F(0)] * 4
    return {
        "profile_rate_jets": jets,
        "third_or_higher_nonzero": any(value != 0 for value in jets[3:]),
        "frozen_rate_misses_jets": jets != frozen,
        "inverse_identity": series_mul(a, series_mul(inv, inv, 4), 4) == [F(1), F(0), F(0), F(0), F(0)],
    }


def review():
    order = topological_order()
    cycle = False
    try:
        topological_order(edges=EDGES + (("full_path_contract", "same_background_tail"),))
    except ValueError:
        cycle = True
    ode = projected_ode(
        (2, 1, 1), (1, -1, 0),
        ((0, -2, 1), (3, 0, -1), (1, 2, 0)),
        (3, -2, 5), (1, -2, 0), 7, -2, F(3, 5),
    )
    pure_normal = projected_ode(
        (0, 0, 1), (0, 0, 0), ((0, 0, 0),) * 3,
        (0, 0, 2), (1, 0, 0), 5, 1, 1,
    )
    ambient = ambient_lift(
        (2, 1, 1), (1, -1, 0),
        ((0, -2, 1), (3, 0, -1), (1, 2, 0)), F(3, 5), F(4, 7),
    )
    cut = cutoff_review()
    growth = varying_growth_review()
    checks = {
        "provider_route_is_acyclic": len(order) == len(NODES),
        "downstream_to_background_cycle_is_rejected": cycle,
        "projected_equation_preserves_moving_constraint": ode["constraint_derivative"] == 0,
        "projected_and_original_equations_agree_exactly": all(value == 0 for value in ode["original_residual"]),
        "pure_normal_source_has_zero_projected_velocity_source": all(value == 0 for value in pure_normal["projected_forcing"]),
        "pure_normal_source_pressure_is_nonzero": pure_normal["pressure_over_i"] != 0,
        "ambient_lift_has_exact_profile_normal_response": ambient["normal_row"] == ambient["expected_normal_row"],
        "ambient_lift_agrees_on_tangent_vectors": ambient["tangent_agreement"],
        "potential_cutoff_derivative_is_retained": cut["potential_derivative_commutator"] != 0,
        "pressure_cutoff_derivative_is_retained": cut["pressure_derivative_commutator"] != 0,
        "varying_growth_inverse_identity_is_exact": growth["inverse_identity"],
        "varying_growth_has_nonzero_higher_jets": growth["third_or_higher_nonzero"] and growth["frozen_rate_misses_jets"],
    }
    def canon(value):
        if isinstance(value, F): return str(value)
        if isinstance(value, dict): return {key: canon(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)): return [canon(item) for item in value]
        return value
    return {
        "route": order,
        "controls": canon({"projected_ode": ode, "pure_normal_source": pure_normal,
                           "ambient_lift": ambient, "cutoffs": cut, "varying_growth": growth}),
        "checks": checks,
    }
