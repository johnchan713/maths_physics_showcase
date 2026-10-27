"""Exact finite controls for the fixed-profile/projector contract."""
from fractions import Fraction as Q
from math import comb, factorial


NODES = (
    "current_same_tail_full_path",
    "fixed_pressure_profile",
    "moving_pressure_compensator",
    "reserved_moving_support",
    "rho_squared_projectors",
    "same_choice_signed_recompute",
    "finite_support_closure",
    "fixed_profile_projector_contract",
)
EDGES = (
    ("current_same_tail_full_path", "fixed_profile_projector_contract"),
    ("fixed_pressure_profile", "moving_pressure_compensator"),
    ("fixed_pressure_profile", "rho_squared_projectors"),
    ("moving_pressure_compensator", "rho_squared_projectors"),
    ("reserved_moving_support", "rho_squared_projectors"),
    ("rho_squared_projectors", "same_choice_signed_recompute"),
    ("same_choice_signed_recompute", "fixed_profile_projector_contract"),
    ("finite_support_closure", "fixed_profile_projector_contract"),
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


def beta(a, b):
    """Exact beta integral for positive integer arguments."""
    return Q(factorial(a - 1) * factorial(b - 1), factorial(a + b - 1))


def diagnostic_profile():
    """Finite-regularity polynomial diagnostic on [1,2], not the C-infinity profile."""
    mass = 30 * beta(3, 3)
    square_mass = 30**2 * beta(5, 5)
    return {"mass": mass, "A2": square_mass, "length": Q(1)}


def dilation_moment(s, e):
    """Exact factors after R=sqrt(s)x for rational-square s."""
    if s not in (Q(1, 4), Q(1), Q(4)) or e not in (1, 2):
        raise ValueError("Certified rational-square scale and e=1,2 required")
    root = {Q(1, 4): Q(1, 2), Q(1): Q(1), Q(4): Q(2)}[s]
    # R^e * R^-e * s^-1/2 * dR, with dR=sqrt(s) dx.
    return root**e * root**(-e) * root**(-1) * root


def support_row(s, lower=Q(1), upper=Q(2)):
    if s not in (Q(1, 4), Q(1), Q(4)):
        raise ValueError("Certified rational-square scale required")
    root = {Q(1, 4): Q(1, 2), Q(1): Q(1), Q(4): Q(2)}[s]
    return root * lower, root * upper


def strict_support_inclusion(s, support_lower, support_upper, lower=Q(1), upper=Q(2)):
    """Check a fixed compact support inside, rather than equal to, its hull."""
    lower,upper,support_lower,support_upper=map(Q,(lower,upper,support_lower,support_upper))
    if not lower < support_lower < support_upper < upper:
        raise ValueError("The fixed profile support must lie strictly inside its declared open hull")
    hull=support_row(s,lower,upper); actual=support_row(s,support_lower,support_upper)
    return {"hull":hull,"actual":actual,
            "strict":hull[0] < actual[0] < actual[1] < hull[1]}


def fixed_profile_jet_bound(n, profile_majorants, coordinate_majorants):
    """A finite positive composition/product majorant at one fixed order.

    This is deliberately profile dependent. It is an algebraic diagnostic
    for finiteness, not a numerical enclosure of the selected C-infinity
    profile.
    """
    if type(n) is not int or n < 0:
        raise ValueError("Nonnegative integer order required")
    if len(profile_majorants) != n + 1 or len(coordinate_majorants) != n + 1:
        raise ValueError("Every fixed order through n is required")
    values = []
    for item in tuple(profile_majorants) + tuple(coordinate_majorants):
        if not isinstance(item, (int, Q)) or item < 0:
            raise ValueError("Nonnegative exact majorants required")
        values.append(Q(item))
    profile = values[: n + 1]
    coordinate = values[n + 1 :]
    # Leibniz plus a coarse finite Bell-polynomial envelope.
    return sum(
        Q(comb(n, k) * factorial(k + 1))
        * (profile[k] + 1)
        * (coordinate[n - k] + 1) ** (k + 1)
        for k in range(n + 1)
    )


def arbitrary_family_obstruction(e, frequency):
    """Exact moment cancellation and linear Z-derivative growth in (4)."""
    if e not in (1, 2) or type(frequency) is not int or frequency <= 0:
        raise ValueError("e=1,2 and a positive integer frequency required")
    # Phi=30 y^2(1-y)^2, x=1+y. Compute weighted moments exactly.
    def moment(power):
        return 30 * sum(Q(comb(power, j)) * beta(j + 3, 3) for j in range(power + 1))

    ae = moment(e)
    mean = moment(e + 1) / ae
    weighted_derivative_moment = Q(frequency, 2 * ae) * (moment(e + 1) - mean * ae)
    # At x=5/4 the coefficient is nonzero and exactly proportional to N.
    x = Q(5, 4)
    y = x - 1
    phi = 30 * y**2 * (1 - y) ** 2
    point_derivative = Q(frequency, 2 * ae) * phi * (x - mean)
    return {
        "Ae": ae,
        "weighted_mean": mean,
        "weighted_derivative_moment": weighted_derivative_moment,
        "point_derivative": point_derivative,
    }


def review():
    profile = diagnostic_profile()
    scales = (Q(1, 4), Q(1), Q(4))
    checks = {
        "dependency_route_is_acyclic": set(topological_order()) == set(NODES),
        "diagnostic_profile_has_unit_mass": profile["mass"] == 1,
        "diagnostic_profile_A2_is_positive": profile["A2"] == Q(10, 7),
        "cauchy_schwarz_lower_bound": profile["A2"] * profile["length"] >= profile["mass"] ** 2,
        "reciprocal_A2_bound": 1 / profile["A2"] <= profile["length"],
    }
    for e in (1, 2):
        for s in scales:
            checks[f"e{e}_s{s}_exact_scaled_moment"] = dilation_moment(s, e) == 1
        one = arbitrary_family_obstruction(e, 1)
        seven = arbitrary_family_obstruction(e, 7)
        checks[f"e{e}_arbitrary_family_moment_preserved"] = one["weighted_derivative_moment"] == 0
        checks[f"e{e}_arbitrary_family_derivative_grows_linearly"] = (
            one["point_derivative"] != 0
            and seven["point_derivative"] == 7 * one["point_derivative"]
        )
    hulls = [support_row(s) for s in scales]
    checks["moving_support_enclosure_union_is_l_over_2_to_2u"] = min(a for a, _ in hulls) == Q(1, 2) and max(b for _, b in hulls) == 4
    inclusions=[strict_support_inclusion(s,Q(5,4),Q(7,4)) for s in scales]
    checks["compact_profile_support_is_strictly_inside_scaled_enclosure"] = all(row["strict"] for row in inclusions)
    checks["support_enclosure_is_not_identified_with_actual_support"] = all(row["hull"] != row["actual"] for row in inclusions)
    checks["fixed_finite_profile_data_give_finite_positive_bound"] = fixed_profile_jet_bound(
        4, [1, 3, 7, 15, 31], [1, 2, 4, 8, 16]
    ) > 0
    order = topological_order()
    checks["same_choice_recompute_follows_explicit_projector"] = order.index("rho_squared_projectors") < order.index("same_choice_signed_recompute")
    checks["contract_follows_same_choice_recompute"] = order.index("same_choice_signed_recompute") < order.index("fixed_profile_projector_contract")
    return {
        "checks": checks,
        "profile_diagnostic": {key: str(value) for key, value in profile.items()},
        "support_rows": [[str(a), str(b)] for a, b in hulls],
        "strict_support_rows": [{key:[str(x) for x in value] if isinstance(value,tuple) else value for key,value in row.items()} for row in inclusions],
        "scope": "exact finite rational controls for normalization, support and quantifier guards; the continuum conclusion uses one fixed C-infinity profile and finite-order composition",
    }
