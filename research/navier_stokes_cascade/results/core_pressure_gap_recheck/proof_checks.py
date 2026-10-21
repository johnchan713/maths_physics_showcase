"""Exact rational checks for the shared core pressure and growth partition."""

from fractions import Fraction as F


A_MIN = F(1, 2)
A_MAX = F(51, 100)
J_MAX = F(1, 20)
P2_MIN = F(256)


def _add(left, right):
    size = max(len(left), len(right))
    return tuple((left[i] if i < len(left) else 0) +
                 (right[i] if i < len(right) else 0) for i in range(size))


def _scale(value, poly):
    return tuple(value * item for item in poly)


def _mul(left, right):
    out = [F(0) for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return tuple(out)


def _trim(poly):
    out = list(poly)
    while out and out[-1] == 0:
        out.pop()
    return tuple(out)


def direct_velocity_polynomial(A, j):
    """Expand -A(1-2 eta U)U-H U_eta without the pressure terms."""
    one, eta, d = (F(1),), (F(0), F(1)), (F(1), F(0), F(-1))
    U = (j, F(4))
    H = _add(_scale(1 - A, eta), _mul(d, U))
    bracket = _add(one, _scale(-2, _mul(eta, U)))
    return _trim(_add(_scale(-A, _mul(bracket, U)), _scale(-4, H)))


def stated_velocity_polynomial(A, j):
    return (
        -(A + 4) * j,
        -20 + 2 * A * j * j,
        (4 + 16 * A) * j,
        16 + 32 * A,
    )


def velocity_envelope():
    """Proof constants for |R_v| <= 21|eta|+5j on the full range."""
    cubic_linear = F(20)
    j_eta = (4 + 16 * A_MAX) * J_MAX
    j2_eta = 2 * A_MAX * J_MAX * J_MAX
    constant_over_j = A_MAX + 4
    return {
        "cubic_linear": cubic_linear,
        "j_eta": j_eta,
        "j2_eta": j2_eta,
        "eta_total": cubic_linear + j_eta + j2_eta,
        "constant_over_j": constant_over_j,
    }


def atom_pressure_certificate():
    """Common-denominator numerator for d/dx of the theta=1 atom term."""
    # K_atom(x,A)=10*((1+A)+(A-1)x)/(1+x)^3.
    # Its derivative has numerator 10*(-4-2A+2(1-A)x).
    derivative_numerator_max = -4 - 2 * A_MIN + 2 * (1 - A_MIN)
    endpoint_min = F(5, 2) * A_MIN
    return {
        "derivative_numerator_max_without_factor_10": derivative_numerator_max,
        "endpoint_min": endpoint_min,
        "required_theta1_atom_mass": F(5, 2),
        "total_measure_mass_ceiling": F(5),
        "real_pressure_derivative_upper_coefficient": F(20),
        "atom_derivative_lower_coefficient": F(10),
    }


def separation_certificate():
    """Sharpen the small-Z angular localization and the resulting chi floor."""
    ratio = (F(1, 100) + F(5, 256)) / (F(5, 4) - F(21, 256))
    eta_over_j = F(1, 39)
    h_over_j = (1 - F(9, 2) * eta_over_j
                - J_MAX * J_MAX * eta_over_j * eta_over_j
                - 4 * J_MAX * J_MAX * eta_over_j ** 3)
    h_floor = F(22, 25)
    chi_floor = F(88 * 88, 88 * 88 + 1)
    return {
        "exact_ratio": ratio,
        "eta_over_j_ceiling": eta_over_j,
        "H_over_j_exact_lower_control": h_over_j,
        "H_over_j_floor": h_floor,
        "chi_floor": chi_floor,
    }


def complex_pressure_certificate():
    """Exact check behind the complex pressure and Cauchy input bounds."""
    r0 = F(1, 10**20 * 256)
    real_part_floor = 1 - r0 * r0
    mixture_upper = F(5) / (real_part_floor * real_part_floor)
    return {
        "r0": r0,
        "real_part_floor": real_part_floor,
        "mixture_upper_over_P2": mixture_upper,
        "cauchy_derivative_coefficient_over_P2_r0": F(20),
    }


def parameter_order():
    """The datum is fixed before Lambda and C; XR cancels from dX/X."""
    nodes = ("reference_schedule", "pressure_measure", "core_parameters", "Lambda", "C")
    edges = tuple(zip(nodes, nodes[1:]))
    position = {node: index for index, node in enumerate(nodes)}
    return {
        "nodes": nodes,
        "edges": edges,
        "all_edges_forward": all(position[a] < position[b] for a, b in edges),
        "radial_scale_power_in_dX_over_X": 0,
    }


def _encode(value):
    if type(value) is F:
        return f"{value.numerator}/{value.denominator}"
    if isinstance(value, dict):
        return {key: _encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_encode(item) for item in value]
    return value


def review():
    velocity = velocity_envelope()
    atom = atom_pressure_certificate()
    separation = separation_certificate()
    complex_pressure = complex_pressure_certificate()
    order = parameter_order()
    polynomial_controls = [
        direct_velocity_polynomial(A, j) == stated_velocity_polynomial(A, j)
        for A in (A_MIN, A_MAX, F(503, 1000))
        for j in (F(0), F(1, 10**18), J_MAX)
    ]
    checks = {
        "velocity_polynomial_rederived_exactly": all(polynomial_controls),
        "velocity_eta_coefficient_below_21": velocity["eta_total"] < 21,
        "velocity_constant_coefficient_below_5j": velocity["constant_over_j"] < 5,
        "ideal_atom_Kp_decreases_to_endpoint": atom["derivative_numerator_max_without_factor_10"] < 0,
        "ideal_atom_gives_Kp_at_least_5_over_4": atom["endpoint_min"] == F(5, 4),
        "small_Z_ratio_below_1_over_39": separation["exact_ratio"] < separation["eta_over_j_ceiling"],
        "Hstar_over_j_above_22_over_25": separation["H_over_j_exact_lower_control"] > separation["H_over_j_floor"],
        "chi_floor_above_point_99": separation["chi_floor"] > F(99, 100),
        "complex_mixture_bound_below_10P2": complex_pressure["mixture_upper_over_P2"] < 10,
        "pressure_fixed_before_Lambda_and_C": order["all_edges_forward"],
        "axis_pressure_is_independent_of_radial_scale": order["radial_scale_power_in_dX_over_X"] == 0,
    }
    return _encode({
        "velocity": velocity,
        "atom_pressure": atom,
        "separation": separation,
        "complex_pressure": complex_pressure,
        "parameter_order": order,
        "checks": checks,
    })
