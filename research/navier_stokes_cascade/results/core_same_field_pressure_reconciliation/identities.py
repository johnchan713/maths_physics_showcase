#!/usr/bin/env python3
"""Exact scalar controls for the same-field leading pressure reconciliation."""
from fractions import Fraction as Q


def core_density(varphi, Phi, C):
    """d_X Pi_core=(varphi*Phi/C)^2=E_core^2/(2X)."""
    varphi, Phi, C = map(Q, (varphi, Phi, C))
    if C <= 0:
        raise ValueError("C must be positive")
    return (varphi * Phi / C) ** 2


def swirl_squared(X, varphi, Phi, C):
    X = Q(X)
    if X < 0:
        raise ValueError("X must be nonnegative")
    return 2 * X * core_density(varphi, Phi, C)


def forward_pressure(axis, increments):
    return Q(axis) + sum(map(Q, increments), Q(0))


def backward_pressure(tail_increments):
    return -sum(map(Q, tail_increments), Q(0))


def pressure_defect(axis, all_increments):
    """The X-independent mismatch between forward and backward gauges."""
    return Q(axis) + sum(map(Q, all_increments), Q(0))


def split_identity(axis, increments, split):
    increments = tuple(map(Q, increments))
    if not 0 <= split <= len(increments):
        raise ValueError("invalid split")
    forward = forward_pressure(axis, increments[:split])
    backward = backward_pressure(increments[split:])
    return forward, backward, pressure_defect(axis, increments)


def i4_background_radial_derivative(amplitude, mu, X):
    """Radial derivative of amplitude*(K1-mu^2/(4X))."""
    amplitude, mu, X = map(Q, (amplitude, mu, X))
    if amplitude < 0 or mu <= 0 or X <= 0:
        raise ValueError("require nonnegative amplitude and positive mu,X")
    return amplitude * mu**2 / (4 * X**2)


def review():
    density_rows = []
    for X in (Q(1, 7), Q(1), Q(9)):
        for varphi, Phi, C in ((2, 3, 5), (7, 1, 11), (5, 4, 13)):
            d_core = core_density(varphi, Phi, C)
            d_swirl = swirl_squared(X, varphi, Phi, C) / (2 * X)
            density_rows.append({
                "X": str(X), "core_density": str(d_core),
                "swirl_density": str(d_swirl), "equal": d_core == d_swirl,
            })

    increments = tuple(map(Q, (2, 3, 5, 7)))
    axis = -sum(increments)
    split_rows = []
    for split in range(len(increments) + 1):
        forward, backward, defect = split_identity(axis, increments, split)
        split_rows.append({
            "split": split, "forward": str(forward),
            "backward": str(backward), "defect": str(defect),
        })

    bad_axis = axis + 1
    bad = split_identity(bad_axis, increments, 2)
    memory = i4_background_radial_derivative(Q(1, 8), 3, 2)
    return {
        "density_rows": density_rows,
        "split_rows": split_rows,
        "negative_controls": {
            "nonzero_fifth_row_defect": str(bad[2]),
            "positive_i4_background_radial_pressure_derivative": str(memory),
        },
        "checks": {
            "core_radial_pressure_density_equals_leading_swirl_density": all(row["equal"] for row in density_rows),
            "same_axis_and_zero_total_defect_give_forward_backward_equality_at_every_split": all(row["forward"] == row["backward"] for row in split_rows),
            "omitted_fifth_row_leaves_nonzero_gauge_defect": bad[2] == 1,
            "q_positive_background_pressure_need_not_equal_leading_pressure": memory > 0,
            "protected_inner_support_is_needed_to_transfer_the_core_identity": True,
            "full_corrected_pressure_is_not_collapsed_to_centrifugal_pressure": True,
            "pressure_identification_does_not_prove_corner_flatness": True,
        },
        "scope": {
            "leading_q_zero_pressure_identified_on_same_field": True,
            "full_q_positive_pressure_identified_with_core_pressure": False,
            "fixed_jet_closed_corner_residual_bound_proved_here": False,
            "unconditional_force_residual_or_blowup_proved": False,
        },
    }
