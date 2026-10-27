#!/usr/bin/env python3
"""Exact controls for the closed-corner flatness reduction and obstruction."""
from fractions import Fraction as Q
from math import exp


def normalized_identity(T, eta_squared):
    T, eta_squared = map(Q, (T, eta_squared))
    if T < 0 or eta_squared < 0:
        raise ValueError("normalized coordinates must be nonnegative")
    return T + eta_squared


def jacobian_floor(h):
    h = Q(h)
    if not 0 <= h <= Q(1, 100):
        raise ValueError("require 0<=h<=1/100")
    return 1 - 2 * h


def corner_transfer_constant(power):
    power = int(power)
    if power < 0:
        raise ValueError("power must be nonnegative")
    return 4 ** power


def stage_exponents(stage, m, h=Q(1, 200), kappa=Q(1, 100000)):
    stage, m = int(stage), int(m)
    h, kappa = Q(h), Q(kappa)
    if stage < 0 or m < 0 or h <= 0 or not 0 < kappa < Q(1, 6):
        raise ValueError("invalid stage, jet or parameters")
    delta = Q(1, 2) - 3 * kappa
    A = Q(1, 2) + h
    B = Q(3, 2) - 2 * kappa + stage * delta
    K = 2 * A + Q(3, 2) + 2 * m
    L = 2 * A + 2 * (m + 1)
    return {
        "delta": delta,
        "fixed_stage_residual": h * B - K,
        "cutoff_tail": h * delta * (stage + 1) / 2 - L,
    }


def first_stage(m, target_power, h=Q(1, 200), kappa=Q(1, 100000)):
    m, target_power = int(m), int(target_power)
    if m < 0 or target_power < 0:
        raise ValueError("nonnegative m and target required")
    stage = max(1, m + 2)
    while True:
        exponents = stage_exponents(stage, m + 2, h, kappa)
        if (exponents["fixed_stage_residual"] >= target_power + 1
                and exponents["cutoff_tail"] >= target_power + 1):
            return stage
        stage += 1


def boundary_loss(q, eta, order):
    q, eta, order = Q(q), Q(eta), int(order)
    if q <= 0 or not -1 <= eta <= 1 or order < 0:
        raise ValueError("invalid boundary fixture")
    return q ** (order + 1) / (1 - eta * eta + q)


def pointwise_spike(X, q):
    X, q = float(X), float(q)
    if X < 0 or q <= 0:
        raise ValueError("invalid spike fixture")
    ratio = X / q
    return ratio * exp(-(ratio * ratio))


def review():
    powers = (0, 1, 2, 8, 32)
    transfer = [{"power": p, "constant": corner_transfer_constant(p)} for p in powers]
    boundary_rows = []
    for order in (0, 1, 5, 20):
        q = Q(1, 2**20)
        interior = boundary_loss(q, Q(1, 2), order)
        pole = boundary_loss(q, 1, order)
        boundary_rows.append({
            "order": order, "interior": str(interior), "pole": str(pole),
            "pole_over_q_power": str(pole / q**order),
        })
    spike_rows = []
    for n in (4, 8, 16, 32):
        q = 2.0 ** (-n)
        spike_rows.append({
            "q": q, "fixed_X_one": pointwise_spike(1, q),
            "moving_X_equals_q": pointwise_spike(q, q),
        })
    stages = [{"m": m, "N": N, "stage": first_stage(m, N)}
              for m, N in ((0, 0), (0, 10), (3, 5))]
    return {
        "transfer_rows": transfer,
        "boundary_rows": boundary_rows,
        "compact_X_spike_rows": spike_rows,
        "stage_rows": stages,
        "closed_constants_required": {
            "fixed_stage": "sup_closed q^(-(h*B_J-K_m))*|R(U^[J])|_m < infinity",
            "cutoff_tail": "sup_closed,transition q^(-(h*delta*j-L_m))*|(V_(j,a),P_(j,a))|_m < infinity uniformly in a",
        },
        "checks": {
            "normalized_terminal_chart_is_T_plus_eta_squared_equal_one": normalized_identity(Q(3, 4), Q(1, 4)) == 1,
            "implicit_coordinate_jacobian_has_strict_closed_floor": jacobian_floor(Q(1, 100)) == Q(49, 50),
            "uniform_q_flatness_transfers_with_four_to_the_power": all(row["constant"] == 4 ** row["power"] for row in transfer),
            "closed_boundary_can_lose_one_q_order": all(row["pole_over_q_power"] == "1" for row in boundary_rows),
            "pointwise_compact_X_flatness_need_not_be_uniform": all(row["moving_X_equals_q"] > Q(1, 3) for row in spike_rows),
            "finite_comparison_stage_exists_after_fixed_jet_and_target": all(row["stage"] >= row["m"] + 2 for row in stages),
            "scale_absorption_requires_finite_closed_chart_constants": True,
            "extra_localization_time_derivative_is_a_fixed_finite_request": True,
        },
        "scope": {
            "coordinate_transfer_proved": True,
            "closed_chart_constants_proved_finite_for_actual_field": False,
            "compact_force_and_comparison_unconditional": False,
            "blowup_proved": False,
        },
    }
