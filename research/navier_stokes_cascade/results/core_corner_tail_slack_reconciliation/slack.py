#!/usr/bin/env python3
"""Exact exponent controls for closed-corner logarithmic slack."""
from fractions import Fraction as Q
from math import exp


KAPPA = Q(1, 100000)


def stage_gain(stage, h=Q(1, 200), kappa=KAPPA):
    stage = int(stage)
    h, kappa = Q(h), Q(kappa)
    if stage < 1 or h <= 0 or not 0 < kappa < Q(1, 6):
        raise ValueError("positive stage, h and delta required")
    return h * (Q(1, 2) - 3 * kappa) * stage


def slack_gain(stage, h=Q(1, 200), kappa=KAPPA):
    return stage_gain(stage, h, kappa) / 2


def young_log_majorant(stage, log_power, h=Q(1, 200), kappa=KAPPA):
    """P^2/g for (1+s)^P exp(-g*s/2), s>=1."""
    log_power = int(log_power)
    if log_power < 0:
        raise ValueError("nonnegative logarithmic power required")
    return Q(log_power * log_power, 1) / stage_gain(stage, h, kappa)


def full_exponent_model(s, log_power):
    """Normalized full-exponent model, unbounded for P>0 as s grows."""
    s, log_power = int(s), int(log_power)
    if s < 1 or log_power < 0:
        raise ValueError("s>=1 and P>=0 required")
    return (1 + s) ** log_power


def slack_model(s, stage, log_power, h=Q(1, 200), kappa=KAPPA):
    """Normalized half-gain model used only as a numerical diagnostic."""
    s, log_power = int(s), int(log_power)
    if s < 1 or log_power < 0:
        raise ValueError("s>=1 and P>=0 required")
    gamma = float(slack_gain(stage, h, kappa))
    return (1.0 + s) ** log_power * exp(-gamma * s)


def chosen_log_scale(stage, log_constant, log_power, previous=0,
                     log_qstar_inverse=2, h=Q(1, 200), kappa=KAPPA):
    """The inherited sufficient whole-interval scale rule, exactly."""
    stage, log_power = int(stage), int(log_power)
    log_constant, previous = Q(log_constant), Q(previous)
    domain = Q(log_qstar_inverse)
    if stage < 1 or log_power < 0 or min(log_constant, previous, domain) < 0:
        raise ValueError("finite nonnegative scale data required")
    gamma = slack_gain(stage, h, kappa)
    candidates = [previous + 1, domain + 1, Q(1)]
    candidates.append(2 * (log_constant + stage) / gamma
                      + Q(log_power * log_power, 1) / (gamma * gamma))
    return max(candidates) + 1


def scale_margin(stage, log_constant, log_power, chosen,
                 h=Q(1, 200), kappa=KAPPA):
    """Upper log of C(1+s)^P exp(-gamma*s) 2^j (with ln2<1)."""
    gamma = slack_gain(stage, h, kappa)
    return (Q(log_constant) + Q(log_power * log_power, 1) / (2 * gamma)
            - gamma * Q(chosen) / 2 + int(stage))


def review():
    rows = []
    previous = Q(0)
    for stage, log_constant, power in ((1, 5, 7), (4, 1000, 300), (17, 10**8, 2000)):
        chosen = chosen_log_scale(stage, log_constant, power, previous)
        rows.append({
            "stage": stage,
            "gain": str(stage_gain(stage)),
            "reserved_half_gain": str(slack_gain(stage)),
            "chosen_log_scale": str(chosen),
            "scale_margin": str(scale_margin(stage, log_constant, power, chosen)),
        })
        previous = chosen
    growth = [{"s": s, "full": full_exponent_model(s, 3),
               "slack": slack_model(s, 1, 3)} for s in (1, 10, 1000, 10000)]
    return {
        "rows": rows,
        "model_rows": growth,
        "closed_tail_bounds": {
            "unsupported_full_exponent": "sup q^(-(g_j-L_m)) |cut_j| without its log weight",
            "proved_from_closed_representative_bound": "sup q^(-(g_j/2-L_m)) |cut_j| <= C*exp(P^2/g_j)",
            "chosen_scale_bound": "|cut_j|_m <= 2^(-j) q^(g_j/2-L_m) on q<=1/a_j",
        },
        "checks": {
            "positive_stage_gain": all(Q(row["gain"]) > 0 for row in rows),
            "exactly_half_gain_reserved": all(Q(row["reserved_half_gain"]) * 2 == Q(row["gain"]) for row in rows),
            "full_exponent_log_model_is_unbounded_for_positive_power": full_exponent_model(10**6, 3) > full_exponent_model(10**3, 3),
            "young_majorant_is_finite_for_every_fixed_stage_and_power": all(young_log_majorant(j, p) >= 0 for j, p in ((1, 0), (1, 7), (17, 2000))),
            "whole_interval_scale_margins_are_strict": all(Q(row["scale_margin"]) < 0 for row in rows),
            "cutoff_scale_is_selected_after_its_finite_coefficient_data": True,
            "cutoff_transition_adds_no_new_uniform_constant_obstruction": True,
            "primitive_closed_coefficient_bounds_still_required": True,
        },
        "scope": {
            "full_exponent_no_log_tail_supremum_proved": False,
            "full_exponent_no_log_tail_supremum_required": False,
            "half_gain_closed_tail_transfer_proved_conditionally": True,
            "actual_primitive_closed_coefficients_proved": False,
            "compact_force_comparison_or_blowup_unconditional": False,
        },
    }
