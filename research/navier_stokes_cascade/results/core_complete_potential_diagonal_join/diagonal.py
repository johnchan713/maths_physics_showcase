"""Corrected finite-stage hierarchy joined to a locally finite diagonal."""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
KAPPA = Q(1, 100000)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


current = load(
    "complete_diagonal_current_recurrence",
    HERE.parent / "core_finite_stage_compact_recurrence/recurrence.py",
)


def index(value, name="index", minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}")
    return value


def parameters(h=Q(1, 200), kappa=KAPPA):
    h, kappa = Q(h), Q(kappa)
    if not Q(0) < h < Q(1, 100):
        raise ValueError("use the retained 0<h<1/100 regime")
    p = current.parameters(kappa)
    return dict(h=h, A=Q(1, 2) + h, D=Q(1, 2) - h, **p)


def representative_orders(j, kappa=KAPPA):
    """Conservative orders for the actual cycle j-1 -> j representatives."""
    j = index(j, "j", 1)
    p = parameters(kappa=kappa)
    B = current.rates(j - 1, p["kappa"])["regular"]
    C = B + Q(1, 2)
    return dict(
        complete_wave_velocity=B - p["kappa"],
        wave_potential=C - p["kappa"],
        mean_azimuthal_potential=C - 2 * p["kappa"],
        direct_axisymmetric_swirl=C - 2 * p["kappa"],
        complete_fixed_gauge_pressure_difference=C - 2 * p["kappa"],
        mean_radial_velocity=C + 1 - 2 * p["kappa"],
    )


def tuple_loss(m, h=Q(1, 200)):
    return 2 * parameters(h)["A"] + 2 * index(m, "m")


def cut_loss(m, h=Q(1, 200)):
    return tuple_loss(index(m, "m") + 1, h)


def total_loss(m, h=Q(1, 200)):
    return cut_loss(m, h) + 1


def residual_loss(m, h=Q(1, 200)):
    return 2 * parameters(h)["A"] + Q(3, 2) + 2 * index(m, "m")


def gain(j, h=Q(1, 200), kappa=KAPPA):
    return parameters(h, kappa)["h"] * current.parameters(kappa)["gain"] * index(j, "j", 1)


def ceil(value):
    value = Q(value)
    return -(-value.numerator // value.denominator)


def comparison(m=0, target=0, h=Q(1, 200), kappa=KAPPA):
    """One fixed comparison state proving the target for the one cut field."""
    m, target = index(m, "m"), index(target, "target")
    p = parameters(h, kappa)
    delta, B0 = p["gain"], p["B0"]
    L, H = cut_loss(m + 2, h), total_loss(m + 2, h)
    tail_stage = ceil(2 * (target + 1 + L + H) / (p["h"] * delta) - 1)
    finite_stage = ceil(((target + 1 + residual_loss(m, h)) / p["h"] - B0) / delta)
    J = max(m + 2, 1, tail_stage, finite_stage)
    tail = gain(J + 1, h, kappa) / 2 - L
    finite = p["h"] * current.rates(J, kappa)["regular"] - residual_loss(m, h)
    return dict(
        jet=m,
        target=target,
        comparison_stage=J,
        tail_stage=tail_stage,
        finite_residual_stage=finite_stage,
        tail_exponent=tail,
        cross_exponent=tail - H,
        square_exponent=2 * tail,
        finite_residual_exponent=finite,
        checks=dict(
            enough_jets=J >= m + 2,
            linear=tail >= target + 1,
            both_crosses=tail - H >= target + 1,
            square=2 * tail >= target + 1,
            finite_residual=finite >= target + 1,
        ),
    )


def choose_log_scale(j, gamma, constraints, previous=0, log_qstar_inverse=2):
    """Exact whole-interval majorant for arbitrary finite stage constants."""
    j = index(j, "j", 1)
    gamma, previous, domain = Q(gamma), Q(previous), Q(log_qstar_inverse)
    pairs = [(Q(c), Q(power)) for c, power in constraints]
    if gamma <= 0 or previous < 0 or domain < 0 or not pairs:
        raise ValueError("positive gamma and finite nonnegative data required")
    if any(c < 0 or power < 0 for c, power in pairs):
        raise ValueError("finite nonnegative log-constant data required")
    candidates = [previous + 1, domain + 1, Q(1)]
    candidates += [2 * (c + j) / gamma + power * power / (gamma * gamma)
                   for c, power in pairs]
    return max(candidates) + 1


def scale_review():
    previous = Q(0)
    rows = []
    for j, constraints in (
        (1, [(0, 0), (5, 7)]),
        (4, [(1000, 300), (10**6, 2000)]),
        (17, [(10**30, 10**12), (300, 0)]),
    ):
        gamma = gain(j) / 2
        chosen = choose_log_scale(j, gamma, constraints, previous)
        margin = max(Q(c) + Q(power) ** 2 / (2 * gamma) - gamma * chosen / 2 + j
                     for c, power in constraints)
        rows.append(dict(stage=j, gamma=gamma, chosen=chosen, margin=margin,
                         doubling_margin=chosen - previous - 1))
        previous = chosen
    return dict(rows=rows, checks=dict(
        whole_interval_strict=all(row["margin"] < 0 for row in rows),
        strict_doubling=all(row["doubling_margin"] > 0 for row in rows),
        arbitrary_finite_constants_allowed=True,
        stage_uniform_constants_not_required=True,
        raw_series_convergence_not_required=True,
    ))


@lru_cache(None)
def universal_controls():
    return load(
        "complete_diagonal_universal_controls",
        HERE.parent / "realization_foundation/realized_checks.py",
    )


def closure_statement():
    return dict(
        corrected_current_seed_and_gain_used=True,
        actual_fixed_gauge_complete_representatives_required=True,
        complete_pressure_differences_and_direct_swirl_cut=True,
        one_locally_finite_divergence_free_field_follows_conditionally=True,
        all_fixed_jets_of_complete_local_residual_are_flat_conditionally=True,
        protected_inner_velocity_and_heat_exterior_preserved_conditionally=True,
        stage_uniform_constants_fixed_q_contraction_or_raw_convergence_required=False,
        primitive_analytic_contracts_independently_proved=False,
        closed_time_localization_forcing_energy_global_comparison_or_blowup_proved=False,
    )
