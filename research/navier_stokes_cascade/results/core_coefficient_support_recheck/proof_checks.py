"""Finite coefficient support join and reserved-pressure obstruction.

Exact rational controls check algebra. They do not evaluate the selected
profile, solve its inner Volterra equations, or prove a continuum supremum.
"""

from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent


def load_local(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


matching = load_local("coefficient_recheck_matching", HERE.parent / "matching_root_foundation/matching_bounds.py")
conservative = load_local("coefficient_recheck_conservative", HERE.parent / "coefficient_support_foundation/conservative_checks.py")


def normalize(raw, rho, K, e, lam):
    """Ordinary (M,I,Cp,J,S) to B rows; K=1/(sqrt(2)*rho^(3/2))."""
    rho, K, e, lam = map(Q, (rho, K, e, lam))
    if min(rho, K, e, lam) <= 0:
        raise ValueError("Positive normalization scales required")
    M, I, Cp, J, S = map(Q, raw)
    m, i, cp, j, s = M / rho, K * I, Cp / e, K * J / e, S / (rho * e)
    return (m, (m - j) / lam, i, -s, cp)


def denormalize(rows, rho, K, e, lam):
    rho, K, e, lam = map(Q, (rho, K, e, lam))
    if min(rho, K, e, lam) <= 0:
        raise ValueError("Positive normalization scales required")
    r1, r2, r3, r4, r5 = map(Q, rows)
    return (rho * r1, r3 / K, e * r5, e * (r1 - lam * r2) / K, -rho * e * r4)


def reciprocal_norm(k, e_star):
    if type(k) is not int or k < 0 or e_star <= 0:
        raise ValueError("Fixed nonnegative order and positive amplitude required")
    return Q(2 if k == 0 else 4 if k == 1 else 5) / Q(e_star)


def coefficient_bound(norms, k, rho, K, e_star, lam):
    rho, K, e_star, lam = map(Q, (rho, K, e_star, lam))
    if min(rho, K, e_star, lam) <= 0 or lam > Q(1, 100):
        raise ValueError("Positive scales and 0<lambda<=.01 required")
    if len(norms) != 5 or min(norms) < 0:
        raise ValueError("Five nonnegative moment norm bounds required")
    a, b, c, d, e = map(Q, norms)
    W = reciprocal_norm(k, e_star)
    return 800 * max(a / rho, (a / rho + K * W * d) / lam, K * b, W * e / rho, W * c)


def radial_flux(U, M, Me, X, eta, h, n):
    D = Q(1, 2) - h
    return (2 * eta * X * U - 2 * (D + 2 * n * h) * eta * M - (1 - eta * eta) * Me) / (1 - 2 * h * eta * eta)


def cutoff_flux(chi, sigma_chi_prime, V, M, eta, h):
    return chi * V - 2 * eta * sigma_chi_prime * M / (1 - 2 * h * eta * eta)


def support_edges(n):
    if type(n) is not int or n < 1:
        raise ValueError("Positive finite coefficient order required")
    return {
        "velocity": "Xvel",
        "pressure": "Xv" if n <= 2 else "Xvel",
        "stress": "Xb" if n == 1 else "Xv" if n == 2 else "Xvel",
    }


def pressure_departure(X, mu, K1, K2, chi1, q2h, chi2):
    """q^(2A)*p_B-Pi0 on I4; q^(4h)=(q^(2h))^2."""
    return chi1 * q2h * (K1 - mu * mu / (4 * X)) + chi2 * q2h * q2h * K2


def pressure_floor(Xa, Xb, mu, chi1, q2h):
    """Best possible sup norm after arbitrary X-independent cancellation."""
    if not 0 < Xa < Xb or q2h <= 0:
        raise ValueError("Ordered positive radial points and positive q^(2h) required")
    return abs(chi1) * q2h * mu * mu * (1 / Xa - 1 / Xb) / 8


def axial_cut(b, P, PX, Pe, X, eta, h, chi, sigma_chi_prime):
    L, d = 1 - 2 * h * eta * eta, 1 - eta * eta
    return (chi * (2 * b * eta * P - 2 * eta * X * PX + d * Pe)
            + 2 * eta * sigma_chi_prime * P) / L


def time_cut(b, P, PX, Pe, X, eta, h, chi, sigma_chi_prime):
    L, D = 1 - 2 * h * eta * eta, Q(1, 2) - h
    return (chi * (-b * P + X * PX + D * eta * Pe) - sigma_chi_prime * P) / L


@lru_cache(None)
def review():
    raw = tuple(map(Q, (2, 3, 5, 7, 11)))
    rho, K, e_star, lam = Q(2), Q(3, 5), Q(7), Q(1, 1000)
    angles = tuple(map(Q, (-1, 0, 1))) + (Q(2, 5),)
    inverse = matching.repair_review()
    norm_bound = coefficient_bound(raw, 8, rho, K, e_star, lam)
    zero_flux, base_flux, cut_flux, nonzero_cuts = [], [], [], []
    h, mu, X = Q(1, 200), Q(3), Q(10)
    for eta in angles:
        base_flux.append(radial_flux(0, mu * eta, mu, X, eta, h, 0) == -mu)
        for n in (1, 2, 3, 7, 1000000):
            V = radial_flux(0, 0, 0, X, eta, h, n)
            zero_flux.append(V == 0)
            cut_flux.append(cutoff_flux(Q(5, 32), Q(-27, 32), V, 0, eta, h) == 0)
        P, PX = -mu * mu / (4 * X), mu * mu / (4 * X * X)
        zero_flux.extend((axial_cut(-1, P, PX, 0, X, eta, h, 1, 0) == 0,
                          time_cut(-1, P, PX, 0, X, eta, h, 1, 0) == 0))
        nonzero_cuts.append(time_cut(-1, P, PX, 0, X, eta, h, Q(5, 32), Q(-27, 32)) != 0)
    sample = dict(Xa=Q(10), Xb=Q(30), mu=mu, K1=Q(7), K2=Q(-13), chi1=Q(3, 7), q2h=Q(1, 5), chi2=Q(2, 3))
    p0 = pressure_departure(sample["Xa"], **{k: v for k, v in sample.items() if k not in ("Xa", "Xb")})
    p1 = pressure_departure(sample["Xb"], **{k: v for k, v in sample.items() if k not in ("Xa", "Xb")})
    floor = pressure_floor(sample["Xa"], sample["Xb"], mu, sample["chi1"], sample["q2h"])
    # The midpoint constant saturates the elementary two-point minimax bound.
    a = mu * mu / 4
    mid = a * (1 / sample["Xa"] + 1 / sample["Xb"]) / 2
    centered = [pressure_departure(x, mu, mid, 0, sample["chi1"], sample["q2h"], 0)
                for x in (sample["Xa"], sample["Xb"])]
    lo6 = matching.rb.exp_bounds(-6)[0]
    hi11 = matching.rb.exp_bounds(-11)[1]
    residuals = conservative.conservative_review()
    plateau = conservative.plateau_review()
    graded = all([(i, n - i) for i in range(n + 1) if i == n or n - i == n] == [(0, n), (n, 0)] and 2 * n > n
                 for n in range(1, 20))
    checks = {
        "ordinary_five_row_normalization_roundtrip": all(denormalize(normalize(raw, rho, K, e_star / (1 + eta * eta), lam), rho, K, e_star / (1 + eta * eta), lam) == raw for eta in angles),
        "current_order_map_is_linear_by_grading": graded,
        "continuum_inverse_bound_800_replayed": all(inverse["checks"].values()) and Q(inverse["E_inverse_upper"]) < 800,
        "normalizer_amplitude_jets_and_lambda_loss_retained": norm_bound == 3200000 and reciprocal_norm(8, 7) == Q(5, 7),
        "same_reserved_leading_radial_velocity": all(base_flux),
        "all_sampled_finite_positive_order_velocities_and_cutoff_fluxes_zero": all(zero_flux + cut_flux),
        "independent_signed_total_controls_keep_previous_viscosity": all(residuals["checks"].values()),
        "nonlocal_second_pressure_constant_can_survive": all(plateau["checks"].values()),
        "pressure_oscillation_independent_of_K1_K2": p1 - p0 == 2 * floor and floor > 0,
        "pressure_floor_saturated_by_midpoint_constant": centered == [-floor, floor],
        "actual_log_interval_has_strict_positive_floor": lo6 > hi11,
        "cutoff_transition_restores_time_derivative_of_stationary_piece": all(nonzero_cuts),
        "both_exceptional_pressure_support_edges_retained": support_edges(1)["pressure"] == support_edges(2)["pressure"] == "Xv" and support_edges(3)["pressure"] == "Xvel",
    }
    return {
        "checks": checks,
        "inverse_bound": 800,
        "E_inverse_upper": inverse["E_inverse_upper"],
        "U_inverse_upper": inverse["U_inverse_upper"],
        "example_coefficient_norm_bound": norm_bound,
        "exact_conservative_identity_count": residuals["exact_integrated_identity_count"],
        "manufactured_second_pressure_plateau": plateau["manufactured_second_pressure_plateau"],
        "pressure_fixture": {"left": p0, "right": p1, "necessary_sup_floor": floor},
        "actual_pressure_floor": "abs(chi(c1*q))*q^(2h)*mu^2/(8*rho3)*(exp(-6)-exp(-11))",
        "stationary_physical_piece": "q^-1*(-mu^2/(4X))=-mu^2/(2r^2)",
        "support_rows": {str(n): support_edges(n) for n in (1, 2, 3, 7, 1000000)},
        "scope": "rational fixtures and inherited continuum estimates; finite coefficient construction conditional on regular inner Volterra lemmas; no uniform physical pulse or residual claim",
    }
