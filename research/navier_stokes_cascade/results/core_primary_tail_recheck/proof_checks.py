"""One primary-pulse tail and quantitative dyadic/edge flatness constants."""

from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from math import factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("primary_tail_inherited_arithmetic", HERE.parent / "explicit_tilt_pulse_tail/tail_bounds.py")
inherited = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = inherited
spec.loader.exec_module(inherited)

KEYS = {"h", "cost", "rb", "omega", "lpsi", "cmesh", "cpert", "dbox", "dchart", "epsilon", "radius", "qstar", "length", "slow_cost", "slow_tolerance"}


def exact(x):
    if type(x) not in (int, F):
        raise TypeError("Exact integer/Fraction input required")
    return F(x)


def select_tail(inputs):
    if set(inputs) != KEYS:
        raise ValueError("All actual finite budget inputs required; no fixture defaults")
    p = {k: exact(v) for k, v in inputs.items()}
    if min(p.values()) <= 0 or not 0 < p["h"] < F(1, 100):
        raise ValueError("Positive finite budgets and construction h required")
    if p["epsilon"] > 1 or p["qstar"] > 1 or p["radius"] > F(1, 6057072) or p["slow_cost"] < 1 or p["slow_tolerance"] > 1:
        raise ValueError("Construction radius, domain, deficit and slow-budget ranges required")
    tilt = inherited.tilt_certificate(p["epsilon"])
    gamma = p["omega"] ** 2 * tilt["chosen_d0"] / (9 * tilt["kappa"])
    B0 = min(p["rb"] / 8, p["omega"] ** 2 * tilt["chosen_d0"] / (16 * p["lpsi"]), gamma / (128 * p["cpert"]))
    m = inherited.ceil_log2(p["cost"] / B0)
    slow_power = F(3 * factorial(8) * 4**8) * p["slow_cost"] / (p["slow_tolerance"] * p["h"]**8)
    thresholds = {
        "mixed_value": 2 + inherited.ceil_fraction(F(m) / (2 * p["h"])),
        "mesh_and_chart": inherited.ceil_root(max(F(1), 16 * p["cmesh"] / gamma, p["dbox"] / p["dchart"]), 6),
        "pulse_length": inherited.ceil_root(p["length"] / (2 * p["radius"]), 2),
        "background": 3 + inherited.ceil_log2(1 / p["qstar"]),
        "padded_slow_coordinates": inherited.ceil_root(slow_power, 4),
        "geometry_floor": 16,
    }
    ell = max(thresholds.values())
    checks = band_checks(p, ell, gamma, m, slow_power)
    return {"ell_tail": ell, "thresholds": thresholds, "gamma": gamma,
            "pointwise_budget": B0, "binary_exponent_required": m,
            "slow_fourth_power_required": slow_power,
            "binary_scale_exponent": inherited.floor_fraction(2 * p["h"] * (ell - 2)),
            "checks": checks}


def band_checks(p, ell, gamma, m, slow_power):
    if type(ell) is not int or ell < 2:
        raise ValueError("Positive integer band at least two required")
    return {
        "binary_value_gate": inherited.floor_fraction(2 * p["h"] * (ell - 2)) >= m,
        "mesh_half_budget": F(ell**6) >= 16 * p["cmesh"] / gamma,
        "independent_chart_gate": F(ell**6) >= p["dbox"] / p["dchart"],
        "whole_pulse_length": 2 * p["radius"] * ell**2 >= p["length"],
        "background_strict_half_margin": ell >= 3 + inherited.ceil_log2(1 / p["qstar"]),
        "whole_half_line_slow_gate": F(ell**4) >= slow_power,
        "geometry_floor": ell >= 16,
    }


def fixture(h=F(1, 200)):
    """Manufactured finite budgets, never actual selected norm enclosures."""
    return dict(h=h, cost=F(1000), rb=F(1, 9), omega=F(1, 8), lpsi=F(13),
                cmesh=F(17), cpert=F(19), dbox=F(6), dchart=F(1, 23),
                epsilon=F(1, 64), radius=F(1, 6057072), qstar=F(1, 128),
                length=F(100000), slow_cost=F(23), slow_tolerance=F(1, 11))


def flat_constant(a, B, loss, target):
    a, loss = exact(a), exact(loss)
    if a <= 0 or type(B) is not int or type(target) is not int or min(B, target) < 0:
        raise ValueError("Positive curvature and nonnegative fixed orders required")
    D = F(B) + max(F(0), loss + target)
    return {"D": D, "log_constant": D * D / (4 * a)}


def edge_constant(B, ca, cb):
    ca, cb = exact(ca), exact(cb)
    if type(B) is not int or B < 0 or min(ca, cb) <= 0:
        raise ValueError("Fixed nonnegative edge power and positive exact edge constants required")
    j = (B + 1) // 2
    # Keep half of BOTH original exponential edge factors after absorption.
    return F(3) + factorial(j) * ((2 / ca)**j + (2 / cb)**j)


def dyadic(power):
    return F(2**power) if power >= 0 else F(1, 2**(-power))


def active_bands(q):
    q = exact(q)
    if q <= 0:
        raise ValueError("Positive physical q required")
    e = q.numerator.bit_length() - q.denominator.bit_length()
    if dyadic(e) > q:
        e -= 1
    center = -e
    return [ell for ell in range(center - 4, center + 5) if F(1, 4) <= q / dyadic(-ell) <= 4]


def physical_losses(h, a, b):
    return inherited.physical_exponents(exact(h), a, b)


@lru_cache(None)
def review():
    p = fixture()
    selected = select_tail(p)
    base = selected["ell_tail"]
    later = [band_checks(p, base + delta, selected["gamma"], selected["binary_exponent_required"], selected["slow_fourth_power_required"]) for delta in (0, 1, 2, 17, 10000)]
    tiny = select_tail(fixture(F(1, 2**200)))
    sums = []
    for q in (dyadic(-20), F(3, 7) * dyadic(-20), F(1, 3) * dyadic(-30)):
        bands = active_bands(q)
        for K in (0, 1, 3, 9):
            sums.append(len(bands) <= 5 and sum((dyadic(-ell)**K for ell in bands), F(0)) <= 5 * 4**K * q**K)
    squares = []
    for B in (0, 3, 11):
        for loss in (F(-2), F(0), F(13, 7)):
            cert = flat_constant(F(1, 101), B, loss, 7)
            for ell in map(F, (1, 2, 16, 1000)):
                gap = cert["log_constant"] - (cert["D"] * ell - F(1, 101) * ell**2)
                squares.append(gap == F(1, 101) * (ell - cert["D"] / (2 * F(1, 101)))**2 and gap >= 0)
    exponents = [physical_losses(F(1, 200), a, b) for a in range(7) for b in range(7 - a)]
    return {
        "selector": selected,
        "tiny_h_selector_summary": {"h": F(1, 2**200), "band_integer_bits": tiny["ell_tail"].bit_length(),
                                     "scale_exponent_bits": tiny["binary_scale_exponent"].bit_length()},
        "flatness_example": flat_constant(F(1, 101), 11, F(13, 7), 7),
        "edge_example": edge_constant(8, F(1, 2 * 100**2), F(2)),
        "physical_exponent_rows": exponents,
        "checks": {
            "all_starting_gates_hold": all(selected["checks"].values()),
            "all_checked_later_bands_hold": all(all(x.values()) for x in later),
            "tiny_h_without_materializing_binary_scale": all(tiny["checks"].values()) and tiny["binary_scale_exponent"] > 10**20,
            "all_enlarged_band_sums_have_fixed_overlap_bound": all(sums),
            "same_domain_square_completion_is_exact": all(squares),
            "primary_pressure_velocity_potential_losses_retained": all(x["pressure"] - x["velocity"] == x["sx"] and x["velocity"] - x["potential"] == x["sx"] for x in exponents),
            "narrow_inner_edge_cost_is_not_erased": edge_constant(8, F(1, 20000), F(2)) > edge_constant(8, F(1), F(2)),
        },
    }
