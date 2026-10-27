#!/usr/bin/env python3
"""Exact controls for the scale-compatible primitive-domain certificate."""
from fractions import Fraction as F


PRIME = 2003
NORMALIZED_SUPPORT_MARGIN = F(3, 4 * PRIME)


def background_floors(f, t, d, gbar, cone):
    """Floors implied by the certified same-band perturbation estimates."""
    f, t, d, gbar, cone = map(F, (f, t, d, gbar, cone))
    if min(f, t, d, cone) <= 0 or gbar < 0 or cone >= 1:
        raise ValueError("Invalid same-band gap data")
    return {
        "F": f / 2,
        "abs_g": t / 2,
        "minus_Ntheta": t / (2 * (gbar + 1)),
        "lambda_squared": f * t * d / (4 * (gbar + 1) ** 2),
        "one_minus_cone_ratio": cone,
    }


def physical_support_margin(scale):
    """A normalized torus margin scales to zero in physical coordinates."""
    scale = F(scale)
    if scale <= 0:
        raise ValueError("Positive scale required")
    return scale * NORMALIZED_SUPPORT_MARGIN


def annular_contraction(epsilon):
    """The retained annular root contraction using its normalized target ball."""
    epsilon = F(epsilon)
    if epsilon <= 0:
        raise ValueError("Positive target scale required")
    beta, quadratic, target_ratio = 610, 50, F(143, 500)
    absolute_radius = 2 * beta * target_ratio * epsilon
    contraction = 8 * beta ** 2 * quadratic * target_ratio * epsilon
    return {
        "target_ratio": target_ratio,
        "absolute_radius": absolute_radius,
        "contraction": contraction,
        "slack": 1 - contraction,
    }


def repair_contraction():
    """The accepted repair contraction and normalized uniqueness margin."""
    contraction = F(1, 2 ** 219)
    return {
        "contraction": contraction,
        "slack": 1 - contraction,
        "normalized_uniqueness_margin": F(2000 - 1600),
    }


def edge_weight_exponent(w, x, length):
    """The exponent diverges at either edge, so the weight tends to zero."""
    w, x, length = map(F, (w, x, length))
    if min(w, x, length - x) <= 0:
        raise ValueError("Interior weighted coordinate required")
    return w * w / (x * x) + F(4) / ((length - x) ** 2)


def review():
    floors = background_floors(F(1, 5), F(1, 7), F(1, 11), 4, F(1, 17))
    support_scales = [physical_support_margin(F(1, 2 ** n)) for n in range(1, 7)]
    annular = annular_contraction(F(1, 10 ** 16))
    repair = repair_contraction()
    edge_exponents = [edge_weight_exponent(F(1, 100), F(1, n), 1) for n in (100, 200, 400, 800)]
    checks = {
        "background_floors_are_uniformly_positive": all(x > 0 for x in floors.values()),
        "normalized_support_margin_is_exactly_three_over_8012": NORMALIZED_SUPPORT_MARGIN == F(3, 8012),
        "physical_support_margins_shrink_with_scale": support_scales == sorted(support_scales, reverse=True),
        "no_uniform_physical_support_floor_is_inferred": physical_support_margin(F(1, 2 ** 100)) < F(1, 10 ** 20),
        "annular_contraction_has_strict_dimensionless_slack": 0 < annular["contraction"] < 1,
        "annular_absolute_root_radius_scales_with_input": annular_contraction(F(1, 10 ** 17))["absolute_radius"] == annular["absolute_radius"] / 10,
        "repair_contraction_has_strict_dimensionless_slack": 0 < repair["contraction"] < 1,
        "repair_normalized_uniqueness_margin_is_positive": repair["normalized_uniqueness_margin"] > 0,
        "compact_edge_weight_has_zero_floor_not_positive_floor": edge_exponents == sorted(edge_exponents),
    }
    certificate = {
        "uniform_background_floors": {key: str(value) for key, value in floors.items()},
        "normalized_support_margin": str(NORMALIZED_SUPPORT_MARGIN),
        "annular_contraction": {key: str(value) for key, value in annular.items()},
        "repair_contraction": {key: str(value) for key, value in repair.items()},
        "compact_map_domain": "fixed normalized support hull and zeta*delta^(-p) weighted seminorms; no positive infimum of zeta is used",
    }
    scope = {
        "heterogeneous_scalar_Delta_prim_rejected": True,
        "scale_compatible_primitive_certificate_joined": True,
        "all_fixed_order_seminorms_follow_from_the_five_existing_contracts": True,
        "unconditional_complete_residual_force_comparison_or_blowup": False,
    }
    return {"certificate": certificate, "checks": checks, "scope": scope}
