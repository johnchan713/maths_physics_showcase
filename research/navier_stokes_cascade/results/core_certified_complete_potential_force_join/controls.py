"""Exact controls for the certified complete-potential force join."""
from fractions import Fraction as Q
from localization import (
    domain_review, residual_review, extension_review, energy_review,
    closure_statement,
)


def encode(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def review():
    sections = dict(
        domain=domain_review(), residual=residual_review(),
        extension=extension_review(), energy=energy_review(),
    )
    checks = {
        name: bool(section["checks"]) and all(section["checks"].values())
        for name, section in sections.items()
    }
    scope = closure_statement()
    checks["scope"] = (
        scope["same_certified_complete_diagonal_field_localized"]
        and scope["certified_fixed_stage_domains_and_single_cutoff_sequence_retained"]
        and not scope["uniform_stage_constant_growth_required"]
        and scope["compact_smooth_force_and_zero_initial_data_follow_conditionally"]
        and not scope["primitive_corner_and_punctured_bounds_independently_proved"]
        and not scope["future_velocity_or_pressure_constructed"]
        and not scope["global_competitor_comparison_or_blowup_proved"]
    )
    return dict(sections=encode(sections), scope=scope, checks=checks,
                exact_controls_are_not_continuum_proofs=True)
