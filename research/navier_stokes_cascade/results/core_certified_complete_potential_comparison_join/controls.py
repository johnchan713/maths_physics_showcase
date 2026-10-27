"""Exact controls for the certified complete-potential comparison join."""
from fractions import Fraction as Q
from comparison import (
    interface_review, pressure_review, commutator_review, radius_review,
    limit_review, closure_statement,
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
        interface=interface_review(), pressure=pressure_review(),
        commutator=commutator_review(), radius=radius_review(),
        limit=limit_review(),
    )
    checks = {
        name: bool(section["checks"]) and all(section["checks"].values())
        for name, section in sections.items()
    }
    scope = closure_statement()
    checks["scope"] = (
        scope["same_certified_complete_potential_field_force_and_fixed_gauge_pressure_used"]
        and scope["certified_fixed_stage_domains_and_single_cutoff_sequence_retained"]
        and scope["positive_core_contract_retained_as_dependency"]
        and not scope["uniform_stage_constant_growth_required"]
        and scope["pressure_gradient_recovered_without_competitor_pressure_decay"]
        and scope["global_smooth_bounded_energy_competitor_exclusion_follows_conditionally_from_inherited_core_growth"]
        and not scope["primitive_field_force_and_positive_core_contracts_independently_proved"]
        and not scope["unconditional_admissible_force_or_blowup_proved"]
    )
    return dict(sections=encode(sections), scope=scope, checks=checks,
                exact_controls_are_not_continuum_proofs=True)
