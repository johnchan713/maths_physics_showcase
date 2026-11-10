"""Bind Q154 endpoint jets to exact localization, compact force, and energy."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q154-physical-localization-compact-force-energy-join"
PROVIDERS = (
    "core_current_head_q8_localization_force_join",
    "core_certified_complete_potential_force_join",
    "force_foundation",
    "current_localized_force",
    "localized_force_extension",
    "core_scale_compatible_corner_force_replay",
)

LOCALIZATION_TERMS = (
    "c R_1(u0,p0)",
    "c_t u0",
    "d_t",
    "-(Delta c)u0",
    "-2 sum_i c_i partial_i u0",
    "-Delta d",
    "p0 grad c",
    "c(c-1)(u0 dot grad)u0",
    "c(u0 dot grad c)u0",
    "c(u0 dot grad)d",
    "c(d dot grad)u0",
    "(d dot grad)d",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def localization_ledger(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q155 identity required")
    return {
        "terms": LOCALIZATION_TERMS,
        "count": len(LOCALIZATION_TERMS),
        "vanishing_thirteenth": "(d dot grad c)u0=0 because d=grad(c) cross A",
        "viscosity": Q(1),
        "complete_pressure": "p0 grad c",
        "force_divergence_retained": True,
    }


def derivative_requests(m, identity=IDENTITY):
    if identity != IDENTITY or m < 0:
        raise ValueError("Q155 identity and nonnegative order required")
    return {
        "force": {"potential": m + 2, "direct_swirl": m + 1, "pressure": m},
        "cauchy": {"potential": m + 3, "direct_swirl": m + 2,
                   "pressure": m + 1, "residual": m + 1},
        "Q154_spare": {"potential": 1, "direct_swirl": 1, "pressure": 1},
    }


def corner_transfer(N, identity=IDENTITY):
    if identity != IDENTITY or N < 0:
        raise ValueError("Q155 identity and nonnegative power required")
    return {
        "q_power": N + 1,
        "rho_power": N + 1,
        "constant_factor": 4 ** (N + 1),
        "cartesian_flat_for_arbitrary_N": True,
    }


def borel_term_power(j, time_order, identity=IDENTITY):
    if identity != IDENTITY or j < 0 or time_order < 0:
        raise ValueError("Q155 identity and nonnegative orders required")
    return time_order - j


def borel_admitted(j, spatial_order, time_order, identity=IDENTITY):
    if identity != IDENTITY or min(j, spatial_order, time_order) < 0:
        raise ValueError("Q155 identity and nonnegative orders required")
    return spatial_order + time_order <= j // 2 and time_order < j


def energy_budget(force_l1, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q155 identity required")
    F = Q(force_l1)
    if F < 0:
        raise ValueError("Nonnegative force integral required")
    return {
        "fixed_time": "E(t)+2 integral_0^t D <= F(t)^2",
        "sup_energy": F * F,
        "integrated_dissipation": F * F / 2,
        "sum_of_separate_extrema_budget_claimed": False,
    }


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q155 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    required = {
        "force_foundation": (
            "same_current_single_realization_is_localized",
            "all_twelve_pressure_viscous_time_and_nonlinear_localization_terms_retained",
            "axisymmetric_potential_cutoff_and_direct_swirl_preserve_divergence",
            "both_background_pressure_tails_and_full_fixed_gauge_memory_retained",
            "one_Borel_sequence_controls_every_fixed_mixed_force_jet",
            "whole_space_preterminal_energy_and_total_dissipation_bounds_joined",
            "zero_initial_datum_joined_to_same_field",
        ),
        "current_localized_force": (
            "all_twelve_curl_pressure_viscous_and_nonlinear_terms_retained",
            "complete_current_terminal_force_jets_smooth_compact_and_flat_at_origin_conditionally_derived",
            "conditional_one_smooth_compact_global_force_matching_exact_preterminal_residual_supplied",
            "conditional_current_whole_space_energy_and_integrated_dissipation_bound_supplied",
            "force_divergence_in_exact_pressure_equation_retained",
        ),
        "localized_force_extension": (
            "complete_curl_pressure_viscosity_and_transport_localization_terms_retained",
            "conditional_all_order_compact_smooth_force_extension_supplied",
            "conditional_whole_space_energy_and_integrated_dissipation_supplied",
        ),
    }
    for name, keys in required.items():
        flags = data[name].get("flags", {})
        if not all(flags.get(key) for key in keys):
            raise ValueError("Missing provider contract: " + name)
    if data["core_current_head_q8_localization_force_join"]["status"] != \
            "current-head-Q8-twelve-term-localization-and-compact-force-joined-conditionally":
        raise ValueError("Q8 exact localization theorem required")
    if data["core_certified_complete_potential_force_join"]["status"] != \
            "certified-domain-conditional-localization-force-and-energy-join":
        raise ValueError("Certified force theorem required")
    return {
        "Q154_supplies_Q8_punctured_endpoint_premise": True,
        "Q153_supplies_corner_complete_residual_flatness": True,
        "same_field_cutoffs_support_pressure_gauge": True,
        "compact_smooth_force_and_preterminal_energy": "conditional",
        "future_velocity_or_pressure": False,
    }


def record():
    return {
        "status": "Q155-physical-localization-compact-force-and-energy-bound-conditionally",
        "correction": "Q8's punctured endpoint premise is supplied by Q154, while Q153's complete residual gives the Cartesian-flat corner trace; the sharper direct force ledger is A:m+2, B:m+1, p:m with one further Cauchy time derivative",
        "established": "the same Q154 physical field admits exact divergence-free axisymmetric localization, zero initial velocity, a twelve-term unit-viscosity residual with complete pressure, one compact smooth global force matching every terminal force jet, and the fixed-time preterminal energy/dissipation estimate",
        "scope": "conditional analytical composition of named continuum providers; the force trace is regularized but the terminal velocity and pressure at the origin are not, and no future solution is constructed",
        "next": "recover the localized fixed-gauge pressure and perform the expanding-ball comparison with any admissible smooth finite-energy solution having the same zero datum and force, without assuming competitor pressure, derivative or boundary-flux bounds",
        "compact_force_and_energy_certified_conditionally": True,
        "future_solution_pressure_comparison_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
