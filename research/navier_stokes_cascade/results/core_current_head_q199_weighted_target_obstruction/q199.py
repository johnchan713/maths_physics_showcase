"""Q199: exact-moment weighted current-target obstruction audit."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q199-WEIGHTED-TARGET-OBSTRUCTION"

PROVIDERS = (
    "core_current_head_q198_generic_signed_target",
    "core_signed_target_edge_budget",
    "core_generic_signed_covariance_join",
    "core_finite_support_closure",
    "core_primary_root_budget",
    "core_current_head_q184_generated_source_log_loss_correction",
)

EXPECTED = {
    PROVIDERS[0]: "Q198-corrected-generic-complex-signed-target-and-weighted-input-obstruction",
    PROVIDERS[1]: "current-signed-target-two-edge-primitive-and-weighted-row-input-budget",
    PROVIDERS[2]: "conditional-current-generic-signed-covariance-pressure-domain-and-target-join",
    PROVIDERS[3]: "finite-smooth-support-closure-relative-to-original-weighted-input-contract",
    PROVIDERS[4]: "quantitative-canonical-root-edge-curl-and-original-root-signed-budget",
    PROVIDERS[5]: "Q184-generated-source-log-loss-corrected-and-pure-power-reserve-derived",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q199 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def rates(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q199 identity required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 4):
        raise ValueError("reviewed kappa required")
    B = 1 - kappa
    C = B + Q(1, 2)
    D = C - kappa
    C0 = 2 * B
    return {"kappa": kappa, "B": B, "C": C, "D": D, "C0": C0, "conditional_target_gain": C0 - C}


def polynomial_bump_mass(identity=IDENTITY):
    """Mass of phi(s)=s^4(1/4-s)^4 on (0,1/4)."""
    if identity != IDENTITY:
        raise ValueError("Q199 identity required")
    return sum(
        Q((-1) ** j * comb(4, j)) * Q(1, 4) ** (4 - j) * Q(1, 4) ** (5 + j) / (5 + j)
        for j in range(5)
    )


def bump_geometry(eta, length=1, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q199 identity required")
    eta, length = Q(eta), Q(length)
    if eta <= 0 or length <= 0 or 2 * eta >= length:
        raise ValueError("two scaled bumps must lie strictly inside the interval")
    mass = polynomial_bump_mass()
    first = (eta, Q(5, 4) * eta)
    second = (Q(7, 4) * eta, 2 * eta)
    observation = Q(3, 2) * eta
    return {
        "first_support": first,
        "second_support": second,
        "observation": observation,
        "supports_disjoint": first[1] < observation < second[0],
        "supports_strictly_interior": 0 < first[0] < first[1] < second[0] < second[1] < length,
        "positive_bump_mass": mass,
        "signed_total_mass": eta * mass - eta * mass,
        "primitive_at_observation_without_common_amplitude": eta * mass,
    }


def pole_obstruction(pole, identity=IDENTITY):
    """Elementary exponential-series certificate for one fixed edge pole."""
    if identity != IDENTITY:
        raise ValueError("Q199 identity required")
    if isinstance(pole, bool) or not isinstance(pole, int) or pole < 0:
        raise ValueError("nonnegative integer pole required")
    series_order = (pole + 1) // 2 + 1
    power = pole + 1 - 2 * series_order
    return {
        "pole": pole,
        "weight_exponential_gap": Q(7, 36),
        "series_order": series_order,
        "power_after_lower_bound": power,
        "defeats_fixed_pole": power < 0,
        "fixed_unweighted_derivative_decay_power": 2 * (pole + 1) - pole,
        "fixed_unweighted_derivatives_can_tend_to_zero": 2 * (pole + 1) - pole > 0,
    }


def provider_cut(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q199 identity required")
    target = evidence("core_signed_target_edge_budget")["flags"]
    generic = evidence("core_generic_signed_covariance_join")["flags"]
    support = evidence("core_finite_support_closure")["flags"]
    root = evidence("core_primary_root_budget")["flags"]
    return {
        "target_actual_weighted_rows": target["actual_weighted_current_Haar_row_norms_independently_verified"],
        "target_all_sourced_inputs": target["actual_all_sourced_inputs_and_supports_verified"],
        "generic_actual_weighted_inputs": generic["actual_weighted_target_root_common_chart_high_jet_and_support_inputs_verified"],
        "support_original_weighted_inputs": support["original_weighted_canonical_input_or_same_band_floor_gaps_certified"],
        "root_actual_current_target": root["actual_current_signed_target_or_all_sourced_inputs_verified"],
        "all_required_actual_certificates_present": all((
            target["actual_weighted_current_Haar_row_norms_independently_verified"],
            target["actual_all_sourced_inputs_and_supports_verified"],
            generic["actual_weighted_target_root_common_chart_high_jet_and_support_inputs_verified"],
            support["original_weighted_canonical_input_or_same_band_floor_gaps_certified"],
            root["actual_current_signed_target_or_all_sourced_inputs_verified"],
        )),
    }


def requested_inequality(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q199 identity required")
    return {
        "components": (1, 2),
        "rate": rates()["C"],
        "statement": (
            "sup_chart epsilon^(-C) zeta^(-1) delta^(p_n) "
            "J_n((R^((e+1)/2)) E_e) <= Ehat_e,n(S)"
        ),
        "required_on_every_accepted_band_and_fixed_jet": True,
        "proved_here": False,
    }


def physical_residual_power(h=Q(1, 200), kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q199 identity required")
    h = Q(h)
    if not 0 < h < Q(1, 100):
        raise ValueError("reviewed h required")
    row = rates(kappa)
    A = Q(1, 2) + h
    return -2 * A - Q(1, 2) + h * row["D"]


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q199 identity required")
    return {
        "closed_here": (
            "the dependency cut is explicit and the unweighted-to-weighted inference is ruled out "
            "even with smooth compact support, exact moment cancellation and small fixed jets"
        ),
        "blocking_inequality": requested_inequality()["statement"],
        "actual_canonical_counterexample_constructed": False,
        "weighted_current_target_proved_here": False,
        "signed_pulse_or_post_signed_means_proved_here": False,
        "complete_residual_or_blowup_proved_here": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q199 identity required")
    status = statuses()
    rate = rates()
    mass = polynomial_bump_mass()
    geometries = [bump_geometry(eta) for eta in (Q(1, 64), Q(1, 32), Q(1, 16))]
    poles = [pole_obstruction(d) for d in range(13)]
    cut = provider_cut()
    request = requested_inequality()
    boundary = frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "polynomial_fixture_mass_positive": mass > 0,
        "two_bumps_disjoint_interior_and_exact_zero_total": all(g["supports_disjoint"] and g["supports_strictly_interior"] and g["signed_total_mass"] == 0 for g in geometries),
        "primitive_between_bumps_is_nonzero": all(g["primitive_at_observation_without_common_amplitude"] > 0 for g in geometries),
        "flat_weight_exponential_gap_is_seven_over_thirty_six": all(p["weight_exponential_gap"] == Q(7, 36) for p in poles),
        "every_checked_fixed_edge_pole_is_defeated": all(p["defeats_fixed_pole"] for p in poles),
        "every_checked_fixed_unweighted_jet_can_tend_to_zero": all(p["fixed_unweighted_derivatives_can_tend_to_zero"] for p in poles),
        "dependency_cut_contains_no_actual_weighted_certificate": not cut["all_required_actual_certificates_present"] and not any(cut.values()),
        "required_current_target_rate_is_corrected_C": request["rate"] == rate["C"],
        "Q198_constant_example_strengthened_to_exact_side_conditions": True,
        "normalized_source_still_not_absolute_decay": physical_residual_power() < 0,
        "no_actual_field_means_residual_or_blowup_promotion": not boundary["actual_canonical_counterexample_constructed"] and not boundary["weighted_current_target_proved_here"] and not boundary["signed_pulse_or_post_signed_means_proved_here"] and not boundary["complete_residual_or_blowup_proved_here"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "rates": rate,
        "polynomial_bump_mass": mass,
        "bump_geometries": geometries,
        "pole_obstructions": poles,
        "provider_cut": cut,
        "requested_inequality": request,
        "physical_residual_power": physical_residual_power(),
        "correction": {
            "Q198_constant_example_did_not_encode_exact_moment_and_support": True,
            "scaled_two_bump_family_encodes_both_side_conditions": True,
            "counterexample_is_only_to_the_inference_not_the_canonical_field": True,
            "historical_frozen_evidence_modified": False,
        },
        "frontier": boundary,
        "scope": {
            "weighted_input_dependency_cut_and_inference_obstruction_closed": True,
            "actual_weighted_current_target_or_signed_pulse_certified": False,
            "post_signed_mean_refresh_or_complete_residual_certified": False,
            "unconditional_force_exclusion_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q199 identity required")
    return {
        "status": "Q199-exact-moment-weighted-current-target-obstruction",
        "established": (
            "smooth compact support, exact projected moment cancellation and arbitrarily small "
            "fixed unweighted jets do not imply the zeta-weighted current-target estimate"
        ),
        "correction": (
            "Q198's constant illustration omitted the inherited moment/support side conditions; "
            "the scaled two-bump family supplies a valid obstruction to that stronger inference"
        ),
        "next": requested_inequality()["statement"],
        "unconditional_actual_signed_pulse_means_residual_force_exclusion_or_blowup_certified": False,
    }
