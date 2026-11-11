"""Q178 physical-primary-pulse source audit and covariance obstruction."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q178-PHYSICAL-PRIMARY-SOURCE-AUDIT"
PROVIDERS = (
    "core_current_head_q177_weighted_support_source_reconciliation",
    "core_current_head_q148_coefficient_support_pressure",
    "core_current_head_q149_primary_pulse_join",
    "physical_pulse_transfer",
    "physical_pulse_foundation",
)
EXPECTED = {
    "core_current_head_q177_weighted_support_source_reconciliation": "Q177-selected-input-weighted-stress-and-reserved-support-route-reconciled",
    "core_current_head_q148_coefficient_support_pressure": "Q147-final-profile-finite-coefficient-supports-and-reserved-pressure-obstruction-joined",
    "core_current_head_q149_primary_pulse_join": "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined",
    "physical_pulse_transfer": "actual-background-to-physical-primary-pulse-transfer-audited",
    "physical_pulse_foundation": "actual-primary-physical-pulses-all-fixed-jets-and-common-domain-tails-relative-to-exact-core-and-finite-lemmas",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q178 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def physical_losses(a, b, h=Q(1, 200)):
    if type(a) is not int or type(b) is not int or a < 0 or b < 0:
        raise ValueError("nonnegative integer derivative orders required")
    h = Q(h)
    if not 0 < h < Q(1, 100):
        raise ValueError("inherited h range required")
    sx = Q(1, 2) + h / 2
    st = 1 + 3 * h / 2
    shift = a * sx + b * st
    return {
        "spatial_loss": sx,
        "time_loss": st,
        "velocity_growth": shift + sx,
        "pressure_growth": shift + 2 * sx,
        "potential_growth": shift,
    }


def borel_tail(N, qh):
    if type(N) is not int or N < 1:
        raise ValueError("positive integer tail index required")
    qh = Q(qh)
    if not 0 <= qh <= 1:
        raise ValueError("0<=q^h<=1 required")
    exact = (qh / 2) ** N / (1 - qh / 2)
    return {"exact": exact, "uniform_upper": Q(2) ** (1 - N), "bounded": exact <= Q(2) ** (1 - N)}


def covariance_growth(ell, h=Q(1, 200), kappa=Q(1, 100000)):
    if type(ell) is not int or ell < 0:
        raise ValueError("nonnegative integer dyadic band required")
    h, kappa = map(Q, (h, kappa))
    if not 0 < h < Q(1, 100) or kappa <= 0:
        raise ValueError("positive inherited h and kappa required")
    q_exponent = -1 - 2 * h * kappa
    binary_growth_exponent = -q_exponent * ell
    return {
        "Q_exponent": q_exponent,
        "binary_growth_exponent": binary_growth_exponent,
        "grows_with_band": q_exponent < 0,
        "absolutely_decaying": q_exponent > 0,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q178 identity required")
    status = statuses()
    q177 = evidence("core_current_head_q177_weighted_support_source_reconciliation")
    q148 = evidence("core_current_head_q148_coefficient_support_pressure")["symbolic"]
    q149 = evidence("core_current_head_q149_primary_pulse_join")["symbolic"]
    transfer = evidence("physical_pulse_transfer")["flags"]
    pulse = evidence("physical_pulse_foundation")["flags"]
    zero = physical_losses(0, 0)
    mixed = physical_losses(2, 1)
    tail = borel_tail(5, Q(3, 4))
    covariance = covariance_growth(10)
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "Q177_requests_this_exact_source_audit": "physical_pulse_foundation" in q177["review"]["next"] and "Q149" in q177["review"]["next"],
        "Q148_retains_reserved_pressure_and_support_frontier": "Pi1" in q148["reserved_mean"] and "Pi2" in q148["reserved_mean"],
        "Q149_corrects_the_rule_defined_schedule_frontier": "follows by construction" in q149["correction"],
        "same_profile_uniform_primary_family_is_joined": "same finalized leading/background profile" in q149["established"] and transfer["actual_background_inputs_to_primary_pulse_theorem_discharged"],
        "all_actual_used_bands_and_fixed_derivatives_are_covered": transfer["common_all_later_band_primary_pulse_values_and_fixed_jets_derived"] and pulse["one_starting_band_for_all_later_scales_and_every_fixed_physical_derivative_order"],
        "complete_velocity_pressure_potential_and_supports_are_retained": transfer["complete_primary_curl_velocity_and_wave_pressure_physically_transferred"] and transfer["exact_closed_support_and_one_sided_boundary_conditions_retained"] and pulse["uniform_primary_complete_curl_pressure_and_fixed_potential_bounds_derived"],
        "moving_pressure_and_background_pressure_memory_are_retained": pulse["exact_physical_graph_frame_and_moving_pressure_terms_retained"] and pulse["complete_background_radial_flow_and_I4_pressure_tails_retained"],
        "zero_and_mixed_physical_exponents_are_exact": zero["velocity_growth"] == Q(201, 400) and zero["pressure_growth"] == Q(201, 200) and mixed["potential_growth"] == Q(161, 80),
        "rule_defined_tail_has_one_uniform_bound": tail["bounded"],
        "sourced_inverse_input_class_boundary_is_retained": transfer["sourced_inverse_requires_input_class_and_smooth_endpoint_extension"] and not transfer["all_nonlinear_generated_source_class_estimates_verified"] and not pulse["generated_nonlinear_source_bounds_and_supports_all_verified"],
        "uncut_amplitude_is_not_a_cutoff_velocity_floor": pulse["uncut_radial_ODE_pulse_has_positive_lower_bound_relative_to_P"] and not pulse["positive_pointwise_lower_bound_for_oscillating_cutoff_velocity"],
        "primary_covariance_remainder_grows_absolutely": covariance["grows_with_band"] and not covariance["absolutely_decaying"],
        "complete_covariance_residual_and_blowup_are_not_promoted": not pulse["unconditional_complete_covariance_mean_corrections_verified"] and not pulse["unconditional_complete_residual_improvement_verified"] and not pulse["blowup_verified"],
        "formal_or_numerical_certificate_is_absent": not pulse["formal_certificate"] and not pulse["actual_giant_parameters_fields_or_starting_band_numerically_evaluated"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "physical_loss_fixtures": {
            "zero": {key: str(value) for key, value in zero.items()},
            "mixed_2_1": {key: str(value) for key, value in mixed.items()},
        },
        "borel_tail_fixture": {key: str(value) if isinstance(value, Q) else value for key, value in tail.items()},
        "covariance_fixture": {key: str(value) if isinstance(value, Q) else value for key, value in covariance.items()},
        "established": {
            "conditional_uniform_actual_primary_pulse_values_fixed_derivatives_and_supports": True,
            "one_rule_defined_schedule_and_one_later_tail": True,
            "full_graph_frame_curl_axial_velocity_and_pressure_memory_retained": True,
        },
        "correction": {
            "uncut_radial_amplitude_lower_bound_is_not_a_pointwise_cutoff_velocity_lower_bound": True,
            "finite_primary_bounds_do_not_verify_generated_source_membership": True,
            "negative_Q_exponent_cannot_be_called_absolute_decay": True,
            "historical_frozen_evidence_modified": False,
        },
        "scope": {
            "analytical_provider_theorems_replaced_by_tests": False,
            "all_generated_nonlinear_source_inputs_verified": False,
            "complete_covariance_mean_residual_or_blowup": False,
        },
        "next": "source-audit the complete covariance/mean construction beginning with Q150, retaining radial and axial velocity, every tensor and curl term, Pi1/Pi2, moving pressure and cutoff derivatives",
    }


def record():
    return {
        "status": "Q178-physical-primary-pulse-source-boundary-and-covariance-obstruction-audited",
        "established": "relative to the retained exact-core and finite lemmas, one rule-defined schedule gives uniform complete primary curl velocity, wave pressure and fixed-potential value, fixed-derivative and closed-support estimates on every used physical scale",
        "correction": "the uncut radial ODE lower bound is not a pointwise floor for the oscillating cutoff velocity, and the sourced inverse remains conditional on explicit source bounds and smooth endpoint support",
        "obstruction": "the uncancelled primary covariance remainder is Q^(-1-2h*kappa); its negative Q exponent makes it grow on later bands, so complete covariance/mean identities are required before any nonlinear hierarchy or residual claim",
        "historical_frozen_records_modified": False,
        "unconditional_generated_pulses_complete_covariance_residual_or_blowup_certified": False,
    }
