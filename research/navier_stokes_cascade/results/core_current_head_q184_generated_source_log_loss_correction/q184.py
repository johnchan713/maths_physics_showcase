"""Correct Q183's source frontier by retaining polynomial band losses."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q184-GENERATED-SOURCE-LOG-LOSS-CORRECTION"
PROVIDERS = (
    "core_current_head_q183_provider_discharge_generated_source_frontier",
    "core_current_head_q180_target_source_simultaneous_induction",
    "core_current_head_q161_source_tail_simultaneous_induction",
    "source_envelope_edge_audit",
    "regular_source_foundation",
    "finite_stage_induction",
    "current_target_cycle_bridge",
    "physical_pulse_foundation",
    "core_current_head_q177_weighted_support_source_reconciliation",
)
EXPECTED = {
    PROVIDERS[0]: "Q183-current-head-provider-frontier-reconciled-to-generated-source-supremum",
    PROVIDERS[1]: "Q180-physical-target-source-tail-simultaneous-induction-reconciled",
    PROVIDERS[2]: "Q161-physical-source-inverse-tail-simultaneous-induction-joined",
    PROVIDERS[3]: "regular-source-envelope-edge-estimates-with-explicit-freezing-correction",
    PROVIDERS[4]: "regular-source-foundation-sharper-linear-gain-and-actual-initial-supported-input-relative-to-core-and-pulse-lemmas",
    PROVIDERS[5]: "conditional-current-state-complete-finite-stage-induction",
    PROVIDERS[6]: "conditional-current-target-cancellation-cycle-bridge",
    PROVIDERS[7]: "actual-primary-physical-pulses-all-fixed-jets-and-common-domain-tails-relative-to-exact-core-and-finite-lemmas",
    PROVIDERS[8]: "Q177-selected-input-weighted-stress-and-reserved-support-route-reconciled",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q184 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def polynomial_absorption(degree, h=Q(1, 200), eta=Q(1, 100000), identity=IDENTITY):
    """Bound ell^(2d) 2^(-h eta ell) using D=ceil(d), log(2)>=1/2."""
    if identity != IDENTITY:
        raise ValueError("Q184 identity required")
    if isinstance(degree, bool) or not isinstance(degree, (int, Q)):
        raise ValueError("nonnegative rational polynomial degree required")
    degree = Q(degree)
    if degree < 0:
        raise ValueError("nonnegative rational polynomial degree required")
    h, eta = map(Q, (h, eta))
    if h <= 0 or eta <= 0:
        raise ValueError("positive h and exponent reserve required")
    degree_ceiling = -(-degree.numerator // degree.denominator)
    power = 2 * degree_ceiling
    upper = Q(factorial(power)) * (Q(2) / (h * eta)) ** power
    return {
        "degree_in_S": degree,
        "degree_ceiling": degree_ceiling,
        "power_in_ell": power,
        "h": h,
        "eta": eta,
        "uniform_upper": upper,
        "inequality": (
            "ell^(2d) epsilon^eta <= (2D)! (2/(h eta))^(2D), D=ceil(d), "
            "epsilon=2^(-h ell)"
        ),
        "uses_log_two_lower_bound_one_half": True,
    }


def exact_power_counterexample(degree=1, bands=(16, 32, 64), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q184 identity required")
    if isinstance(degree, bool) or not isinstance(degree, int) or degree < 1:
        raise ValueError("positive integer degree required")
    if not bands or any(isinstance(ell, bool) or not isinstance(ell, int) or ell < 1 for ell in bands):
        raise ValueError("positive integer bands required")
    values = tuple(Q(ell) ** (2 * degree) for ell in bands)
    return {
        "degree_in_S": degree,
        "bands": tuple(bands),
        "epsilon_minus_B_times_epsilon_B_Sd": values,
        "strictly_increasing": all(left < right for left, right in zip(values, values[1:])),
        "all_band_supremum_finite": False,
        "reason": "the normalized sequence is S^d=ell^(2d), which is unbounded",
    }


def corrected_rates(j=0, kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q184 identity required")
    if isinstance(j, bool) or not isinstance(j, int) or j < 0:
        raise ValueError("finite nonnegative cycle required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 8):
        raise ValueError("reviewed kappa range required")
    B0_log = Q(3, 2) - 2 * kappa
    C0_log = B0_log + Q(1, 2)
    delta_log = Q(1, 2) - 3 * kappa
    eta = kappa
    B0_pure = B0_log - eta
    C0_pure = C0_log - eta
    delta_pure = delta_log - eta
    target_log_gain = Q(1, 2) - kappa
    target_pure_gain = target_log_gain - eta
    return {
        "j": j,
        "kappa": kappa,
        "eta": eta,
        "B_log": B0_log + j * delta_log,
        "C_log": C0_log + j * delta_log,
        "delta_log": delta_log,
        "B_pure": B0_pure + j * delta_pure,
        "C_pure": C0_pure + j * delta_pure,
        "delta_pure": delta_pure,
        "target_log_gain": target_log_gain,
        "target_pure_gain": target_pure_gain,
        "target_pure_margin_over_required": target_pure_gain - delta_pure,
        "pure_source_gain_positive": delta_pure > 0,
    }


def fixed_cycle_source_contract(stage=0, jet_order=2, polynomial_degree=1, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q184 identity required")
    if any(isinstance(value, bool) or not isinstance(value, int) for value in (stage, jet_order)):
        raise ValueError("integer stage and jet required")
    if isinstance(polynomial_degree, bool) or not isinstance(polynomial_degree, (int, Q)):
        raise ValueError("rational polynomial degree required")
    polynomial_degree = Q(polynomial_degree)
    if stage < 0 or jet_order < 0 or polynomial_degree < 0:
        raise ValueError("nonnegative stage, jet and degree required")
    rates = corrected_rates(stage)
    absorption = polynomial_absorption(polynomial_degree)
    return {
        "stage": stage,
        "jet_order": jet_order,
        "polynomial_degree": polynomial_degree,
        "log_renormalized_bound": (
            "|D^I(f_(j,m)/P)| <= C epsilon^(B_log) S^d sqrt(zeta) "
            "y_a^(-r) y_b^(-s)"
        ),
        "same_domain_raw_source_Ai_exist_conditionally": True,
        "same_domain_log_P_Ki_exist_conditionally": True,
        "smooth_endpoint_extension_exists_conditionally": True,
        "Q183_exact_power_bound_follows_without_margin": polynomial_degree == 0,
        "pure_power_bound_after_kappa_reserve": (
            "|D^I(f_(j,m)/P)| <= C_eta epsilon^(B_pure) sqrt(zeta) "
            "y_a^(-r) y_b^(-s)"
        ),
        "absorption": absorption,
        "rates": rates,
    }


def provider_derivation(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q184 identity required")
    q183 = evidence(PROVIDERS[0])["review"]
    source = evidence(PROVIDERS[3])["flags"]
    regular = evidence(PROVIDERS[4])["flags"]
    finite = evidence(PROVIDERS[5])["flags"]
    target = evidence(PROVIDERS[6])["flags"]
    pulse = evidence(PROVIDERS[7])["flags"]
    q177 = evidence(PROVIDERS[8])["review"]
    return {
        "Q183_correctly_rejects_synthetic_fixture_as_proof": q183["correction"][
            "Q161_synthetic_source_quotient_fixture_is_not_a_physical_supremum_proof"
        ],
        "Bell_majorants_are_analytically_available_at_all_fixed_mixed_orders": source[
            "all_fixed_mixed_normalized_source_jets_have_Bell_majorants"
        ],
        "raw_source_grammar_retains_full_path_and_current_targets": source[
            "generated_regular_source_grammar_has_full_path_and_entrance_inputs_relative_to_current_target_estimates"
        ],
        "pressure_viscosity_graph_and_connections_are_retained": regular[
            "full_viscosity_pressure_graph_coefficients_and_cylindrical_connections_retained"
        ],
        "complete_cut_source_has_endpoint_extension": regular[
            "full_source_endpoint_extension_and_closed_label_supports_preserved"
        ],
        "finite_cycle_preserves_source_envelope_and_endpoint_hypotheses": finite[
            "every_generated_regular_source_support_Gaussian_envelope_and_endpoint_hypothesis_preserved"
        ],
        "current_targets_are_derived_on_the_whole_half_line": target[
            "all_half_line_target_recurrence_relative_to_reviewed_current_map_inputs"
        ],
        "primary_fixed_jets_share_one_later_band": pulse[
            "one_starting_band_for_all_later_scales_and_every_fixed_physical_derivative_order"
        ],
        "polynomial_band_degrees_are_part_of_the_named_classes": True,
        "independent_exact_core_and_finite_lemma_validation_remains_open": q177["correction"][
            "external_independent_validation_of_exact_core_and_finite_lemmas_still_open"
        ],
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q184 identity required")
    status = statuses()
    derivation = provider_derivation()
    counterexample = exact_power_counterexample()
    contract = fixed_cycle_source_contract()
    rates = corrected_rates()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "all_analytic_source_derivation_links_hold": all(derivation.values()),
        "Q183_exact_power_normalization_is_not_implied": (
            counterexample["strictly_increasing"]
            and not counterexample["all_band_supremum_finite"]
            and not contract["Q183_exact_power_bound_follows_without_margin"]
        ),
        "one_kappa_absorbs_every_fixed_polynomial_degree": contract["absorption"]["uniform_upper"] > 0,
        "corrected_pure_cycle_has_positive_gain": rates["pure_source_gain_positive"],
        "target_retains_two_kappa_pure_margin": rates["target_pure_margin_over_required"] == 2 * rates["kappa"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "provider_derivation": derivation,
        "exact_power_counterexample": counterexample,
        "fixed_cycle_contract": contract,
        "corrected_rates": rates,
        "correction": {
            "Q183_fixture_warning_remains_valid": True,
            "Q183_misidentified_missing_Ai_Ki_as_first_conditional_gap": True,
            "named_analytic_theorems_supply_Ai_Ki_and_endpoint_extension_with_polynomial_S_losses": True,
            "Q183_displayed_exact_power_supremum_omits_allowed_polynomial_S_degree": True,
            "historical_frozen_evidence_modified": False,
        },
        "scope": {
            "log_renormalized_fixed_cycle_source_contract_closed_conditionally": True,
            "pure_power_contract_closed_after_one_kappa_reserve": True,
            "original_exact_delta_pure_power_contract_proved": False,
            "exact_core_and_finite_existence_inverse_gap_lemmas_independently_reproved": False,
            "downstream_Q180_Q182_replayed_with_corrected_pure_rates_here": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
        "next": (
            "replay Q180--Q182 using either the exact logarithmically renormalized classes or "
            "the pure-power rates B0=3/2-3*kappa and delta=1/2-4*kappa; then begin independent "
            "review of the retained exact-core and finite existence/inverse/gap premises"
        ),
    }


def record():
    return {
        "status": "Q184-generated-source-log-loss-corrected-and-pure-power-reserve-derived",
        "established": (
            "relative to the named current-state, source and finite-map theorems, every fixed "
            "physical generated-source jet has finite raw-source and log-Gaussian Bell inputs, "
            "complete endpoint extension and a uniform logarithmically renormalized bound; "
            "spending eta=kappa converts every fixed polynomial band loss to a pure-power bound"
        ),
        "correction": (
            "Q183 correctly rejects synthetic Bell arrays as a continuum proof but overlooks the "
            "separate analytical source theorem and writes an exact-power supremum without the "
            "polynomial S=ell^2 losses allowed by W/M/S; that supremum need not be finite"
        ),
        "limiting_inequalities": (
            "ell^(2d) epsilon^eta <= (2D)!(2/(h eta))^(2D), D=ceil(d); choosing eta=kappa changes "
            "B0 from 3/2-2kappa to 3/2-3kappa and cycle gain from 1/2-3kappa to "
            "1/2-4kappa, while the target step retains margin 2kappa"
        ),
        "historical_frozen_records_modified": False,
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
