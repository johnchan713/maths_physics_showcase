"""Q190: bind Q189 through Q177--Q185 and correct the live frontier."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q190-DOWNSTREAM-ROUTE-FRONTIER-AUDIT"

PROVIDERS = (
    "core_current_head_q189_same_profile_envelope_gap_audit",
    "core_current_head_q164_core_fixed_point_pressure_audit",
    "core_current_head_q165_finite_map_frontier_audit",
    "core_current_head_q166_complete_residual_frontier_audit",
    "core_current_head_q177_weighted_support_source_reconciliation",
    "core_current_head_q178_physical_primary_source_audit",
    "core_current_head_q179_covariance_mean_cycle_frontier",
    "core_current_head_q180_target_source_simultaneous_induction",
    "core_current_head_q181_diagonal_endpoint_force_replay",
    "core_current_head_q182_comparison_core_exclusion_replay",
    "core_current_head_q183_provider_discharge_generated_source_frontier",
    "core_current_head_q184_generated_source_log_loss_correction",
    "core_current_head_q185_corrected_rate_diagonal_comparison_replay",
    "physical_pulse_foundation",
)

EXPECTED = {
    PROVIDERS[0]: "Q189-same-profile-envelope-and-gap-join-closed-relative-to-region-theorems",
    PROVIDERS[1]: "Q164-all-index-core-contraction-and-leading-pressure-audited",
    PROVIDERS[2]: "Q165-current-head-finite-map-and-weighted-source-chain-audited",
    PROVIDERS[3]: "Q166-current-head-complete-six-group-local-residual-audited-conditionally",
    PROVIDERS[4]: "Q177-selected-input-weighted-stress-and-reserved-support-route-reconciled",
    PROVIDERS[5]: "Q178-physical-primary-pulse-source-boundary-and-covariance-obstruction-audited",
    PROVIDERS[6]: "Q179-initial-and-first-feedback-covariance-means-closed-P3-universal-target-obstructed",
    PROVIDERS[7]: "Q180-physical-target-source-tail-simultaneous-induction-reconciled",
    PROVIDERS[8]: "Q181-physical-diagonal-endpoint-localized-force-replayed-conditionally",
    PROVIDERS[9]: "Q182-physical-pressure-comparison-protected-core-replayed-conditionally",
    PROVIDERS[10]: "Q183-current-head-provider-frontier-reconciled-to-generated-source-supremum",
    PROVIDERS[11]: "Q184-generated-source-log-loss-corrected-and-pure-power-reserve-derived",
    PROVIDERS[12]: "Q185-corrected-rate-diagonal-localization-comparison-replayed-conditionally",
    PROVIDERS[13]: "actual-primary-physical-pulses-all-fixed-jets-and-common-domain-tails-relative-to-exact-core-and-finite-lemmas",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q190 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def _ceil(value):
    value = Q(value)
    return -(-value.numerator // value.denominator)


def comparison_stage(target, loss_L=0, loss_H=0, route="log", identity=IDENTITY):
    """Reconstruct the least Q185 fixed comparison stage exactly."""
    if identity != IDENTITY:
        raise ValueError("Q190 identity required")
    values = (target, loss_L, loss_H)
    if any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in values):
        raise ValueError("nonnegative integer request and losses required")
    h, kappa = Q(1, 200), Q(1, 100000)
    if route == "log":
        delta = Q(1, 2) - 3 * kappa
    elif route == "pure":
        delta = Q(1, 2) - 4 * kappa
    else:
        raise ValueError("route must be log or pure")
    required = Q(target + 1 + loss_L + loss_H)
    J0 = max(0, _ceil(2 * required / (h * delta)) - 1)
    J = J0 + 1
    threshold = h * delta * (J0 + 1) / 2 - loss_L - loss_H
    replay = h * delta * (J + 1) / 2 - loss_L - loss_H
    return {
        "route": route,
        "h": h,
        "kappa": kappa,
        "delta": delta,
        "target": target,
        "loss_L": loss_L,
        "loss_H": loss_H,
        "J0": J0,
        "replay_stage": J,
        "threshold": threshold,
        "retained_margin": h * delta / 2,
        "threshold_holds": threshold >= target + 1,
        "replay_holds_with_margin": replay >= Q(target + 1) + h * delta / 2,
    }


def quantifier_order(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q190 identity required")
    log = comparison_stage(7, 11, 13, "log")
    pure = comparison_stage(7, 11, 13, "pure")
    return {
        "order": (
            "fix Cartesian jet m and target N; choose a finite comparison stage J; "
            "collect the finitely many constants through J; choose one later band "
            "beyond the Q2, primary-pulse and finite-stage thresholds; then send q to zero"
        ),
        "comparison_stage_depends_on_q": False,
        "stage_uniform_constant_required": False,
        "log_route": log,
        "pure_route": pure,
        "log_margin": Q(49997, 40000000),
        "pure_margin": Q(49996, 40000000),
        "A0_A1_each_own_margin": Q(1, 40000000),
        "all_strict": (
            log["threshold_holds"]
            and log["replay_holds_with_margin"]
            and pure["threshold_holds"]
            and pure["replay_holds_with_margin"]
            and log["retained_margin"] == Q(49997, 40000000)
            and pure["retained_margin"] == Q(49996, 40000000)
        ),
    }


def downstream_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q190 identity required")
    q189 = evidence(PROVIDERS[0])["review"]
    q164 = evidence(PROVIDERS[1])["symbolic"]
    q165 = evidence(PROVIDERS[2])["symbolic"]
    q177 = evidence(PROVIDERS[4])["review"]
    q178 = evidence(PROVIDERS[5])["review"]
    q180 = evidence(PROVIDERS[7])["review"]
    q181 = evidence(PROVIDERS[8])["review"]
    q182 = evidence(PROVIDERS[9])["review"]
    q184 = evidence(PROVIDERS[11])["review"]
    q185 = evidence(PROVIDERS[12])["review"]
    pulse = evidence(PROVIDERS[13])["flags"]
    return {
        "same_profile_envelope_and_gaps": q189["checks"]["same_selected_full_pressure_profile"],
        "exact_core_and_leading_pressure": q164["U1_core_contraction_and_leading_pressure_closed"],
        "finite_maps_and_weighted_sources": q165["U2_finite_map_and_weighted_source_chain_closed"],
        "weighted_stress_coefficients_and_reserved_supports": all(q177["established"].values()),
        "uniform_primary_all_scale_fixed_jet_bounds": (
            q178["established"]["conditional_uniform_actual_primary_pulse_values_fixed_derivatives_and_supports"]
            and pulse["uniform_primary_complete_curl_pressure_and_fixed_potential_bounds_derived"]
            and pulse["one_starting_band_for_all_later_scales_and_every_fixed_physical_derivative_order"]
        ),
        "complete_covariance_and_generated_source_cycle": (
            q180["established"]["physical_P3_target_and_both_moments_for_each_finite_cycle"]
            and q180["established"]["physical_P4_generated_source_quotient_and_endpoint_extension"]
            and q180["established"]["physical_P5_fixed_stage_tail_for_each_fixed_jet_and_power"]
        ),
        "log_source_contract_and_corrected_pure_contract": (
            q184["scope"]["log_renormalized_fixed_cycle_source_contract_closed_conditionally"]
            and q184["scope"]["pure_power_contract_closed_after_one_kappa_reserve"]
        ),
        "complete_residual_endpoint_localization_force": (
            q181["checks"]["complete_six_group_residual_is_flat_for_every_fixed_jet"]
            and q181["checks"]["localization_has_exact_twelve_groups"]
            and q185["checks"]["complete_residual_endpoint_and_localization_memory_survives"]
        ),
        "fixed_gauge_comparison_and_protected_core": (
            q182["scope"]["P7_same_force_competitor_exclusion_closed_conditionally"]
            and q185["checks"]["comparison_and_core_are_rate_independent_after_rebuild"]
        ),
        "one_selected_profile_and_log_route_object": True,
    }


def correction_and_frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q190 identity required")
    q185 = evidence(PROVIDERS[12])["review"]
    return {
        "Q185_named_exact_core_and_finite_maps_as_next": not q185["scope"][
            "exact_core_and_finite_existence_inverse_gap_lemmas_independently_reproved"
        ],
        "Q164_Q165_and_Q189_now_discharge_that_provider_level_frontier": True,
        "missing_scalar_selector_or_diagonal_inequality_through_Q185": False,
        "missing_named_analytical_provider_through_Q185": False,
        "all_provider_theorems_rederived_from_initial_equations_in_Q190": False,
        "first_fresh_independent_rederivation": (
            "the physical_pulse_foundation all-band complete-curl, pressure and potential "
            "bounds, including the frame/graph coefficient estimates and smooth closed support"
        ),
        "blocking_supremum": (
            "for each fixed a,b, bound the normalized suprema of grad^a dt^b(u_p,p_p,A_p) "
            "over every accepted band by finite constants with losses "
            "Q^{-((a+1)sx+b st)}, Q^{-((a+2)sx+b st)}, Q^{-(a sx+b st)} respectively, "
            "and prove smooth zero extension on the same supports"
        ),
        "historical_first_band_or_giant_profile_numerically_certified": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q190 identity required")
    status = statuses()
    chain = downstream_contract()
    quantifiers = quantifier_order()
    frontier = correction_and_frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "Q189_object_reaches_every_downstream_provider": all(chain.values()),
        "both_corrected_diagonal_routes_have_exact_positive_margin": quantifiers["all_strict"],
        "Q185_frontier_is_stale_at_current_head": (
            frontier["Q185_named_exact_core_and_finite_maps_as_next"]
            and frontier["Q164_Q165_and_Q189_now_discharge_that_provider_level_frontier"]
        ),
        "no_missing_provider_is_not_an_unconditional_proof": (
            not frontier["missing_named_analytical_provider_through_Q185"]
            and not frontier["all_provider_theorems_rederived_from_initial_equations_in_Q190"]
        ),
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "downstream_contract": chain,
        "quantifier_order": quantifiers,
        "correction_and_frontier": frontier,
        "scope": {
            "Q189_through_Q185_provider_route_closed_conditionally": True,
            "same_object_complete_pressure_log_route_retained": True,
            "all_primitive_continuum_providers_independently_reproved_here": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q190 identity required")
    return {
        "status": "Q190-Q189-through-Q185-downstream-route-closed-conditionally-and-frontier-corrected",
        "established": (
            "Q189's same-profile full-pressure envelope/gap object composes through the Q177--Q185 "
            "primary-pulse, covariance/mean, generated-source, residual, localization and comparison "
            "providers with the exact fixed-request quantifier order and positive log/pure margins"
        ),
        "correction": (
            "Q185's next exact-core/finite-map frontier is stale after Q164--Q165 and Q187--Q189; "
            "there is no remaining scalar selector or named-provider gap through Q185, but this "
            "composition is not an independent reproof of every continuum provider"
        ),
        "next": correction_and_frontier()["first_fresh_independent_rederivation"],
        "unconditional_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
