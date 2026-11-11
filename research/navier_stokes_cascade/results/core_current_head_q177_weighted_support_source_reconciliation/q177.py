"""Q177 selected-input, weighted-stress and reserved-support source audit."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q177-WEIGHTED-SUPPORT-SOURCE-RECONCILIATION"
PROVIDERS = (
    "core_current_head_q176_modulation_repair_dependency_audit",
    "prescribed_continuation_foundation",
    "compact_jet_envelope",
    "loop_modulation_bounds",
    "repair_state_bounds",
    "regional_cone_join",
    "annular_exterior_foundation",
    "global_pressure_moment_foundation",
    "weighted_stress_foundation",
    "coefficient_support_foundation",
)
EXPECTED = {
    "core_current_head_q176_modulation_repair_dependency_audit": "Q176-compact-modulation-exact-repair-audited-with-refined-dependency-order",
    "prescribed_continuation_foundation": "pre-modulation-prescribed-continuation-through-shape-all-fixed-jets",
    "compact_jet_envelope": "actual-compact-input-envelope-bounded",
    "loop_modulation_bounds": "fixed-phase-loop-and-modulation-errors-bounded",
    "repair_state_bounds": "repair-state-bound-and-conditional-frequency-acceptance",
    "regional_cone_join": "regional-input-gaps-and-leading-value-cone-joined",
    "annular_exterior_foundation": "annular-exterior-all-order-leading-regularity-and-global-moment-foundation-relative-to-exact-roots",
    "global_pressure_moment_foundation": "current-same-root-global-pressure-moments-and-leading-exterior-stresses-joined-relative-to-retained-inputs",
    "weighted_stress_foundation": "actual-weighted-first-stress-and-background-foundation-joined-relative-to-exact-core-and-finite-lemmas",
    "coefficient_support_foundation": "current-linear-coefficient-jets-reserved-velocity-and-exceptional-pressure-supports-joined-relative-to-leading-inputs",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q177 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def frequency_acceptance_gate(G=3):
    if type(G) is not int or G < 1:
        raise ValueError("positive integer G required")
    # The repair proof bounds every logarithmic threshold by 16 G^4 + 193 G.
    dominant_slack = 16 * G ** 4 - 193 * G
    cone_slack = G ** 30 - 512 * 4 ** 10
    return {
        "G": G,
        "log_threshold_slack": dominant_slack,
        "log_threshold_below_32_G4": dominant_slack > 0,
        "cone_absorption_slack": cone_slack,
        "cone_absorption": cone_slack > 0,
    }


def support_geometry(i3_right=Q(51, 20), velocity_edge=3, i4_left=6, i4_right=11):
    values = tuple(map(Q, (i3_right, velocity_edge, i4_left, i4_right)))
    if any(value <= 0 for value in values):
        raise ValueError("positive support coordinates required")
    i3, vel, left, right = values
    return {
        "i3_right": i3,
        "velocity_edge": vel,
        "i4_left": left,
        "i4_right": right,
        "strict_separation": i3 < vel < left < right,
        "reserved_gap": left - vel,
    }


def linear_inverse_gate(U_inverse=Q(2), E_inverse=Q(393781, 500), target=800):
    U_inverse, E_inverse, target = map(Q, (U_inverse, E_inverse, target))
    if min(U_inverse, E_inverse, target) <= 0:
        raise ValueError("positive inverse data required")
    return {
        "U_inverse_upper": U_inverse,
        "E_inverse_upper": E_inverse,
        "target": target,
        "max_inverse_upper": max(U_inverse, E_inverse),
        "strict_inverse_gate": max(U_inverse, E_inverse) < target,
    }


def outer_quotient_jets(Xb):
    Xb = Q(Xb)
    if Xb <= 0:
        raise ValueError("positive outer radius required")
    return {
        "cubic_coefficient": Xb / 8,
        "third_delta_derivative": 3 * Xb / 4,
        "second_eta_derivative": 4 * Xb,
    }


def reserved_pressure(mu, X, K1, K2):
    mu, X, K1, K2 = map(Q, (mu, X, K1, K2))
    if X <= 0:
        raise ValueError("positive radius required")
    return {
        "Pi1": K1 - mu * mu / (4 * X),
        "Pi2": K2,
        "higher_pressure": Q(0),
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q177 identity required")
    status = statuses()
    q176 = evidence("core_current_head_q176_modulation_repair_dependency_audit")
    continuation = evidence("prescribed_continuation_foundation")["flags"]
    envelope = evidence("compact_jet_envelope")["flags"]
    loop = evidence("loop_modulation_bounds")["flags"]
    repair = evidence("repair_state_bounds")["flags"]
    regional = evidence("regional_cone_join")["flags"]
    exterior = evidence("annular_exterior_foundation")["flags"]
    pressure = evidence("global_pressure_moment_foundation")["flags"]
    weighted = evidence("weighted_stress_foundation")["flags"]
    coefficients = evidence("coefficient_support_foundation")["flags"]
    frequency = frequency_acceptance_gate(4)
    supports = support_geometry()
    inverse = linear_inverse_gate()
    quotient = outer_quotient_jets(13)
    pressure_fixture = reserved_pressure(2, 5, 7, 11)
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "Q176_left_selected_suprema_as_unrecomputed": not q176["review"]["scope"]["actual_selected_continuum_input_suprema_independently_recomputed_here"],
        "prescribed_continuation_has_all_fixed_jets_relative_to_exact_core": continuation["every_fixed_jet_of_prescribed_continuation_through_shape_bounded_relative_to_exact_core"],
        "compact_envelope_is_analytically_bounded": envelope["actual_compact_jet_envelope_verified"],
        "loop_and_repair_constants_are_supplied_in_sequence": loop["actual_modulation_error_constants_bounded"] and repair["actual_repair_error_constants_bounded"] and repair["positive_repair_neighborhood_bounded"],
        "selected_frequency_is_conditionally_accepted": repair["frequency_acceptance_given_inherited_estimates"] and frequency["log_threshold_below_32_G4"],
        "same_frequency_regional_cone_is_joined": regional["compact_input_cone_hypotheses_discharged"] and regional["same_finite_frequency_leading_value_cone_joined"] and frequency["cone_absorption"],
        "whole_compact_and_exterior_inputs_are_composed": exterior["whole_compact_smooth_input_to_same_finite_N_maps_supplied_relative_to_exact_core_and_finite_existence_bounds"],
        "all_five_global_moments_and_same_gauge_pressure_are_joined": pressure["all_five_global_leading_moment_identities_joined"] and pressure["complete_core_backward_pressure_same_gauge_trace_discharged"] and pressure["full_pressure_memory_and_eta_derivative_retained"],
        "weighted_first_stress_is_smooth_at_eta_zero": weighted["T1_over_abs_T0_smooth_on_closed_annulus_including_eta_zero"] and weighted["actual_Gamma_eta_zero_cubic_quotient_Xb_over_eight_derived"] and quotient["cubic_coefficient"] == Q(13, 8),
        "weighted_background_edges_are_joined_at_all_fixed_orders": weighted["all_fixed_weighted_leading_and_first_stress_edges_independently_joined"] and weighted["one_admissible_diagonal_schedule_controls_all_fixed_relative_background_jets"],
        "positive_order_inverse_and_supports_are_joined": coefficients["actual_positive_order_linear_map_and_continuum_inverse_800_supplied"] and coefficients["every_finite_coefficient_support_and_reserved_mean_velocity_joined"] and inverse["strict_inverse_gate"] and supports["strict_separation"],
        "reserved_pressure_and_leading_radial_memory_are_not_erased": coefficients["first_two_pressure_tails_and_nonzero_leading_radial_memory_retained"] and not coefficients["reserved_mean_pressure_equal_to_leading_pressure"] and pressure_fixture["Pi1"] != 0 and pressure_fixture["Pi2"] != 0,
        "external_or_formal_independent_certificate_remains_absent": not continuation["formal_certificate"] and not weighted["formal_certificate"] and not coefficients["formal_certificate"],
        "actual_pulse_and_complete_mean_residual_are_not_promoted": not weighted["unconditional_uniform_actual_physical_pulses_verified"] and not coefficients["unconditional_complete_covariance_mean_residual_or_blowup_verified"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "frequency_acceptance": frequency,
        "support_geometry": {k: str(v) if isinstance(v, Q) else v for k, v in supports.items()},
        "linear_inverse": {k: str(v) if isinstance(v, Q) else v for k, v in inverse.items()},
        "outer_quotient_fixture": {k: str(v) for k, v in quotient.items()},
        "reserved_pressure_fixture": {k: str(v) for k, v in pressure_fixture.items()},
        "correction": {
            "Q176_missing_selected_supremum_route_was_too_broad": True,
            "in_repository_conditional_analytic_route_already_supplies_selected_inputs": True,
            "external_independent_validation_of_exact_core_and_finite_lemmas_still_open": True,
            "historical_frozen_evidence_modified": False,
        },
        "established": {
            "conditional_selected_input_frequency_and_cone_chain": True,
            "same_root_global_pressure_and_five_moments": True,
            "all_fixed_weighted_first_stress_and_background_edges": True,
            "finite_coefficient_and_reserved_velocity_supports_with_pressure_memory": True,
        },
        "scope": {
            "provider_theorems_replaced_by_tests_or_numerical_diagnostics": False,
            "actual_huge_fields_or_scales_materialized": False,
            "unconditional_uniform_actual_pulses": False,
            "complete_covariance_mean_residual_or_blowup": False,
        },
        "next": "source-audit physical_pulse_foundation and the Q149 same-profile join, retaining the sourced-inverse input-class boundary; then audit complete covariance means with radial and axial velocity and Pi1/Pi2",
    }


def record():
    return {
        "status": "Q177-selected-input-weighted-stress-and-reserved-support-route-reconciled",
        "established": "the repository already supplies a conditional analytical selected-input/frequency/cone chain, exact same-gauge five-moment pressure identities, all-fixed weighted first-stress edges, and separated finite coefficient/reserved-velocity supports retaining Pi1/Pi2",
        "correction": "Q176's missing selected-profile supremum frontier was too broad as an in-repository route claim; the missing qualification is independent validation of the retained exact-core and finite existence/inverse/gap lemmas",
        "obstruction": "no test or diagnostic independently certifies those continuum providers, and this checkpoint does not audit the actual primary pulse theorem, sourced nonlinear inverse inputs, complete covariance means or the full residual",
        "historical_frozen_records_modified": False,
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
