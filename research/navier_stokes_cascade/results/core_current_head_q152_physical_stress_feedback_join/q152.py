"""Close Q151's physical stress gate from the inherited same-family assembly."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q151-inherited-physical-same-family-stress-feedback-join"
PROVIDERS = (
    "core_current_head_q149_primary_pulse_join",
    "actual_band_pulse_bounds",
    "physical_pulse_foundation",
    "complete_covariance_initial_means",
    "core_actual_scalar_gain_join",
    "core_fast_row_gain_join",
    "core_haar_row_stress_join",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def rates(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q152 identity required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 4):
        raise ValueError("Inherited kappa range required")
    B = 1 - kappa
    C = B + Q(1, 2)
    D = C - kappa
    R = 2 - 2 * kappa
    return {
        "B": B,
        "C": C,
        "D": D,
        "principal_complete_curl_remainder": R,
        "required_stress_after_epsilon_inverse": D - 1,
        "available_complete_curl_after_epsilon_inverse": R - 1,
        "complete_curl_margin": (R - 1) - (D - 1),
        "feedback": min(C, D),
    }


def scale_identity(h=Q(1, 200), identity=IDENTITY):
    """Check the physical square-scale identity symbolically in exponents."""
    if identity != IDENTITY:
        raise ValueError("Q152 identity required")
    h = Q(h)
    if not 0 < h < Q(1, 100):
        raise ValueError("Inherited h range required")
    A = Q(1, 2) + h
    q_power = -(A + Q(1, 2))
    band_power = -2 * A + h + A + Q(1, 2)
    return {
        "A": A,
        "a_scale": "Q_beta^(-2A) epsilon_beta",
        "target": "T_beta=(Q_beta/q)^(A+1/2) T0",
        "stress": "Sigma0=q^(-A-1/2) T0",
        "band_power": band_power,
        "q_power": q_power,
        "exact": band_power == 0,
    }


def principal_assembly(labels, sigma0, identity=IDENTITY):
    """Finite exact model of the same-family H y=T and squared cover identity."""
    if identity != IDENTITY or not labels:
        raise ValueError("Q152 identity and nonempty physical labels required")
    actual = Q(0)
    cover = Q(0)
    root_defect = Q(0)
    for item in labels:
        eta2, a_scale, H, y, target = map(Q, item)
        cover += eta2
        actual += eta2 * a_scale * H * y
        root_defect += eta2 * a_scale * (H * y - target)
    sigma0 = Q(sigma0)
    return {
        "actual": actual,
        "cover": cover,
        "root_defect": root_defect,
        "defect": actual - sigma0,
        "exact": actual == sigma0 and cover == 1 and root_defect == 0,
    }


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q152 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    q149 = data["core_current_head_q149_primary_pulse_join"]["symbolic"]
    if q149["status"] != "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined":
        raise ValueError("Q149 physical family required")
    required = {
        "actual_band_pulse_bounds": (
            "conditional_primary_covariance_positive_inverse_and_edge_bounds_supplied",
            "conditional_real_primary_curl_covariance_remainder_gain_supplied",
            "conditional_signed_inverse_at_fixed_primary_supplied",
        ),
        "physical_pulse_foundation": (
            "actual_background_inputs_joined_to_primary_pulse_frame_and_covariance",
            "one_starting_band_for_all_later_scales_and_every_fixed_physical_derivative_order",
            "uniform_primary_complete_curl_pressure_and_fixed_potential_bounds_derived",
            "exact_closed_supports_and_one_sided_terminal_jets_preserved",
        ),
        "complete_covariance_initial_means": (
            "actual_primary_full_six_entry_covariance_initial_mean_inputs_joined",
            "both_background_pressure_tails_retained_via_same_base_subtraction",
            "complete_axial_velocity_and_moving_pressure_moment_retained",
        ),
        "core_actual_scalar_gain_join": (
            "postfast_scalar_C_composed_without_new_incoming_scalar_callback",
            "fresh_slow_C_size_follows_from_same_actual_five_map",
            "whole_moving_compensator_product_and_actual_fast_A1_fields_retained",
        ),
        "core_fast_row_gain_join": (
            "postwave_zero_Haar_fast_rows_D_and_fresh_actual_fast_D_derived_conditionally",
            "conditional_mean_feedback_D_and_next_regular_source_D_follow",
            "all_pressure_viscosity_cylindrical_connections_and_mean_crosses_retained",
            "original_exact_additive_tails_excluded_from_sourced_inverse",
        ),
        "core_haar_row_stress_join": (
            "complete_Haar_row_envelope_with_all_mean_correlations_and_pressure_derived",
            "principal_stress_gap_and_fine_flux_C_separated_explicitly",
            "proper_poles_profiles_and_exact_tail_ledger_unchanged",
        ),
    }
    for name, keys in required.items():
        flags = data[name].get("flags", {})
        if not all(flags.get(key) for key in keys):
            raise ValueError("Missing provider contract: " + name)
    return {
        "same_physical_family": True,
        "canonical_and_physical_entrances_identified": False,
        "principal_stress_identity": "H_beta y_beta=T_beta, a_scale T_beta=Sigma0, sum eta_beta^2=1",
        "complete_curl_selected_remainder": "M_(2-2kappa)",
        "postfast_scalars": "C",
        "fresh_fast": "D retaining the exact A1 remainder",
        "fresh_slow": "C",
        "mean_feedback": "D",
        "next_regular_source": "D with pole 9p; additive tails excluded",
    }


def stress_gate(kappa=Q(1, 100000), identity=IDENTITY):
    r = rates(kappa, identity)
    return {
        "principal_discrepancy": "L_e=bar(W0_principal)_(r,e)-Sigma0_e=0",
        "required": "||epsilon^-1 L_e||_n <= epsilon^(D-1) Q_n(S) zeta delta^-p",
        "principal_gate_closed_exactly": True,
        "complete_curl_alternative_order": r["available_complete_curl_after_epsilon_inverse"],
        "required_order": r["required_stress_after_epsilon_inverse"],
        "complete_curl_margin": r["complete_curl_margin"],
        "polynomial_family_closeness_used": False,
    }


def record():
    return {
        "status": "Q152-physical-principal-stress-gate-closed-and-feedback-D-joined",
        "correction": "Q151 correctly rejected transfer by physical/canonical frame closeness, but overlooked the inherited physical family's own H_beta y_beta=T_beta root equation, physical square scaling and squared partition identity",
        "established": "on Q149's inherited physical family, the localized principal selected stress equals Sigma0 exactly; the complete-curl selected remainder has order 2-2kappa, so the stress gate and conditional mean-feedback/source order D follow with full pressure, axial covariance, viscosity and A1 memory",
        "scope": "same-family conditional finite-stage join; it does not identify the physical and canonical entrances or prove the generated-source endpoint class, additive physical tails, a complete residual theorem or blow-up",
        "next": "prove every generated regular harmonic has the physical full-path Gaussian/input bounds and smooth endpoint extension at the exposed high jets, and bound every excluded additive Gaussian/spatial/scale/A0/A1 tail after physical differentiation before invoking the complete residual diagonal",
        "generated_source_endpoint_and_complete_tail_join_certified": False,
        "canonical_and_physical_primary_families_identified": False,
        "complete_residual_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
