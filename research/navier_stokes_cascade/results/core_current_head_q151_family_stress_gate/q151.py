"""Correct Q150's frontier and isolate the physical/canonical stress gate."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q150-initial-physical-primary-to-canonical-later-cycle-gate"
PROVIDERS = (
    "core_actual_scalar_gain_join",
    "core_fast_row_gain_join",
    "core_haar_row_stress_join",
    "core_localized_stress_assembly",
    "canonical_primary_foundation",
    "core_canonical_full_path_ode_join",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def provider_flags():
    return {name: evidence(name).get("flags", {}) for name in PROVIDERS}


def conditional_join(q150_record, identity):
    """Separate the closed canonical feedback algebra from family transfer."""
    if identity != IDENTITY:
        raise ValueError("Q151 family identity required")
    if q150_record.get("status") != "Q149-uniform-primary-initial-complete-covariance-mean-correction-joined":
        raise ValueError("Q150 initial covariance/mean join required")
    if q150_record.get("historical_frozen_records_modified"):
        raise ValueError("Historical evidence must remain immutable")
    flags = provider_flags()
    required_true = {
        "core_actual_scalar_gain_join": (
            "postfast_scalar_C_composed_without_new_incoming_scalar_callback",
            "fresh_slow_C_size_follows_from_same_actual_five_map",
            "complete_axial_velocity_covariance_and_weighted_pressure_moment_retained",
            "whole_moving_compensator_product_and_actual_fast_A1_fields_retained",
        ),
        "core_fast_row_gain_join": (
            "postwave_zero_Haar_fast_rows_D_and_fresh_actual_fast_D_derived_conditionally",
            "conditional_mean_feedback_D_and_next_regular_source_D_follow",
            "complete_sourced_covariance_and_proper_operator_poles_retained",
            "all_pressure_viscosity_cylindrical_connections_and_mean_crosses_retained",
            "large_finite_A1_source_and_raw_kernel_dependencies_explicit",
        ),
        "core_haar_row_stress_join": (
            "same_field_canonical_identity_suffices_conditionally_for_Haar_C",
            "two_different_primary_families_need_not_be_identified_for_canonical_route",
        ),
        "core_localized_stress_assembly": (
            "canonical_positive_roots_and_physical_squared_partition_match_same_stress",
            "full_regular_postwave_tangential_mean_rows_D_derived_conditionally",
            "precise_postwave_stage_full_pressure_and_high_jet_A1_recipe_retained",
        ),
        "canonical_primary_foundation": (
            "canonical_principal_roots_rebuilt_not_reused_from_inherited_frame",
            "canonical_manuscript_reference_frame_and_growing_entrance_used",
            "same_operator_distinct_entrances_have_Gaussian_relative_difference_bound",
        ),
        "core_canonical_full_path_ode_join": (
            "same_rule_defined_background_tail_and_original_physical_representatives_used",
            "whole_path_raw_solution_bounds_and_gaussian_ratio_cancellation_joined",
            "full_moving_normal_pressure_and_potential_retained",
        ),
    }
    for name, keys in required_true.items():
        if not all(flags[name].get(k) for k in keys):
            raise ValueError("Missing provider contract: " + name)
    required_false = {
        "core_fast_row_gain_join": (
            "original_actual_support_or_primary_frame_contract_independently_certified",
            "high_jet_actual_input_or_raw_kernel_polynomials_numerically_enclosed",
        ),
        "core_haar_row_stress_join": (
            "actual_whole_support_principal_stress_identity_or_frame_transfer_certified",
        ),
        "core_localized_stress_assembly": (
            "two_different_homogeneous_families_identified_by_closeness",
            "independent_actual_column_root_profile_or_support_inputs_certified",
        ),
        "canonical_primary_foundation": (
            "canonical_and_inherited_physical_primaries_or_roots_exactly_equal",
        ),
        "core_canonical_full_path_ode_join": (
            "canonical_and_inherited_primary_families_identified",
        ),
    }
    for name, keys in required_false.items():
        if any(flags[name].get(k) for k in keys):
            raise ValueError("Invalid family promotion: " + name)
    return {
        "canonical_feedback_inequality_D": "proved conditionally for the consistently rebuilt canonical family",
        "canonical_postfast_scalars": "P,Jtheta,Jz at C with full axial covariance and moving compact pressure",
        "canonical_fresh_fast": "tangential D and radial D+1, retaining the exact A1 remainder",
        "canonical_fresh_slow": "tangential C and radial C+1",
        "canonical_next_regular_source": "D with pole 9p and additive tails excluded",
        "q150_inherited_physical_family_transferred_to_canonical_cycle": False,
        "reason": "the inherited physical entrance/frame/roots are not the rebuilt canonical entrance/frame/roots",
        "complete_residual_or_blowup_certified": False,
    }


def rates(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q151 identity required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 4):
        raise ValueError("Inherited kappa range required")
    B = 1 - kappa
    C = B + Q(1, 2)
    D = C - kappa
    return {
        "B": B,
        "C": C,
        "D": D,
        "D_minus_B": D - B,
        "D_minus_1": D - 1,
        "canonical_fast_feedback": D,
        "canonical_slow_feedback": C,
        "canonical_combined_feedback": min(D, C),
    }


def stress_gate(h, kappa=Q(1, 100000), polynomial_power=2, identity=IDENTITY):
    """Exact physical-family principal stress gate and closeness obstruction."""
    if identity != IDENTITY:
        raise ValueError("Q151 identity required")
    h, polynomial_power = Q(h), int(polynomial_power)
    r = rates(kappa, identity)
    if not (0 < h < Q(1, 100) and polynomial_power >= 0):
        raise ValueError("Inherited h and nonnegative polynomial power required")
    exponent = h * r["D_minus_1"]
    return {
        "principal_discrepancy": "L_e=bar(W0_principal)_(r,e)-Sigma0_e",
        "required_inequality": "||epsilon^-1 L_e||_n <= epsilon^(D-1) Q_n(S) zeta delta^-p",
        "D_minus_1": r["D_minus_1"],
        "exponential_rate_in_ell": exponent,
        "positive_rate": exponent > 0,
        "closeness_model": f"1/ell^{polynomial_power}",
        "required_ratio": f"exp({exponent}*ell)/ell^{polynomial_power}",
        "polynomial_closeness_suffices": False,
        "exact_same_family_stress_identity_would_close_gate": True,
    }


def a1_request(n, kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY or type(n) is not int or n < 0:
        raise ValueError("Q151 identity and nonnegative integer order required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 4):
        raise ValueError("Inherited kappa range required")
    B = 1 - kappa
    H = 2 - 3 * kappa
    q = (H - B) // kappa + 1
    J = n + 2 * q + 3
    return {
        "compact_flatness_iterations": q,
        "old_prefast_source_order": J,
        "order_zero_fresh_fast_old_source": 200006 if kappa == Q(1, 100000) and n == 0 else None,
        "finite_for_each_fixed_output_order": True,
        "same_family_raw_kernel_and_source_recipes_required": True,
        "numerical_materialization_required_for_continuum_implication": False,
        "actual_same_family_enclosure_certified_here": False,
    }


def record():
    return {
        "status": "Q150-frontier-corrected-to-physical-canonical-principal-stress-gate",
        "correction": "the fresh mean/current-wave feedback estimate at D already follows conditionally for the consistently rebuilt canonical family; it is not available for Q150's distinct inherited physical family",
        "established": "canonical exact row subtraction plus scalar C and fast-row D estimates close the finite regular feedback/source algebra with complete covariance, axial velocity, full pressure and exact A1 memory",
        "scope": "a conditional canonical finite-stage estimate; no transfer to the inherited physical primary and no complete residual estimate",
        "obstruction": "either prove the exact same-family physical principal-stress inequality ||epsilon^-1 L_e||<=epsilon^(D-1)Q_n(S)zeta delta^-p, or carry the rebuilt canonical family through independently validated uniform inputs; O(1/L) frame/root closeness is insufficient",
        "next": "choose one route without merging families: audit the inherited physical stress identity on its full localized support, or bind the canonical columns/roots/high-jet compact inputs and supports to the finalized background before the complete tail/residual join",
        "inherited_and_canonical_primary_families_identified": False,
        "complete_residual_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
