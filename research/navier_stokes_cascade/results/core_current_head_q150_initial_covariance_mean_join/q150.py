"""Join Q149 to the audited initial complete covariance/mean correction."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
IDENTITY = "Q149-physical-primary-to-initial-complete-covariance-mean"
PROVIDERS = (
    "complete_covariance_initial_means",
    "core_covariance_mean_recheck",
    "core_dynamic_mean_recheck",
    "core_multiharmonic_mean_residual_budget",
)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


multi = load("q150_multiharmonic", "core_multiharmonic_mean_residual_budget/budgets.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def provider_flags():
    return {name: evidence(name).get("flags", {}) for name in PROVIDERS}


def initial_join(q149_record, identity):
    """Bind Q149's actual physical primaries to the initial mean theorem."""
    if identity != IDENTITY:
        raise ValueError("Q150 same-primary identity required")
    if q149_record.get("status") != "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined":
        raise ValueError("Q149 uniform physical primary family required")
    if q149_record.get("historical_frozen_records_modified"):
        raise ValueError("Historical evidence must remain immutable")
    flags = provider_flags()
    required = {
        "complete_covariance_initial_means": (
            "actual_primary_full_six_entry_covariance_initial_mean_inputs_joined",
            "initial_fast_update_has_exact_zero_Haar_velocity",
            "initial_recomputed_three_compatibility_defects_in_S2",
            "initial_angular_mean_annular_residual_in_M2_minus_3kappa_modulo_flat",
            "complete_axial_velocity_and_moving_pressure_moment_retained",
            "both_background_pressure_tails_retained_via_same_base_subtraction",
        ),
        "core_covariance_mean_recheck": (
            "complete_covariance_axial_velocity_and_curl_squares_retained",
            "both_pressure_tails_and_dynamic_mean_terms_retained",
            "initial_nonlinear_defects_recomputed_after_linear_five_map",
            "joint_pressure_axial_dynamic_cancellation_gate_derived",
        ),
        "core_dynamic_mean_recheck": (
            "all_six_covariance_entries_kept_in_moment_or_target_flux",
            "nonzero_old_Haar_means_and_complete_covariance_changes_retained",
            "whole_moving_pressure_product_and_extra_derivative_cost_retained",
            "old_mean_and_same_phase_wave_cross_rate_controls_derived",
        ),
        "core_multiharmonic_mean_residual_budget": (
            "exact_wave_cancellation_and_both_fresh_mean_wave_feedback_orientations_retained",
            "zero_mode_full_pressure_and_all_three_scalars_reconstructed",
            "actual_mean_feedback_order_D_is_precise_remaining_inequality",
            "class_level_phase_advection_witness_shows_coarse_size_alone_insufficient",
        ),
    }
    for name, keys in required.items():
        if not all(flags[name].get(key) for key in keys):
            raise ValueError("Missing provider contract: " + name)
    return {
        "stage": "same-physical-primary-initial-complete-covariance-mean-join",
        "correction": "Q149's negative-power covariance term blocks absolute decay but is canceled in the already audited exact initial covariance/mean identities",
        "covariance": "all six W_ab entries, including axial velocity and complete curl remainders",
        "pressure": "Pi1=K1-mu^2/(4X), Pi2=K2, compact p_m, moving c_rho*P product and both background tails via the same exact base subtraction",
        "initial_fast_means": "exactly zero after the zero-Haar temporal inverse and compact curl",
        "initial_recomputed_defects": "P,Jtheta,Jz belong to S_2 and are generally nonzero",
        "initial_angular_mean_residual": "M_(2-3*kappa) plus separately retained flat errors",
        "uniform_domain": "every sufficiently late actual physical primary band from Q149, for every fixed coefficient derivative order",
        "later_signed_dynamic_cycle_certified": False,
        "complete_residual_or_blowup_certified": False,
    }


def initial_rates(h, kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q150 identity required")
    h, kappa = Q(h), Q(kappa)
    if not (0 < h < Q(1, 100) and 0 < kappa < Q(1, 4)):
        raise ValueError("Inherited h and kappa ranges required")
    desired = Q(3, 2) - 2 * kappa
    initial = 2 - 3 * kappa
    return {
        "epsilon_initial_residual_order": initial,
        "Q_initial_residual_exponent": h * initial,
        "later_desired_order_D": desired,
        "initial_margin_over_D": initial - desired,
        "positive_initial_margin": initial > desired,
        "transfer_to_later_nonzero_mean_state": False,
    }


def frontier(kappa=Q(1, 100000), identity=IDENTITY):
    """The precise first open later-cycle inequality, not a new premise."""
    if identity != IDENTITY:
        raise ValueError("Q150 identity required")
    rates = multi.rates(Q(kappa))
    return {
        "B": rates["B"],
        "C": rates["C"],
        "D": rates["desired"],
        "coarse_feedback_order": rates["guaranteed"],
        "missing_gain": rates["desired"] - rates["guaranteed"],
        "inequality": "||Pperp[N(d,wcur)+N(wcur,d)]||_(n,half-weight,9p) <= epsilon^D Q_n(S)",
        "for_every_fixed_n_on_original_accepted_band": True,
        "coarse_mean_size_alone_suffices": False,
        "sufficient_unproved_route": "post-wave zero-Haar fast rows at D and fresh postfast P,Jtheta,Jz at C, including the extra Jtheta derivative",
        "class_witness_is_actual_prescribed_map_output": False,
    }


def record():
    return {
        "status": "Q149-uniform-primary-initial-complete-covariance-mean-correction-joined",
        "correction": "the initial complete covariance/mean step was already audited; Q149's negative absolute power is canceled there by exact six-entry covariance, compact pressure and zero-Haar identities",
        "established": "on Q149's uniform physical primary family, the initial fast/slow mean update retains axial velocity, both pressure tails and moving pressure, gives recomputed P,Jtheta,Jz in S2, and bounds the angular-mean residual by M_(2-3kappa) plus flat errors",
        "scope": "the initial zero-old-mean stage only; later signed waves have nonzero old means and changing covariances, so the initial cancellation cannot be reused",
        "obstruction": "the first open later-cycle estimate is the actual fresh mean/current-wave feedback inequality at D=3/2-2kappa; the available coarse order B=1-kappa misses D by 1/2-kappa",
        "next": "prove that actual feedback inequality, or its stated sharper post-wave row/scalar sufficient bounds, while retaining complete zero-mode pressure, axial moments, supports, viscosity, cutoffs and all interactions",
        "later_signed_dynamic_cycle_or_complete_residual_certified": False,
        "blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
