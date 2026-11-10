"""Correct Q148's schedule frontier and join uniform primary pulses."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
IDENTITY = "Q148-final-profile-rule-schedule-and-physical-primary-pulse"
PULSE_PROVIDERS = (
    "core_leading_background_borel_join",
    "core_same_band_gap_join",
    "physical_pulse_transfer",
)
FRONTIER_PROVIDERS = (
    "core_canonical_full_path_ode_join",
    "core_fixed_profile_projector_join",
    "core_compact_slow_inverse_domain_join",
)
PROVIDERS = PULSE_PROVIDERS + FRONTIER_PROVIDERS


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


transfer = load("q149_transfer", "physical_pulse_transfer/transfer.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def schedule_row(n, m, nh, log2_C_upper, ell, identity):
    """Exact binary form of the rule-defined diagonal inequality."""
    if identity != IDENTITY or type(n) is not int or type(m) is not int:
        raise ValueError("Q149 identity and integer orders required")
    nh, log2_C_upper, ell = map(Q, (nh, log2_C_upper, ell))
    if not (n >= 1 and 0 <= m <= n and nh > 0 and log2_C_upper >= 0 and ell > 0):
        raise ValueError("Positive schedule data with 0<=m<=n required")
    margin = nh * ell - log2_C_upper - n
    return {
        "binary_margin": margin,
        "strict_padded_rule_holds": margin >= nh,
        "consequence": "C_nm*q^(nh)<=2^-n*2^-nh<2^-n for q<=2^-ell",
        "scale_is_chosen_after_finite_C_nm": True,
        "preexisting_numeric_scale_certificate_required": False,
    }


def provider_flags():
    return {name: evidence(name).get("flags", {}) for name in PROVIDERS}


def close(q148_record, identity):
    """Join the same finalized profile to the audited primary pulse family."""
    if identity != IDENTITY:
        raise ValueError("Q149 same-object identity required")
    if q148_record.get("status") != "Q147-final-profile-finite-coefficient-supports-and-reserved-pressure-obstruction-joined":
        raise ValueError("Q148 finalized-profile coefficient support required")
    if q148_record.get("historical_frozen_records_modified"):
        raise ValueError("Historical evidence must remain immutable")
    flags = provider_flags()
    required = {
        "core_leading_background_borel_join": (
            "same_completed_leading_profile_used_by_background_recursion",
            "rule_defined_borel_sum_identified_with_actual_background",
            "full_pressure_and_curl_cutoff_commutator_retained",
        ),
        "core_same_band_gap_join": (
            "same_rule_defined_background_and_completed_leading_profile_used",
            "same_tail_has_positive_background_spectral_and_same_tilt_cone_gaps",
            "reference_relative_full_phase_normal_floor_joined_on_that_tail",
        ),
        "physical_pulse_transfer": (
            "actual_background_inputs_to_primary_pulse_theorem_discharged",
            "common_all_later_band_primary_pulse_values_and_fixed_jets_derived",
            "complete_primary_curl_velocity_and_wave_pressure_physically_transferred",
            "exact_closed_support_and_one_sided_boundary_conditions_retained",
        ),
    }
    for name, keys in required.items():
        if not all(flags[name].get(key) for key in keys):
            raise ValueError("Missing provider contract: " + name)
    return {
        "stage": "same-final-profile-uniform-physical-primary-pulse-join",
        "schedule_correction": "Q148 equation (7) follows from the rule-defined scale choice after each finite coefficient seminorm exists; it is not an external unresolved inequality",
        "scale_domain": "every actual used dyadic Q=2^-ell with ell>=ell0 on the common derivative-independent later tail",
        "spatial_loss": "s_x=1/2+h/2",
        "time_loss": "s_t=1+3h/2",
        "velocity": "|nabla^a dt^b u_primary|<=C_ab Q^(-(a+1)s_x-b*s_t) S^B sqrt(zeta) delta^-B P",
        "pressure": "|nabla^a dt^b p_primary|<=C_ab Q^(-(a+2)s_x-b*s_t) S^B sqrt(zeta) delta^-B P",
        "potential": "|nabla^a dt^b A_primary|<=C_ab Q^(-a*s_x-b*s_t) S^B sqrt(zeta) delta^-B P",
        "supports": {
            "transverse": "abs(xi)<=r0/2",
            "pulse": "abs(v/Ls-1/2)<=2/5",
            "radial": "fixed closed active shell",
            "slow": "fixed eta cutoff support",
        },
        "full_pressure_memory": "Pi1=K1-mu^2/(4X), Pi2=K2, moving-normal pressure and cutoff derivatives retained",
        "uniform_primary_pulse_family_certified": True,
        "canonical_and_physical_primary_families_identified": False,
        "complete_nonlinear_generated_pulse_hierarchy_certified": False,
        "complete_covariance_residual_or_blowup_certified": False,
    }


def primary_powers(h, a, b, identity):
    if identity != IDENTITY:
        raise ValueError("Q149 identity required")
    return transfer.primary_powers(Q(h), a, b)


def covariance_obstruction(h, identity):
    """The known primary covariance remainder grows as Q tends to zero."""
    if identity != IDENTITY:
        raise ValueError("Q149 identity required")
    h = Q(h)
    if not 0 < h < Q(1, 100):
        raise ValueError("Inherited h range required")
    exponent = -1 - 2 * h * transfer.KAPPA
    return {
        "Q_exponent": exponent,
        "absolutely_decaying": exponent > 0,
        "necessary_next_step": "exact complete covariance/mean cancellation with axial velocity, all tensor entries, curl remainders and both pressure tails",
    }


def record():
    return {
        "status": "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined",
        "correction": "the exact-real rule-defined Borel prescription chooses each scale after all finite coefficient seminorms exist, so Q148's diagonal inequality follows by construction rather than remaining an external numerical-scale premise",
        "established": "on the same finalized leading/background profile, complete primary curl velocity, wave pressure and fixed potential have uniform value, fixed physical derivative and closed-support bounds on every actual used dyadic scale",
        "scope": "uniform independently audited physical primary family only; it is not identified with the separate canonical family, sourced inverses remain conditional on their input class and smooth endpoint support, and no complete nonlinear pulse hierarchy is claimed",
        "obstruction": "the known primary covariance remainder has Q exponent -1-2h*kappa<0 and cannot be absorbed by absolute decay",
        "next": "derive the complete wave covariance and mean corrections, retaining axial velocity, all tensor entries, curl remainders, Pi1/Pi2 and moving pressure before the nonlinear source and full residual ledger",
        "complete_covariance_residual_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
