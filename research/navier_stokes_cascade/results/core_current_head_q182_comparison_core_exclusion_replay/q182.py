"""Replay pressure comparison and protected-core exclusion after Q181."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
IDENTITY = "Q182-COMPARISON-PROTECTED-CORE-REPLAY"
PROVIDERS = (
    "core_current_head_q181_diagonal_endpoint_force_replay",
    "core_current_head_q163_comparison_core_exclusion_join",
    "core_current_head_q169_pressure_comparison_frontier_audit",
    "core_current_head_q170_protected_core_limit_frontier_audit",
)
EXPECTED = {
    PROVIDERS[0]: "Q181-physical-diagonal-endpoint-localized-force-replayed-conditionally",
    PROVIDERS[1]: "Q163-same-object-comparison-and-core-exclusion-joined-conditionally",
    PROVIDERS[2]: "Q169-current-head-fixed-gauge-pressure-comparison-audited-conditionally",
    PROVIDERS[3]: "Q170-current-head-protected-core-and-limit-order-audited-conditionally",
}


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q163 = load("q182_q163", "core_current_head_q163_comparison_core_exclusion_join/q163.py")
q169 = load("q182_q169", "core_current_head_q169_pressure_comparison_frontier_audit/q169.py")
q170 = load("q182_q170", "core_current_head_q170_protected_core_limit_frontier_audit/q170.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q182 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def identity_chain(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q182 identity required")
    q181 = evidence(PROVIDERS[0])
    q163_record = evidence(PROVIDERS[1])
    q169_record = evidence(PROVIDERS[2])
    q170_record = evidence(PROVIDERS[3])
    q162_path = "results/core_current_head_q162_diagonal_endpoint_force_join/evidence.json"
    q168_path = "results/core_current_head_q168_localization_force_frontier_audit/evidence.json"
    q169_path = "results/core_current_head_q169_pressure_comparison_frontier_audit/evidence.json"
    return {
        "Q181_and_Q163_share_Q162": (
            q181["dependency_evidence_sha256"][q162_path]
            == q163_record["dependency_evidence_sha256"][q162_path]
            == sha(PROJECT / q162_path)
        ),
        "Q181_and_Q169_share_Q168": (
            q181["dependency_evidence_sha256"][q168_path]
            == q169_record["dependency_evidence_sha256"][q168_path]
            == sha(PROJECT / q168_path)
        ),
        "Q170_uses_frozen_Q169": (
            q170_record["dependency_evidence_sha256"][q169_path]
            == sha(PROJECT / q169_path)
        ),
        "same_field_force_pressure_gauge_and_localization": q169_record[
            "providers"
        ]["same_Q162_field_force_pressure_gauge_and_localization"],
    }


def pressure_comparison_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q182 identity required")
    pressure = q169.pressure_contract()
    commutator = q169.commutator_contract()
    radius = q169.radius_contract()
    return {
        "individual_pressures": pressure["individual_pressures"],
        "difference_equation": pressure["difference_equation"],
        "force_cancellations": pressure["force_cancellations"],
        "constructed_pressure_gauge": pressure["constructed_pressure_gauge"],
        "pressure_recovery": pressure["recovery"],
        "recovery_logic": pressure["recovery_logic"],
        "finite_energy_essential": pressure["finite_energy_essential"],
        "finite_energy_counterexample": pressure["finite_energy_counterexample"],
        "competitor_pressure_decay_assumed": pressure[
            "competitor_pressure_decay_assumed"
        ],
        "Poisson_alone_removes_harmonic_gradient": pressure[
            "Poisson_equation_alone_removes_harmonic_gradient"
        ],
        "commutator_exponent": commutator["exponent"],
        "commutator_radial_coefficient": commutator["radial_coefficient"],
        "commutator_norm_decay": commutator["norm_decay"],
        "cutoff_power": commutator["cutoff_power"],
        "seventh_cutoff_power_sufficient": commutator[
            "seventh_power_sufficient"
        ],
        "compact_mixed_flux_decay": radius["compact_mixed_flux_decay"],
        "young_decays": tuple(item["young_decay"] for item in radius["terms"]),
        "sobolev_offset_decays": tuple(
            item["sobolev_offset_decay"] for item in radius["terms"]
        ),
        "cutoff_laplacian_decay": radius["cutoff_laplacian_decay"],
        "worst_radius_decay": radius["worst_decay"],
        "localized_inequality": radius["localized_inequality"],
        "gronwall": radius["gronwall"],
        "conclusion": radius["conclusion"],
        "competitor_derivative_bound_assumed": radius[
            "competitor_derivative_bound_assumed"
        ],
        "competitor_dissipation_bound_assumed": radius[
            "competitor_dissipation_bound_assumed"
        ],
        "competitor_boundary_flux_bound_assumed": radius[
            "competitor_boundary_flux_bound_assumed"
        ],
        "terminal_uniform_constants_assumed": radius[
            "terminal_uniform_constants_assumed"
        ],
    }


def protected_core_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q182 identity required")
    path = q170.core_path_contract(Q(1, 1000))
    threshold = q170.threshold_contract()
    limit = q170.limit_contract()
    return {
        "h": path["h"],
        "A": path["A"],
        "X_a": path["X_a"],
        "X_in": path["X_in"],
        "support_hole_margin": path["support_hole_margin"],
        "path": path["path"],
        "leading_coefficient": path["leading_coefficient"],
        "background_remainder": path["background_remainder"],
        "background_terms_vanish_exactly": path[
            "background_terms_vanish_exactly"
        ],
        "annular_diagonal_terms_vanish_exactly": path[
            "annular_diagonal_terms_vanish_exactly"
        ],
        "canonical_physical_primary_identification_used": path[
            "canonical_physical_primary_identification_used"
        ],
        "pressure_gauge_changed": path["pressure_gauge_changed"],
        "q_core": threshold["q_core"],
        "tau_star": threshold["tau_star"],
        "threshold_strictly_positive": threshold["strictly_positive"],
        "normalized_error_within_half": threshold[
            "normalized_error_within_half"
        ],
        "time_plateau_bound": threshold["time_plateau_bound"],
        "radial_plateau_bound": threshold["radial_plateau_bound"],
        "lower_bound": threshold["lower_bound"],
        "strict_time_order": limit["strict_time_order"],
        "first_limit": limit["first"],
        "second_limit": limit["second"],
        "third_limit": limit["third"],
        "uniform_comparison_constant_as_T_to_one_required": limit[
            "uniform_comparison_constant_as_T_to_one_required"
        ],
        "conditional_conclusion": limit["conditional_conclusion"],
        "future_velocity_or_pressure_constructed": limit[
            "future_velocity_or_pressure_constructed"
        ],
        "unconditional_blowup_claimed": limit["unconditional_blowup_claimed"],
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q182 identity required")
    status = statuses()
    chain = identity_chain()
    comparison = pressure_comparison_contract()
    core = protected_core_contract()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "all_same_object_links_hold": all(chain.values()),
        "common_force_is_retained_then_cancelled_once": (
            "-div f" in comparison["individual_pressures"]
            and comparison["force_cancellations"] == 1
        ),
        "pressure_gradient_uses_momentum_and_finite_energy": (
            comparison["finite_energy_essential"]
            and not comparison["Poisson_alone_removes_harmonic_gradient"]
            and not comparison["competitor_pressure_decay_assumed"]
        ),
        "expanding_ball_remainder_is_R_minus_two": (
            comparison["worst_radius_decay"] == 2
            and "R^-2" in comparison["localized_inequality"]
        ),
        "eighth_cutoff_power_is_required_by_this_factorization": (
            comparison["cutoff_power"] == 8
            and not comparison["seventh_cutoff_power_sufficient"]
        ),
        "protected_path_separates_annular_vanishing_and_background_smallness": (
            core["annular_diagonal_terms_vanish_exactly"]
            and not core["background_terms_vanish_exactly"]
            and core["threshold_strictly_positive"]
            and core["normalized_error_within_half"]
        ),
        "limit_order_is_fixed_T_then_radius_then_tau": (
            core["strict_time_order"]
            and "R to infinity" in core["second_limit"]
            and "tau" in core["third_limit"]
            and not core["uniform_comparison_constant_as_T_to_one_required"]
        ),
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "identity_chain": chain,
        "pressure_comparison": comparison,
        "protected_core": core,
        "correction": {
            "individual_pressure_equations_retain_force_divergence": True,
            "common_force_cancels_exactly_once_in_difference": True,
            "Poisson_equation_alone_does_not_remove_harmonic_gradient": True,
            "annular_corrections_vanish_but_background_is_only_small": True,
            "no_comparison_constant_uniform_as_T_to_one_is_required": True,
            "historical_frozen_evidence_modified": False,
        },
        "scope": {
            "P7_same_force_competitor_exclusion_closed_conditionally": True,
            "globally_smooth_competitor_class_requires_fixed_interval_L2": True,
            "all_named_continuum_providers_independently_reproved": False,
            "future_velocity_or_pressure_constructed": False,
            "unconditional_admissible_force_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
        "next": (
            "perform a current-head provider-discharge audit of Q172-Q182, avoiding the stale "
            "Q171 frontier, and identify the first named continuum estimate still only routed "
            "rather than independently proved"
        ),
    }


def record():
    return {
        "status": "Q182-physical-pressure-comparison-protected-core-replayed-conditionally",
        "established": (
            "on Q181's exact localized field, fixed-energy pressure recovery and the expanding-ball "
            "estimate give preterminal uniqueness for every fixed T<1; the same protected physical "
            "core path then conditionally excludes a globally smooth same-force zero-datum "
            "fixed-interval finite-energy competitor"
        ),
        "correction": (
            "the force divergence is retained in both individual pressure equations and cancels "
            "once; momentum plus finite energy removes the harmonic pressure gradient; annular "
            "increments vanish on the core path while the shared Borel background is only small"
        ),
        "limiting_estimates": (
            "the expanding-ball remainder is O(R^-2), the core exponent is A=101/200, and limits "
            "are taken at fixed T_tau<1, then R to infinity, then tau to zero"
        ),
        "historical_frozen_records_modified": False,
        "unconditional_admissible_force_breakdown_or_blowup_certified": False,
    }
