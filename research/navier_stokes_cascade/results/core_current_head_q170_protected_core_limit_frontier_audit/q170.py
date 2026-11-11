"""Audit the protected physical core path and terminal limit order at Q169."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import importlib.util
import json
import sys


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
IDENTITY = "Q170-current-head-protected-core-limit-frontier-audit"
PROVIDERS = (
    "core_current_head_q169_pressure_comparison_frontier_audit",
    "core_current_head_q157_physical_core_exclusion",
    "core_current_head_q163_comparison_core_exclusion_join",
    "core_same_field_positive_core_reconciliation",
    "core_current_head_q149_primary_pulse_join",
)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q157 = load(
    "q170_q157",
    "core_current_head_q157_physical_core_exclusion/q157.py",
)
q163 = load(
    "q170_q163",
    "core_current_head_q163_comparison_core_exclusion_join/q163.py",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q170 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "Q169-current-head-fixed-gauge-pressure-comparison-audited-conditionally",
        PROVIDERS[1]: "Q157-same-physical-core-path-and-competitor-exclusion-conditionally",
        PROVIDERS[2]: "Q163-same-object-comparison-and-core-exclusion-joined-conditionally",
        PROVIDERS[3]: "same-complete-diagonal-positive-core-bound-reconciled",
        PROVIDERS[4]: "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined",
    }
    for name, status in expected.items():
        if data[name].get("status") != status:
            raise ValueError("Q170 provider status mismatch: " + name)

    paths = {
        "q168": "results/core_current_head_q168_localization_force_frontier_audit/evidence.json",
        "q167": "results/core_current_head_q167_corner_endpoint_frontier_audit/evidence.json",
        "q166": "results/core_current_head_q166_complete_residual_frontier_audit/evidence.json",
        "q165": "results/core_current_head_q165_finite_map_frontier_audit/evidence.json",
        "q164": "results/core_current_head_q164_core_fixed_point_pressure_audit/evidence.json",
        "q163": "results/core_current_head_q163_comparison_core_exclusion_join/evidence.json",
        "q157": "results/core_current_head_q157_physical_core_exclusion/evidence.json",
        "q156": "results/core_current_head_q156_physical_pressure_comparison/evidence.json",
        "q155": "results/core_current_head_q155_physical_force_energy_join/evidence.json",
        "q153": "results/core_current_head_q153_physical_source_tail_diagonal_join/evidence.json",
        "q149": "results/core_current_head_q149_primary_pulse_join/evidence.json",
        "core": "results/core_same_field_positive_core_reconciliation/evidence.json",
    }
    q169, q157_record, q163_record, core, q149 = (data[name] for name in PROVIDERS)
    q168 = evidence("core_current_head_q168_localization_force_frontier_audit")
    q167 = evidence("core_current_head_q167_corner_endpoint_frontier_audit")
    q166 = evidence("core_current_head_q166_complete_residual_frontier_audit")
    q165 = evidence("core_current_head_q165_finite_map_frontier_audit")
    q164 = evidence("core_current_head_q164_core_fixed_point_pressure_audit")
    chain = {
        "Q169_to_Q168": q169["dependency_evidence_sha256"][paths["q168"]] == sha(PROJECT / paths["q168"]),
        "Q168_to_Q167": q168["dependency_evidence_sha256"][paths["q167"]] == sha(PROJECT / paths["q167"]),
        "Q167_to_Q166": q167["dependency_evidence_sha256"][paths["q166"]] == sha(PROJECT / paths["q166"]),
        "Q166_to_Q165": q166["dependency_evidence_sha256"][paths["q165"]] == sha(PROJECT / paths["q165"]),
        "Q165_to_Q164": q165["dependency_evidence_sha256"][paths["q164"]] == sha(PROJECT / paths["q164"]),
        "Q164_to_Q163": q164["dependency_evidence_sha256"][paths["q163"]] == sha(PROJECT / paths["q163"]),
        "Q163_to_Q157": q163_record["dependency_evidence_sha256"][paths["q157"]] == sha(PROJECT / paths["q157"]),
        "Q157_to_Q153": q157_record["dependency_evidence_sha256"][paths["q153"]] == sha(PROJECT / paths["q153"]),
        "Q157_and_Q168_share_Q155": (
            q157_record["dependency_evidence_sha256"][paths["q155"]]
            == q168["dependency_evidence_sha256"][paths["q155"]]
        ),
        "Q157_and_Q169_share_Q156": (
            q157_record["dependency_evidence_sha256"][paths["q156"]]
            == q169["dependency_evidence_sha256"][paths["q156"]]
        ),
        "Q157_uses_frozen_positive_core": (
            q157_record["dependency_evidence_sha256"][paths["core"]]
            == sha(PROJECT / paths["core"])
        ),
        "Q157_uses_same_Q149_background": (
            q157_record["dependency_evidence_sha256"][paths["q149"]]
            == sha(PROJECT / paths["q149"])
        ),
    }
    conclusions = core["review"]["analytic_conclusions"]
    return {
        "identity_chain": chain,
        "all_identity_links_hold": all(chain.values()),
        "same_finalized_background": q157_record["providers"][
            "same_finalized_background"
        ],
        "physical_diagonal_has_unchanged_annular_support_hole": q157_record[
            "providers"
        ]["physical_diagonal_has_unchanged_annular_support_hole"],
        "localization_plateau_bound_to_same_field": q157_record["providers"][
            "Q155_localization_plateau"
        ],
        "Q169_preterminal_comparison": q169["symbolic"][
            "U6_pressure_comparison_closed_conditionally"
        ],
        "positive_core_provider_conclusions": conclusions,
        "Q149_same_background_phrase_present": (
            "same finalized leading/background profile" in q149["symbolic"]["established"]
        ),
    }


def core_path_contract(tau=Q(1, 1000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q170 identity required")
    tau = Q(tau)
    row = q157.path_coordinates(tau)
    params = q157.fixed_parameters()
    binding = q157.same_physical_binding()
    return {
        "h": params["h"],
        "A": params["A"],
        "X_a": params["X_a"],
        "X_in": params["X_in"],
        "support_hole_margin": params["X_a"] - params["X_in"],
        "two_h_reciprocal": params["two_h_reciprocal"],
        "coordinates": row,
        "path": "x_tau=(sqrt(2 X_in tau),0,0), t=1-tau",
        "leading_coefficient": (
            "e0=E0(X_in,0)=sqrt(2 X_in)*varphi_star*Phi/C>0"
        ),
        "background_remainder": binding["background_remainder"],
        "annular_corrections": binding["support_hole"],
        "background_terms_vanish_exactly": False,
        "annular_diagonal_terms_vanish_exactly": True,
        "canonical_physical_primary_identification_used": binding[
            "canonical_physical_primary_identification_used"
        ],
        "pressure_gauge_changed": binding["pressure_gauge_changed"],
    }


def threshold_contract(e0=1, background_constant=2, q0=Q(1, 8),
                       tau0=Q(1, 16), radius_squared=Q(1, 4),
                       identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q170 identity required")
    e0, background_constant, q0 = map(Q, (e0, background_constant, q0))
    tau0, radius_squared = map(Q, (tau0, radius_squared))
    row = q157.combined_threshold(
        e0, background_constant, q0, tau0, radius_squared
    )
    params = q157.fixed_parameters()
    tau = row["tau_star"]
    ratio = min(Q(1), e0 / (2 * background_constant))
    normalized_error_bound = background_constant * ratio
    path_radius_squared = 2 * params["X_in"] * tau
    return {
        "q_core": row["q_core"],
        "plateau": row["plateau"],
        "tau_star": tau,
        "strictly_positive": tau > 0,
        "normalized_error_bound_on_interval": normalized_error_bound,
        "half_leading_value": e0 / 2,
        "normalized_error_within_half": normalized_error_bound <= e0 / 2,
        "time_plateau_bound": tau <= tau0 / 2,
        "radial_plateau_bound": tau <= radius_squared / (8 * params["X_in"]),
        "path_radius_squared": path_radius_squared,
        "path_inside_quarter_radius_squared": path_radius_squared <= radius_squared / 4,
        "lower_bound": (
            "|u_theta(x_tau,1-tau)|>=(e0/2) tau^(-101/200) "
            "for 0<tau<=tau_star"
        ),
    }


def limit_contract(tau=None, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q170 identity required")
    row = q163.limit_order(tau)
    contradiction = q163.contradiction_contract()
    return {
        "tau": row["tau"],
        "path_time": row["t"],
        "comparison_time": row["T_tau"],
        "strict_time_order": row["strict_time_order"],
        "first": row["first"],
        "second": row["second"],
        "third": row["third"],
        "uniform_comparison_constant_as_T_to_one_required": row[
            "uniform_comparison_constant_as_T_to_one_required"
        ],
        "path_tends_to_terminal_origin": True,
        "growth_exponent": q157.fixed_parameters()["A"],
        "conditional_conclusion": contradiction["conditional_conclusion"],
        "global_smooth_competitor_property": contradiction[
            "global_smooth_competitor_property"
        ],
        "future_velocity_or_pressure_constructed": contradiction[
            "future_velocity_or_pressure_constructed"
        ],
        "unconditional_blowup_claimed": contradiction["unconditional_blowup_claimed"],
    }


def remaining_gate(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q170 identity required")
    return {
        "U7_protected_core_and_limit_order": "closed-conditionally-by-Q170",
        "next": (
            "perform a consolidated provider-discharge audit of Q164-Q170, reconcile the "
            "older comparison-foundation unresolved flag with the later independent audits, "
            "and prove or precisely obstruct the first remaining continuum premise"
        ),
        "Q163_stale_next_core_audit_already_completed_by_Q164": True,
        "all_named_continuum_providers_independently_discharged_here": False,
        "unconditional_admissible_force_or_blowup_certified": False,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q170 identity required")
    return (
        {"id": "U1", "status": "closed-by-Q164-core-pressure-audit"},
        {"id": "U2", "status": "closed-by-Q165-finite-map-audit"},
        {"id": "U3", "status": "closed-conditionally-by-Q166-residual-audit"},
        {"id": "U4", "status": "closed-conditionally-by-Q167-corner-endpoint-audit"},
        {"id": "U5", "status": "closed-conditionally-by-Q168-localization-force-energy-audit"},
        {"id": "U6", "status": "closed-conditionally-by-Q169-pressure-comparison-audit"},
        {"id": "U7", "status": "closed-conditionally-by-Q170-protected-core-limit-audit"},
        {"id": "U8", "status": "next-consolidated-provider-discharge-audit"},
    )


def record():
    return {
        "status": "Q170-current-head-protected-core-and-limit-order-audited-conditionally",
        "correction": (
            "annular diagonal increments vanish exactly in the support hole, whereas the "
            "positive-order shared Borel background is only O(tau^(2h)); the exponent-100 "
            "threshold controls that remainder without identifying canonical and physical pulses"
        ),
        "established": (
            "the Q169 localized physical field has a protected path with A=101/200 and a "
            "strictly positive threshold on which |u_theta|>=(e0/2)tau^(-A); completing the "
            "radius limit separately at T_tau=1-tau/2 conditionally excludes a globally smooth "
            "same-force zero-datum fixed-interval finite-energy competitor"
        ),
        "scope": (
            "conditional analytical join of named continuum providers; no future field is "
            "constructed, dependency hashes and tests are not continuum proofs, and the result "
            "is not an unconditional admissible-force or Navier-Stokes blow-up theorem"
        ),
        "U7_protected_core_limit_closed_conditionally": True,
        "conditional_same_force_competitor_exclusion": True,
        "unconditional_admissible_force_or_blowup_certified": False,
        "future_velocity_or_pressure_constructed": False,
        "historical_frozen_records_modified": False,
        "next": remaining_gate()["next"],
    }
