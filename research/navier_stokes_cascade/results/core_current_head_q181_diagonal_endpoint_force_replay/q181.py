"""Replay the physical diagonal, endpoint, localization and force chain after Q180."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
IDENTITY = "Q181-PHYSICAL-DIAGONAL-ENDPOINT-FORCE-REPLAY"
PROVIDERS = (
    "core_current_head_q180_target_source_simultaneous_induction",
    "core_current_head_q162_diagonal_endpoint_force_join",
    "core_current_head_q166_complete_residual_frontier_audit",
    "core_current_head_q167_corner_endpoint_frontier_audit",
    "core_current_head_q168_localization_force_frontier_audit",
)
EXPECTED = {
    PROVIDERS[0]: "Q180-physical-target-source-tail-simultaneous-induction-reconciled",
    PROVIDERS[1]: "Q162-physical-diagonal-endpoint-localized-force-chain-joined",
    PROVIDERS[2]: "Q166-current-head-complete-six-group-local-residual-audited-conditionally",
    PROVIDERS[3]: "Q167-current-head-cartesian-corner-and-punctured-endpoint-audited-conditionally",
    PROVIDERS[4]: "Q168-current-head-localization-force-energy-audited-conditionally",
}


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q162 = load("q181_q162", "core_current_head_q162_diagonal_endpoint_force_join/q162.py")
q166 = load("q181_q166", "core_current_head_q166_complete_residual_frontier_audit/q166.py")
q167 = load("q181_q167", "core_current_head_q167_corner_endpoint_frontier_audit/q167.py")
q168 = load("q181_q168", "core_current_head_q168_localization_force_frontier_audit/q168.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q181 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def identity_chain(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q181 identity required")
    q180 = evidence(PROVIDERS[0])
    q162_record = evidence(PROVIDERS[1])
    q167_record = evidence(PROVIDERS[3])
    q168_record = evidence(PROVIDERS[4])
    q161_path = "results/core_current_head_q161_source_tail_simultaneous_induction/evidence.json"
    q162_path = "results/core_current_head_q162_diagonal_endpoint_force_join/evidence.json"
    q167_path = "results/core_current_head_q167_corner_endpoint_frontier_audit/evidence.json"
    return {
        "Q180_and_Q162_share_Q161": (
            q180["dependency_evidence_sha256"][q161_path]
            == q162_record["dependency_evidence_sha256"][q161_path]
            == sha(PROJECT / q161_path)
        ),
        "Q167_uses_frozen_Q162": (
            q167_record["dependency_evidence_sha256"][q162_path]
            == sha(PROJECT / q162_path)
        ),
        "Q168_uses_frozen_Q167": (
            q168_record["dependency_evidence_sha256"][q167_path]
            == sha(PROJECT / q167_path)
        ),
        "same_physical_field_cutoff_support_and_pressure_gauge": q162_record[
            "providers"
        ]["same_physical_field_cutoff_sequence_support_and_pressure_gauge"],
    }


def diagonal_residual_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q181 identity required")
    diagonal = q162.diagonal_contract()
    margin = q166.universal_margin()
    residual = q166.residual_contract()
    return {
        "stage_zero_primitive_ceiling": diagonal["physical_stage_zero_primitive_ceiling"],
        "all_fixed_graphs_finite_and_acyclic": diagonal["all_fixed_graphs_finite_and_acyclic"],
        "cutoff_sequence_chosen_once": diagonal["cutoff_sequence_chosen_once"],
        "locally_finite_divergence_free_field": diagonal["one_locally_finite_divergence_free_field"],
        "raw_series_convergence_used": diagonal["raw_series_convergence_used"],
        "complete_residual_groups": diagonal["complete_residual_groups"],
        "all_six_groups_retain_complete_memory": (
            residual["all_six_groups_retained_in_order"]
            and residual["complete_pressure_covariance_axial_and_cutoff_memory_retained"]
        ),
        "fixed_jet_target_then_finite_stage": residual[
            "fixed_jet_target_then_finite_stage_quantifier_retained"
        ],
        "stage_uniform_constant_or_J_of_q_used": residual[
            "stage_uniform_constant_or_J_of_q_used"
        ],
        "regular_and_cross_margin": margin["regular_retained_margin"],
        "A0_own_log_margin": margin["A0_own_log_retained_margin"],
        "A1_own_log_margin": margin["A1_own_log_retained_margin"],
        "distinct_log_reserves_not_added": (
            not margin["disjoint_reserves_added_to_strengthen_common_exponent"]
            and not margin["compact_reserves_may_pay_regular_or_cutoff_logs"]
        ),
        "all_fixed_cartesian_residual_jets_flat": diagonal["all_fixed_residual_jets_flat"],
    }


def endpoint_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q181 identity required")
    endpoint = q167.endpoint_contract()
    corner = q167.corner_contract()
    return {
        "jacobian_floor": corner["jacobian_floor"],
        "arbitrary_fixed_cartesian_power_closes": corner[
            "arbitrary_fixed_cartesian_power_closes"
        ],
        "physical_support_distance_or_edge_floor_used": corner[
            "physical_support_distance_or_edge_weight_floor_used"
        ],
        "endpoint_inequality": endpoint["endpoint_inequality"],
        "endpoint_inequality_supplied_conditionally": endpoint[
            "endpoint_inequality_supplied_conditionally"
        ],
        "closed_band_window": endpoint["full_closed_band_window"],
        "localization_request": endpoint["localization_request"],
        "full_pressure_Leibniz_convolution_retained": endpoint[
            "full_pressure_Leibniz_convolution_retained"
        ],
        "pressure_S2_excess": endpoint["pressure_S2_exceeds_two_A2_by"],
        "compatible_one_sided_fixed_order_limits": endpoint[
            "compatible_one_sided_fixed_order_limits"
        ],
        "exterior_time_Taylor_radius_positive": endpoint[
            "exterior_time_Taylor_radius_positive"
        ],
        "singular_origin_included": endpoint["singular_origin_included"],
    }


def localization_force_energy_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q181 identity required")
    localization = q168.localization_contract()
    force = q168.force_contract()
    energy = q168.energy_contract()
    return {
        "field": localization["field"],
        "divergence_free": localization["divergence_free"],
        "zero_initial_velocity": localization["zero_initial_velocity"],
        "protected_late_core_unchanged": localization["protected_late_core_unchanged"],
        "term_count": localization["term_count"],
        "terms": localization["terms"],
        "vanishing_thirteenth": localization["vanishing_thirteenth"],
        "unit_viscosity": localization["unit_viscosity"],
        "complete_pressure": localization["complete_pressure"],
        "force_divergence_retained": localization["force_divergence_retained"],
        "direct_force_orders": localization["direct_force_orders"],
        "Cauchy_orders": localization["Cauchy_orders"],
        "supplied_orders": localization["Q167_supplied_orders"],
        "one_unused_field_derivative": localization["unused_field_derivative"],
        "one_compact_smooth_force": force["all_samples_have_geometric_tails"],
        "raw_Taylor_convergence_used": force["raw_Taylor_convergence_used"],
        "uniform_in_j_norm_assumed": force["uniform_in_j_norm_assumed"],
        "energy_identity": energy["identity"],
        "zero_safe_energy_step": energy["zero_safe_step"],
        "factor_two_retained": energy["factor_two_retained"],
        "combined_endpoint_extrema_budget_claimed": energy[
            "combined_endpoint_extrema_budget_claimed"
        ],
        "future_velocity_or_pressure_constructed": force[
            "future_velocity_or_pressure_constructed"
        ],
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q181 identity required")
    status = statuses()
    chain = identity_chain()
    diagonal = diagonal_residual_contract()
    endpoint = endpoint_contract()
    localization = localization_force_energy_contract()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "all_same_object_links_hold": all(chain.values()),
        "complete_six_group_residual_is_flat_for_every_fixed_jet": (
            diagonal["all_six_groups_retain_complete_memory"]
            and diagonal["all_fixed_cartesian_residual_jets_flat"]
        ),
        "residual_margin_is_strict": diagonal["regular_and_cross_margin"] == Q(49997, 40000000),
        "compact_log_margins_are_separate_and_strict": (
            diagonal["A0_own_log_margin"] == Q(1, 40000000)
            and diagonal["A1_own_log_margin"] == Q(1, 40000000)
            and diagonal["distinct_log_reserves_not_added"]
        ),
        "endpoint_is_separate_from_residual_flatness": (
            endpoint["endpoint_inequality_supplied_conditionally"]
            and endpoint["compatible_one_sided_fixed_order_limits"]
            and not endpoint["singular_origin_included"]
        ),
        "localization_has_exact_twelve_groups": (
            localization["term_count"] == 12
            and "=0" in localization["vanishing_thirteenth"]
        ),
        "force_extension_and_energy_are_zero_safe": (
            localization["one_compact_smooth_force"]
            and localization["factor_two_retained"]
            and not localization["combined_endpoint_extrema_budget_claimed"]
        ),
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "identity_chain": chain,
        "diagonal_residual": diagonal,
        "endpoint": endpoint,
        "localization_force_energy": localization,
        "correction": {
            "Q162_P6_is_now_bound_to_Q180_coupled_physical_induction": True,
            "residual_flatness_does_not_supply_endpoint_field_bounds": True,
            "zero_energy_division_is_replaced_by_epsilon_regularization": True,
            "smooth_force_gluing_does_not_construct_future_velocity_or_pressure": True,
            "historical_frozen_evidence_modified": False,
        },
        "scope": {
            "P6_closed_relative_to_named_continuum_providers": True,
            "terminal_origin_regularized": False,
            "global_residual_improvement_proved": False,
            "future_velocity_or_pressure_constructed": False,
            "comparison_core_exclusion_or_blowup_replayed": False,
            "tests_are_continuum_proof": False,
        },
        "next": (
            "replay fixed-gauge pressure recovery, expanding-ball comparison and the protected-core "
            "lower bound on this exact localized field, taking fixed T<1, then R to infinity, "
            "then tau to zero"
        ),
    }


def record():
    return {
        "status": "Q181-physical-diagonal-endpoint-localized-force-replayed-conditionally",
        "established": (
            "Q180's coupled physical induction feeds one complete-potential diagonal with all six "
            "residual groups flat in every fixed Cartesian jet, separate punctured endpoint bounds, "
            "exact twelve-term divergence-free localization, one compact smooth force and the "
            "zero-safe preterminal energy estimate"
        ),
        "correction": (
            "local residual flatness and endpoint field bounds are separate obligations; the latter "
            "comes from Q154/Q167, while the zero-datum energy step uses epsilon regularization "
            "rather than division by sqrt(E)"
        ),
        "strict_margins": (
            "regular/cross reserve h*delta/2=49997/40000000; independent A0 and A1 logarithmic "
            "reserves h*kappa/2=1/40000000"
        ),
        "historical_frozen_records_modified": False,
        "unconditional_complete_residual_force_exclusion_or_blowup_certified": False,
    }
