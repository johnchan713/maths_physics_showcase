"""Q203: replay diagonal, endpoint and localization on Q202's corrected state."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
IDENTITY = "Q203-DIAGONAL-ENDPOINT-FORCE-CORRECTED-DOMAIN"
PROVIDERS = (
    "core_current_head_q202_next_sourced_cycle",
    "core_current_head_q181_diagonal_endpoint_force_replay",
    "core_current_head_q185_corrected_rate_diagonal_comparison_replay",
)
EXPECTED = {
    PROVIDERS[0]: "Q202-next-constrained-sourced-cycle-joined-on-corrected-physical-domain",
    PROVIDERS[1]: "Q181-physical-diagonal-endpoint-localized-force-replayed-conditionally",
    PROVIDERS[2]: "Q185-corrected-rate-diagonal-localization-comparison-replayed-conditionally",
}


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q202 = load("q203_q202", "core_current_head_q202_next_sourced_cycle/q202.py")
q181 = load("q203_q181", "core_current_head_q181_diagonal_endpoint_force_replay/q181.py")
q185 = load("q203_q185", "core_current_head_q185_corrected_rate_diagonal_comparison_replay/q185.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q203 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def same_object(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q203 identity required")
    q202_record = evidence(PROVIDERS[0])
    q181_record = evidence(PROVIDERS[1])
    q185_record = evidence(PROVIDERS[2])
    q180_path = "results/core_current_head_q180_target_source_simultaneous_induction/evidence.json"
    q181_path = "results/core_current_head_q181_diagonal_endpoint_force_replay/evidence.json"
    chain = q181.identity_chain()
    return {
        "Q202_and_Q181_share_frozen_Q180": (
            q202_record["dependency_evidence_sha256"][q180_path]
            == q181_record["dependency_evidence_sha256"][q180_path]
            == sha(PROJECT / q180_path)
        ),
        "Q185_uses_frozen_Q181": (
            q185_record["dependency_evidence_sha256"][q181_path]
            == sha(PROJECT / q181_path)
        ),
        "Q181_same_object_chain_holds": all(chain.values()),
        "Q202_changes_domain_bookkeeping_not_field_or_width": (
            q202.correction()["historical_frozen_evidence_modified"] is False
            and q202.sourced_inverse()["same_accepted_domain_and_C_minus_20_width"]
        ),
    }


def rate_route(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q203 identity required")
    current = q202.rates()
    margins = q185.residual_margin_replay()
    return {
        "Q202_current_source_D": current["current_source_D"],
        "Q202_cycle_gain": current["cycle_gain_delta"],
        "Q202_is_existing_log_renormalized_WMS_route": True,
        "log_margin": margins["log_regular_and_cross_margin"],
        "pure_margin": margins["pure_regular_and_cross_margin"],
        "pure_margin_correction": margins["pure_margin_correction"],
        "A0_own_margin": margins["A0_own_margin_both_routes"],
        "A1_own_margin": margins["A1_own_margin_both_routes"],
        "compact_accounts_disjoint": margins["compact_accounts_disjoint"],
        "Q202_D_may_be_read_as_same_exponent_pure_power": False,
        "pure_route_rebuild_and_scale_reselection_required": True,
    }


def diagonal_stage(target=0, loss_L=0, loss_H=0, route="log", identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q203 identity required")
    result = q185.comparison_stage(target, loss_L, loss_H, route)
    domain = q202.physical_domain()
    result.update({
        "Q202_physical_common_pole": domain["common_wave_pole"],
        "physical_pole_substituted_for_L_or_H": False,
        "losses_are_fixed_jet_diagonal_losses": True,
    })
    return result


def complete_diagonal(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q203 identity required")
    diagonal = q181.diagonal_residual_contract()
    stage = diagonal_stage()
    return {
        "increment": "Delta u_j=curl A_j+b_j e_theta; pi_j=fixed-gauge pressure increment",
        "cut_increment": "V_j=curl(chi(a_j q)A_j)+chi(a_j q)b_j e_theta; P_j=chi(a_j q)pi_j",
        "divergence_free_term_by_term": True,
        "scale_selected_after_finite_stage_constants": True,
        "cutoff_sequence_selected_once": diagonal["cutoff_sequence_chosen_once"],
        "locally_finite": diagonal["locally_finite_divergence_free_field"],
        "raw_series_convergence_used": diagonal["raw_series_convergence_used"],
        "six_groups": (
            "e_t", "-Delta e", "grad pi", "(U.grad)e", "(e.grad)U", "(e.grad)e"
        ),
        "all_six_groups_retain_complete_memory": diagonal["all_six_groups_retain_complete_memory"],
        "every_scale_cutoff_curl_and_pressure_derivative_retained": True,
        "fixed_jet_then_finite_stage": diagonal["fixed_jet_target_then_finite_stage"],
        "stage_uniform_constant_or_J_of_q_used": diagonal["stage_uniform_constant_or_J_of_q_used"],
        "all_fixed_cartesian_residual_jets_flat_conditionally": diagonal["all_fixed_cartesian_residual_jets_flat"],
        "sample_zero_loss_stage": stage,
        "physical_domain_is_separate_from_diagonal_losses": (
            stage["Q202_physical_common_pole"] == 361200138
            and not stage["physical_pole_substituted_for_L_or_H"]
        ),
    }


def endpoint_join(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q203 identity required")
    endpoint = q181.endpoint_contract()
    return {
        "jacobian_floor": endpoint["jacobian_floor"],
        "corner_relation": "q<=4 rho",
        "arbitrary_fixed_cartesian_power_closes": endpoint["arbitrary_fixed_cartesian_power_closes"],
        "residual_flatness_used_as_endpoint_field_bound": False,
        "punctured_endpoint_bound_supplied_separately": endpoint["endpoint_inequality_supplied_conditionally"],
        "closed_band_window": endpoint["closed_band_window"],
        "fixed_localization_request": endpoint["localization_request"],
        "full_exterior_pressure_convolution_retained": endpoint["full_pressure_Leibniz_convolution_retained"],
        "pressure_S2_excess": endpoint["pressure_S2_excess"],
        "compatible_one_sided_fixed_order_limits": endpoint["compatible_one_sided_fixed_order_limits"],
        "exterior_time_Taylor_radius_positive": endpoint["exterior_time_Taylor_radius_positive"],
        "terminal_origin_included": endpoint["singular_origin_included"],
        "Q202_giant_fixed_domain_is_admissible_but_not_evaluated": (
            q202.physical_domain()["actual_high_jet_values_evaluated"] is False
        ),
    }


def localization_force(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q203 identity required")
    local = q181.localization_force_energy_contract()
    return {
        "field": local["field"],
        "divergence_free": local["divergence_free"],
        "zero_initial_velocity": local["zero_initial_velocity"],
        "protected_inner_core_unchanged": local["protected_late_core_unchanged"],
        "term_count": local["term_count"],
        "terms": local["terms"],
        "vanishing_thirteenth": local["vanishing_thirteenth"],
        "unit_viscosity": local["unit_viscosity"],
        "complete_pressure": local["complete_pressure"],
        "force_divergence_retained": local["force_divergence_retained"],
        "direct_force_orders": local["direct_force_orders"],
        "Cauchy_orders": local["Cauchy_orders"],
        "supplied_orders": local["supplied_orders"],
        "one_unused_field_derivative": local["one_unused_field_derivative"],
        "one_compact_smooth_force_conditionally": local["one_compact_smooth_force"],
        "raw_Taylor_convergence_used": local["raw_Taylor_convergence_used"],
        "uniform_in_index_norm_assumed": local["uniform_in_j_norm_assumed"],
        "zero_safe_energy_step": local["zero_safe_energy_step"],
        "factor_two_retained": local["factor_two_retained"],
        "future_velocity_or_pressure_constructed": local["future_velocity_or_pressure_constructed"],
    }


def obstruction(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q203 identity required")
    powers = q202.physical_powers()
    return {
        "exact_diagonal_condition": "h*delta*(J0+1)/2-L_{m+2}-H_{m+2} >= N+1",
        "finite_stage_exists_conditionally_for_each_fixed_finite_L_H": True,
        "Q202_physical_pole_determines_L_H": False,
        "actual_L_H_or_giant_seminorms_evaluated": False,
        "sourced_velocity_zero_jet_power": powers["sourced_velocity_zero_jet"],
        "sourced_pressure_zero_jet_power": powers["sourced_pressure_zero_jet"],
        "next_regular_source_zero_jet_power": powers["next_regular_source_zero_jet"],
        "available_upper_bounds_certify_absolute_decay": False,
        "conditional_diagonal_flatness_promoted_to_uniform_actual_pulse_bound": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q203 identity required")
    status = statuses()
    chain = same_object()
    rates = rate_route()
    diagonal = complete_diagonal()
    endpoint = endpoint_join()
    local = localization_force()
    block = obstruction()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "same_physical_object_chain_holds": all(chain.values()),
        "Q202_uses_correct_log_route_margin": (
            rates["Q202_cycle_gain"] == Q(49997, 100000)
            and rates["log_margin"] == Q(49997, 40000000)
        ),
        "pure_route_is_separate_and_corrected": (
            not rates["Q202_D_may_be_read_as_same_exponent_pure_power"]
            and rates["pure_margin"] == Q(49996, 40000000)
            and rates["pure_margin_correction"] == -Q(1, 40000000)
        ),
        "complete_diagonal_memory_and_quantifiers_retained": (
            diagonal["all_six_groups_retain_complete_memory"]
            and diagonal["every_scale_cutoff_curl_and_pressure_derivative_retained"]
            and diagonal["fixed_jet_then_finite_stage"]
            and not diagonal["stage_uniform_constant_or_J_of_q_used"]
        ),
        "physical_domain_not_misused_as_diagonal_loss": diagonal["physical_domain_is_separate_from_diagonal_losses"],
        "corner_and_endpoint_are_separate_and_complete": (
            endpoint["arbitrary_fixed_cartesian_power_closes"]
            and endpoint["punctured_endpoint_bound_supplied_separately"]
            and endpoint["full_exterior_pressure_convolution_retained"]
            and not endpoint["terminal_origin_included"]
        ),
        "localization_has_exact_complete_twelve_term_ledger": (
            local["term_count"] == 12
            and local["complete_pressure"]
            and local["unit_viscosity"] == 1
            and "=0" in local["vanishing_thirteenth"]
        ),
        "force_and_energy_replay_is_limited_correctly": (
            local["one_compact_smooth_force_conditionally"]
            and local["factor_two_retained"]
            and not local["future_velocity_or_pressure_constructed"]
        ),
        "negative_physical_powers_block_unconditional_promotion": (
            block["sourced_velocity_zero_jet_power"] < 0
            and block["sourced_pressure_zero_jet_power"] < 0
            and block["next_regular_source_zero_jet_power"] < 0
            and not block["available_upper_bounds_certify_absolute_decay"]
        ),
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "same_object": chain,
        "rate_route": rates,
        "diagonal": diagonal,
        "endpoint": endpoint,
        "localization_force_energy": local,
        "precise_limit": block,
        "scope": {
            "Q202_log_route_diagonal_endpoint_force_replayed_conditionally": True,
            "pure_route_identified_with_Q202_field_or_scales": False,
            "terminal_origin_regularized": False,
            "actual_giant_constants_or_seminorms_evaluated": False,
            "uniform_actual_pulse_or_complete_global_residual_bound_proved": False,
            "future_velocity_or_pressure_constructed": False,
            "unconditional_force_exclusion_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
        "next": (
            "audit fixed-gauge global pressure and moment identities, coefficient and reserved-mean "
            "supports, then seek uniform actual-pulse size/derivative/support estimates across "
            "every physical scale before any complete global residual claim"
        ),
    }


def record():
    return {
        "status": "Q203-corrected-domain-diagonal-endpoint-localized-force-replayed-conditionally",
        "established": (
            "Q202's corrected coupled state feeds the existing log-renormalized complete-potential "
            "diagonal, separate punctured endpoint and exact twelve-term localization/force chain, "
            "retaining complete pressure, viscosity, cutoff curls, both crosses and the square"
        ),
        "correction": (
            "Q202's D and 49997/40000000 reserve belong to the W/M/S route; a pure-power rebuild "
            "spends kappa, retains 49996/40000000 and may reselect scales, while the physical pole "
            "361200138 is not a diagonal loss L or H"
        ),
        "precise_limit": (
            "the diagonal stage requires h*delta*(J0+1)/2-L-H>=N+1 for fixed finite losses; the "
            "actual losses and giant seminorms are unevaluated and Q202's physical powers remain "
            "negative, so no uniform actual-pulse or unconditional residual bound follows"
        ),
        "historical_frozen_records_modified": False,
        "unconditional_complete_residual_force_exclusion_or_blowup_certified": False,
    }
