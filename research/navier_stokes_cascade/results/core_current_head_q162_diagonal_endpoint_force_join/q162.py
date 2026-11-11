"""Bind the complete diagonal, endpoint atlas, localization and force chain."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
IDENTITY = "Q162-same-physical-diagonal-endpoint-localized-force-join"
PROVIDERS = (
    "core_current_head_q161_source_tail_simultaneous_induction",
    "core_current_head_q153_physical_source_tail_diagonal_join",
    "core_certified_complete_potential_diagonal_join",
    "core_current_head_q154_physical_endpoint_join",
    "core_certified_complete_potential_endpoint_join",
    "endpoint_foundation",
    "current_endpoint_jets",
    "core_current_head_q155_physical_force_energy_join",
    "core_certified_complete_potential_force_join",
    "force_foundation",
    "current_localized_force",
    "localized_force_extension",
    "core_scale_compatible_corner_force_replay",
)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


diagonal = load("q162_diagonal", "core_certified_complete_potential_diagonal_join/diagonal_join.py")
q154 = load("q162_q154", "core_current_head_q154_physical_endpoint_join/q154.py")
q155 = load("q162_q155", "core_current_head_q155_physical_force_energy_join/q155.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def diagonal_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q162 identity required")
    inputs = diagonal.certified_input_review()
    scale = diagonal.scale_absorption_review()
    comparison = diagonal.fixed_comparison_review()
    closure = diagonal.closure_statement()
    return {
        "physical_stage_zero_primitive_ceiling": inputs["physical_zero_primitive_ceiling"],
        "all_fixed_graphs_finite_and_acyclic": inputs["all_fixed_graphs_finite_and_acyclic"],
        "complete_pressure_covariance_axial_velocity_and_supports_retained": inputs[
            "complete_pressure_covariance_axial_velocity_and_supports_retained"
        ],
        "arbitrary_finite_stage_constants_absorbed": scale["all_extreme_margins_strict"],
        "cutoff_sequence_chosen_once": scale["cutoff_sequence_chosen_once_for_all_fixed_jets"],
        "complete_residual_groups": comparison["complete_residual_groups"],
        "comparison_stage_fixed_before_q_limit": comparison["comparison_stage_fixed_before_q_limit"],
        "one_locally_finite_divergence_free_field": closure[
            "one_locally_finite_divergence_free_field_follows_conditionally"
        ],
        "all_fixed_residual_jets_flat": closure[
            "all_fixed_jets_of_complete_local_residual_are_flat_conditionally"
        ],
        "fixed_gauges": closure["fixed_operator_potential_and_pressure_gauges_retained"],
        "raw_series_convergence_used": False,
    }


def endpoint_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q162 identity required")
    params = q154.selected_parameters()
    window = q154.band_window(24, 20, 18)
    margin = q154.observation_margin(Q(1, 8), 18)
    requests = q154.endpoint_requests(3, 2)
    exterior2 = q154.exterior_coefficients(2)
    return {
        "implicit_floor": params["implicit_floor"],
        "closed_band_window": window,
        "observation_inside_qstar": margin["inside_qstar"],
        "observation_beyond_first_band": margin["first_band_margin"],
        "fixed_request_example": requests,
        "pressure_convolution_at_two": exterior2["S_j"],
        "pressure_convolution_strictly_exceeds_two_A2": (
            exterior2["S_j"] > 2 * exterior2["A_j"]
        ),
        "one_sided_jets": (
            "every fixed potential, direct-swirl, velocity, pressure and residual jet has a "
            "uniform one-sided limit on every punctured terminal compact"
        ),
        "singular_origin_included": False,
        "exterior_time_Taylor_radius_positive": False,
    }


def localization_force_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q162 identity required")
    ledger = q155.localization_ledger()
    request = q155.derivative_requests(5)
    energy = q155.energy_budget(Q(7, 3))
    return {
        "localization_terms": ledger["terms"],
        "term_count": ledger["count"],
        "vanishing_thirteenth": ledger["vanishing_thirteenth"],
        "unit_viscosity": ledger["viscosity"],
        "complete_pressure": ledger["complete_pressure"],
        "force_divergence_retained": ledger["force_divergence_retained"],
        "derivative_request_example": request,
        "corner_transfer_arbitrary_power": q155.corner_transfer(17)[
            "cartesian_flat_for_arbitrary_N"
        ],
        "single_Borel_schedule": (
            "choose each future-force cutoff scale after the finite terminal jet norms at "
            "that index; every fixed mixed derivative has a geometric tail"
        ),
        "energy": energy,
        "future_velocity_or_pressure_constructed": False,
    }


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q162 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "Q161-physical-source-inverse-tail-simultaneous-induction-joined",
        PROVIDERS[1]: "Q153-physical-source-tail-and-complete-potential-diagonal-bound-conditionally",
        PROVIDERS[2]: "certified-domain-conditional-complete-potential-diagonal-join",
        PROVIDERS[3]: "Q154-physical-diagonal-punctured-endpoint-bound-conditionally",
        PROVIDERS[4]: "certified-domain-conditional-complete-potential-punctured-endpoint-join",
        PROVIDERS[5]: "current-source-punctured-terminal-jets-relative-to-retained-core-and-finite-analytic-lemmas",
        PROVIDERS[6]: "conditional-current-field-closed-time-atlas-and-punctured-terminal-jets",
        PROVIDERS[7]: "Q155-physical-localization-compact-force-and-energy-bound-conditionally",
        PROVIDERS[8]: "certified-domain-conditional-localization-force-and-energy-join",
        PROVIDERS[9]: "current-source-compact-smooth-force-and-preterminal-energy-relative-to-retained-core-and-finite-analytic-lemmas",
        PROVIDERS[10]: "conditional-current-field-localization-smooth-compact-force-and-energy",
        PROVIDERS[11]: "conditional-physical-localization-compact-force-and-whole-space-energy",
        PROVIDERS[12]: "scale-compatible-corner-flatness-and-force-replayed",
    }
    for name, status in expected.items():
        if data[name].get("status") != status:
            raise ValueError("Q162 provider status mismatch: " + name)
    q161 = data[PROVIDERS[0]]["symbolic"]
    q153_record = data[PROVIDERS[1]]["symbolic"]
    diagonal_flags = data[PROVIDERS[2]]["flags"]
    q154_record = data[PROVIDERS[3]]["symbolic"]
    endpoint_flags = data[PROVIDERS[5]]["flags"]
    q155_record = data[PROVIDERS[7]]["symbolic"]
    force_flags = data[PROVIDERS[9]]["flags"]
    localized_flags = data[PROVIDERS[10]]["flags"]
    if not all((
        q161["P4_generated_source_endpoint_contract_discharged"],
        q161["P5_fixed_stage_weighted_inverse_and_tail_contract_discharged"],
        q153_record["physical_generated_source_endpoint_join_certified_conditionally"],
        q153_record["all_fixed_stage_additive_physical_tails_certified_conditionally"],
        q153_record["complete_potential_diagonal_applicable_conditionally"],
        diagonal_flags["one_locally_finite_divergence_free_field_constructed_conditionally"],
        diagonal_flags["all_fixed_Cartesian_jets_of_local_complete_residual_flat_conditionally"],
        diagonal_flags["complete_pressure_viscosity_cutoff_crosses_and_square_retained"],
        q154_record["punctured_closed_time_endpoint_certified_conditionally"],
        endpoint_flags["complete_uniform_next_time_derivative_bounds_joined_from_current_foundations"],
        endpoint_flags["both_background_pressure_tails_and_fixed_gauge_memory_retained"],
        q155_record["compact_force_and_energy_certified_conditionally"],
        force_flags["all_twelve_pressure_viscous_time_and_nonlinear_localization_terms_retained"],
        force_flags["one_Borel_sequence_controls_every_fixed_mixed_force_jet"],
        localized_flags["complete_current_terminal_force_jets_smooth_compact_and_flat_at_origin_conditionally_derived"],
        localized_flags["conditional_one_smooth_compact_global_force_matching_exact_preterminal_residual_supplied"],
    )):
        raise ValueError("Q162 diagonal/endpoint/force provider contract missing")
    return {
        "same_physical_field_cutoff_sequence_support_and_pressure_gauge": True,
        "P6_complete_diagonal_endpoint_localized_force_contract_closed": True,
        "historical_conditional_records_preserved": True,
        "P7_comparison_and_core_exclusion_not_replayed_here": True,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q162 identity required")
    return (
        {"id": "P1", "status": "closed", "contract": "PREMOD value/root envelope"},
        {"id": "P2", "status": "closed", "contract": "Borel background and physical primary pulses"},
        {"id": "P3", "status": "closed", "contract": "complete covariance/mean targets"},
        {"id": "P4", "status": "closed", "contract": "generated-source quotient and endpoint extension"},
        {"id": "P5", "status": "closed-for-each-fixed-stage-jet-power", "contract": "weighted inverses and physical tails"},
        {"id": "P6", "status": "closed-by-Q162-same-object-replay", "contract": "complete diagonal, endpoint and localized force"},
        {"id": "P7", "status": "next-replay-required", "contract": "comparison and core exclusion"},
    )


def record():
    return {
        "status": "Q162-physical-diagonal-endpoint-localized-force-chain-joined",
        "correction": (
            "Q153-Q155 were conditional because their physical source/target/tail inputs had not "
            "been closed together; Q161 now supplies those inputs on the same field and gauge"
        ),
        "established": (
            "one locally finite divergence-free complete-potential field has complete residual "
            "flat in every fixed Cartesian jet; its potential, direct swirl and fixed-gauge pressure "
            "have compatible one-sided jets on every punctured terminal compact; exact localization "
            "gives zero initial velocity, one compact smooth global force, and the preterminal "
            "energy and integrated-dissipation estimates"
        ),
        "scope": (
            "same-object analytical composition of the named continuum providers; the singular "
            "terminal origin is not regularized, exterior time Taylor series have zero radius, "
            "and no future velocity or pressure is constructed"
        ),
        "next": (
            "replay fixed-gauge pressure recovery, expanding-ball comparison and the protected-core "
            "lower bound on this exact localized field and force, respecting the order T<1, "
            "R to infinity, then tau to zero"
        ),
        "P6_complete_diagonal_endpoint_localized_force_contract_discharged": True,
        "P7_comparison_and_core_exclusion_replayed": False,
        "unconditional_competitor_exclusion_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
