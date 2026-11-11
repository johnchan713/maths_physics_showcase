"""Bind the existing Q2--Q5 finite-map chain to Q164's audited core."""
from decimal import Decimal
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
IDENTITY = "Q165-current-head-finite-map-and-weighted-source-frontier-audit"
PROVIDERS = (
    "core_current_head_q164_core_fixed_point_pressure_audit",
    "core_current_head_q1_exact_map_join",
    "core_current_head_q2_regional_gap_join",
    "core_current_head_q3_full_frame_join",
    "core_current_head_q4_fixed_profile_join",
    "core_current_head_q5_weighted_inverse_source_join",
)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q2 = load("q165_q2", "core_current_head_q2_regional_gap_join/regional_gaps.py")
q3 = load("q165_q3", "core_current_head_q3_full_frame_join/full_frame.py")
q4 = load("q165_q4", "core_current_head_q4_fixed_profile_join/fixed_profiles.py")
q5 = load("q165_q5", "core_current_head_q5_weighted_inverse_source_join/weighted_inverse.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def regional_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q165 identity required")
    radius = q2.perturbation_radius(
        Q(1, 5), 3, Q(1, 7), 4, Q(1, 11), Q(1, 13), Q(1, 17), 2
    )
    edge = q2.edge_contract()
    cover = q2.regional_cover()
    floor = q2.corrected_outer_ramp_floor()
    return {
        "region_count": len(cover["regions"]),
        "exact_shared_endpoints": cover["exact_shared_endpoints"],
        "corrected_outer_ramp_floor": str(floor),
        "floor_exceeds_one_half": floor > Decimal("0.5"),
        "absolute_flat_edge_stress_infimum": edge["absolute_open_annulus_stress_infimum"],
        "positive_absolute_flat_edge_floor_claimed": edge["positive_absolute_stress_floor"],
        "weighted_flat_edge_bound": edge["valid_lower_bound"],
        "smooth_edge_directions": edge["normalized_directions_extend_smoothly"],
        "same_band_radius": radius["radius"],
        "same_band_rows": radius["rows"],
        "all_same_band_rows_positive": all(value > 0 for value in radius["rows"].values()),
        "uncertified_historical_existing_tilt_gate_used": False,
    }


def frame_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q165 identity required")
    budget = q3.frame_budget(
        Q(961, 103424), Q(1, 38500), Q(1, 101), Q(1, 11), 3,
        2, Q(5, 3), Q(7, 4), 13, 17,
    )
    normal_source = q3.projected_pressure(
        (0, 0, 1), (0, 0, 0), (0, 0, 2), (1, 0, 0), 5, 1
    )
    moving = q3.projected_pressure(
        (2, 1, 1), (1, -1, 0), (0, 0, 0), (1, -2, 0), 7, -2
    )
    compactness = q3.compactness_contract()
    return {
        "budget": budget,
        "positive_normal_frame_and_diffusion_denominators": (
            budget["normal_reciprocal"] > 0
            and budget["dstar"] > 0
            and budget["tangent_constant"] > 0
        ),
        "picard_factorial_ratio_at_most_one_half": budget["picard_next_ratio"] <= Q(1, 2),
        "pure_normal_source_pressure": normal_source,
        "moving_normal_pressure": moving,
        "unprojected_source_pressure_retained": normal_source["source_normal_over_i"] != 0,
        "moving_normal_pressure_retained": moving["moving_normal_over_i"] != 0,
        "old_numeric_Cphi_gate_required": compactness["old_numeric_Cphi_gate_required"],
        "unspecified_historical_first_band_certified": compactness["unspecified_historical_first_band_certified"],
        "potential_and_pressure_cutoff_commutators_retained": True,
    }


def projector_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q165 identity required")
    rows = [
        q4.projector_contract(1, 2, Q(5, 4), Q(7, 4), scale, Q(10, 7), exponent)
        for exponent in (1, 2) for scale in (Q(1, 4), Q(1), Q(4))
    ]
    memory = q4.moving_pressure_memory(Q(1, 5), Q(7, 3), Q(11, 4), Q(13, 6), Q(17, 9))
    return {
        "fixed_order_majorant": q4.fixed_jet_majorant(
            4, [1, 3, 7, 15, 31], [1, 2, 4, 8, 16]
        ),
        "projector_rows": rows,
        "all_six_moments_exact": all(row["moment"] == 1 for row in rows),
        "actual_support_strictly_inside_hull": all(row["strict_support"] for row in rows),
        "moving_pressure_memory": memory,
        "crho_Z_times_P_retained": memory["crho_Z_times_P"] != 0,
        "same_choice_recompute_after_projector_selection": True,
        "uniformity_over_arbitrary_profiles_or_all_orders_claimed": False,
    }


def weighted_inverse_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q165 identity required")
    requests = q5.request_graph(0, 399995, 200001)
    reserve = q5.log_reserve(Q(1, 200), Q(1, 100000), 37)
    equality = q5.pure_power_equality(Q(1, 100000))
    domain = q5.fixed_stage_domain(7, 13, 19, 37)
    return {
        "two_edge_kernel": q5.two_edge_constant(Q(4, 49), Q(401, 100), Q(1, 4900)),
        "request_graph": requests,
        "log_reserve": reserve,
        "pure_power_equality": equality,
        "fixed_stage_domain": domain,
        "retained_log_reserve_positive": reserve["retained"] > 0,
        "old_pure_power_margin_is_exactly_zero": (
            equality["equality"]
            and equality["source_margin"] == 0
            and equality["viscosity_margin"] == 0
        ),
        "physical_characteristic_shift_retained": True,
        "A1_axial_velocity_retained": True,
        "complete_moving_pressure_and_Pi1_Pi2_retained": True,
        "stage_or_derivative_uniform_constant_required": False,
    }


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q165 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "Q164-all-index-core-contraction-and-leading-pressure-audited",
        PROVIDERS[1]: "current-head-Q1-exact-core-roots-and-leading-pressure-bound-on-one-profile",
        PROVIDERS[2]: "current-head-Q2-regional-continuum-factored-edges-and-same-band-gap-bound",
        PROVIDERS[3]: "current-head-Q3-full-frame-whole-path-and-phase-domain-bound",
        PROVIDERS[4]: "current-head-Q4-fixed-profile-compensator-and-same-choice-projectors-bound",
        PROVIDERS[5]: "current-head-Q5-weighted-inverse-and-log-source-domain-bound",
    }
    for name, status in expected.items():
        if data[name].get("status") != status:
            raise ValueError("Q165 provider status mismatch: " + name)
    q164 = data[PROVIDERS[0]]["symbolic"]
    q1 = data[PROVIDERS[1]]["review"]
    q2_scope = data[PROVIDERS[2]]["review"]["scope"]
    q3_scope = data[PROVIDERS[3]]["review"]["scope"]
    q4_scope = data[PROVIDERS[4]]["review"]["scope"]
    q5_scope = data[PROVIDERS[5]]["review"]["scope"]
    required = (
        q164["U1_core_contraction_and_leading_pressure_closed"],
        q1["scope"]["Q1_provider_chain_bound_on_one_current_head_object"],
        q2_scope["Q2_regional_continuum_and_factored_edge_chain_bound_on_that_object"],
        not q2_scope["positive_absolute_stress_infimum_at_flat_edges_claimed"],
        q3_scope["Q3_full_frame_primitive_path_and_phase_domain_bound_on_same_later_tail"],
        not q3_scope["old_numeric_Cphi_sufficient_gate_used"],
        q4_scope["Q4_fixed_profile_compensator_projectors_bound_on_same_later_tail"],
        q4_scope["same_choice_projector_dependent_objects_recomputed"],
        q5_scope["Q5_weighted_inverse_and_log_source_contract_bound_on_same_later_tail"],
        q5_scope["finite_stage_log_source_induction_instantiated"],
        not q5_scope["raw_pure_power_supremum_promoted_as_exact_current_source_domain"],
        not q5_scope["six_group_complete_residual_replayed"],
    )
    if not all(required):
        raise ValueError("Q165 finite-map provider contract missing")
    return {
        "Q164_U2_frontier_was_stale_relative_to_live_Q2_Q5_chain": True,
        "Q2_regional_and_factored_edge_contract_bound": True,
        "Q3_full_frame_and_whole_path_contract_bound": True,
        "Q4_fixed_projector_and_compensator_contract_bound": True,
        "Q5_weighted_inverse_and_log_source_contract_bound": True,
        "Q6_complete_residual_replayed_here": False,
    }


def remaining_gate(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q165 identity required")
    return {
        "U1_core_contraction_and_leading_pressure": "closed-by-Q164",
        "U2_finite_regional_frame_projector_inverse_source_chain": "closed-by-current-head-Q2-Q5-providers-and-Q165-audit",
        "next": (
            "audit Q6's complete six-group local residual replay on the same Q164-Q165 "
            "object, retaining time, unit viscosity, full pressure, both nonlinear "
            "cross orientations, the new square, every potential/pressure cutoff "
            "derivative, repeated covariance means and fresh axial velocity"
        ),
        "Q6_current_head_residual_audited_here": False,
        "corner_endpoint_localization_force_energy_comparison_promoted_here": False,
        "dependency_hashes_or_passing_tests_suffice": False,
        "unconditional_theorem_gate_closed": False,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q165 identity required")
    return (
        *({"id": f"P{i}", "status": "closed-relative-to-named-providers"}
          for i in range(1, 7)),
        {"id": "P7", "status": "closed-as-conditional-consequence-by-Q163"},
        {"id": "U1", "status": "closed-by-Q164-core-pressure-audit"},
        {"id": "U2", "status": "closed-by-Q165-finite-map-audit"},
        {"id": "U3", "status": "next-Q6-complete-residual-audit"},
    )


def record():
    return {
        "status": "Q165-current-head-finite-map-and-weighted-source-chain-audited",
        "correction": (
            "Q164's U2 frontier was stale: current-head Q2--Q5 already bind the "
            "regional gaps, full frame/path domain, fixed same-choice projectors and "
            "log-renormalized weighted inverse/source chain; Q165 joins them directly "
            "to Q164's audited core"
        ),
        "established": (
            "the valid weighted flat-edge bound and positive same-band radius, finite "
            "full-frame constrained ODE with complete pressure, exact fixed-profile "
            "projector moments/supports, and a strict log-renormalized source reserve "
            "after rejecting the zero-margin pure-power route"
        ),
        "scope": (
            "same-object analytical join relative to the named continuum providers; "
            "exact fixtures and evidence hashes are regressions, not proofs of actual "
            "continuum infima or selected-field values"
        ),
        "next": remaining_gate()["next"],
        "U2_finite_map_and_weighted_source_chain_closed": True,
        "Q6_complete_residual_audited_here": False,
        "unconditional_pulse_residual_force_exclusion_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
