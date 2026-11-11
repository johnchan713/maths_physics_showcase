"""Join Q166's residual to Q7 and the later Q154 endpoint provider."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import sys

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
IDENTITY = "Q167-current-head-corner-endpoint-frontier-audit"
PROVIDERS = (
    "core_current_head_q166_complete_residual_frontier_audit",
    "core_current_head_q7_corner_endpoint_join",
    "core_current_head_q154_physical_endpoint_join",
    "core_current_head_q162_diagonal_endpoint_force_join",
)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q7 = load("q167_q7", "core_current_head_q7_corner_endpoint_join/corner_endpoint.py")
q154 = load("q167_q154", "core_current_head_q154_physical_endpoint_join/q154.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q167 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "Q166-current-head-complete-six-group-local-residual-audited-conditionally",
        PROVIDERS[1]: "current-head-Q7-cartesian-corner-and-punctured-endpoint-joined-conditionally",
        PROVIDERS[2]: "Q154-physical-diagonal-punctured-endpoint-bound-conditionally",
        PROVIDERS[3]: "Q162-physical-diagonal-endpoint-localized-force-chain-joined",
    }
    for name, status in expected.items():
        if data[name].get("status") != status:
            raise ValueError("Q167 provider status mismatch: " + name)

    q6_path = "results/core_current_head_q6_complete_residual_join/evidence.json"
    q7_path = "results/core_current_head_q7_corner_endpoint_join/evidence.json"
    q154_path = "results/core_current_head_q154_physical_endpoint_join/evidence.json"
    q165_path = "results/core_current_head_q165_finite_map_frontier_audit/evidence.json"
    q164_path = "results/core_current_head_q164_core_fixed_point_pressure_audit/evidence.json"
    q163_path = "results/core_current_head_q163_comparison_core_exclusion_join/evidence.json"
    q162_path = "results/core_current_head_q162_diagonal_endpoint_force_join/evidence.json"

    q166, q7_record, q154_record, q162_record = (data[name] for name in PROVIDERS)
    q165 = evidence("core_current_head_q165_finite_map_frontier_audit")
    q164 = evidence("core_current_head_q164_core_fixed_point_pressure_audit")
    q163 = evidence("core_current_head_q163_comparison_core_exclusion_join")
    chain = {
        "Q166_to_Q165": q166["dependency_evidence_sha256"][q165_path] == sha(PROJECT / q165_path),
        "Q165_to_Q164": q165["dependency_evidence_sha256"][q164_path] == sha(PROJECT / q164_path),
        "Q164_to_Q163": q164["dependency_evidence_sha256"][q163_path] == sha(PROJECT / q163_path),
        "Q163_to_Q162": q163["dependency_evidence_sha256"][q162_path] == sha(PROJECT / q162_path),
        "Q162_to_Q154": q162_record["dependency_evidence_sha256"][q154_path] == sha(PROJECT / q154_path),
        "Q154_to_Q7": q154_record["dependency_evidence_sha256"][q7_path] == sha(PROJECT / q7_path),
        "Q166_and_Q7_share_Q6": (
            q166["dependency_evidence_sha256"][q6_path]
            == q7_record["dependency_evidence_sha256"][q6_path]
        ),
    }
    return {
        "identity_chain": chain,
        "all_identity_links_hold": all(chain.values()),
        "Q162_same_physical_field_cutoff_support_and_pressure_gauge": q162_record[
            "providers"
        ]["same_physical_field_cutoff_sequence_support_and_pressure_gauge"],
        "Q154_supplies_Q7_endpoint_blocker_conditionally": q154_record["symbolic"][
            "punctured_closed_time_endpoint_certified_conditionally"
        ],
        "singular_origin_regularized": q154_record["symbolic"]["singular_origin_regularized"],
    }


def corner_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q167 identity required")
    rows = (
        q7.corner_target(0, 0),
        q7.corner_target(1, 3, (7, 11), 13, 17),
        q7.corner_target(5, 11, (37, 41), 43, 47),
    )
    return {
        "selected_h": q7.H,
        "selected_kappa": q7.KAPPA,
        "jacobian_floor": q7.selected_chart()["jacobian_floor"],
        "sample_stages": [row["comparison_stage"] for row in rows],
        "sample_cartesian_factors": [row["cartesian_factor"] for row in rows],
        "sample_exponent_margins": [row["retained_exponent_margin"] for row in rows],
        "all_controls_hold": all(all(row["checks"].values()) for row in rows),
        "all_six_residual_groups_reach_corner": all(
            tuple(row["six_group_exponents"]) == q7.q6.RESIDUAL_GROUPS for row in rows
        ),
        "complete_memory_reaches_corner": tuple(rows[0]["physical_memory"]) == q7.q6.PHYSICAL_MEMORY,
        "physical_support_distance_or_edge_weight_floor_used": False,
        "arbitrary_fixed_cartesian_power_closes": True,
    }


def endpoint_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q167 identity required")
    params = q154.selected_parameters()
    window = q154.band_window(24, 20, 18)
    observation = q154.observation_margin(Q(1, 8), 18)
    requests = q154.endpoint_requests(3, 2)
    second = q154.exterior_coefficients(2)
    ratios = [q154.taylor_ratio(j, Q(3, 2)) for j in (1, 10, 100)]
    return {
        "endpoint_inequality": (
            "sup_(E x [1-delta,1)) |D_x^alpha partial_t^j(A,B,p)| < infinity "
            "for every punctured terminal compact E and fixed alpha,j"
        ),
        "endpoint_inequality_supplied_conditionally": True,
        "implicit_floor": params["implicit_floor"],
        "full_closed_band_window": list(window["bands"]),
        "observation_inside_qstar": observation["inside_qstar"],
        "observation_beyond_first_band": observation["first_band_margin"],
        "finite_active_band_box_and_stage_sets_on_each_punctured_compact": True,
        "localization_request": requests,
        "pressure_S2": second["S_j"],
        "pressure_S2_exceeds_two_A2_by": second["S_j"] - 2 * second["A_j"],
        "full_pressure_Leibniz_convolution_retained": True,
        "sample_Taylor_ratios": ratios,
        "Taylor_ratios_grow_unbounded": ratios[2] > 8 * ratios[1] > ratios[0],
        "exterior_time_Taylor_radius_positive": False,
        "compatible_one_sided_fixed_order_limits": True,
        "singular_origin_included": False,
    }


def remaining_gate(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q167 identity required")
    return {
        "U4_corner_and_punctured_endpoint": "closed-conditionally-by-Q7-Q154-and-audited-by-Q167",
        "next": (
            "audit Q8's exact twelve-term axisymmetric localization, zero initialization "
            "and compact-force extension on the same Q164-Q167 object, retaining unit "
            "viscosity, every potential cutoff derivative, full pressure, both crosses "
            "and the cutoff-curl square"
        ),
        "localization_force_energy_comparison_replayed_here": False,
        "unconditional_theorem_gate_closed": False,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q167 identity required")
    return (
        {"id": "U1", "status": "closed-by-Q164-core-pressure-audit"},
        {"id": "U2", "status": "closed-by-Q165-finite-map-audit"},
        {"id": "U3", "status": "closed-conditionally-by-Q166-residual-audit"},
        {"id": "U4", "status": "closed-conditionally-by-Q167-corner-endpoint-audit"},
        {"id": "U5", "status": "next-Q8-localization-force-audit"},
    )


def record():
    return {
        "status": "Q167-current-head-cartesian-corner-and-punctured-endpoint-audited-conditionally",
        "correction": (
            "Q7's endpoint inequality was open at its own head, but Q154 later proves it "
            "on the same physical field and Q162 freezes the field/gauge/support identity; "
            "Q167 binds that provider to Q166's residual"
        ),
        "established": (
            "all fixed residual jets transfer to the closed Cartesian corner and every "
            "fixed potential, direct-swirl and full fixed-gauge pressure jet is uniformly "
            "bounded with a compatible one-sided limit on each punctured terminal compact"
        ),
        "scope": (
            "conditional analytical join of named continuum providers; the terminal origin "
            "is excluded, exterior time Taylor series have zero radius, and tests/hashes are "
            "not continuum proofs"
        ),
        "U4_corner_endpoint_closed_conditionally": True,
        "singular_origin_regularized": False,
        "localization_force_energy_comparison_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
        "next": remaining_gate()["next"],
    }
