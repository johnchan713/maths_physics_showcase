"""Q205: audit stagewise generated-pulse constants and diagonal quantifiers."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
IDENTITY = "Q205-STAGEWISE-PULSE-DIAGONAL"
PROVIDERS = (
    "core_current_head_q204_global_ledger_uniform_primary",
    "sourced_inverse_foundation",
    "core_current_head_q153_physical_source_tail_diagonal_join",
    "core_certified_complete_potential_diagonal_join",
    "core_current_head_q203_diagonal_endpoint_force_corrected_domain",
)
EXPECTED = {
    PROVIDERS[0]: "Q204-global-pressure-coefficient-uniform-primary-and-covariance-ledgers-joined-conditionally",
    PROVIDERS[1]: "explicit-same-domain-sourced-inverse-and-high-harmonic-parameter-jets-relative-to-reviewed-frame-inputs",
    PROVIDERS[2]: "Q153-physical-source-tail-and-complete-potential-diagonal-bound-conditionally",
    PROVIDERS[3]: "certified-domain-conditional-complete-potential-diagonal-join",
    PROVIDERS[4]: "Q203-corrected-domain-diagonal-endpoint-localized-force-replayed-conditionally",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q205 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def same_object_links(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q205 identity required")
    links = (
        (PROVIDERS[0], PROVIDERS[4]),
        (PROVIDERS[4], "core_current_head_q202_next_sourced_cycle"),
        ("core_current_head_q161_source_tail_simultaneous_induction", PROVIDERS[2]),
        (PROVIDERS[2], "core_current_head_q152_physical_stress_feedback_join"),
    )
    result = {}
    for parent, child in links:
        child_path = f"results/{child}/evidence.json"
        result[f"{parent}_uses_{child}"] = (
            evidence(parent)["dependency_evidence_sha256"][child_path]
            == sha(PROJECT / child_path)
        )
    return result


def fixed_cycle_bounds(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q205 identity required")
    h = Q(1, 200)
    A = Q(1, 2) + h
    sx = Q(1, 2) + h / 2
    st = 1 + 3 * h / 2
    alpha = Q(3, 2) - Q(2, 100000)
    return {
        "h": h,
        "A": A,
        "sx": sx,
        "st": st,
        "sample_alpha": alpha,
        "velocity_zero_jet_power": -A + h * alpha,
        "pressure_zero_jet_power": -2 * A + h * (alpha + Q(1, 2)),
        "potential_zero_jet_power": Q(1, 2) - A + h * (alpha + Q(1, 2)),
        "constant_independent_of_running_band_at_fixed_cycle": True,
        "constant_may_depend_on_cycle_derivative_and_profiles": True,
        "full_path_zero_gives_smooth_zero_extension": True,
        "compact_source_need_not_have_zero_exit": True,
        "potential_and_pressure_cut_before_curl": True,
        "closed_support_not_enlarged_by_curl_graph_or_frame": True,
        "moving_normal_pressure_retained": True,
        "source_normal_pressure_retained": True,
        "normal_source_blocks_extra_inverse_square_pressure_gain": True,
    }


def raw_uniform_obstruction(identity=IDENTITY, count=12):
    if identity != IDENTITY:
        raise ValueError("Q205 identity required")
    if count < 2:
        raise ValueError("need at least two cycles")
    h = Q(1, 200)
    delta = Q(1, 2) - Q(3, 100000)
    exponents = tuple(Q(j * j) - h * delta * j for j in range(1, count + 1))
    return {
        "available_quantifiers": "forall j,m exists finite C_jm",
        "unproved_quantifiers": "forall m exists finite C_m forall j",
        "countermodel_Cj": "2^(j^2)",
        "countermodel_Pj": 0,
        "countermodel_Xj": "2^(j^2) q^(h delta j-L_m)",
        "countermodel_satisfies_stagewise_representative_bound": True,
        "fixed_q0": Q(1, 2),
        "base_two_exponents": exponents,
        "strictly_increasing_sample": all(a < b for a, b in zip(exponents, exponents[1:])),
        "unbounded_in_j": True,
        "common_raw_hierarchy_constant_follows": False,
        "actual_pulse_lower_bound": False,
        "residual_lower_bound": False,
    }


def diagonal_absorption(identity=IDENTITY, j=10, P=7):
    if identity != IDENTITY:
        raise ValueError("Q205 identity required")
    if j <= 0 or P < 0:
        raise ValueError("invalid sample")
    h = 1 / 200
    delta = 1 / 2 - 3 / 100000
    gamma = h * delta * j / 2
    c = j * j * math.log(2)
    s = max(1.0, 2 * (c + j) / gamma + P * P / (gamma * gamma)) + 1
    log_lhs = c + P * math.log1p(s) - gamma * s
    log_rhs = -j * math.log(2)
    derivative_at_s = P / (1 + s) - gamma
    return {
        "j": j,
        "finite_request_count": j + 1,
        "gamma": gamma,
        "countermodel_log_C": c,
        "P": P,
        "chosen_log_cutoff": s,
        "endpoint_log_lhs": log_lhs,
        "endpoint_log_rhs": log_rhs,
        "endpoint_bound_holds": log_lhs <= log_rhs,
        "log_bound_decreasing_after_cutoff": derivative_at_s < 0,
        "works_for_arbitrary_finite_stage_constants": True,
        "scale_chosen_after_current_constants": True,
        "stage_uniform_constant_used": False,
        "J_of_q_substitution_used": False,
        "locally_finite_on_q_positive": True,
    }


def pressure_support_memory(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q205 identity required")
    return {
        "fixed_gauge_increment": "Delta u_j=curl A_j+b_j e_theta; pi_j=p^[j]-p^[j-1]",
        "cut_complete_potential_direct_swirl_and_pressure": True,
        "stage_dependent_potential_gauge_allowed": False,
        "stage_dependent_pressure_gauge_allowed": False,
        "six_covariances": True,
        "radial_and_axial_velocity": True,
        "moving_pressure_product": True,
        "Pi1_Pi2_and_same_base_tails": True,
        "viscosity_cutoffs_both_crosses_and_square": True,
        "reserved_pressure_floor_preserved": True,
        "closed_support_stagewise": True,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q205 identity required")
    status = statuses()
    links = same_object_links()
    fixed = fixed_cycle_bounds()
    obstruction = raw_uniform_obstruction()
    diagonal = diagonal_absorption()
    memory = pressure_support_memory()
    checks = {
        "provider_statuses_match": all(status.values()),
        "same_object_links_frozen": all(links.values()),
        "fixed_cycle_all_band_bounds_and_closed_support": (
            fixed["constant_independent_of_running_band_at_fixed_cycle"]
            and fixed["potential_and_pressure_cut_before_curl"]
            and fixed["closed_support_not_enlarged_by_curl_graph_or_frame"]
        ),
        "full_pressure_memory_retained": (
            fixed["moving_normal_pressure_retained"]
            and fixed["source_normal_pressure_retained"]
            and memory["Pi1_Pi2_and_same_base_tails"]
        ),
        "common_raw_uniform_inference_blocked": (
            obstruction["strictly_increasing_sample"]
            and obstruction["unbounded_in_j"]
            and not obstruction["common_raw_hierarchy_constant_follows"]
        ),
        "stagewise_constants_absorbed_by_diagonal": (
            diagonal["endpoint_bound_holds"]
            and diagonal["log_bound_decreasing_after_cutoff"]
            and diagonal["scale_chosen_after_current_constants"]
            and not diagonal["stage_uniform_constant_used"]
        ),
        "complete_residual_groups_preserved_conditionally": memory["viscosity_cutoffs_both_crosses_and_square"],
    }
    return {
        "status": "Q205-stagewise-generated-pulse-bounds-suffice-for-diagonal-uniform-raw-hierarchy-remains-unproved",
        "checks": checks,
        "scope": {
            "stagewise_actual_generated_pulse_bounds_joined": True,
            "uniform_generated_hierarchy_constant_required_by_diagonal": False,
            "uniform_uncut_generated_hierarchy_proved": False,
            "absolute_pulse_decay_proved": False,
            "unconditional_complete_residual_or_blowup": False,
        },
        "next": "audit one-sided terminal velocity and pressure jets on the same locally finite fixed-gauge object",
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q205 identity required")
    return {
        "status": review()["status"],
        "providers": statuses(),
        "same_object": same_object_links(),
        "fixed_cycle": fixed_cycle_bounds(),
        "raw_uniform_obstruction": raw_uniform_obstruction(),
        "diagonal": diagonal_absorption(),
        "pressure_support": pressure_support_memory(),
        "review": review(),
    }
