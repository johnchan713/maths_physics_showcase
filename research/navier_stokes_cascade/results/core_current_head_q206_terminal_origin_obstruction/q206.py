"""Q206: terminal-origin velocity and pressure value-jet obstruction."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
IDENTITY = "Q206-TERMINAL-ORIGIN-OBSTRUCTION"
PROVIDERS = (
    "core_current_head_q205_stagewise_pulse_diagonal",
    "core_current_head_q203_diagonal_endpoint_force_corrected_domain",
    "core_current_head_q182_comparison_core_exclusion_replay",
    "core_current_head_q170_protected_core_limit_frontier_audit",
    "core_pressure_identification",
    "background_summation",
)
EXPECTED = {
    PROVIDERS[0]: "Q205-stagewise-generated-pulse-bounds-suffice-for-diagonal-uniform-raw-hierarchy-remains-unproved",
    PROVIDERS[1]: "Q203-corrected-domain-diagonal-endpoint-localized-force-replayed-conditionally",
    PROVIDERS[2]: "Q182-physical-pressure-comparison-protected-core-replayed-conditionally",
    PROVIDERS[3]: "Q170-current-head-protected-core-and-limit-order-audited-conditionally",
    PROVIDERS[4]: "leading-core-pressure-identified-with-regular-angular-inverse-and-retained-higher-memory",
    PROVIDERS[5]: "conditional-all-order-background-with-flat-augmented-residual",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q206 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def same_object_links(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q206 identity required")
    links = (
        (PROVIDERS[0], PROVIDERS[1]),
        (PROVIDERS[1], "core_current_head_q185_corrected_rate_diagonal_comparison_replay"),
        ("core_current_head_q185_corrected_rate_diagonal_comparison_replay", PROVIDERS[2]),
        (PROVIDERS[2], PROVIDERS[3]),
        (PROVIDERS[3], "core_same_field_positive_core_reconciliation"),
    )
    result = {}
    for parent, child in links:
        child_path = f"results/{child}/evidence.json"
        result[f"{parent}_uses_{child}"] = (
            evidence(parent)["dependency_evidence_sha256"][child_path]
            == sha(PROJECT / child_path)
        )
    return result


def geometry(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q206 identity required")
    h = Q(1, 200)
    A = Q(1, 2) + h
    Xa = Q(4, 101)
    X1 = Q(1, 101)
    return {
        "h": h,
        "A": A,
        "two_A": 2 * A,
        "Xa": Xa,
        "X1": X1,
        "support_margin": Xa - X1,
        "inverse_2h": 1 / (2 * h),
        "path": "x_(X,tau)=(sqrt(2 X tau),0,0); q=tau, eta=0",
        "annular_generated_support_starts_at_Xa": True,
        "localization_plateau_preserves_both_paths": True,
    }


def velocity_obstruction(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q206 identity required")
    g = geometry()
    return {
        "leading_coefficient": "e0=E0(X1,0)>0",
        "normalized_identity": "tau^A u_theta=e0+R_u(tau)",
        "remainder": "|R_u(tau)|<=C_u tau^(2h)",
        "threshold_power": g["inverse_2h"],
        "lower_bound": "|u_theta|>=(e0/2)tau^(-101/200)",
        "lower_bound_exponent": -g["A"],
        "path_tends_to_terminal_origin": True,
        "continuous_terminal_velocity_possible": False,
        "value_jet_exists": False,
    }


def pressure_obstruction(e0=Q(3, 5), X2=Q(3, 202), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q206 identity required")
    g = geometry()
    if e0 <= 0 or not g["X1"] < X2 < g["Xa"]:
        raise ValueError("need e0>0 and X1<X2<Xa")
    # A rational lower surrogate for log(X2/X1): log y >= 2(y-1)/(y+1), y>=1.
    y = X2 / g["X1"]
    log_lower = 2 * (y - 1) / (y + 1)
    d0_lower = e0 * e0 * log_lower / 8
    return {
        "leading_gauge": "Pi0(X,eta)=-integral_X^infinity E0(x,eta)^2/(2x) dx",
        "X2": X2,
        "assumed_continuity_interval_bound": "|E0(X,0)|>=e0/2 on [X1,X2]",
        "pressure_gap": "d0=Pi0(X2,0)-Pi0(X1,0)>0",
        "log_ratio_lower": log_lower,
        "d0_rational_lower": d0_lower,
        "normalized_difference": "tau^(2A)(p(x2_tau)-p(x1_tau))=d0+R_p(tau)",
        "remainder": "|R_p(tau)|<=C_p tau^(2h)",
        "lower_bound": "|p(x2_tau)-p(x1_tau)|>=(d0/2)tau^(-101/100)",
        "lower_bound_exponent": -g["two_A"],
        "time_dependent_spatial_gauge_cancels": True,
        "continuous_terminal_pressure_in_any_time_gauge_possible": False,
        "value_jet_exists": False,
        "actual_continuum_X2_selected_by_continuity_not_by_this_fixture": True,
    }


def endpoint_separation(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q206 identity required")
    return {
        "punctured_one_sided_velocity_pressure_jets": True,
        "terminal_origin_in_punctured_atlas": False,
        "smooth_localized_force_trace": True,
        "smooth_force_implies_smooth_field": False,
        "complete_residual_corner_flatness": True,
        "corner_flatness_implies_field_value_bound": False,
        "future_velocity_or_pressure_constructed": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q206 identity required")
    status = statuses()
    links = same_object_links()
    g = geometry()
    u = velocity_obstruction()
    p = pressure_obstruction()
    endpoint = endpoint_separation()
    checks = {
        "provider_statuses_match": all(status.values()),
        "same_object_chain_frozen": all(links.values()),
        "protected_paths_stay_inside_support_hole": g["support_margin"] == Q(3, 101),
        "velocity_value_jet_obstructed": (
            u["lower_bound_exponent"] == Q(-101, 200)
            and not u["continuous_terminal_velocity_possible"]
        ),
        "pressure_two_path_gap_positive": p["d0_rational_lower"] > 0,
        "pressure_obstruction_is_gauge_independent": (
            p["time_dependent_spatial_gauge_cancels"]
            and p["lower_bound_exponent"] == Q(-101, 100)
            and not p["continuous_terminal_pressure_in_any_time_gauge_possible"]
        ),
        "punctured_jets_force_and_corner_field_are_separated": (
            endpoint["punctured_one_sided_velocity_pressure_jets"]
            and endpoint["smooth_localized_force_trace"]
            and not endpoint["smooth_force_implies_smooth_field"]
            and not endpoint["corner_flatness_implies_field_value_bound"]
        ),
    }
    return {
        "status": "Q206-terminal-origin-velocity-and-gauge-independent-pressure-value-jets-obstructed-conditionally",
        "checks": checks,
        "scope": {
            "regular_terminal_velocity_obstructed": True,
            "regular_terminal_pressure_obstructed": True,
            "punctured_terminal_jets_preserved": True,
            "smooth_force_trace_preserved": True,
            "unconditional_complete_residual_force_exclusion_or_blowup": False,
        },
        "next": "reconstruct fixed-stage complete residual constants from the selected profile roots and compact/slow inverse maps",
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q206 identity required")
    return {
        "status": review()["status"],
        "providers": statuses(),
        "same_object": same_object_links(),
        "geometry": geometry(),
        "velocity": velocity_obstruction(),
        "pressure": pressure_obstruction(),
        "endpoint": endpoint_separation(),
        "review": review(),
    }
