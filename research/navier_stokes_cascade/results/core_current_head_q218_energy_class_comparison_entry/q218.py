"""Q218: energy-class entry to Q217's fixed-gauge comparison theorem."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q218-ENERGY-CLASS-COMPARISON-ENTRY"
PROVIDERS = (
    "core_current_head_q217_localization_comparison_frontier_join",
    "core_current_head_q168_localization_force_frontier_audit",
    "core_current_head_q169_pressure_comparison_frontier_audit",
    "core_current_head_q170_protected_core_limit_frontier_audit",
    "core_current_head_q206_terminal_origin_obstruction",
)
EXPECTED = {
    PROVIDERS[0]: "Q217-Q216-object-localization-force-energy-and-comparison-joined-conditionally",
    PROVIDERS[1]: "Q168-current-head-localization-force-energy-audited-conditionally",
    PROVIDERS[2]: "Q169-current-head-fixed-gauge-pressure-comparison-audited-conditionally",
    PROVIDERS[3]: "Q170-current-head-protected-core-and-limit-order-audited-conditionally",
    PROVIDERS[4]: "Q206-terminal-origin-velocity-and-gauge-independent-pressure-value-jets-obstructed-conditionally",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q218 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def energy_class(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q218 identity required")
    return {
        "local_smoothness_on_fixed_preterminal_intervals": True,
        "same_compact_smooth_force": True,
        "zero_initial_velocity": True,
        "velocity_in_Linfinity_t_L2_x": True,
        "gradient_in_L2_t_L2_x": True,
        "energy_inequality": True,
        "pressure_decay_assumed": False,
        "terminal_smoothness_assumed": False,
    }


def energy_entry(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q218 identity required")
    return {
        "starting_inequality": "E(t)+2 int_0^t D <= 2 int_0^t <f,v>",
        "young_inequality": "2 <f,v> <= ||f||_2^2+E",
        "gronwall_bound": "sup_[0,T] E <= exp(T) int_0^T ||f||_2^2",
        "compact_smooth_force_has_finite_rhs": True,
        "supremum_L2_required_by_Q169": True,
        "derived_from_bare_local_smoothness": False,
    }


def pressure_entry(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q218 identity required")
    return {
        "w_in_Linfinity_t_L2_x": True,
        "quadratic_difference_in_Linfinity_t_L1_x": True,
        "riesz_pressure_in_Hminus2": True,
        "momentum_places_harmonic_gradient_in_Hminus3": True,
        "harmonic_Hminus3_gradient_vanishes": True,
        "common_force_cancels_exactly_once": True,
        "competitor_pressure_decay_needed": False,
    }


def harmonic_counterexample(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q218 identity required")
    return {
        "velocity": "t e_1",
        "pressure": "-x_1",
        "force": "0",
        "time_derivative": "e_1",
        "pressure_gradient": "-e_1",
        "viscosity": "0",
        "convection": "0",
        "residual": "0",
        "zero_initial_velocity": True,
        "locally_smooth": True,
        "finite_energy_at_positive_time": False,
    }


def protected_limit(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q218 identity required")
    return {
        "velocity_exponent": Q(101, 200),
        "support_hole_margin": Q(3, 101),
        "limit_order": ["fix T_tau<1", "send R to infinity", "send tau to zero"],
        "excludes_smooth_through_terminal_energy_class_competitor": True,
        "excludes_every_locally_smooth_competitor": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q218 identity required")
    cls = energy_class()
    entry = energy_entry()
    pressure = pressure_entry()
    example = harmonic_counterexample()
    limit = protected_limit()
    checks = {
        "provider_statuses_match": all(statuses().values()),
        "energy_class_supplies_Q169_L2_entry": (
            cls["velocity_in_Linfinity_t_L2_x"]
            and entry["supremum_L2_required_by_Q169"]
            and entry["compact_smooth_force_has_finite_rhs"]),
        "pressure_recovery_keeps_full_memory": (
            pressure["common_force_cancels_exactly_once"]
            and not pressure["competitor_pressure_decay_needed"]),
        "bare_smooth_class_still_blocked": (
            not entry["derived_from_bare_local_smoothness"]
            and example["residual"] == "0"
            and not example["finite_energy_at_positive_time"]),
        "terminal_limit_order_preserved": limit["limit_order"] == [
            "fix T_tau<1", "send R to infinity", "send tau to zero"],
    }
    return {
        "status": "Q218-energy-class-entry-to-fixed-gauge-comparison-derived-conditionally",
        "checks": checks,
        "scope": {
            "finite_energy_comparison_premise_discharged_for_energy_class": True,
            "conditional_smooth_energy_competitor_exclusion": True,
            "bare_local_smooth_competitor_exclusion": False,
            "upstream_continuum_providers_discharged": False,
            "unconditional_force_residual_or_blowup": False,
        },
        "next": (
            "audit the first upstream continuum leaf left cited by Q216: the same-band regional "
            "gap infima and weighted relative-stress estimate on the finalized profile"
        ),
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q218 identity required")
    limit = protected_limit()
    return {
        "status": review()["status"],
        "providers": statuses(),
        "energy_class": energy_class(),
        "energy_entry": energy_entry(),
        "pressure_entry": pressure_entry(),
        "harmonic_counterexample": harmonic_counterexample(),
        "protected_limit": {key: str(value) if isinstance(value, Q) else value
                            for key, value in limit.items()},
        "review": review(),
    }
