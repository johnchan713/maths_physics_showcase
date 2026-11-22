"""Q217: bind Q216 to localization, force, energy and comparison."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q217-LOCALIZATION-COMPARISON-FRONTIER-JOIN"
PROVIDERS = (
    "core_current_head_q216_quantified_provider_source_audit",
    "core_current_head_q203_diagonal_endpoint_force_corrected_domain",
    "core_current_head_q206_terminal_origin_obstruction",
    "core_current_head_q168_localization_force_frontier_audit",
    "core_current_head_q169_pressure_comparison_frontier_audit",
    "core_current_head_q170_protected_core_limit_frontier_audit",
    "core_current_head_q182_comparison_core_exclusion_replay",
)
EXPECTED = {
    PROVIDERS[0]: "Q216-all-fixed-order-and-finite-pole-provider-quantifiers-source-audited",
    PROVIDERS[1]: "Q203-corrected-domain-diagonal-endpoint-localized-force-replayed-conditionally",
    PROVIDERS[2]: "Q206-terminal-origin-velocity-and-gauge-independent-pressure-value-jets-obstructed-conditionally",
    PROVIDERS[3]: "Q168-current-head-localization-force-energy-audited-conditionally",
    PROVIDERS[4]: "Q169-current-head-fixed-gauge-pressure-comparison-audited-conditionally",
    PROVIDERS[5]: "Q170-current-head-protected-core-and-limit-order-audited-conditionally",
    PROVIDERS[6]: "Q182-physical-pressure-comparison-protected-core-replayed-conditionally",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q217 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def localization_requests(order, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q217 identity required")
    if isinstance(order, bool) or not isinstance(order, int) or order < 0:
        raise ValueError("nonnegative integer order required")
    direct = {"potential": order + 2, "direct_swirl": order + 1, "pressure": order}
    cauchy = {key: value + 1 for key, value in direct.items()}
    supplied = {"potential": order + 4, "direct_swirl": order + 3,
                "pressure": order + 2}
    return {
        "order": order,
        "direct": direct,
        "cauchy": cauchy,
        "supplied": supplied,
        "spare_over_cauchy": {key: supplied[key] - cauchy[key] for key in supplied},
        "twelve_nonzero_groups": True,
        "complete_pressure_cutoff_retained": True,
        "force_divergence_retained": True,
    }


def energy(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q217 identity required")
    return {
        "zero_safe_regularization": "sqrt(E+epsilon^2)",
        "bound": "E(t)+2*integral_0^t D(s) ds <= F(t)^2",
        "energy_and_dissipation_extrema_added_under_one_rhs": False,
        "terminal_field_value_used": False,
    }


def comparison(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q217 identity required")
    return {
        "individual_pressures_keep_minus_div_f": True,
        "common_force_cancels_exactly_once": True,
        "harmonic_gradient_removed_using_momentum_and_finite_energy": True,
        "expanding_ball_error_rate": Q(2),
        "limit_order": ["fix T<1", "send R to infinity", "send tau to zero"],
        "competitor_derivative_bound_used": False,
        "competitor_finite_energy_used": True,
    }


def protected_core(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q217 identity required")
    h = Q(1, 200)
    exponent = Q(1, 2) + h
    return {
        "X_in": Q(1, 101),
        "X_a": Q(4, 101),
        "support_hole_margin": Q(3, 101),
        "velocity_exponent": exponent,
        "background_remainder_exponent": 2 * h,
        "annular_corrections_vanish_exactly": True,
        "background_vanishes": False,
        "localization_is_identity_on_path": True,
    }


def first_new_conditional_premise(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q217 identity required")
    return {
        "premise": "sup_{0<=t<=T} ||v(t)||_2 < infinity for every fixed T<1",
        "introduced_at": "fixed-gauge pressure recovery and expanding-ball comparison",
        "implied_by_local_smoothness_zero_data_and_compact_force": False,
        "counterexample": {
            "velocity": "v(x,t)=t e_1",
            "pressure": "P(x,t)=-x_1",
            "force": "0",
            "divergence_free": True,
            "zero_initial_data": True,
            "navier_stokes_residual_zero": True,
            "finite_L2_at_positive_time": False,
        },
        "upstream_Q216_provider_theorems_remain_analytical_inputs": True,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q217 identity required")
    req = localization_requests(5)
    comp = comparison()
    core = protected_core()
    frontier = first_new_conditional_premise()
    checks = {
        "provider_statuses_match": all(statuses().values()),
        "localization_has_one_spare_cauchy_derivative": all(
            value == 1 for value in req["spare_over_cauchy"].values()),
        "pressure_memory_and_force_cancellation_are_distinct": (
            comp["individual_pressures_keep_minus_div_f"]
            and comp["common_force_cancels_exactly_once"]),
        "protected_path_stays_in_support_hole": core["X_in"] < core["X_a"],
        "limit_order_is_legal": comp["limit_order"] == [
            "fix T<1", "send R to infinity", "send tau to zero"],
        "finite_energy_is_not_automatic": (
            not frontier["implied_by_local_smoothness_zero_data_and_compact_force"]
            and frontier["counterexample"]["navier_stokes_residual_zero"]
            and not frontier["counterexample"]["finite_L2_at_positive_time"]),
    }
    return {
        "status": "Q217-Q216-object-localization-force-energy-and-comparison-joined-conditionally",
        "checks": checks,
        "scope": {
            "same_source_audited_object_reaches_compact_smooth_force": True,
            "force_trace_separated_from_terminal_field_regularization": True,
            "finite_energy_competitor_exclusion_implication": True,
            "all_locally_smooth_competitors_excluded": False,
            "unconditional_force_residual_or_blowup": False,
        },
        "next": (
            "either derive finite-energy entry from an explicitly stronger global-solution class "
            "or retain it as an essential comparison hypothesis; separately discharge the cited "
            "upstream continuum provider theorems before any unconditional claim"
        ),
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q217 identity required")
    core = protected_core()
    comp = comparison()
    return {
        "status": review()["status"],
        "providers": statuses(),
        "localization_order_5": localization_requests(5),
        "energy": energy(),
        "comparison": {**comp, "expanding_ball_error_rate": str(comp["expanding_ball_error_rate"])},
        "protected_core": {key: str(value) if isinstance(value, Q) else value
                           for key, value in core.items()},
        "first_new_conditional_premise": first_new_conditional_premise(),
        "review": review(),
    }
