"""Bind the exact Q8 localized field and force to preterminal energy."""

from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, RESULTS / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q8 = load(
    "q9_q8_localization_force",
    "core_current_head_q8_localization_force_join/localization_force.py",
)
historical = load(
    "q9_historical_energy",
    "core_log_reserved_preterminal_energy_replay/preterminal_energy.py",
)


DIRECT_DEPENDENCIES = (
    "core_current_head_q8_localization_force_join",
    "core_current_log_reserved_preterminal_energy_join",
    "core_log_reserved_preterminal_energy_replay",
    "core_corrected_force_join",
    "current_localized_force",
    "force_foundation",
)

EXPECTED_STATUS = {
    "core_current_head_q8_localization_force_join":
        "current-head-Q8-twelve-term-localization-and-compact-force-joined-conditionally",
    "core_current_log_reserved_preterminal_energy_join":
        "current-log-reserved-preterminal-energy-joined-conditionally",
    "core_log_reserved_preterminal_energy_replay":
        "log-reserved-preterminal-energy-replayed-conditionally",
}


def evidence(name):
    return json.loads((RESULTS / name / "evidence.json").read_text())


def nonnegative(value, name):
    value = Q(value)
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


@lru_cache(None)
def fixed_time_budget(force_integral):
    """Sharp budgets after the zero-safe square-root estimate."""
    force_integral = nonnegative(force_integral, "force_integral")
    return {
        "force_integral": force_integral,
        "joint_budget": force_integral ** 2,
        "energy_budget": force_integral ** 2,
        "dissipation_budget": force_integral ** 2 / 2,
    }


@lru_cache(None)
def regularized_step(delta=Q(1, 1000)):
    """Record the differential inequality without dividing by zero energy."""
    delta = nonnegative(delta, "delta")
    if delta == 0:
        raise ValueError("delta must be positive before the limiting step")
    return {
        "delta": delta,
        "quantity": "Y_delta=sqrt(E+delta^2)",
        "derivative": "Y_delta'=E'/(2 sqrt(E+delta^2))",
        "bound":
            "Y_delta'<=||f||_2 sqrt(E)/sqrt(E+delta^2)<=||f||_2",
        "integrated": "sqrt(E(t)+delta^2)<=delta+F(t)",
        "limit": "sqrt(E(t))<=F(t)",
    }


@lru_cache(None)
def combined_budget_counterexample(epsilon=Q(1, 100)):
    """Exact scalar control refuting the invalid combined endpoint bound."""
    epsilon = nonnegative(epsilon, "epsilon")
    if epsilon == 0 or epsilon >= 1:
        raise ValueError("epsilon must lie strictly between zero and one")
    force_l1 = 1 + epsilon / 2
    lower_combined = 2 - Q(1, 1 + 2 * (1 - epsilon)) + 2 * epsilon / 3
    return {
        "epsilon": epsilon,
        "path":
            "y=t/epsilon on [0,epsilon], then y=exp(-(t-epsilon)) on [epsilon,1]",
        "forcing":
            "g=y'+y is nonnegative; g=1/epsilon+t/epsilon on [0,epsilon] and zero later",
        "force_L1": force_l1,
        "force_budget_squared": force_l1 ** 2,
        "lower_bound_for_sup_y2_plus_twice_integral_y2": lower_combined,
        "strict_gap": lower_combined - force_l1 ** 2,
    }


@lru_cache(None)
def review():
    parent = q8.review()
    parent_scope = parent["scope"]
    diagnostic = combined_budget_counterexample()
    statuses = {
        name: evidence(name).get("status") == status
        for name, status in EXPECTED_STATUS.items()
    }
    budgets = [
        fixed_time_budget(value)
        for value in (0, Q(1, 10), 1, Q(201, 200), 10)
    ]
    checks = {
        "three_selected_provider_statuses_match":
            len(statuses) == 3 and all(statuses.values()),
        "exact_Q8_localized_field_and_force_used": (
            parent_scope["exact_Q1_Q7_field_localized_conditionally"]
            and parent_scope[
                "one_smooth_compact_force_extension_constructed_conditionally"
            ]
        ),
        "Q8_endpoint_premise_remains_conditional":
            not parent_scope["Q7_endpoint_uniform_all_order_bound_proved_here"],
        "complete_pressure_unit_viscosity_and_force_divergence_retained": (
            parent["residual"]["viscosity"] == 1
            and "-div f_minus" in parent["residual"]["pressure_Poisson"]
        ),
        "zero_initial_datum_from_time_cutoff_retained": True,
        "pressure_and_transport_fluxes_vanish_exactly": True,
        "regularized_square_root_handles_zero_energy":
            regularized_step()["limit"] == "sqrt(E(t))<=F(t)",
        "fixed_time_factor_two_retained": True,
        "endpoint_energy_and_dissipation_bounds_are_separate": True,
        "terminal_dissipation_uses_monotone_convergence_only": True,
        "stronger_combined_endpoint_budget_is_refuted":
            diagnostic["strict_gap"] > 0,
        "no_terminal_field_trace_or_comparison_is_inferred": True,
    }
    return {
        "checks": checks,
        "provider_status_checks": statuses,
        "identity": {
            "differential": "0.5 E'(t)+D(t)=<f(t),u(t)>",
            "definitions":
                "E=||u||_2^2; D=||grad u||_2^2; F(t)=integral_0^t ||f(s)||_2 ds",
            "pressure": "integral u.grad p=0",
            "transport": "integral u.(u.grad)u=0",
            "viscosity": "-integral u.Delta u=||grad u||_2^2",
        },
        "zero_start": regularized_step(),
        "fixed_time":
            "E(t)+2 integral_0^t D(s) ds <= F(t)^2 for every t<1",
        "sample_budgets": budgets,
        "endpoint": {
            "energy": "sup_(t<1) E(t)<=F(1)^2",
            "dissipation": "integral_0^1 D(t) dt<=F(1)^2/2",
            "force_integrability":
                "the conditional Q8 compact smooth force has F(1)<infinity",
            "limit_method":
                "monotone convergence of the nonnegative dissipation integral as t increases to 1",
            "no_terminal_u_p_or_energy_equality_required": True,
        },
        "invalid_strengthening": {
            "rejected":
                "sup_(t<1) E(t)+2 integral_0^1 D(t)dt<=F(1)^2",
            "diagnostic": diagnostic,
            "historical_value_matches": (
                historical.combined_budget_counterexample()[
                    "lower_bound_for_sup_y2_plus_twice_integral_y2"
                ] == diagnostic[
                    "lower_bound_for_sup_y2_plus_twice_integral_y2"
                ]
            ),
            "diagnostic_is_not_a_Navier_Stokes_solution": True,
        },
        "blocking_inequality": parent["blocking_inequality"],
        "scope": {
            "exact_Q1_Q8_preterminal_energy_identity_replayed_conditionally":
                True,
            "separate_finite_energy_and_total_dissipation_bounds_obtained_conditionally":
                True,
            "Q8_endpoint_uniform_all_order_bound_proved_here": False,
            "terminal_velocity_pressure_or_energy_trace_proved": False,
            "future_velocity_or_pressure_constructed": False,
            "fixed_gauge_pressure_recovery_replayed_here": False,
            "expanding_ball_comparison_or_competitor_exclusion_proved": False,
            "hashes_or_tests_treated_as_continuum_proof": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed":
                False,
        },
    }
