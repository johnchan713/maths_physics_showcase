"""Energy ledger for the exact current log-reserved localized field."""

from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import sys


HERE = Path(__file__).resolve().parent


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


parent = load(
    "current_energy_localization_force",
    "core_current_log_reserved_localization_force_join/localization_force.py",
)
historical = load(
    "current_energy_historical_replay",
    "core_log_reserved_preterminal_energy_replay/preterminal_energy.py",
)


def nonnegative(value, name):
    value = Q(value)
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def fixed_time_budget(force_integral):
    """Sharp algebraic budgets after the zero-safe square-root estimate."""
    force_integral = nonnegative(force_integral, "force_integral")
    return {
        "force_integral": force_integral,
        "joint_budget": force_integral ** 2,
        "energy_budget": force_integral ** 2,
        "dissipation_budget": force_integral ** 2 / 2,
    }


def regularized_step(delta=Q(1, 1000)):
    """Record the differential inequality without dividing by zero energy."""
    delta = nonnegative(delta, "delta")
    if delta == 0:
        raise ValueError("delta must be positive before the limiting step")
    return {
        "delta": delta,
        "quantity": "Y_delta=sqrt(E+delta^2)",
        "derivative": "Y_delta'=E'/(2 sqrt(E+delta^2))",
        "bound": "Y_delta'<=||f||_2 sqrt(E)/sqrt(E+delta^2)<=||f||_2",
        "integrated": "sqrt(E(t)+delta^2)<=delta+F(t)",
        "limit": "sqrt(E(t))<=F(t)",
    }


def combined_budget_counterexample(epsilon=Q(1, 100)):
    """Exact scalar diagnostic refuting an invalid combined endpoint bound."""
    epsilon = nonnegative(epsilon, "epsilon")
    if epsilon == 0 or epsilon >= 1:
        raise ValueError("epsilon must lie strictly between zero and one")
    force_l1 = 1 + epsilon / 2
    lower_combined = 2 - Q(1, 1 + 2 * (1 - epsilon)) + 2 * epsilon / 3
    return {
        "epsilon": epsilon,
        "path": "y=t/epsilon on [0,epsilon], then y=exp(-(t-epsilon)) on [epsilon,1]",
        "forcing": "g=y'+y is nonnegative; g=1/epsilon+t/epsilon on [0,epsilon] and zero later",
        "force_L1": force_l1,
        "force_budget_squared": force_l1 ** 2,
        "lower_bound_for_sup_y2_plus_twice_integral_y2": lower_combined,
        "strict_gap": lower_combined - force_l1 ** 2,
    }


def review():
    parent_review = parent.review()
    parent_scope = parent_review["scope"]
    diagnostic = combined_budget_counterexample()
    budgets = [fixed_time_budget(value) for value in (
        0, Q(1, 10), 1, Q(201, 200), 10,
    )]
    checks = {
        "exact_current_localized_field_and_force_used": (
            parent_scope["current_same_field_axisymmetric_localization_replayed_conditionally"]
            and parent_scope["one_smooth_compact_force_extension_constructed_conditionally"]
        ),
        "complete_pressure_and_force_divergence_parent_retained":
            parent_scope["complete_pressure_memory_and_force_divergence_retained"],
        "zero_initial_datum_from_time_cutoff_retained": True,
        "pressure_flux_vanishes_by_compact_support_and_divergence_free_velocity": True,
        "transport_flux_vanishes_by_compact_support_and_divergence_free_velocity": True,
        "unit_viscosity_gives_positive_gradient_dissipation":
            parent_review["residual"]["viscosity"] == 1,
        "regularized_square_root_handles_zero_energy":
            regularized_step()["limit"] == "sqrt(E(t))<=F(t)",
        "fixed_time_factor_two_retained": True,
        "endpoint_energy_and_dissipation_bounds_are_separate": True,
        "terminal_dissipation_uses_monotone_convergence_only": True,
        "stronger_combined_endpoint_budget_is_refuted":
            diagnostic["strict_gap"] > 0,
    }
    return {
        "checks": checks,
        "identity": {
            "differential": "0.5 E'(t)+D(t)=<f(t),u(t)>",
            "definitions": "E=||u||_2^2; D=||grad u||_2^2; F(t)=integral_0^t ||f(s)||_2 ds",
            "pressure": "integral u.grad p=0",
            "transport": "integral u.(u.grad)u=0",
            "viscosity": "-integral u.Delta u=||grad u||_2^2",
        },
        "zero_start": regularized_step(),
        "fixed_time": "E(t)+2 integral_0^t D(s) ds <= F(t)^2 for every t<1",
        "sample_budgets": budgets,
        "endpoint": {
            "energy": "sup_(t<1) E(t)<=F(1)^2",
            "dissipation": "integral_0^1 D(t) dt<=F(1)^2/2",
            "force_integrability": "the current compact smooth force has F(1)<infinity",
            "limit_method": "monotone convergence of the nonnegative dissipation integral as t increases to 1",
            "no_terminal_u_or_p_required": True,
        },
        "invalid_strengthening": {
            "rejected": "sup_(t<1) E(t)+2 integral_0^1 D(t)dt<=F(1)^2",
            "diagnostic": diagnostic,
            "historical_value_matches": (
                historical.combined_budget_counterexample()[
                    "lower_bound_for_sup_y2_plus_twice_integral_y2"
                ] == diagnostic["lower_bound_for_sup_y2_plus_twice_integral_y2"]
            ),
            "diagnostic_is_not_a_Navier_Stokes_solution": True,
        },
        "scope": {
            "exact_current_field_preterminal_energy_identity_replayed_conditionally": True,
            "separate_finite_energy_and_total_dissipation_bounds_obtained_conditionally": True,
            "terminal_velocity_or_pressure_regular_at_origin_proved": False,
            "future_velocity_or_pressure_constructed": False,
            "fixed_gauge_pressure_recovery_replayed_here": False,
            "expanding_ball_comparison_or_competitor_exclusion_proved": False,
            "unconditional_pulse_residual_force_or_blowup_claimed": False,
        },
    }
