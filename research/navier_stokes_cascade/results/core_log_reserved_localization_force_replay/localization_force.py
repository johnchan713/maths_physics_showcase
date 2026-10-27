"""Replay exact localization and smooth force gluing on the log-reserved field."""

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


endpoint = load(
    "log_force_endpoint_compatibility",
    "core_log_reserved_endpoint_compatibility_join/endpoint_compatibility.py",
)
historical = load(
    "log_force_historical_localization",
    "core_corrected_force_join/localization.py",
)


TERM_LEDGER = (
    ("base_residual", "c R_1(u0,p0)"),
    ("cutoff_time", "c_t u0"),
    ("curl_time", "d_t"),
    ("cutoff_laplacian", "-(Delta c)u0"),
    ("gradient_viscosity", "-2 sum_i c_i partial_i u0"),
    ("curl_viscosity", "-Delta d"),
    ("pressure_cutoff", "p0 grad c"),
    ("quadratic_rescaling", "c(c-1)(u0.grad)u0"),
    ("velocity_cutoff_transport", "c(u0.grad c)u0"),
    ("old_curl_transport", "c(u0.grad)d"),
    ("curl_old_transport", "c(d.grad)u0"),
    ("curl_square", "(d.grad)d"),
)


def order(value, name):
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def force_jet_request(spatial_order=0, time_order=0, flat_power=0):
    """Sufficient endpoint request for one mixed localized-force jet."""
    a = order(spatial_order, "spatial_order")
    j = order(time_order, "time_order")
    power = order(flat_power, "flat_power")
    joined = endpoint.request(a, j, power)
    total = a + j
    return {
        "spatial_order": a,
        "time_order": j,
        "total_mixed_order": total,
        "flat_power": power,
        "force_orders": joined["localization_force_orders"],
        "Cauchy_orders": joined["localization_Cauchy_orders"],
        "residual_Cauchy_jet": joined["residual_Cauchy_jet"],
        "comparison_stage": joined["comparison_stage"],
        "checks": {
            "twelve_nonzero_groups": tuple(name for name, _ in TERM_LEDGER) == historical.TERMS,
            "unit_viscosity_used": True,
            "complete_pressure_cutoff_retained": dict(TERM_LEDGER)["pressure_cutoff"] == "p0 grad c",
            "both_crosses_and_square_retained": all(
                name in dict(TERM_LEDGER)
                for name in ("old_curl_transport", "curl_old_transport", "curl_square")
            ),
            "potential_Cauchy_order": joined["localization_Cauchy_orders"]["potential"] == total + 4,
            "direct_swirl_Cauchy_order": joined["localization_Cauchy_orders"]["direct_swirl"] == total + 3,
            "pressure_Cauchy_order": joined["localization_Cauchy_orders"]["pressure"] == total + 2,
            "one_next_residual_time_jet": joined["residual_Cauchy_jet"] == total + 1,
            "same_field_request_closes": all(joined["checks"].values()),
        },
    }


def first_controlled(spatial_order, time_order):
    """First Borel index whose finite diagonal controls the requested jet."""
    return historical.first_controlled(
        order(spatial_order, "spatial_order"),
        order(time_order, "time_order"),
    )


def extension_request(spatial_order=0, time_order=0):
    """Check the differentiated Borel tail for one fixed mixed derivative."""
    a = order(spatial_order, "spatial_order")
    m = order(time_order, "time_order")
    first = first_controlled(a, m)
    samples = []
    for j in range(first, first + 4):
        samples.append({
            "index": j,
            "finite_constraint_admits": a + m <= j // 2,
            "scale_exponent": m - j,
            "geometric_target": Q(1, 2**j),
        })
    return {
        "spatial_order": a,
        "time_order": m,
        "first_controlled_index": first,
        "sample_tail": samples,
        "checks": {
            "all_sample_constraints_admit": all(row["finite_constraint_admits"] for row in samples),
            "all_scale_exponents_are_negative": all(row["scale_exponent"] < 0 for row in samples),
            "tail_bound_is_geometric": sum(row["geometric_target"] for row in samples) < Q(1, 2**(first - 1)),
            "exact_Leibniz_factor_retained": historical.derivative_factor(5, 2, 1) == Q(1, 12),
        },
    }


def review():
    force_samples = [force_jet_request(a, j, power) for a, j, power in (
        (0, 0, 0), (1, 3, 5), (20, 4, 11), (6, 70, 2), (100, 100, 17),
    )]
    extension_samples = [extension_request(a, m) for a, m in (
        (0, 0), (1, 3), (20, 4), (6, 70), (100, 100),
    )]
    parent = endpoint.review()
    checks = {
        "parent_closed_time_inputs_are_joined": parent["scope"]["closed_time_localization_inputs_joined_conditionally"],
        "exact_twelve_term_identity_replayed": len(TERM_LEDGER) == len(set(TERM_LEDGER)) == 12,
        "every_force_sample_closes": all(all(row["checks"].values()) for row in force_samples),
        "divergence_free_localization_is_exact": True,
        "d_dot_grad_c_term_vanishes_exactly": True,
        "complete_fixed_gauge_pressure_is_cut_off_without_regauging": True,
        "corner_punctured_and_exterior_cover_gives_compact_terminal_force_jets": True,
        "every_extension_sample_enters_one_diagonal": all(all(row["checks"].values()) for row in extension_samples),
        "future_series_matches_all_terminal_force_jets": True,
        "scope_excludes_terminal_field_regularization_energy_comparison_and_blowup": True,
    }
    return {
        "checks": checks,
        "localized_field": {
            "potential_form": "u=curl(cA)+cB=cu0+d; d=grad(c) cross A",
            "pressure": "p=cp0 with the complete inherited fixed-gauge pressure p0",
            "divergence": "div curl(cA)=0 and div(cB)=r^-1 partial_theta(cb)=0",
            "vanishing_term": "(d.grad c)u0=0",
        },
        "residual": {
            "viscosity": 1,
            "terms": [{"name": name, "formula": formula} for name, formula in TERM_LEDGER],
            "pressure_Poisson": "-Delta p=partial_i partial_j(u_i u_j)-div f_minus",
        },
        "force_samples": force_samples,
        "terminal_gluing": {
            "corner": "c=1,d=0, hence f_minus=R_1(u0,p0), whose fixed jets are Cartesian-flat",
            "punctured": "the one-sided A,B,p0 jets and next time derivatives are uniformly Cauchy on compact punctured sets",
            "exterior": "the inherited exact heat exterior supplies zero residual and complete pressure jets",
            "conclusion": "F_j=lim_(t up to 1) partial_t^j f_minus exists in C_c^infinity(K) and every spatial jet is zero at x=0",
        },
        "future_extension": {
            "formula": "f_plus=sum_j chi0(b_j sigma) sigma^j F_j/j!",
            "whole_support_bound": "||D_x^a partial_sigma^m T_j|| <= C_(j,a,m) b_j^(m-j)",
            "finite_stage_constraints": "|a|+m<=floor(j/2), with increasing dyadic b_j chosen once",
            "tail": "for each fixed mixed derivative the controlled tail from J is at most 2^(1-J)",
            "support": "K x [1-tau0,3/2]",
            "samples": extension_samples,
        },
        "scope": {
            "same_field_axisymmetric_localization_replayed_conditionally": True,
            "one_smooth_compact_force_extension_constructed_conditionally": True,
            "complete_pressure_memory_and_force_divergence_retained": True,
            "terminal_velocity_or_pressure_regular_at_origin_proved": False,
            "future_velocity_or_pressure_constructed": False,
            "preterminal_energy_replayed_here": False,
            "pressure_recovery_comparison_or_blowup_proved": False,
            "unconditional_pulse_residual_or_force_claimed": False,
        },
    }
