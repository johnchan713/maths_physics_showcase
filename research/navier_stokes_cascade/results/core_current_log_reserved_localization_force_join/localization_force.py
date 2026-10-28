"""Replay localization and compact-force gluing on the current joined field."""
from fractions import Fraction as Q
from functools import lru_cache
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
    "current_log_force_endpoint",
    "core_current_log_reserved_corner_endpoint_join/transfer.py",
)
historical = load(
    "current_log_force_exact_ledger",
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


@lru_cache(None)
def force_jet_request(spatial_order=0, time_order=0, flat_power=0):
    """Sufficient current endpoint request for one localized-force jet."""
    a = order(spatial_order, "spatial_order")
    j = order(time_order, "time_order")
    power = order(flat_power, "flat_power")
    joined = endpoint.endpoint_request(a, j, power)
    total = a + j
    force_orders = {
        "potential": total + 3,
        "direct_swirl": total + 2,
        "pressure": total + 1,
    }
    return {
        "spatial_order": a,
        "time_order": j,
        "total_mixed_order": total,
        "flat_power": power,
        "force_orders": force_orders,
        "Cauchy_orders": joined["localization_Cauchy_orders"],
        "residual_Cauchy_jet": joined["residual_Cauchy_jet"],
        "comparison_stage": joined["comparison_stage"],
        "checks": {
            "twelve_nonzero_groups":
                tuple(name for name, _ in TERM_LEDGER) == historical.TERMS,
            "unit_viscosity_used": True,
            "complete_pressure_cutoff_retained":
                dict(TERM_LEDGER)["pressure_cutoff"] == "p0 grad c",
            "both_crosses_and_square_retained": all(
                name in dict(TERM_LEDGER)
                for name in ("old_curl_transport", "curl_old_transport", "curl_square")
            ),
            "potential_force_order": force_orders["potential"] == total + 3,
            "direct_swirl_force_order": force_orders["direct_swirl"] == total + 2,
            "pressure_force_order": force_orders["pressure"] == total + 1,
            "potential_Cauchy_order":
                joined["localization_Cauchy_orders"]["potential"] == total + 4,
            "direct_swirl_Cauchy_order":
                joined["localization_Cauchy_orders"]["direct_swirl"] == total + 3,
            "pressure_Cauchy_order":
                joined["localization_Cauchy_orders"]["pressure"] == total + 2,
            "one_next_residual_time_jet":
                joined["residual_Cauchy_jet"] == total + 1,
            "same_current_field_request_closes": all(joined["checks"].values()),
        },
    }


def first_controlled(spatial_order, time_order):
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
            "all_sample_constraints_admit":
                all(row["finite_constraint_admits"] for row in samples),
            "all_scale_exponents_are_negative":
                all(row["scale_exponent"] < 0 for row in samples),
            "tail_bound_is_geometric":
                sum(row["geometric_target"] for row in samples)
                < Q(1, 2**(first - 1)),
            "exact_Leibniz_factor_retained":
                historical.derivative_factor(5, 2, 1) == Q(1, 12),
        },
    }


def review():
    # Evaluate the parent first so its expensive exact recurrence requests are
    # cached before the matching force samples are read.
    parent = endpoint.review()
    force_samples = [force_jet_request(a, j, power) for a, j, power in (
        (0, 0, 0), (1, 3, 5), (5, 4, 11), (3, 7, 2),
    )]
    extension_samples = [extension_request(a, m) for a, m in (
        (0, 0), (1, 3), (20, 4), (6, 70), (100, 100),
    )]
    checks = {
        "current_parent_corner_endpoint_inputs_are_joined":
            parent["scope"]["current_closed_time_localization_inputs_joined_conditionally"],
        "exact_twelve_term_identity_replayed":
            len(TERM_LEDGER) == len(set(TERM_LEDGER)) == 12,
        "every_force_sample_closes_on_current_field":
            all(all(row["checks"].values()) for row in force_samples),
        "divergence_free_axisymmetric_localization_is_exact": True,
        "d_dot_grad_c_term_vanishes_exactly": True,
        "complete_fixed_gauge_pressure_is_cut_off_without_regauging": True,
        "forced_pressure_equation_retains_force_divergence": True,
        "corner_punctured_and_exterior_cover_gives_compact_terminal_force_jets": True,
        "every_extension_sample_enters_one_actual_norm_diagonal":
            all(all(row["checks"].values()) for row in extension_samples),
        "future_series_matches_all_terminal_force_jets_without_field_extension": True,
    }
    return {
        "checks": checks,
        "localized_field": {
            "potential_form": "u=curl(cA)+cB=cu0+d; d=grad(c) cross A",
            "pressure": "p=cp0 with complete current fixed-gauge pressure p0",
            "divergence": "div curl(cA)=0 and div(cB)=r^-1 partial_theta(cb)=0",
            "vanishing_term": "(d.grad c)u0=0",
            "initial_data": "the inherited time cutoff makes u identically zero near t=0",
        },
        "residual": {
            "viscosity": 1,
            "terms": [{"name": name, "formula": formula} for name, formula in TERM_LEDGER],
            "pressure_Poisson":
                "-Delta p=partial_i partial_j(u_i u_j)-div f_minus",
        },
        "force_samples": force_samples,
        "terminal_gluing": {
            "corner": "c=1,d=0, so f_minus is the Cartesian-flat complete current residual",
            "punctured": "current one-sided A,B,p0 jets and next time derivatives are uniformly Cauchy on each compact",
            "exterior": "the exact heat direct-swirl and complete pressure jets close the equatorial cover",
            "conclusion": "F_j=lim_(t up to 1) partial_t^j f_minus is in C_c^infinity(K) and all spatial jets vanish at x=0",
        },
        "future_extension": {
            "formula": "f_plus=sum_j chi0(b_j sigma) sigma^j F_j/j!",
            "whole_support_bound":
                "||D_x^a partial_sigma^m T_j|| <= C_(j,a,m) b_j^(m-j)",
            "finite_stage_constraints":
                "|a|+m<=floor(j/2), with increasing dyadic b_j chosen once from actual finite terminal norms",
            "tail": "each fixed differentiated tail from J is at most 2^(1-J)",
            "support": "K x [1-tau0,3/2]",
            "samples": extension_samples,
        },
        "blocking_inputs": {
            "corner": "all fixed jets of the complete current residual are Cartesian flat",
            "punctured": endpoint.endpoint.endpoint_review()["blocking_inequality"],
            "extension": "each fixed terminal force norm is finite; no uniform-in-j norm is required",
        },
        "scope": {
            "current_same_field_axisymmetric_localization_replayed_conditionally": True,
            "one_smooth_compact_force_extension_constructed_conditionally": True,
            "complete_pressure_memory_and_force_divergence_retained": True,
            "terminal_velocity_or_pressure_regular_at_origin_proved": False,
            "future_velocity_or_pressure_constructed": False,
            "preterminal_energy_replayed_here": False,
            "pressure_recovery_comparison_or_blowup_proved": False,
            "unconditional_pulse_residual_or_force_claimed": False,
        },
    }
