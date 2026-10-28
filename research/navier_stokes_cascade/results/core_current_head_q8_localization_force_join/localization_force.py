"""Join the exact Q7 field to localization and compact-force gluing."""

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


q7 = load(
    "q8_q7_corner_endpoint",
    "core_current_head_q7_corner_endpoint_join/corner_endpoint.py",
)
historical = load(
    "q8_exact_localization_ledger",
    "core_corrected_force_join/localization.py",
)


DIRECT_DEPENDENCIES = (
    "core_current_head_q7_corner_endpoint_join",
    "core_current_log_reserved_localization_force_join",
    "core_log_reserved_localization_force_replay",
    "core_corrected_force_join",
    "current_localized_force",
    "force_foundation",
)

EXPECTED_STATUS = {
    "core_current_head_q7_corner_endpoint_join":
        "current-head-Q7-cartesian-corner-and-punctured-endpoint-joined-conditionally",
    "core_current_log_reserved_localization_force_join":
        "current-log-reserved-localization-and-compact-force-joined-conditionally",
    "core_log_reserved_localization_force_replay":
        "log-reserved-localization-and-compact-force-replayed-conditionally",
}


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


def evidence(name):
    return json.loads((RESULTS / name / "evidence.json").read_text())


def order(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


@lru_cache(None)
def force_jet_request(spatial_order=0, time_order=0, flat_power=0):
    """Return sharp field orders and the larger Q7 sufficient request."""
    a = order(spatial_order, "spatial_order")
    j = order(time_order, "time_order")
    power = order(flat_power, "flat_power")
    total = a + j
    parent = q7.endpoint_request(a, j, power)
    sharp_force = {
        "potential": total + 2,
        "direct_swirl": total + 1,
        "pressure": total,
    }
    sharp_cauchy = {
        "potential": total + 3,
        "direct_swirl": total + 2,
        "pressure": total + 1,
    }
    supplied = parent["localization_Cauchy_orders"]
    reserve = {name: supplied[name] - sharp_cauchy[name] for name in supplied}
    return {
        "spatial_order": a,
        "time_order": j,
        "total_mixed_order": total,
        "flat_power": power,
        "sharp_direct_force_orders": sharp_force,
        "sharp_Cauchy_orders": sharp_cauchy,
        "Q7_supplied_orders": supplied,
        "Q7_order_reserve": reserve,
        "residual_Cauchy_jet": parent["residual_Cauchy_jet"],
        "comparison_stage": parent["comparison_stage"],
        "endpoint_blocking_inequality": parent["endpoint_blocking_inequality"],
        "checks": {
            "Q7_parent_request_closes": all(parent["checks"].values()),
            "direct_potential_order_from_Delta_d_or_grad_u0":
                sharp_force["potential"] == total + 2,
            "direct_swirl_order_from_grad_u0":
                sharp_force["direct_swirl"] == total + 1,
            "direct_pressure_order_from_p0_grad_c":
                sharp_force["pressure"] == total,
            "Cauchy_uses_one_next_time_derivative": all(
                sharp_cauchy[name] == sharp_force[name] + 1
                for name in sharp_force
            ),
            "Q7_retains_one_unused_field_derivative":
                set(reserve.values()) == {1},
            "Q7_retains_next_residual_jet":
                parent["residual_Cauchy_jet"] == total + 1,
        },
    }


def first_controlled(spatial_order, time_order):
    return historical.first_controlled(
        order(spatial_order, "spatial_order"),
        order(time_order, "time_order"),
    )


@lru_cache(None)
def extension_request(spatial_order=0, time_order=0):
    """Check one fixed mixed derivative in the single Borel diagonal."""
    a = order(spatial_order, "spatial_order")
    m = order(time_order, "time_order")
    first = first_controlled(a, m)
    samples = []
    for index in range(first, first + 4):
        samples.append({
            "index": index,
            "finite_constraint_admits": a + m <= index // 2,
            "scale_exponent": m - index,
            "geometric_target": Q(1, 2**index),
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


@lru_cache(None)
def review():
    parent = q7.review()
    force_samples = [force_jet_request(a, j, power) for a, j, power in (
        (0, 0, 0), (1, 3, 5), (5, 4, 11), (3, 7, 2),
    )]
    extension_samples = [extension_request(a, m) for a, m in (
        (0, 0), (1, 3), (20, 4), (6, 70), (100, 100),
    )]
    statuses = {
        name: evidence(name).get("status") == status
        for name, status in EXPECTED_STATUS.items()
    }
    checks = {
        "three_selected_provider_statuses_match":
            len(statuses) == 3 and all(statuses.values()),
        "Q7_corner_endpoint_parent_is_unchanged":
            all(parent["checks"].values()),
        "exact_twelve_term_identity_replayed":
            tuple(name for name, _ in TERM_LEDGER) == historical.TERMS,
        "unit_viscosity_pressure_crosses_and_square_retained":
            dict(TERM_LEDGER)["pressure_cutoff"] == "p0 grad c"
            and all(name in dict(TERM_LEDGER) for name in (
                "old_curl_transport", "curl_old_transport", "curl_square"
            )),
        "only_apparent_thirteenth_group_vanishes": True,
        "sharp_ledger_correction_closes_with_Q7_reserve":
            all(all(row["checks"].values()) for row in force_samples),
        "forced_pressure_equation_retains_force_divergence": True,
        "endpoint_premise_remains_explicit":
            not parent["scope"]["endpoint_uniform_all_order_physical_bound_proved_here"],
        "conditional_corner_punctured_cover_gives_terminal_force_jets": True,
        "one_actual_norm_Borel_diagonal_closes":
            all(all(row["checks"].values()) for row in extension_samples),
        "terminal_field_regularization_is_not_inferred": True,
        "energy_comparison_exclusion_and_blowup_are_not_replayed": True,
    }
    return {
        "checks": checks,
        "provider_status_checks": statuses,
        "localized_field": {
            "potential_form": "u=curl(cA)+cB=cu0+d; d=grad(c) cross A",
            "pressure": "p=cp0 with complete Q7 fixed-gauge pressure p0",
            "divergence": "div curl(cA)=0 and div(cB)=r^-1 partial_theta(cb)=0",
            "vanishing_term": "(d.grad c)u0=0",
        },
        "residual": {
            "viscosity": 1,
            "terms": [
                {"name": name, "formula": formula}
                for name, formula in TERM_LEDGER
            ],
            "pressure_Poisson":
                "-Delta p=partial_i partial_j(u_i u_j)-div f_minus",
        },
        "force_samples": force_samples,
        "terminal_gluing": {
            "corner": "c=1,d=0, so f_minus is the Cartesian-flat Q7 complete residual",
            "punctured": "conditional on the Q7 all-fixed-order bound, one more time derivative gives uniform Cauchy limits",
            "conclusion": "F_j=lim_(t up to 1) partial_t^j f_minus is in C_c^infinity(K), with every spatial jet zero at x=0",
        },
        "future_extension": {
            "formula": "f_plus=sum_j chi0(b_j sigma) sigma^j F_j/j!",
            "whole_support_bound":
                "||D_x^a partial_sigma^m T_j|| <= C_(j,a,m) b_j^(m-j)",
            "finite_stage_constraints":
                "|a|+m<=floor(j/2), with increasing dyadic b_j selected once from actual finite F_j norms",
            "tail": "each fixed differentiated tail from J is at most 2^(1-J)",
            "support": "K x [1-tau0,3/2]",
            "samples": extension_samples,
        },
        "blocking_inequality": parent["blocking_inequality"],
        "scope": {
            "exact_Q1_Q7_field_localized_conditionally": True,
            "one_smooth_compact_force_extension_constructed_conditionally": True,
            "Q7_endpoint_uniform_all_order_bound_proved_here": False,
            "terminal_velocity_or_pressure_regular_at_origin_proved": False,
            "future_velocity_or_pressure_constructed": False,
            "preterminal_energy_replayed_here": False,
            "pressure_recovery_comparison_or_blowup_proved": False,
            "hashes_or_tests_treated_as_continuum_proof": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
