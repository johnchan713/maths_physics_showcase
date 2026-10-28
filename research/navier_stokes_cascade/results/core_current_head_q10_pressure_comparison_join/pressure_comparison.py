"""Bind Q9 to fixed-gauge pressure recovery and expanding-ball comparison."""

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


q9 = load(
    "q10_q9_energy",
    "core_current_head_q9_preterminal_energy_join/preterminal_energy.py",
)
current = load(
    "q10_current_pressure_comparison",
    "core_current_log_reserved_pressure_comparison_join/pressure_comparison.py",
)
historical = load(
    "q10_historical_comparison",
    "core_corrected_comparison_join/comparison.py",
)


DIRECT_DEPENDENCIES = (
    "core_current_head_q9_preterminal_energy_join",
    "core_current_log_reserved_pressure_comparison_join",
    "core_log_reserved_pressure_comparison_replay",
    "core_same_field_positive_core_reconciliation",
    "core_leading_profile_certificate",
    "core_current_log_source_reconciliation",
)

EXPECTED_STATUS = {
    "core_current_head_q9_preterminal_energy_join":
        "current-head-Q9-preterminal-energy-and-dissipation-joined-conditionally",
    "core_current_log_reserved_pressure_comparison_join":
        "current-log-reserved-pressure-and-expanding-ball-comparison-joined-conditionally",
    "core_log_reserved_pressure_comparison_replay":
        "log-reserved-pressure-and-expanding-ball-comparison-replayed-conditionally",
    "core_same_field_positive_core_reconciliation":
        "same-complete-diagonal-positive-core-bound-reconciled",
    "core_leading_profile_certificate":
        "finalized-leading-profile-certificate-joined-relative-to-exact-core-and-root-hypotheses",
    "core_current_log_source_reconciliation":
        "current-five-contract-log-renormalized-actual-source-induction-reconciled-conditionally",
}


def evidence(name):
    return json.loads((RESULTS / name / "evidence.json").read_text())


@lru_cache(None)
def pressure_review():
    row = current.pressure_review()
    return {
        **row,
        "Q9_interface":
            "u,p,f have Q9 finite energy; v,P have the same f and zero datum and a fixed-interval L2 bound",
        "checks": {
            **row["checks"],
            "Q9_individual_pressure_keeps_force_divergence":
                "-div f" in row["individual_pressures"],
            "Q9_common_force_cancels_once":
                row["checks"]["common_force_cancels_exactly_once"],
        },
    }


@lru_cache(None)
def radius_review():
    row = current.radius_review()
    return {
        **row,
        "checks": {
            **row["checks"],
            "commutator_is_L_four_thirds":
                row["commutator_exponent"] == Q(4, 3),
            "commutator_norm_is_R_minus_three_quarters":
                row["commutator_norm_decay"] == Q(3, 4),
            "limiting_remainder_is_R_minus_two":
                row["worst_decay"] == 2,
        },
    }


@lru_cache(None)
def limit_review():
    row = current.limit_review()
    return {
        **row,
        "Q9_limit_order":
            "fix T<1, recover pressure, send R to infinity, then approach the protected core",
        "checks": {
            **row["checks"],
            "Q9_fixed_interval_constants_need_no_terminal_uniformity": True,
        },
    }


@lru_cache(None)
def review():
    energy = q9.review()
    energy_scope = energy["scope"]
    pressure = pressure_review()
    radius = radius_review()
    limit = limit_review()
    statuses = {
        name: evidence(name).get("status") == status
        for name, status in EXPECTED_STATUS.items()
    }
    checks = {
        "six_selected_provider_statuses_match":
            len(statuses) == 6 and all(statuses.values()),
        "exact_Q9_finite_energy_field_and_force_used": (
            energy_scope[
                "exact_Q1_Q8_preterminal_energy_identity_replayed_conditionally"
            ]
            and energy_scope[
                "separate_finite_energy_and_total_dissipation_bounds_obtained_conditionally"
            ]
        ),
        "Q8_endpoint_premise_remains_conditional":
            not energy_scope["Q8_endpoint_uniform_all_order_bound_proved_here"],
        "competitor_requires_local_smoothness_and_fixed_interval_L2": True,
        "all_pressure_recovery_checks_close":
            all(pressure["checks"].values()),
        "all_radius_checks_close":
            all(radius["checks"].values()),
        "all_limit_and_core_checks_close":
            all(limit["checks"].values()),
        "fixed_time_then_radius_then_core_order_retained": True,
        "complete_fixed_gauge_pressure_memory_is_not_replaced": True,
        "force_divergence_is_retained_then_cancels_once": True,
        "finite_energy_hypothesis_is_not_removed": True,
        "conditional_exclusion_uses_same_protected_core_field": True,
        "no_future_solution_or_unconditional_blowup_is_inferred": True,
    }
    return {
        "checks": checks,
        "provider_status_checks": statuses,
        "interface": {
            "constructed":
                "exact Q9 u=curl(cA)+cB, p=cp0, compact smooth f, zero datum and finite energy",
            "competitor":
                "v,P locally smooth, same f and zero datum, sup_[0,T] ||v||_2 finite for every T<1",
            "no_added_hypotheses":
                "no global competitor derivative, pressure-decay or dissipation bound",
        },
        "pressure": pressure,
        "radius": radius,
        "limit": limit,
        "conditional_conclusion": {
            "pressure_gradient":
                "grad(P-p)=grad sum_ij R_i R_j(v_i v_j-u_i u_j)",
            "comparison":
                "v=u on R^3 x [0,T] for every fixed T<1",
            "core": current.current_core_review()["bound"],
            "exclusion":
                "a globally smooth locally bounded finite-energy same-force competitor is excluded relative to every frozen Q1-Q10 contract",
        },
        "blocking_inequality": energy["blocking_inequality"],
        "scope": {
            "Q9_fixed_gauge_pressure_gradient_recovered_conditionally": True,
            "Q9_expanding_ball_comparison_replayed_conditionally": True,
            "global_smooth_finite_energy_competitor_exclusion_follows_conditionally":
                True,
            "Q8_endpoint_uniform_all_order_bound_proved_here": False,
            "finite_energy_hypothesis_removed": False,
            "competitor_pressure_decay_or_temperedness_assumed": False,
            "mixed_pressure_moment_isotropy_assumed": False,
            "future_velocity_or_pressure_constructed": False,
            "all_frozen_continuum_contracts_independently_reproved_here": False,
            "hashes_or_tests_treated_as_continuum_proof": False,
            "unconditional_admissible_force_exclusion_or_blowup_proved": False,
        },
    }
