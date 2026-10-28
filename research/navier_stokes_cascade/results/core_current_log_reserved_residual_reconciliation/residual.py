"""Instantiate the three-account six-group residual replay on the current chain."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


current_source = load(
    "current_residual_source_reconciliation",
    "core_current_log_source_reconciliation/reconcile.py",
)
frozen_residual = load(
    "current_residual_frozen_replay",
    "core_log_reserved_complete_residual_replay/replay.py",
)


RESIDUAL_GROUPS = frozen_residual.RESIDUAL_GROUPS
MEMORY = current_source.MEMORY


def order(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


@lru_cache(None)
def current_stage_request(stage, physical_jet, requested_power,
                          h=Q(1, 200), kappa=Q(1, 100000)):
    """Compute the residual-stage current request without materializing pole lists."""
    stage = order(stage, "stage")
    physical_jet = order(physical_jet, "physical_jet")
    requested_power = order(requested_power, "requested_power")
    recurrence = current_source.source.recurrence
    state = recurrence.state_ceiling(stage, physical_jet, kappa)
    physical = recurrence.physical_request(
        stage, physical_jet, requested_power, h, kappa
    )
    raw = physical["raw_sources"]
    return {
        "stage": stage,
        "physical_jet": physical_jet,
        "requested_power": requested_power,
        "first_real_interface_ceiling": state["first_real_ceiling"],
        "raw_sources_before_extra_integration": raw,
        "raw_sources_with_extra_integration": {
            name: value + 2 for name, value in raw.items()
        },
        "extra_source_derivatives_per_compact_map": 2,
        "finite_formula_graph": state["stage_decreases_on_every_recursive_edge"],
        "base_strict_margins": {
            name: str(value) for name, value in physical["strict_margins"].items()
        },
    }


def reserve_accounts(row, compact_log_degrees=(0, 0),
                     regular_log_degree=0, cutoff_log_degree=0):
    if (not isinstance(compact_log_degrees, tuple)
            or len(compact_log_degrees) != 2):
        raise ValueError("compact_log_degrees must be an (A0,A1) tuple")
    degrees = {
        "A0": order(compact_log_degrees[0], "A0 log degree"),
        "A1": order(compact_log_degrees[1], "A1 log degree"),
        "regular": order(regular_log_degree, "regular log degree"),
        "cutoff": order(cutoff_log_degree, "cutoff log degree"),
    }
    compact = row["compact_accounts"]
    result = {}
    for name in ("A0", "A1"):
        sigma = compact[name]["spent_on_own_source_log"]
        result[name] = {
            "owner": name + "_actual_source",
            "reserve_class": "compact_source",
            "reserve": str(sigma),
            "log_degree": degrees[name],
            "absorption_exponent_bound": str(
                frozen_residual.logarithmic_absorption(degrees[name], sigma)
            ),
            "may_pay_other_account": False,
        }
    for name, sigma in (
        ("regular", row["regular_log_spend"]),
        ("cutoff", row["cutoff_log_spend"]),
    ):
        result[name] = {
            "owner": name + "_residual",
            "reserve_class": (
                "regular_residual" if name == "regular" else "cutoff_tail"
            ),
            "reserve": str(sigma),
            "log_degree": degrees[name],
            "absorption_exponent_bound": str(
                frozen_residual.logarithmic_absorption(degrees[name], sigma)
            ),
            "may_pay_other_account": False,
        }
    return result


def comparison(m=0, target=0, compact_log_degrees=(0, 0),
               regular_log_degree=0, cutoff_log_degree=0,
               h=Q(1, 200), kappa=Q(1, 100000)):
    m, target = order(m, "m"), order(target, "target")
    row = frozen_residual.comparison(m, target, h, kappa)
    stage = row["log_reserved_comparison_stage"]
    request = current_stage_request(stage, m, target, h, kappa)
    accounts = reserve_accounts(
        row, compact_log_degrees, regular_log_degree, cutoff_log_degree
    )
    groups = row["residual_group_exponents"]
    checks = {
        "frozen_residual_controls_hold": all(row["checks"].values()),
        "current_source_reconciliation_controls_hold":
            all(current_source.review()["checks"].values()),
        "same_log_reserved_stage_used_for_source_and_residual":
            request["stage"] == stage,
        "extra_compact_derivative_cost_carried_into_residual_stage":
            all(request["raw_sources_with_extra_integration"][name]
                == request["raw_sources_before_extra_integration"][name] + 2
                for name in ("A0", "A1")),
        "all_three_reserve_classes_are_disjoint":
            len({entry["owner"] for entry in accounts.values()}) == 4
            and {entry["reserve_class"] for entry in accounts.values()}
            == {"compact_source", "regular_residual", "cutoff_tail"}
            and all(not entry["may_pay_other_account"] for entry in accounts.values()),
        "A0_and_A1_each_retain_positive_margin":
            all(row["compact_accounts"][name]["retained_strict_margin"] > 0
                for name in ("A0", "A1")),
        "regular_stage_retains_positive_margin":
            row["regular_retained_margin"] > 0,
        "cutoff_tail_spends_only_its_half_gain":
            row["cutoff_log_spend"] * 2
            == row["first_omitted_representative_gain"],
        "all_six_residual_groups_retained": tuple(groups) == RESIDUAL_GROUPS,
        "both_cross_orientations_have_the_same_proved_exponent":
            groups["old_new"] == groups["new_old"],
        "all_six_groups_and_fixed_state_beat_target":
            row["fixed_state_exponent_after_logs"] >= target + 1
            and all(value >= target + 1 for value in groups.values()),
        "complete_pressure_memory_retained":
            tuple(current_source.review()["memory"]) == MEMORY,
    }
    return {
        "jet": m,
        "target": target,
        "base_comparison_stage": row["base_comparison_stage"],
        "log_reserved_comparison_stage": stage,
        "current_stage_request": request,
        "reserve_accounts": accounts,
        "regular_stage_gain": str(row["regular_stage_gain"]),
        "regular_retained_margin": str(row["regular_retained_margin"]),
        "cutoff_retained_tail_exponent": str(row["tail_exponent_after_cutoff_logs"]),
        "fixed_state_exponent_after_logs": str(row["fixed_state_exponent_after_logs"]),
        "residual_group_exponents": {name: str(value) for name, value in groups.items()},
        "complete_minimum_exponent": str(row["complete_minimum_exponent"]),
        "residual_identity": row["identity"],
        "checks": checks,
    }


def review():
    samples = [
        comparison(0, 0, (0, 0), 0, 0),
        comparison(1, 3, (7, 11), 13, 17),
        comparison(5, 11, (37, 41), 43, 47),
    ]
    checks = {
        "current_source_precedes_current_residual":
            all(row["checks"]["current_source_reconciliation_controls_hold"]
                for row in samples),
        "frozen_three_account_residual_arithmetic_retained":
            all(row["checks"]["frozen_residual_controls_hold"] for row in samples),
        "extra_derivative_cost_present_at_each_residual_stage":
            all(row["checks"]["extra_compact_derivative_cost_carried_into_residual_stage"]
                for row in samples),
        "compact_regular_and_cutoff_reserves_not_reused":
            all(row["checks"]["all_three_reserve_classes_are_disjoint"]
                for row in samples),
        "all_six_groups_and_both_cross_orientations_retained":
            all(row["checks"]["all_six_residual_groups_retained"]
                and row["checks"]["both_cross_orientations_have_the_same_proved_exponent"]
                for row in samples),
        "sample_fixed_jet_targets_close_with_positive_reserve":
            all(all(row["checks"].values()) for row in samples),
        "same_current_pressure_covariance_and_axial_memory_retained":
            all(row["checks"]["complete_pressure_memory_retained"] for row in samples),
        "corner_force_energy_comparison_not_silently_replayed": True,
    }
    return {
        "checks": checks,
        "sample_comparisons": samples,
        "residual_groups": list(RESIDUAL_GROUPS),
        "memory": list(MEMORY),
        "scope": {
            "six_group_log_reserved_local_residual_instantiated_on_current_chain": True,
            "same_current_locally_finite_field_and_cutoff_sequence_retained": True,
            "stage_uniform_constants_or_fixed_q_contraction_used": False,
            "corner_endpoint_force_energy_comparison_replayed_here": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
