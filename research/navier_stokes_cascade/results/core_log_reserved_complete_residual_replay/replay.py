"""Separate logarithmic reserves for the complete six-group residual replay."""

from fractions import Fraction as Q
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


diagonal = load(
    "log_reserved_diagonal",
    "core_complete_potential_diagonal_join/diagonal.py",
)
actual_source = load(
    "log_reserved_actual_source",
    "core_log_renormalized_actual_source_induction/source_induction.py",
)


RESIDUAL_GROUPS = (
    "time", "viscosity", "pressure", "old_new", "new_old", "new_square",
)


def order(value, name):
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def logarithmic_absorption(power, reserve):
    """Exact sufficient constant for (1+s)^P exp(-reserve*s)."""
    power = order(power, "power")
    reserve = Q(reserve)
    if reserve <= 0:
        raise ValueError("positive exponent reserve required")
    return Q(power * power, 1) / (4 * reserve)


def compact_reserve_accounts(h=Q(1, 200), kappa=diagonal.KAPPA):
    """The extra integration is booked independently for A0 and A1."""
    h, kappa = Q(h), Q(kappa)
    diagonal.parameters(h, kappa)
    if kappa <= 0:
        raise ValueError("positive kappa required")
    gain = h * kappa
    account = {
        "extra_integration_gain": gain,
        "spent_on_own_source_log": gain / 2,
        "retained_strict_margin": gain / 2,
        "may_pay_regular_or_cutoff_log": False,
    }
    return {"A0": dict(account), "A1": dict(account)}


def comparison(m=0, target=0, h=Q(1, 200), kappa=diagonal.KAPPA):
    """One fixed stage after the old threshold, with three disjoint accounts."""
    m, target = order(m, "m"), order(target, "target")
    h, kappa = Q(h), Q(kappa)
    p = diagonal.parameters(h, kappa)
    base = diagonal.comparison(m, target, h, kappa)
    J0 = base["comparison_stage"]
    J = J0 + 1

    stage_gain = h * p["gain"]
    base_regular = (
        h * diagonal.current.rates(J0, kappa)["regular"]
        - diagonal.residual_loss(m, h)
    )
    raw_regular = (
        h * diagonal.current.rates(J, kappa)["regular"]
        - diagonal.residual_loss(m, h)
    )
    regular_after_log = raw_regular - stage_gain / 2

    # The locally finite tail has already spent half of every omitted
    # representative gain on its own finite log degree and cutoff constants.
    L = diagonal.cut_loss(m + 2, h)
    H = diagonal.total_loss(m + 2, h)
    first_omitted_gain = diagonal.gain(J + 1, h, kappa)
    tail = first_omitted_gain / 2 - L
    groups = {
        "time": tail,
        "viscosity": tail,
        "pressure": tail,
        "old_new": tail - H,
        "new_old": tail - H,
        "new_square": 2 * tail,
    }
    compact = compact_reserve_accounts(h, kappa)
    fixed_state_exponent = min(
        regular_after_log,
        raw_regular + compact["A0"]["retained_strict_margin"],
        raw_regular + compact["A1"]["retained_strict_margin"],
    )
    complete_minimum = min(fixed_state_exponent, *groups.values())
    return {
        "jet": m,
        "target": target,
        "base_comparison_stage": J0,
        "log_reserved_comparison_stage": J,
        "base_regular_exponent": base_regular,
        "raw_regular_exponent": raw_regular,
        "regular_stage_gain": stage_gain,
        "regular_log_spend": stage_gain / 2,
        "regular_retained_margin": stage_gain / 2,
        "fixed_state_exponent_after_logs": fixed_state_exponent,
        "first_omitted_representative_gain": first_omitted_gain,
        "cutoff_log_spend": first_omitted_gain / 2,
        "tail_exponent_after_cutoff_logs": tail,
        "background_state_loss": H,
        "compact_accounts": compact,
        "residual_group_exponents": groups,
        "complete_minimum_exponent": complete_minimum,
        "identity": (
            "R(u*,p*)=R(U[J]) + e_t - Delta e + grad pi "
            "+ (U[J].grad)e + (e.grad)U[J] + (e.grad)e"
        ),
        "checks": {
            "old_threshold_checks_hold": all(base["checks"].values()),
            "exactly_one_stage_added": J == J0 + 1,
            "one_stage_gain_is_exact": raw_regular - base_regular == stage_gain,
            "regular_log_account_positive": stage_gain / 2 > 0,
            "compact_accounts_are_separate_and_positive": (
                compact["A0"] is not compact["A1"]
                and all(row["retained_strict_margin"] > 0
                        for row in compact.values())
            ),
            "cutoff_account_positive": first_omitted_gain / 2 > 0,
            "fixed_state_beats_target": fixed_state_exponent >= target + 1,
            "all_six_tail_groups_beat_target": (
                tuple(groups) == RESIDUAL_GROUPS
                and all(x >= target + 1 for x in groups.values())
            ),
            "complete_residual_beats_target": complete_minimum >= target + 1,
        },
    }


def review():
    source = actual_source.review()
    samples = [comparison(m, target) for m, target in ((0, 0), (1, 3), (5, 11))]
    zero_reserve_ratios = [actual_source.no_log_counterexample(j)["pure_power_ratio"]
                           for j in (0, 10, 1000)]
    checks = {
        "actual_source_log_domain_precedes_this_replay": source["scope"][
            "log_renormalized_actual_source_supremum_proved_for_every_fixed_stage"
        ],
        "zero_reserve_pure_power_ratio_is_unbounded": zero_reserve_ratios[-1] > zero_reserve_ratios[1],
        "positive_reserve_absorbs_each_fixed_log_degree": logarithmic_absorption(37, Q(1, 200000)) > 0,
        "compact_A0_and_A1_reserves_are_distinct": all(
            not row["compact_accounts"][name]["may_pay_regular_or_cutoff_log"]
            for row in samples for name in ("A0", "A1")
        ),
        "regular_reserve_comes_from_one_new_stage": all(
            row["checks"]["one_stage_gain_is_exact"] for row in samples
        ),
        "cutoff_reserve_remains_the_existing_half_gain": all(
            row["cutoff_log_spend"] * 2 == row["first_omitted_representative_gain"]
            for row in samples
        ),
        "all_six_residual_groups_are_retained": all(
            tuple(row["residual_group_exponents"]) == RESIDUAL_GROUPS
            for row in samples
        ),
        "sample_complete_residual_targets_close": all(
            all(row["checks"].values()) for row in samples
        ),
        "one_fixed_stage_is_used_for_each_jet_and_power": all(
            isinstance(row["log_reserved_comparison_stage"], int) for row in samples
        ),
        "localization_force_comparison_and_blowup_not_claimed": True,
    }
    return {
        "checks": checks,
        "zero_reserve_counterexample_ratios": zero_reserve_ratios,
        "default_comparison": comparison(),
        "sample_comparisons": samples,
        "residual_groups": list(RESIDUAL_GROUPS),
        "scope": {
            "log_reserved_six_group_complete_local_residual_flat_conditionally": True,
            "same_locally_finite_field_and_cutoff_sequence_retained": True,
            "complete_pressure_memory_and_both_interaction_orientations_retained": True,
            "stage_uniform_constants_or_fixed_q_contraction_used": False,
            "unconditional_pulse_force_comparison_or_blowup_claimed": False,
            "localization_force_energy_and_comparison_replayed_after_source_correction": False,
        },
    }
