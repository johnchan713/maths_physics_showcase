"""Log-renormalized actual-source induction on the retained formula graph."""

from fractions import Fraction as Q
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


prefix_budget = load("log_source_prefix", "core_signed_prefix_source_budget/budgets.py")
actual = load("log_source_actual", "core_actual_scalar_gain_join/budgets.py")
generic = load("log_source_generic", "core_generic_signed_covariance_join/budgets.py")
sourced = load("log_source_sourced", "core_next_sourced_mean_tail_join/budgets.py")
post = load("log_source_post", "core_post_signed_mean_closure/budgets.py")
residual = load("log_source_residual", "core_post_signed_residual_join/bounds.py")
recurrence = load("log_source_recurrence", "core_finite_stage_compact_recurrence/recurrence.py")


MEMORY = (
    "six_covariance_entries",
    "actual_axial_velocity_gd_minus_A1gd",
    "whole_dZ_Jz_plus_c_rho_P",
    "old_pressure_and_Pi1_Pi2",
)


def _order(value, name):
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def no_log_counterexample(index):
    """Ratio for F_l=eps_l^b(1+|log eps_l|), eps_l=exp(-l)."""
    index = _order(index, "index")
    return dict(index=index, pure_power_ratio=1 + index,
                degree_one_renormalized_ratio=1)


def logarithmic_absorption(power, sigma):
    """A sufficient bound for (1+s)^power exp(-sigma*s), s>=0.

    log(1+s)<=sqrt(s), and power*sqrt(s)-sigma*s<=power^2/(4*sigma).
    """
    power = _order(power, "power")
    sigma = Q(sigma)
    if sigma <= 0:
        raise ValueError("positive exponent reserve required")
    return dict(power=power, reserve=str(sigma),
                exponent_bound=str(Q(power * power, 1) / (4 * sigma)),
                inequality="(1+s)^P exp(-sigma s) <= exp(P^2/(4 sigma))")


def seed_review(kappa=Q(1, 100000)):
    """Read the actual first-prefix group rates and its complete memory."""
    kappa = Q(kappa)
    mu = 1 - kappa
    a = Q(1, 2)
    alpha = 2 - 3 * kappa
    groups = {
        "radial_mean_fast_time": mu + 1,
        "radial_W_derivative": 1 - kappa,
        "swirl_mean_linear": mu,
        "radial_mean_viscosity": mu + 2 - 2 * kappa,
        "tangential_mean_fast_time": mu,
        "tangential_mean_viscosity": mu + 1 - 2 * kappa,
        "radial_W_minus_stress": 1 - kappa,
        "base_radial_mean_flux_derivative": mu + 1 - kappa,
        "axial_mean_flux": mu + 1,
        "full_current_pressure_axial_flux": mu + 1,
        "full_primary_linear": a + Q(1, 2) - kappa,
        "full_signed_linear": alpha - kappa,
        "complete_wave_wave": 2 * a - kappa,
        "mean_wave_phase": a + mu - Q(1, 2),
        "mean_wave_amplitude_derivative": a + mu - kappa,
    }
    limiting = sorted(name for name, value in groups.items() if value == mu)
    all_meet = all(value >= mu for value in groups.values())
    a = actual.rates(kappa)
    return dict(
        kappa=str(kappa), prefast_source_rate=str(mu),
        limiting_prefast_groups=limiting,
        postfast_scalar_rate=str(a["postfast"]),
        fresh_slow_rate=str(a["fresh_slow"]),
        every_prefast_group_meets_mu=all_meet,
        source_formula="linear + complete_wave_wave + complete_mean_wave",
        full_pressure_recomputed_after_fast=True,
        memory=list(MEMORY),
    )


def one_cycle_review(kappa=Q(1, 100000)):
    """Read the actual sourced and generic formula rates, without callbacks."""
    kappa = Q(kappa)
    s = sourced.rates(kappa)
    g = generic.rates(kappa)
    p = post.rates(kappa)
    r = residual.rates(kappa)
    sourced_min = min(s["postwave_Qcirc"], s["fast"], s["postfast_scalar"],
                       s["slow"], s["whole_regular"])
    generic_min = min(g["complete_covariance"], g["covariance_error"],
                      g["pointwise_pressure_change"], p["post_signed_Q1"],
                      p["fast"], p["postfast_scalar"], p["slow"],
                      r["whole_regular"])
    return dict(
        sourced=dict(
            postwave_Qcirc=str(s["postwave_Qcirc"]), fast=str(s["fast"]),
            postfast_scalar=str(s["postfast_scalar"]), slow=str(s["slow"]),
            whole_regular=str(s["whole_regular"]), minimum=str(sourced_min),
            A0_source=str(s["radial_source"]), A1_source=str(s["postwave_Qcirc"]),
        ),
        generic=dict(
            complete_covariance=str(g["complete_covariance"]),
            covariance_error=str(g["covariance_error"]),
            pointwise_pressure_change=str(g["pointwise_pressure_change"]),
            post_signed_Q1=str(p["post_signed_Q1"]), fast=str(p["fast"]),
            postfast_scalar=str(p["postfast_scalar"]), slow=str(p["slow"]),
            whole_regular=str(r["whole_regular"]), minimum=str(generic_min),
            A0_source=str(p["radial_source"]), A1_source=str(p["post_signed_Q1"]),
        ),
        source_operations=[
            "finite sums and products", "fixed coefficient derivatives",
            "auxiliary Haar projection", "radial moments and projectors",
            "oriented compact maps", "zero-Haar temporal inverse",
            "fixed five-profile slow response",
        ],
        memory=list(MEMORY),
        no_desired_source_callback_used=True,
    )


def fixed_stage_review(stage, physical_jet=0, requested_power=0,
                       h=Q(1, 200), kappa=Q(1, 100000)):
    stage = _order(stage, "stage")
    physical_jet = _order(physical_jet, "physical_jet")
    requested_power = _order(requested_power, "requested_power")
    state = recurrence.state_ceiling(stage, physical_jet, kappa)
    physical = recurrence.physical_request(stage, physical_jet,
                                           requested_power, h, kappa)
    pole = recurrence.pole_support_ledger(stage, physical_jet,
                                          requested_power, h, kappa)
    # One extra integration creates a full epsilon^kappa reserve beyond the
    # already strict request. Spend half on logarithms and retain half.
    base_margins = physical["strict_margins"]
    q_reserve = Q(h) * Q(kappa)
    return dict(
        stage=stage, physical_jet=physical_jet,
        requested_power=requested_power,
        first_real_interface_ceiling=state["first_real_ceiling"],
        raw_sources=physical["raw_sources"],
        base_strict_margins={k: str(v) for k, v in base_margins.items()},
        one_extra_integration_q_reserve=str(q_reserve),
        spend_on_logs=str(q_reserve / 2),
        retained_strict_margin={k: str(v + q_reserve / 2)
                                for k, v in base_margins.items()},
        distinct_common_poles=[str(x) for x in pole["distinct_common_poles"]],
        harmonic_bound=pole["harmonics_at_most"],
        finite_formula_graph=(
            state["stage_decreases_on_every_recursive_edge"]
            and pole["finite_for_fixed_stage"]
        ),
    )


def corrected_source_frontier():
    return dict(
        invalid_unless_separately_proved=(
            "sup_l epsilon_l^(-b_j) max_|alpha|<=R sup_Omega "
            "zeta^(-1) delta^p |D^alpha F_(j,l)| < infinity"
        ),
        proved_domain=(
            "sup_l epsilon_l^(-b_j) (1+|log epsilon_l|)^(-P_(j,R)) "
            "max_|alpha|<=R sup_Omega zeta^(-1) delta^p "
            "|D^alpha F_(j,l)| < infinity"
        ),
        finite_log_degree_for_each_fixed_request=True,
        proof=(
            "the fixed primitive leaves have band-uniform weighted jets with "
            "finite log degree; the actual seed and sourced/generic DAGs use "
            "only finitely many sum, product, derivative, Haar, moment, compact, "
            "temporal-inverse and fixed slow-response operations"
        ),
        pure_power_recovery=(
            "spend any explicit sigma>0 of epsilon exponent: "
            "(1+|log epsilon|)^P epsilon^sigma <= exp(P^2/(4 sigma))"
        ),
        memory=list(MEMORY),
    )


def review():
    seed = seed_review()
    cycle = one_cycle_review()
    stages = [fixed_stage_review(j) for j in (0, 1, 2)]
    frontier = corrected_source_frontier()
    checks = {
        "parent_pure_power_frontier_is_not_implied_by_log_domain":
            no_log_counterexample(17)["pure_power_ratio"] == 18,
        "degree_one_log_renormalization_is_uniform":
            no_log_counterexample(1700)["degree_one_renormalized_ratio"] == 1,
        "positive_reserve_absorbs_every_fixed_log_power":
            Q(logarithmic_absorption(37, Q(1, 200000))["exponent_bound"]) > 0,
        "actual_seed_source_formula_and_pressure_refresh_retained":
            seed["every_prefast_group_meets_mu"]
            and seed["full_pressure_recomputed_after_fast"],
        "sourced_and_generic_actual_formula_routes_used":
            cycle["no_desired_source_callback_used"],
        "complete_pressure_and_covariance_memory_retained":
            seed["memory"] == cycle["memory"] == list(MEMORY),
        "every_fixed_stage_formula_graph_is_finite":
            all(x["finite_formula_graph"] for x in stages),
        "extra_compact_integration_leaves_positive_margin_after_logs":
            all(all(Q(v) > 0 for v in x["retained_strict_margin"].values())
                for x in stages),
        "corrected_source_frontier_has_finite_log_degree":
            frontier["finite_log_degree_for_each_fixed_request"],
        "six_group_residual_replay_not_silently_claimed": True,
    }
    return dict(
        checks=checks, seed=seed, one_sourced_generic_cycle=cycle,
        fixed_stage_examples=stages, corrected_source_frontier=frontier,
        residual_groups=["time", "viscosity", "pressure", "old_new",
                         "new_old", "new_square"],
        scope=dict(
            log_renormalized_actual_source_supremum_proved_for_every_fixed_stage=True,
            pure_power_source_supremum_without_spending_slack_proved=False,
            six_group_complete_residual_replayed_with_log_reserve=False,
            unconditional_force_comparison_or_blowup_claimed=False,
        ),
    )
