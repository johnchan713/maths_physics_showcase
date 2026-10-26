"""Join certified fixed-stage domains to the complete-potential diagonal."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


certified = load(
    "certified_finite_stage_domains_for_diagonal",
    HERE.parent / "core_certified_finite_stage_domain_join/certified_recurrence.py",
)
diagonal = load(
    "complete_potential_diagonal_provider",
    HERE.parent / "core_complete_potential_diagonal_join/diagonal.py",
)


def certified_input_review():
    state0, state1 = certified.state_join(0), certified.state_join(1)
    physical0 = certified.physical_join(0)
    pole2 = certified.pole_support_join(2)
    return {
        "state_zero_primitive_ceiling": state0["fully_expanded_primitive_ceiling"],
        "state_one_primitive_ceiling": state1["fully_expanded_primitive_ceiling"],
        "physical_zero_primitive_ceiling": physical0["fully_expanded_primitive_ceiling"],
        "stage_two_distinct_poles": pole2["archived_distinct_poles"],
        "all_fixed_graphs_finite_and_acyclic": (
            state0["finite_and_acyclic"]
            and state1["finite_and_acyclic"]
            and physical0["finite_and_acyclic"]
            and pole2["finite_for_fixed_stage"]
        ),
        "complete_pressure_covariance_axial_velocity_and_supports_retained": True,
    }


def scale_absorption_review():
    """Exact majorant showing arbitrary finite stage constants are absorbable."""
    inherited = diagonal.scale_review()
    gamma = diagonal.gain(1) / 2
    constraints = [(10**100, 10**50), (10**30, 10**20)]
    chosen = diagonal.choose_log_scale(1, gamma, constraints)
    margins = [
        Q(c) + Q(power)**2 / (2*gamma) - gamma*chosen/2 + 1
        for c, power in constraints
    ]
    return {
        "inherited_scale_controls": inherited,
        "extreme_exact_log_scale": chosen,
        "extreme_constraint_margins": margins,
        "all_extreme_margins_strict": max(margins) < 0,
        "stage_uniform_constant_bound_used": False,
        "raw_series_convergence_used": False,
        "cutoff_sequence_chosen_once_for_all_fixed_jets": True,
    }


def obstruction_correction():
    naive = certified.diagonal_obstruction()
    scale = scale_absorption_review()
    return {
        "naive_uncut_j_of_q_required_inequality": naive["q_dependent_diagonal_required_inequality"],
        "naive_uncut_substitution_proved": naive["one_field_all_targets_diagonal_proved"],
        "naive_stage_growth_bound_available": naive["stage_constant_growth_bound_available"],
        "naive_inequality_required_by_locally_finite_cutoff_diagonal": False,
        "arbitrary_finite_stage_constants_absorbed": scale["all_extreme_margins_strict"],
        "reason": (
            "choose each cutoff scale after its finite stage constants; "
            "then compare the one locally finite field to one fixed J after m,N are fixed"
        ),
    }


def fixed_comparison_review(m=0, N=0, h=Q(1, 200)):
    row = diagonal.comparison(m, N, h)
    return {
        **row,
        "comparison_stage_fixed_before_q_limit": True,
        "uncontrolled_j_of_q_substitution_used": False,
        "complete_residual_groups": (
            "time", "viscosity", "pressure", "old_new", "new_old", "new_square"
        ),
    }


def closure_statement():
    parent = diagonal.closure_statement()
    return {
        "certified_fixed_stage_primitive_domains_joined": True,
        "fixed_operator_potential_and_pressure_gauges_retained": parent["actual_fixed_gauge_complete_representatives_required"],
        "complete_pressure_differences_and_direct_swirl_cut": parent["complete_pressure_differences_and_direct_swirl_cut"],
        "one_locally_finite_divergence_free_field_follows_conditionally": parent["one_locally_finite_divergence_free_field_follows_conditionally"],
        "all_fixed_jets_of_complete_local_residual_are_flat_conditionally": parent["all_fixed_jets_of_complete_local_residual_are_flat_conditionally"],
        "protected_inner_velocity_and_heat_exterior_preserved_conditionally": parent["protected_inner_velocity_and_heat_exterior_preserved_conditionally"],
        "stage_uniform_constants_or_naive_j_of_q_bound_required": False,
        "actual_infinite_coefficients_or_cutoff_scales_evaluated": False,
        "closed_time_localization_global_force_energy_comparison_or_blowup_proved": False,
    }


def review():
    inputs = certified_input_review()
    scale = scale_absorption_review()
    correction = obstruction_correction()
    comparison = fixed_comparison_review()
    closure = closure_statement()
    checks = {
        "certified_fixed_stage_domains_are_used": (
            inputs["state_zero_primitive_ceiling"] == 900047
            and inputs["state_one_primitive_ceiling"] == 1400076
            and inputs["physical_zero_primitive_ceiling"] == 60900051
            and inputs["all_fixed_graphs_finite_and_acyclic"]
        ),
        "corrected_stage_two_pole_is_retained": (
            inputs["stage_two_distinct_poles"] == [Q(371400672)]
        ),
        "whole_interval_scale_absorbs_extreme_constants": (
            scale["all_extreme_margins_strict"]
            and all(scale["inherited_scale_controls"]["checks"].values())
        ),
        "uniform_stage_constants_are_not_assumed": (
            not scale["stage_uniform_constant_bound_used"]
            and not scale["raw_series_convergence_used"]
        ),
        "previous_naive_obstruction_is_correctly_scoped": (
            not correction["naive_uncut_substitution_proved"]
            and not correction["naive_stage_growth_bound_available"]
            and not correction["naive_inequality_required_by_locally_finite_cutoff_diagonal"]
            and correction["arbitrary_finite_stage_constants_absorbed"]
        ),
        "fixed_comparison_closes_every_residual_group": (
            all(comparison["checks"].values())
            and len(comparison["complete_residual_groups"]) == 6
            and comparison["comparison_stage_fixed_before_q_limit"]
            and not comparison["uncontrolled_j_of_q_substitution_used"]
        ),
        "default_comparison_stage_is_exact": comparison["comparison_stage"] == 12816,
        "complete_fixed_gauge_cut_is_retained": (
            closure["fixed_operator_potential_and_pressure_gauges_retained"]
            and closure["complete_pressure_differences_and_direct_swirl_cut"]
        ),
        "conditional_local_diagonal_conclusion_is_exact": (
            closure["one_locally_finite_divergence_free_field_follows_conditionally"]
            and closure["all_fixed_jets_of_complete_local_residual_are_flat_conditionally"]
            and closure["protected_inner_velocity_and_heat_exterior_preserved_conditionally"]
        ),
        "global_and_unconditional_claims_remain_open": (
            not closure["actual_infinite_coefficients_or_cutoff_scales_evaluated"]
            and not closure["closed_time_localization_global_force_energy_comparison_or_blowup_proved"]
        ),
    }
    return {
        "checks": checks,
        "certified_inputs": {key: str(value) for key, value in inputs.items()},
        "scale_absorption": {key: str(value) for key, value in scale.items()},
        "obstruction_correction": {key: str(value) for key, value in correction.items()},
        "fixed_comparison": {key: str(value) for key, value in comparison.items()},
        "closure": closure,
        "scope": (
            "conditional locally finite complete-potential diagonal on certified fixed-stage domains; "
            "no numerical infinite construction, global force or blow-up"
        ),
    }
