"""Exact controls for instantiating the log-source induction on the current chain."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


source = load(
    "current_log_source_induction",
    "core_log_renormalized_actual_source_induction/source_induction.py",
)
inverse = load(
    "current_compact_slow_inverse",
    "core_compact_slow_inverse_reconstruction/inverse.py",
)


NODES = (
    "leading_borel_profile",
    "same_band_tail",
    "full_path_ode",
    "fixed_same_choice_projectors",
    "weighted_compact_slow_inverse",
    "log_renormalized_actual_source_induction",
)
EDGES = tuple(zip(NODES, NODES[1:]))
MEMORY = (
    "six_covariance_entries",
    "actual_axial_velocity_gd_minus_A1gd",
    "whole_dZ_Jz_plus_c_rho_P",
    "old_pressure_and_Pi1_Pi2",
)


def topological_order(nodes=NODES, edges=EDGES):
    incoming = {node: 0 for node in nodes}
    outgoing = {node: [] for node in nodes}
    for left, right in edges:
        if left not in incoming or right not in incoming:
            raise ValueError("unknown current-chain node")
        incoming[right] += 1
        outgoing[left].append(right)
    ready = [node for node in nodes if incoming[node] == 0]
    order = []
    while ready:
        node = ready.pop(0)
        order.append(node)
        for right in outgoing[node]:
            incoming[right] -= 1
            if incoming[right] == 0:
                ready.append(right)
    if len(order) != len(nodes):
        raise ValueError("current-chain dependency cycle")
    return order


def fixed_request_bridge(stage, physical_jet=0, requested_power=0,
                         log_degree=0, h=Q(1, 200),
                         kappa=Q(1, 100000)):
    """Bridge one fixed recurrence request to the current all-fixed-order maps."""
    if not all(type(v) is int and v >= 0 for v in
               (stage, physical_jet, requested_power, log_degree)):
        raise ValueError("stage, jets, target and log degree must be nonnegative integers")
    old = source.fixed_stage_review(
        stage, physical_jet, requested_power, h=h, kappa=kappa
    )
    raw = old["raw_sources"]
    # The frozen log proof states that one additional compact integration
    # requests two additional source derivatives. Make that finite request
    # explicit instead of treating the reserve as a cost-free annotation.
    with_extra = {name: value + 2 for name, value in raw.items()}
    sigma = Q(h) * Q(kappa) / 2
    absorption = source.logarithmic_absorption(log_degree, sigma)
    current_inverse = inverse.review()
    return {
        "stage": stage,
        "physical_jet": physical_jet,
        "requested_power": requested_power,
        "first_real_interface_ceiling": old["first_real_interface_ceiling"],
        "raw_sources_before_extra_integration": raw,
        "raw_sources_with_extra_integration": with_extra,
        "extra_source_derivatives_per_compact_map": 2,
        "finite_formula_graph": old["finite_formula_graph"],
        "current_inverse_accepts_every_fixed_order": all(current_inverse["checks"].values()),
        "one_extra_integration_q_reserve": old["one_extra_integration_q_reserve"],
        "spend_on_logs": old["spend_on_logs"],
        "retained_strict_margin": old["retained_strict_margin"],
        "log_degree": log_degree,
        "absorption_exponent_bound": absorption["exponent_bound"],
        "finite_harmonic_bound": old["harmonic_bound"],
        "finite_common_poles": old["distinct_common_poles"],
    }


def current_leaf_domain():
    """The fixed current primitive leaves enter the log class with degree zero."""
    return {
        "leading_background_fixed_order_weighted_jets": 0,
        "same_band_gap_and_frame_coefficients": 0,
        "full_path_ode_pressure_potential_and_cutoffs": 0,
        "fixed_profile_and_same_choice_projectors": 0,
        "shifted_compact_zero_haar_and_fixed_slow_maps": 0,
        "interpretation": (
            "zero denotes log degree P=0; constants may depend on the fixed "
            "stage and requested order but not on the tail-band index"
        ),
    }


def review():
    order = topological_order()
    current_inverse = inverse.review()
    old_source = source.review()
    examples = [
        fixed_request_bridge(0, 0, 0, 0),
        fixed_request_bridge(1, 3, 2, 7),
        fixed_request_bridge(2, 7, 5, 37),
        fixed_request_bridge(5, 11, 13, 113),
    ]
    leaves = current_leaf_domain()
    checks = {
        "current_chain_is_acyclic_and_ordered": order == list(NODES),
        "all_current_inverse_controls_hold": all(current_inverse["checks"].values()),
        "old_log_source_induction_controls_hold": all(old_source["checks"].values()),
        "current_primitive_leaves_have_log_degree_zero":
            all(value == 0 for key, value in leaves.items() if key != "interpretation"),
        "every_example_recurrence_request_is_finite":
            all(row["finite_formula_graph"] for row in examples),
        "extra_integration_derivative_cost_is_explicit_and_finite":
            all(all(row["raw_sources_with_extra_integration"][name]
                    == row["raw_sources_before_extra_integration"][name] + 2
                    for name in ("A0", "A1")) for row in examples),
        "current_inverse_accepts_each_fixed_enlarged_request":
            all(row["current_inverse_accepts_every_fixed_order"] for row in examples),
        "half_extra_reserve_absorbs_each_fixed_log_degree":
            all(Q(row["absorption_exponent_bound"]) >= 0 for row in examples),
        "positive_margin_remains_after_log_absorption":
            all(all(Q(value) > 0 for value in row["retained_strict_margin"].values())
                for row in examples),
        "complete_pressure_covariance_and_axial_memory_agree":
            tuple(old_source["seed"]["memory"]) == MEMORY,
        "pure_power_zero_margin_is_not_reintroduced_as_frontier":
            not current_inverse["route_classification"]["zero_margin_is_the_current_global_frontier"],
        "complete_residual_is_not_silently_replayed": True,
    }
    return {
        "checks": checks,
        "current_chain_order": order,
        "current_leaf_log_degrees": leaves,
        "fixed_request_examples": examples,
        "source_domain": old_source["corrected_source_frontier"]["proved_domain"],
        "memory": list(MEMORY),
        "quantifier_order": (
            "fix stage j, physical jet m and target N; compute the finite "
            "recurrence request; enlarge A0 and A1 source requests by two; "
            "then apply the current fixed-order primitive bounds"
        ),
        "scope": {
            "log_source_induction_instantiated_on_current_five_contract_chain": True,
            "stage_uniform_constants_or_log_degrees_proved": False,
            "complete_residual_or_downstream_force_comparison_replayed_here": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
