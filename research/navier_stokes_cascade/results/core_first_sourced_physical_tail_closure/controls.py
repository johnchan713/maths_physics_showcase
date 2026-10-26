from fractions import Fraction as Q
from closure import (DEPENDENCIES, PRIMITIVES, fixed_order_induction,
                     normalized_request, physical_request,
                     residual_interface, topological_order)


def review():
    norm, phys = normalized_request(), physical_request()
    interface = residual_interface()
    position = {name: i for i, name in enumerate(topological_order())}
    checks = {
        "changed_state_graph_is_acyclic": all(
            dep in PRIMITIVES or position[dep] < position[name]
            for name, deps in DEPENDENCIES.items() for dep in deps),
        "A1_precedes_fresh_fast_and_A0_is_last": (
            position["postwave_A1_response"] < position["fresh_fast_means"]
            and topological_order()[-1] == "A0_Qcirc_gr_S"),
        "desired_A0_is_not_an_earlier_input": (
            not fixed_order_induction()["desired_A0_output_used_earlier"]
            and all("A0_Qcirc_gr_S" not in deps for name, deps in DEPENDENCIES.items()
                    if name != "A0_Qcirc_gr_S")),
        "normalized_counts_and_safe_boundary_are_exact": (
            norm["integrations"] == {"A0": 50001, "A1": 2}
            and norm["raw_sources"] == {"A0": 100005, "A1": 7}
            and norm["largest_first_real_ceiling"] == 300021),
        "physical_counts_and_safe_boundary_are_exact": (
            phys["integrations"] == {"A0": 30100002, "A1": 30050003}
            and phys["raw_sources"] == {"A0": 60200007, "A1": 60100009}
            and phys["direct_A1_first_real_ceiling"] == 60300013
            and phys["radial_nested_A1_first_real_ceiling"] == 60400023
            and phys["largest_first_real_ceiling"] == 60400023
            and phys["radial_branch_dominates"]),
        "both_compact_physical_margins_are_strict": (
            phys["A0_physical_exponent"]-phys["requested_q_power"] == Q(1, 20000000)
            and phys["A1_physical_exponent"]-phys["requested_q_power"] == Q(1, 20000000)
            and interface["compact_flat_for_each_fixed_m_N"]),
        "regular_physical_exponent_remains_negative": (
            interface["regular_physical_exponent"] == Q(-15025001, 10000000)
            and not interface["current_regular_decay_certified"]),
        "scope_remains_conditional": (
            fixed_order_induction(physical=True)["requested_ceiling_is_finite"]
            and not fixed_order_induction(physical=True)["numerical_high_jet_materialization_required"]
            and not fixed_order_induction(physical=True)["primitive_contracts_independently_certified"]
            and not interface["complete_residual_decay_proved"]),
    }
    return dict(checks=checks, normalized_request=norm,
                physical_illustration=phys, interface=interface,
                induction=fixed_order_induction(physical=True))
