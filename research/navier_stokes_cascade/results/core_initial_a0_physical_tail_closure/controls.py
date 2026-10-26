from fractions import Fraction as Q
from closure import (DEPENDENCIES, PRIMITIVES, closure_statement,
                     fixed_order_induction, normalized_request,
                     physical_interface, physical_request,
                     radial_rate_closure, topological_order)


def review():
    checks = {}
    norm = normalized_request()
    phys = physical_request()
    rate = radial_rate_closure()
    statement = closure_statement()
    position = {name: i for i, name in enumerate(topological_order())}
    checks["extended_first_real_dependency_graph_is_acyclic"] = all(
        dep in PRIMITIVES or position[dep] < position[name]
        for name, deps in DEPENDENCIES.items() for dep in deps)
    checks["desired_A0_output_is_not_a_gr_input"] = (
        not fixed_order_induction()["desired_A0_output_used_as_input"]
        and all("A0_Qcirc_gr_final" not in deps for name, deps in DEPENDENCIES.items()
                if name != "A0_Qcirc_gr_final"))
    checks["normalized_A0_source_and_graph_indices_are_exact"] = (
        norm["A0_integrations"] == 199996
        and norm["A0_source_order"] == 399995
        and norm["largest_same_prefix_order"] == 400009
        and norm["trace"]["maxima"]["next_signed_target_input"] == 400009)
    checks["latest_A1_route_is_retained_inside_the_A0_graph"] = (
        norm["A1_integrations"] == 99995
        and norm["A1_DzPhi_source_order"] == 199993
        and norm["A1_Phi_source_order"] == 199994)
    checks["normalized_radial_tail_has_strict_kappa_margin"] = (
        rate["Qcirc_gr_final"] == Q(1)-Q(1,100000)
        and min(rate["radial_source_groups"].values()) == rate["Qcirc_gr_final"]
        and rate["radial_source_groups"]["linear_centripetal_mean"] == rate["Qcirc_gr_final"]
        and rate["radial_source_groups"]["radial_covariance_flux"] == rate["Qcirc_gr_final"]
        and rate["A0_tail"] == Q(3)-Q(5,100000)
        and rate["full_radial_row"] == Q(3)-Q(5,100000)
        and rate["strict_margin"] == Q(1,100000) and rate["strict"])
    checks["physical_zero_jet_zero_power_request_is_exact"] = (
        phys["A0_integrations"] == 30100002
        and phys["A0_source_order"] == 60200007
        and phys["A1_integrations"] == 30000001
        and phys["A1_DzPhi_source_order"] == 60000005
        and phys["A1_Phi_source_order"] == 60000006
        and phys["largest_same_prefix_order"] == 60200021)
    checks["both_compact_physical_components_have_positive_margin"] = (
        phys["A0_strict_margin"] > 0 and phys["A1_strict_margin"] > 0
        and physical_interface()["each_fixed_jet_and_power_only"])
    checks["all_fixed_order_induction_covers_normalized_and_physical_graphs"] = (
        fixed_order_induction()["requested_order_is_finite"]
        and fixed_order_induction(physical=True)["requested_order_is_finite"]
        and all(v["finite_positive_polynomial"]
                for v in fixed_order_induction(physical=True)["nodes"].values())
        and not fixed_order_induction(physical=True)["numerical_materialization_required"])
    checks["source_inequality_follows_conditionally_not_unconditionally"] = (
        "399995" in statement["conditional_radial_source_inequality"]
        and statement["high_jet_is_not_an_additional_independent_premise"]
        and not statement["actual_primitive_contracts_certified_here"])
    checks["pressure_support_and_residual_scope_is_honest"] = (
        statement["complete_pressure_memory_and_unit_mass_compensation_retained"]
        and statement["whole_support_endpoint_and_seam_closure_reused"]
        and statement["compact_physical_tail_interface_closed_conditionally"]
        and not statement["complete_residual_or_blowup_closed"]
        and not physical_interface()["complete_regular_residual_improvement_proved"])
    return dict(
        checks={name: bool(value) for name, value in checks.items()},
        normalized_request=norm, physical_illustration=phys,
        normalized_rates=rate, physical_interface=physical_interface(),
        induction=fixed_order_induction(), statement=statement,
        parent_commit="ab008e963cb105cfc07aa6a5b3a3b873c5984c0a",
        correction=("the deferred order-400009 graph is the radial A0 graph; its "
                    "order-399995 source bound follows at that fixed order from the "
                    "same acyclic all-fixed-order prefix, conditional on the primitive contracts"),
        scope=("conditional normalized radial A0 closure and fixed-(m,N) compact physical-tail "
               "interface; primitive contracts, the regular residual, later cutoffs and blow-up remain open"),
    )
