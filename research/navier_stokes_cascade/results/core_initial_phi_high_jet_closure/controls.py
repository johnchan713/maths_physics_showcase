from fractions import Fraction as Q
from closure import DEPENDENCIES, PRIMITIVES, closure_statement, fixed_order_induction, rate_closure, topological_order


def review():
    checks = {}
    induction = fixed_order_induction(0)
    graph = induction["request"]
    rates = rate_closure()
    statement = closure_statement()
    position = {name: i for i, name in enumerate(topological_order())}
    checks["first_real_dependency_graph_is_acyclic"] = all(
        dep in PRIMITIVES or position[dep] < position[name]
        for name, deps in DEPENDENCIES.items() for dep in deps)
    checks["desired_A1_output_is_not_a_Phi_input"] = (
        not induction["desired_A1_output_used_as_input"] and
        all("A1_output" not in deps for deps in DEPENDENCIES.values()))
    checks["all_fixed_order_induction_reaches_safe_200018_boundary"] = (
        graph["largest_immediate_order"] == 200018 and induction["requested_order_is_finite"]
        and all(v["finite_positive_polynomial"] for v in induction["nodes"].values()))
    checks["fixed_order_constants_need_not_be_numerically_materialized"] = (
        not induction["numerical_materialization_required"])
    checks["complete_postsigned_Phi_groups_have_order_one"] = (
        rates["Qcirc_Phi"] == 1 and min(rates["Phi_groups"].values()) == 1
        and rates["Phi_groups"]["gamma_square"] > 1)
    checks["explicit_Dz_and_compact_gain_reach_three_minus_five_kappa"] = (
        rates["Dz_Qcirc_Phi"] == 2 and rates["integrations"] == 99995
        and rates["A1_tail"] == Q(3)-Q(5,100000))
    checks["conditional_A1_tail_has_strict_kappa_margin"] = (
        rates["strict"] and rates["strict_margin"] == Q(1,100000))
    checks["Phi_and_DzPhi_inequalities_follow_conditionally_from_primitives"] = (
        "199994" in statement["conditional_Phi_inequality"]
        and "199993" in statement["conditional_DzPhi_inequality"]
        and statement["high_jet_is_not_an_additional_independent_premise"])
    checks["primitive_validity_is_not_silently_promoted"] = (
        not induction["primitive_contract_independently_certified"]
        and not statement["actual_primitive_contracts_certified_here"])
    checks["pressure_support_A0_and_complete_residual_scope_is_honest"] = (
        statement["full_pressure_memory_retained"] and statement["whole_support_and_seam_closure_reused"]
        and not statement["A0_or_complete_residual_closed"])
    return dict(
        checks={name: bool(value) for name, value in checks.items()},
        induction=induction,
        rates=rates,
        statement=statement,
        parent_commit="be781606f848d8d7ce3389f54288a25c6156a756",
        correction=("the order-200018 Phi request is not a separate family of numerical inputs: "
                    "the retained all-fixed-order positive norm DAG implies it at that finite order, "
                    "conditional on the named primitive accepted-tail contracts"),
        scope=("conditional same-state Phi/DzPhi high-jet closure from the finite first-real DAG; "
               "primitive background/phase/root/ODE contracts, A0 and complete residual remain open"),
    )
