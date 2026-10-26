from fractions import Fraction as Q
from requests import pole_domains, qualitative_closure, request_graph, unresolved_inequalities


def review():
    checks = {}
    k = Q(1, 100000)
    g = request_graph(0, k)
    q = pole_domains(0)
    u = unresolved_inequalities(0)
    s = qualitative_closure(0)
    p = g["integrations"]
    checks["A1_count_is_least_strict_for_required_nu"] = (
        g["source_order"] + (p-1)*k <= g["required_order"]
        and g["source_order"] + p*k > g["required_order"]
        and p == 99995)
    checks["raw_DzPhi_and_Phi_orders_are_199993_and_199994"] = (
        g["raw_Dz_Phi_order"], g["raw_Phi_order"]) == (199993, 199994)
    checks["conservative_same_state_pressure_padding_stops_at_200002"] = (
        g["conservative_same_state_order"] == 200002)
    checks["same_choice_target_profile_and_Haar_paddings_are_explicit"] = (
        g["signed_target_order"], g["fixed_profile_order"], g["initial_Haar_order"]) == (200005, 200006, 200005)
    checks["projected_moment_requires_one_extra_sharp_initial_slot"] = (
        g["sharp_initial_order"] == g["initial_Haar_order"] + 1 == 200006)
    checks["initial_field_and_canonical_covariance_safe_orders_are_finite"] = (
        g["initial_field_order"], g["canonical_covariance_safe_order"]) == (200008, 200018)
    checks["immediate_A1_graph_does_not_import_the_deferred_A0_graph"] = (
        g["signed_target_order"] < g["old_combined_A0_A1_target_order_at_output_zero"]
        and g["A0_graph_deferred_to_complete_residual"])
    checks["corrected_A1_tail_keeps_strict_kappa_margin"] = (
        g["achieved_order"] == 3-5*k and g["strict_margin"] == k)
    checks["higher_pole_domains_are_recomputed_from_the_new_state_order"] = (
        q["q"] == 1200030 and q["Kp"] == (2400060, 4800120))
    checks["finite_profile_and_support_closure_needs_no_trace_enumeration"] = (
        s["fixed_C_infinity_profiles_have_finite_requested_seminorms"]
        and s["finite_support_grammar_gives_smooth_zero_extension_and_matching_traces"]
        and not s["numerical_enumeration_of_all_traces_or_factorials_required"])
    checks["qualitative_closure_does_not_promote_missing_weighted_field_rates"] = (
        not s["qualitative_smoothness_implies_weighted_field_rate"]
        and not g["actual_band_uniform_field_inequalities_proved"])
    checks["precise_Phi_DzPhi_and_weighted_prefix_blockers_are_recorded"] = (
        "epsilon P_199994" in u["Phi"] and "epsilon^2 Q_199993" in u["DzPhi"]
        and "200018" in u["weighted_prefix"])
    return dict(
        checks={name: bool(value) for name, value in checks.items()},
        graph=g,
        pole_domains=q,
        unresolved=u,
        qualitative_closure=s,
        parent_commit="50cd291d22ab147bc3a89186b142496fbda5e45f",
        correction=("the immediate corrected first-real final-row A1 request is a finite A1-only graph; "
                    "the older 400009 target index also included the later radial A0 residual tail"),
        scope=("exact finite dependency and support/profile audit for the corrected first-real A1 tail; "
               "no actual weighted high-jet enclosure, A0 estimate, complete residual or blow-up proof"),
    )
