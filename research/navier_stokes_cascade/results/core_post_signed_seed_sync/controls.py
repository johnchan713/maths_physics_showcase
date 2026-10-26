from fractions import Fraction as Q
from ledger import (corrected_seed, interface, normalized_seed_request,
                    physical_seed_request, stage_obligation)


def review():
    norm = normalized_seed_request()
    phys = physical_seed_request()
    seed = corrected_seed()
    later = stage_obligation()
    scope = interface()
    checks = {
        "normalized_post_signed_graph_is_fully_expanded": (
            norm["post_signed_integrations"] == {"A0": 100000, "A1": 50001}
            and norm["post_signed_raw_sources"] == {"A0": 200003, "A1": 100005}
            and norm["first_sourced_target_order"] == 200017
            and norm["first_real_ceiling"] == 500038),
        "physical_post_signed_graph_is_fully_expanded": (
            phys["current_A0"]["raw_source_order"] == 60200007
            and phys["current_A1"]["raw_source_order"] == 60100009
            and phys["first_sourced_target_order"] == 60200021
            and phys["first_real_ceiling"] == 60500042),
        "physical_compact_margins_remain_strict": (
            phys["current_A0"]["strict_margin"] == Q(1, 20000000)
            and phys["current_A1"]["strict_margin"] == Q(1, 20000000)),
        "illustrative_current_common_pole_is_retained": (
            Q(phys["illustrative_domain"]["q"]) == 361200138
            and not phys["illustrative_domain"]["actual_input_bounds_evaluated"]),
        "corrected_seed_satisfies_old_half_line_hypotheses": (
            seed["regular_seed"] == Q(74999, 50000)
            and seed["mean_seed"] == Q(99999, 50000)
            and seed["seed_margin_over_old_minimum"] == Q(1, 25000)
            and seed["mean_relation"] == Q(1, 2)
            and seed["finite_stage_gain"] == Q(49997, 100000)
            and seed["half_line_rate_recurrence_applicable"]),
        "later_physical_integer_rules_have_strict_margins": (
            later["strict_margins"]["A0"] > 0
            and later["strict_margins"]["A1"] > 0),
        "later_stage_compact_sync_is_not_inferred_from_rate_arithmetic": (
            not later["stage_indexed_backward_request_recurrence_supplied"]
            and not later["stage_indexed_common_pole_and_support_recurrence_supplied"]
            and not scope["old_half_line_rate_recurrence_alone_proves_later_compact_inputs"]),
        "scope_remains_conditional_and_nonterminal": (
            scope["seed_compact_graph_synchronized_conditionally"]
            and not scope["complete_finite_stage_residual_decay_proved"]
            and not scope["infinite_realization_or_blowup_proved"]),
    }
    return dict(checks=checks, normalized=norm, physical=phys,
                corrected_seed=seed, next_stage_obligation=later, interface=scope)
