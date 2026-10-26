from fractions import Fraction as Q
from recurrence import (closure_statement, first_stage_for_power, layer,
                        physical_request, pole_support_ledger, rates,
                        state_ceiling)


def review():
    l0, l1 = layer(0), layer(1)
    s0, s1 = state_ceiling(0), state_ceiling(1)
    p0, p1 = physical_request(0), physical_request(1)
    pole = pole_support_ledger(2)
    threshold = first_stage_for_power()
    scope = closure_statement()
    checks = {
        "corrected_stage_rates_and_gain_are_exact": (
            rates(0)["regular"] == Q(74999, 50000)
            and rates(0)["mean"] == Q(99999, 50000)
            and rates(1)["regular"]-rates(0)["regular"] == Q(49997, 100000)),
        "first_repeated_layer_uses_actual_trace_offsets": (
            l0["sourced_integrations"] == {"A0": 100000, "A1": 4}
            and l0["generic_integrations"] == {"A0": 149997, "A1": 50001}
            and l0["sourced_backward_offset"] == 200018
            and l0["generic_backward_offset"] == 300011
            and l0["full_cycle_backward_offset"] == 500029),
        "stage_state_ceiling_recurrence_is_exact_and_decreasing": (
            s0["first_real_ceiling"] == 500038
            and s1["first_real_ceiling"] == 1000067
            and s1["stage_decreases_on_every_recursive_edge"]),
        "seed_physical_request_is_preserved": (
            p0["raw_sources"] == {"A0": 60200007, "A1": 60100009}
            and p0["first_real_ceiling"] == 60500042),
        "next_physical_stage_is_finitely_expanded": (
            p1["raw_sources"] == {"A0": 60200007, "A1": 60000015}
            and p1["first_real_ceiling"] == 60900077
            and min(p1["strict_margins"].values()) == Q(1, 20000000)),
        "distinct_poles_and_supports_are_archived_stagewise": (
            pole["finite_for_fixed_stage"]
            and len(pole["distinct_common_poles"]) == 4
            and pole["harmonics_at_most"] == 16
            and pole["same_closed_supports_preserved_if_inherited_support_contracts_hold"]
            and not pole["actual_constants_or_high_jets_materialized"]),
        "corrected_complete_residual_threshold_is_finite": (
            threshold["stage"] == 1002
            and threshold["achieved_power"] >= 0
            and threshold["preceding_power"] < 0),
        "scope_is_conditional_and_diagonal_remains_open": (
            scope["finite_stage_compact_graph_synchronized_conditionally"]
            and scope["all_fixed_order_and_pole_primitive_contracts_required"]
            and not scope["one_field_all_targets_diagonal_proved"]
            and not scope["localization_forcing_energy_or_blowup_proved"]),
    }
    return dict(checks=checks, first_layer=l0, second_layer=l1,
                stage_zero=s0, stage_one=s1,
                physical_stage_zero=p0, physical_stage_one=p1,
                stage_two_poles=pole, physical_threshold=threshold,
                closure=scope)
