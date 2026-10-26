from fractions import Fraction as Q
from stage_sync import compact_difference_identity, manufactured_change, scope, stages


def review():
    s = stages()
    identity = compact_difference_identity()
    fixture = manufactured_change()
    flags = scope()
    checks = {
        "retained_orders_are_exact": (
            s["first_real"]["regular_order"] == Q(99999, 100000)
            and s["first_sourced"]["regular_order"] == Q(74999, 50000)),
        "physical_exponents_are_exact": (
            s["first_real"]["physical_exponent"] == Q(-30100001, 20000000)
            and s["first_sourced"]["physical_exponent"] == Q(-15025001, 10000000)),
        "sourced_gain_is_positive_but_both_rows_fail": (
            s["sourced_gain"] == Q(49999, 100000)
            and s["physical_gain"] == Q(49999, 20000000)
            and not s["first_real_meets_requested_bound"]
            and not s["first_sourced_meets_requested_bound"]),
        "compact_difference_keeps_covariance_and_pressure": (
            "Bil" in identity["covariance_change"]
            and "Delta_pm" in identity["pressure_change"]
            and "rho Delta_P" in identity["pressure_change"]),
        "initial_DAG_does_not_contain_sourced_state": (
            not identity["initial_graph_contains_first_sourced_wave"]
            and not identity["initial_compact_bound_transfers_to_sourced_state"]),
        "manufactured_changed_compact_inputs_are_nonzero": (
            fixture["radial_changed"] and fixture["axial_changed"]
            and fixture["Qcirc_gr"] == (Q(3, 2), Q(-3, 2))
            and fixture["Qcirc_phi"] == (Q(7), Q(-7))),
        "scope_separates_conditional_parts_from_complete_residual": (
            flags["initial_compact_tail_closed_conditionally"]
            and flags["post_sourced_regular_gain_closed_conditionally"]
            and not flags["synchronized_post_sourced_compact_DAG_closed_here"]
            and not flags["complete_physical_residual_decay_proved"]
            and not flags["blowup_proved"]),
    }
    return dict(checks=checks, stages=s, compact_identity=identity,
                manufactured_change=fixture, scope=flags)
