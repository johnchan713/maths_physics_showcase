"""Exact checks for the corrected complete-potential diagonal join."""
from fractions import Fraction as Q
from diagonal import (
    KAPPA, current, representative_orders, gain, comparison,
    scale_review, cut_loss, tuple_loss, universal_controls,
)


def encode(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def review():
    delta = current.parameters(KAPPA)["gain"]
    first = representative_orders(1)
    second = representative_orders(2)
    margins = {key: value - delta for key, value in first.items()}
    examples = [comparison(m, N, h) for h in (Q(1, 200), Q(9, 1000), Q(1, 10000))
                for m, N in ((0, 0), (1, 3), (5, 11), (20, 2))]
    universal = universal_controls()
    poly = universal.polynomial_review()
    implicit = universal.implicit_review()
    checks = dict(
        corrected_gain=delta == Q(1, 2) - 3 * KAPPA,
        corrected_seed=current.rates(0)["regular"] == Q(3, 2) - 2 * KAPPA,
        every_representative_has_exact_delta_slope=(
            second == {key: value + delta for key, value in first.items()}),
        weakest_representative_margin_is_one=min(margins.values()) == 1,
        one_extra_curl_loss_reserved=cut_loss(3) == tuple_loss(4),
        every_comparison_closes=all(all(row["checks"].values()) for row in examples),
        default_comparison_stage_is_12816=comparison()["comparison_stage"] == 12816,
        whole_interval_scale_majorant=all(scale_review()["checks"].values()),
        complete_cut_velocity_exactly_divergence_free=(
            poly["checks"]["complete_cut_velocity_exactly_divergence_free"]),
        all_six_residual_groups_retained=(
            poly["nonzero_individually_necessary_terms"] ==
            ["new_old", "new_square", "old_new", "pressure", "time", "viscosity"]),
        implicit_q_cutoff_jets_checked=implicit["positive_derivative_cases"] == 54,
    )
    return dict(
        corrected_gain=str(delta),
        first_cycle_orders=encode(first),
        representative_margins=encode(margins),
        comparisons=encode(examples),
        scale_majorant=encode(scale_review()),
        residual_groups=poly["nonzero_individually_necessary_terms"],
        checks=checks,
        diagnostics_are_exact_fixtures_not_actual_pulses=True,
    )
