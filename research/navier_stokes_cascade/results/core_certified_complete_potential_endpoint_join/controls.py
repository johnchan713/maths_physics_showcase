"""Exact controls for the certified complete-potential diagonal endpoint join."""
from fractions import Fraction as Q
from endpoint import (
    diagonal, atlas_review, representative_review, punctured_cover_review,
    exterior_review, heat_A, heat_S, pressure_terminal,
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
    atlas = atlas_review()
    reps = representative_review()
    cover = punctured_cover_review()
    exterior = exterior_review()
    h = Q(1, 200)
    checks = dict(
        certified_diagonal_seed=diagonal.current.rates(0)["regular"] == Q(3, 2) - 2 * Q(1, 100000),
        certified_diagonal_gain=diagonal.current.parameters()["gain"] == Q(1, 2) - 3 * Q(1, 100000),
        enlarged_closed_atlas=all(atlas["checks"].values()),
        complete_representatives=all(reps["checks"].values()),
        punctured_cover=all(cover["checks"].values()),
        complete_exterior_jets=all(exterior["checks"].values()),
        pressure_second_jet_keeps_cross=(heat_S(h, 2) - 2 * heat_A(h, 2) == 2 * heat_A(h, 1) ** 2),
        pressure_time_sign_and_memory=pressure_terminal(h, 0, 1) < 0,
    )
    return dict(
        atlas=encode(atlas),
        representatives=encode(reps),
        punctured_cover=encode(cover),
        exterior=encode(exterior),
        checks=checks,
        exact_controls_are_not_actual_endpoint_fields=True,
    )
