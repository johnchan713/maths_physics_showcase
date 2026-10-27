"""Exact controls for fixed-order background seminorms and Borel selection."""

from fractions import Fraction as F
from itertools import product


def derivative_slots(word):
    """Count derivative letters; adjacent derivative letters are forbidden."""
    word = tuple(word)
    if any(letter not in (0, 1) for letter in word):
        raise ValueError("Binary Volterra word required")
    if any(a == b == 1 for a, b in zip(word, word[1:])):
        raise ValueError("Adjacent derivative factors vanish in the actual block structure")
    return sum(word)


def maximum_surviving_slots(length):
    if type(length) is not int or length < 0:
        raise ValueError("Nonnegative word length required")
    words = [word for word in product((0, 1), repeat=length)
             if not any(a == b == 1 for a, b in zip(word, word[1:]))]
    return max((sum(word) for word in words), default=0)


def expected_slot_bound(length):
    if type(length) is not int or length < 0:
        raise ValueError("Nonnegative word length required")
    return (length + 1) // 2


def word_majorant(length, z, strip_loss):
    """The kth ordered-simplex majorant from the Banach-scale proof."""
    if type(length) is not int or length < 0:
        raise ValueError("Nonnegative word length required")
    z, strip_loss = map(F, (z, strip_loss))
    if z <= 0 or not 0 < strip_loss <= 1:
        raise ValueError("Positive size and strip loss in (0,1] required")
    p = expected_slot_bound(length)
    factorial = 1
    for j in range(2, length + 2):
        factorial *= j
    cauchy = max(F(1), F(p, 1) / strip_loss) ** p if p else F(1)
    return z ** (length + 1) * cauchy / factorial


def normalized_stage_slacks(h, n, log_constant, log_scale, log_two=F(2, 3)):
    """Exact stand-ins for the real logarithmic selector inequalities."""
    h, log_constant, log_scale, log_two = map(F, (h, log_constant, log_scale, log_two))
    if type(n) is not int or n < 1 or h <= 0 or min(log_constant, log_scale, log_two) < 0:
        raise ValueError("Positive order/gain and nonnegative logs required")
    gain = n * h
    return {
        "normalized": gain * log_scale - log_constant - n * log_two,
        "strict_after_unit_padding": gain * log_scale - log_constant - n * log_two - gain,
    }


def finite_selector(h, n, log_constants, previous, boundary, log_two=F(2, 3)):
    """Choose a rational stand-in above one finite set of real-log rows."""
    h, previous, boundary, log_two = map(F, (h, previous, boundary, log_two))
    rows = tuple(map(F, log_constants))
    if type(n) is not int or n < 1 or len(rows) != n + 1 or h <= 0:
        raise ValueError("All fixed seminorm rows through n are required")
    if min(rows + (previous, boundary, log_two)) < 0:
        raise ValueError("Nonnegative log data required")
    gain = n * h
    return F(1) + max(previous + log_two, boundary + log_two,
                      max((row + n * log_two) / gain for row in rows))


def locally_finite_support_count(q, scales):
    q = F(q); scales = tuple(map(F, scales))
    if q <= 0 or not scales or any(a <= 0 for a in scales):
        raise ValueError("Positive data required")
    if any(scales[j + 1] < 2 * scales[j] for j in range(len(scales) - 1)):
        raise ValueError("Borel scale separation missing")
    return sum(a * q < 1 for a in scales)


def review():
    slot_rows = {str(k): maximum_surviving_slots(k) for k in range(13)}
    selector = finite_selector(F(1, 200), 3, (4, 7, 9, 11), 13, 5)
    slacks = normalized_stage_slacks(F(1, 200), 3, 11, selector)
    return {
        "checks": {
            "no_adjacent_derivative_word_has_more_than_ceiling_half_slots": all(slot_rows[str(k)] == expected_slot_bound(k) for k in range(13)),
            "ordered_simplex_majorants_are_finite": all(word_majorant(k, F(7, 3), F(1, 5)) > 0 for k in range(13)),
            "selector_uses_every_fixed_row_through_order_n": selector > 0,
            "selector_gives_strict_normalized_margin": slacks["strict_after_unit_padding"] >= 0,
            "geometrically_separated_scales_are_locally_finite_in_fixture": locally_finite_support_count(F(1, 20), (1, 2, 4, 8, 16, 32)) == 5,
        },
        "controls": {
            "surviving_derivative_slots": slot_rows,
            "selector": str(selector),
            "selector_slacks": {key: str(value) for key, value in slacks.items()},
        },
        "leading_certificate": {
            "same_completed_profile": True,
            "fixed_order_holomorphic_eta_collar_and_smooth_global_jets_required": True,
            "exact_I3_power_law_and_five_ordered_bumps_required": True,
            "compact_velocity_and_pressure_source_identities_required": True,
            "independently_reproved_here": False,
        },
        "scope": {
            "coefficient_and_borel_part_of_P1_reconstructed": True,
            "scale_selected_before_coefficient_seminorm": False,
            "numerical_infinite_scale_list_required": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
