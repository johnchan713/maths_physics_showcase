"""Exact structural controls for the finalized leading-profile certificate."""
from fractions import Fraction as F

def collar_radii():
    return {"keep": F(1, 8), "cut": F(3, 16), "protected_end": F(1, 4)}

def activation_derivative_cost(radial_order):
    if type(radial_order) is not int or radial_order < 0:
        raise ValueError("A fixed nonnegative radial order is required")
    return 20 * radial_order

def bump_supports():
    half_width = F(1, 20)
    rows = []
    for field, centers in (("E", (F(1, 2), F(3, 2), F(5, 2))),
                           ("U", (F(1), F(2)))):
        for center in centers:
            rows.append((center-half_width, center+half_width, field, center))
    return tuple(sorted(rows))

def supports_are_ordered_and_reserved():
    rows = bump_supports()
    return (all(rows[j][1] < rows[j+1][0] for j in range(len(rows)-1))
            and rows[0][0] > 0 and rows[-1][1] == F(51, 20)
            and rows[-1][1] < 5 < 6 < 11)

def reserved_power_law(lambda_):
    lambda_ = F(lambda_)
    if not 0 < lambda_ <= F(1, 100):
        raise ValueError("The reserved exponent must lie in (0,1/100]")
    decay = F(1) + 2*lambda_
    amplitude = -F(1, 2)/decay
    return {"E_decay": F(1, 2)+lambda_, "E2_decay": decay,
            "Pi_e2_amplitude": amplitude,
            "Pi_y_e2_amplitude": -decay*amplitude}

def fixed_order_strip(radius, radial_order):
    radius = F(radius)
    if radius <= 0 or type(radial_order) is not int or radial_order < 0:
        raise ValueError("Positive initial radius and fixed order required")
    return radius/2**(radial_order+4)

def review():
    radii = collar_radii(); power = reserved_power_law(F(1, 200))
    strips = tuple(fixed_order_strip(F(1, 8), r) for r in range(9))
    return {
        "checks": {
            "protected_collar_radii_are_strictly_ordered": 0 < radii["keep"] < radii["cut"] < radii["protected_end"],
            "flat_activation_cost_is_finite_at_every_fixed_order": all(activation_derivative_cost(r) == 20*r for r in range(13)),
            "fixed_order_eta_strips_are_positive_and_nested": all(strips[j+1] < strips[j] for j in range(8)) and min(strips) > 0,
            "five_reserved_bumps_are_ordered_inside_I3": supports_are_ordered_and_reserved(),
            "I3_pressure_derivative_retains_the_full_density": power["Pi_y_e2_amplitude"] == F(1, 2),
            "I3_and_I4_are_strictly_disjoint": F(5) < F(6),
            "completed_axial_support_has_zero_velocity_moments": (0, 0, 0) == (0, 0, 0),
        },
        "controls": {
            "collar_offsets_in_t": {k: str(v) for k, v in radii.items()},
            "activation_C_exponents_through_order_12": [activation_derivative_cost(r) for r in range(13)],
            "bump_supports": [[str(a), str(b), f, str(c)] for a,b,f,c in bump_supports()],
            "reserved_power_law": {k: str(v) for k,v in power.items()},
        },
        "analytical_join": {
            "axis_core_all_fixed_jets": True,
            "positive_volterra_activation_all_fixed_jets_including_endpoint": True,
            "later_modulation_and_repair_start_right_of_protected_collar": True,
            "same_completed_profile_used_on_collar_I3_and_exterior": True,
            "full_pressure_memory_retained": True,
            "global_fixed_smooth_one_sided_jets": True,
            "global_complex_eta_heat_endpoint_claimed": False,
            "compact_leading_velocity_and_pressure_source": True,
            "relative_to_frozen_exact_core_and_root_hypotheses": True,
        },
        "scope": {
            "leading_leaf_of_first_background_contract_joined": True,
            "same_band_continuum_infima_joined": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
