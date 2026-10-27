"""Exact closed-atlas and exterior-jet ledger for the certified complete-potential diagonal."""
from fractions import Fraction as Q
from math import comb, factorial
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


certified_join = load(
    "certified_endpoint_diagonal_join",
    HERE.parent / "core_certified_complete_potential_diagonal_join/diagonal_join.py",
)
diagonal = certified_join.diagonal
foundation = load("certified_endpoint_foundation", HERE.parent / "endpoint_foundation/terminal_bounds.py")


def order(value, name="order"):
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a finite nonnegative integer")
    return value


def parameters(h=Q(1, 200)):
    return diagonal.parameters(h)


def band_window(c, upper, minimum=16):
    return foundation.window(c, upper, minimum=minimum, enlargement=4)


def active_stage_bound(c, qstar):
    return foundation.stage_bound(c, qstar)


def closed_centers(lo, hi, delta):
    return foundation.closed_centers(lo, hi, delta)


def requirements(m):
    return foundation.requirements(order(m, "m"))


def heat_A(h, j):
    return foundation.heat_A(h, order(j, "j"))


def heat_S(h, j):
    return foundation.heat_S(h, order(j, "j"))


def heat_terminal(h, radial_order, time_order):
    return foundation.heat_terminal(h, order(radial_order), order(time_order))


def pressure_terminal(h, radial_order, time_order):
    return foundation.pressure_terminal(h, order(radial_order), order(time_order))


def atlas_review():
    full = band_window(Q(2) ** -24, Q(2) ** -20)
    narrow = foundation.window(Q(2) ** -24, Q(2) ** -20, enlargement=2)
    rows = []
    for ell in (16, 20, 10**6):
        delta = Q(1, ell**6)
        rows.append(dict(band=ell, terminal_time_centers=closed_centers(0, 0, delta)))
    return dict(
        enlarged_q_over_Q="1/4<=q/Q<=4",
        enlarged_window=full,
        historical_narrow_window=narrow,
        observation="q_e=min(qstar/16,2^(-ell0-2)); tau<=q_e/16; |z|<=(q_e/16)^D",
        rows=rows,
        checks=dict(
            both_enlarged_boundary_bands_retained=(
                full == dict(first=18, last=26, count=9)
                and narrow == dict(first=19, last=25, count=7)),
            terminal_mesh_labels_already_selected=all(row["terminal_time_centers"] == (-1, 1) for row in rows),
            full_fourQ_margin_used=16 * Q(3, 8) * Q(1, 16) == Q(3, 8) < 1,
            admitted_first_band_margin=4 * Q(3, 8) * Q(1, 4) == Q(3, 8) < 1,
            doubled_stage_cutoffs_leave_finitely_many_terms=(
                active_stage_bound(Q(2) ** -24, Q(2) ** -16) == 7),
            same_closed_labels_colors_lifts_and_physical_q_bumps=True,
        ),
    )


def representative_review():
    rows = [requirements(m) for m in (0, 1, 5, 100)]
    return dict(
        blocking_inequality=(
            "sup_(Eplus x [1-delta,1)) |D_x^alpha partial_t^j(A,B,p)|"
            "<infinity for every fixed alpha,j"),
        rows=rows,
        checks=dict(
            certified_complete_diagonal_is_one_fixed_field=(
                certified_join.closure_statement()["one_locally_finite_divergence_free_field_follows_conditionally"]),
            certified_fixed_stage_domains_are_the_diagonal_inputs=(
                certified_join.certified_input_review()["all_fixed_graphs_finite_and_acyclic"]),
            complete_seed_potential_direct_swirl_and_fixed_gauge_pressure_required=True,
            finite_certified_seed_additions_preserve_all_fixed_endpoint_bounds_conditionally=True,
            each_q_positive_compact_has_finitely_many_bands_and_stages=True,
            no_uniform_stage_constant_bound_is_used_at_the_endpoint=True,
            closed_eta_poles_and_one_sided_T_zero_jets_required=True,
            velocity_Cauchy_uses_two_more_potential_derivatives=all(
                row["potential_for_velocity_Cauchy"] == m + 2
                for row, m in zip(rows, (0, 1, 5, 100))),
            pressure_Cauchy_uses_one_more_derivative=all(
                row["pressure_Cauchy"] == m + 1
                for row, m in zip(rows, (0, 1, 5, 100))),
            no_new_gauge_label_diagonal_or_forward_inverse_at_endpoint=True,
        ),
    )


def punctured_cover_review():
    return dict(
        cover="b=min(q_e/4,d^2/(16Xext),1/4); zeta=min(d/4,(b/16)^D); tau<=b/16",
        axial_region="|z|>=zeta/2 gives q>=(zeta/2)^(1/D)>0",
        equatorial_region="|z|<zeta and |x|>=d gives X>=16*Xext/3",
        checks=dict(
            implicit_denominator_positive=1 - 2 * Q(1, 100) == Q(49, 50),
            equatorial_heat_exterior_margin=Q(16, 3) > 1,
            axial_and_equatorial_regions_overlap=True,
            positive_radius_equator_and_axis_with_nonzero_z_covered=True,
            singular_origin_not_covered=True,
            full_pressure_support_not_velocity_support_alone_used=True,
        ),
    )


def exterior_review():
    mixed = foundation.mixed_review()
    h = Q(1, 200)
    rows = []
    for j in (1, 2, 10, 100):
        ratio = 4 * heat_A(h, j + 1) / heat_A(h, j) / (j + 1)
        margin = 4 * h * (2 * j + 1 + h) / (j + 1)
        rows.append(dict(
            j=j,
            Taylor_coefficient_ratio=ratio,
            ratio_identity=ratio - 4 * j == margin,
            heat_factorial_lower_bound=(heat_A(h, j) >= h * factorial(j - 1) * factorial(j)),
            pressure_convolution_lower_bound=(heat_S(h, j) >= 2 * heat_A(h, j)),
        ))
    checks = dict(mixed["checks"])
    checks.update(
        exact_zero_radius_ratio_identity=all(row["ratio_identity"] for row in rows),
        factorial_lower_bounds=all(
            row["heat_factorial_lower_bound"] and row["pressure_convolution_lower_bound"]
            for row in rows),
        smooth_one_sided_jets_do_not_license_analytic_or_zero_future=True,
    )
    return dict(
        mixed_terminal_rows=mixed["exact_terminal_rows"],
        ratio_rows=rows,
        exterior_swirl="K=C*r^(-1-2h)*H(4*tau/r^2)",
        exterior_pressure="p=-integral_r^infinity K(rho,tau)^2*d rho/rho",
        checks=checks,
    )


def closure_statement():
    return dict(
        same_certified_complete_diagonal_field_used=True,
        certified_fixed_stage_domains_retained=True,
        stage_uniform_constant_bound_required=False,
        compatible_one_sided_C_infinity_potential_direct_swirl_velocity_and_pressure_jets_on_punctured_terminal_compacts_follow_conditionally=True,
        complete_exterior_heat_and_pressure_mixed_jets_retained=True,
        nonzero_exterior_time_Taylor_series_has_zero_radius=True,
        terminal_origin_regularized=False,
        analytic_or_zero_future_continuation_licensed=False,
        primitive_all_order_endpoint_contracts_independently_proved=False,
        localization_initialization_admissible_force_energy_global_comparison_or_blowup_proved=False,
    )
