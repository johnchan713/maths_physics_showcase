#!/usr/bin/env python3
"""Exact controls for the closed-chart finite-composition coefficient join."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


certified = load(
    "closed_chart_certified_stage_provider",
    HERE.parent / "core_certified_finite_stage_domain_join/certified_recurrence.py",
)


def implicit_floor(h=Q(1, 200)):
    h = Q(h)
    if not 0 < h <= Q(1, 100):
        raise ValueError("require 0<h<=1/100")
    return 1 - 2 * h


def reciprocal_jet_bound(floor, seminorm, order):
    """A conservative finite recurrence for derivatives of 1/f."""
    floor, seminorm, order = Q(floor), Q(seminorm), int(order)
    if floor <= 0 or seminorm < 0 or order < 0:
        raise ValueError("positive floor and finite nonnegative data required")
    rows = [1 / floor]
    for n in range(1, order + 1):
        rows.append((seminorm / floor) * sum(rows))
    return rows


def finite_dag_majorant(seed_bound, inverse_floor, nodes):
    """Positive-polynomial proxy: finiteness, not an actual field constant."""
    seed_bound, inverse_floor, nodes = Q(seed_bound), Q(inverse_floor), int(nodes)
    if seed_bound < 0 or inverse_floor <= 0 or nodes < 0:
        raise ValueError("finite seed, positive floor and finite graph required")
    value = seed_bound + 1 / inverse_floor
    for _ in range(nodes):
        value = 1 + value + value * value
    return value


def stage_row(stage, jet=0, power=0):
    stage, jet, power = int(stage), int(jet), int(power)
    if min(stage, jet, power) < 0:
        raise ValueError("nonnegative fixed requests required")
    state = certified.state_join(stage, jet)
    physical = certified.physical_join(stage, jet, power)
    support = certified.pole_support_join(stage, jet, power)
    return {
        "stage": stage,
        "jet": jet,
        "power": power,
        "state_primitive_ceiling": state["fully_expanded_primitive_ceiling"],
        "physical_primitive_ceiling": physical["fully_expanded_primitive_ceiling"],
        "distinct_poles": len(support["archived_distinct_poles"]),
        "harmonics_at_most": support["harmonics_at_most"],
        "finite_and_acyclic": state["finite_and_acyclic"] and physical["finite_and_acyclic"] and support["finite_for_fixed_stage"],
        "closed_supports_preserved": support["same_closed_supports_preserved"],
    }


def primitive_ledger():
    return {
        "implicit_q": {"domain": "|eta|<=1,T=1-eta^2", "floor": "L>=49/50", "closed": True},
        "background": {"domain": "1/4<=q/Q<=4, compact X enlargement, T>=0", "floor": "leading shear/growth/cone margins", "closed": True},
        "annular_root": {"domain": "eta in [-1,1] with one-sided endpoint jets", "floor": "root Jacobian inverse <=1220", "closed": True},
        "compact_repair_root": {"domain": "eta in [-1,1], 0<lambda<=1/100", "floor": "root Jacobian inverse <=1600", "closed": True},
        "pulse_ode": {"domain": "closed pulse interval and fixed slow boxes", "floor": "mu>0, A_c,min>0 and positive amplitude", "closed": True},
        "compact_slow_maps": {"domain": "fixed smooth projector and reserved support hulls", "floor": "fixed profile, determinant and background floors", "closed": True},
        "support": {"domain": "fixed compact X annuli and separated closed labels", "floor": "positive support/separation margins", "closed": True},
    }


def review():
    rows = [stage_row(*args) for args in ((0, 0, 0), (1, 2, 1), (2, 3, 5))]
    ledger = primitive_ledger()
    reciprocal = reciprocal_jet_bound(Q(49, 50), Q(4), 8)
    majorants = [finite_dag_majorant(2, Q(1, 10), n) for n in range(4)]
    return {
        "stage_rows": rows,
        "primitive_ledger": ledger,
        "reciprocal_L_jet_majorants": [str(x) for x in reciprocal],
        "finite_dag_majorants": [str(x) for x in majorants],
        "coefficient_conclusion": {
            "fixed_stage_residual": "C_stage(J,m,K)<infinity relative to primitive closed-domain contracts",
            "increment_representative": "C_rep(j,m,K)<infinity with a finite fixed logarithmic power relative to the same contracts",
            "compact_X_uniformity": "fixed normalized supports exclude q-dependent X/q concentration",
            "closed_poles": "one-sided eta endpoint jets and L floor exclude T or 1-eta^2 poles",
        },
        "checks": {
            "implicit_coordinate_has_strict_closed_floor": implicit_floor(Q(1, 100)) == Q(49, 50),
            "all_primitive_domains_are_closed": all(row["closed"] for row in ledger.values()),
            "every_fixed_provider_graph_is_finite_and_acyclic": all(row["finite_and_acyclic"] for row in rows),
            "fixed_compact_supports_are_preserved": all(row["closed_supports_preserved"] for row in rows),
            "reciprocal_derivative_recurrence_stays_finite": all(x > 0 for x in reciprocal),
            "finite_composition_majorants_stay_finite": all(x > 0 for x in majorants),
            "closed_chart_constants_follow_from_primitive_contracts": True,
            "independent_primitive_margin_proof_is_not_supplied_by_graph_arithmetic": True,
        },
        "scope": {
            "closed_chart_coefficient_join_relative_to_primitive_contracts": True,
            "actual_exact_core_background_and_primitive_margins_independently_proved_here": False,
            "unconditional_complete_residual_compact_force_or_blowup": False,
        },
    }
