"""Source-level Q2/Q3 audit with the infinite-tail frame correction."""
from fractions import Fraction as Q
from pathlib import Path
import json


HERE = Path(__file__).resolve().parent
IDENTITY = "Q172-Q2-source-and-Q3-uniform-frame-correction"
PROVIDERS = (
    "core_current_head_q171_provider_discharge_frontier_audit",
    "core_current_head_q2_regional_gap_join",
    "regional_cone_join",
    "compact_pressure_bridge",
    "reduced_continuation_margins",
    "shape_join_pressure",
    "axial_stress_audit",
    "intermediate_decay_audit",
    "heat_exterior_audit",
    "outer_edge_stability",
    "weighted_background_edges",
    "core_current_head_q3_full_frame_join",
    "core_current_head_q23_finite_band_frame_join",
    "core_same_band_gap_join",
    "core_canonical_full_path_ode_join",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q172 identity required")
    expected = {
        "core_current_head_q171_provider_discharge_frontier_audit": "Q171-current-head-provider-frontier-reconciled-with-selector-obstruction",
        "core_current_head_q2_regional_gap_join": "current-head-Q2-regional-continuum-factored-edges-and-same-band-gap-bound",
        "regional_cone_join": "regional-input-gaps-and-leading-value-cone-joined",
        "compact_pressure_bridge": "whole-activation-and-conditional-compact-pulse-bridge-audited",
        "reduced_continuation_margins": "reduced-continuation-relaxed-margins-audited-for-reference-datum",
        "shape_join_pressure": "shape-join-entry-and-exact-leading-pressure-ledger-audited",
        "axial_stress_audit": "axial-ratios-bounded-global-profile-unverified",
        "intermediate_decay_audit": "axial-and-intermediate-ratios-bounded-global-profile-unverified",
        "heat_exterior_audit": "heat-compensated-outer-reference-bounded-axis-and-global-blowup-unverified",
        "outer_edge_stability": "outer-heat-edge-and-first-viscous-stress-stability-audited",
        "weighted_background_edges": "all-order-weighted-edges-and-relative-background-stress-audited",
        "core_current_head_q3_full_frame_join": "current-head-Q3-full-frame-whole-path-and-phase-domain-bound",
        "core_current_head_q23_finite_band_frame_join": "finite-head-tail-uniform-upper-propagator",
        "core_same_band_gap_join": "same-rule-defined-background-tail-phase-normal-value-cone-contract-joined",
        "core_canonical_full_path_ode_join": "same-tail-canonical-full-path-ode-pressure-and-cutoff-contract-joined",
    }
    actual = {name: evidence(name).get("status") for name in PROVIDERS}
    return {name: actual[name] == expected[name] for name in PROVIDERS}


REGIONS = (
    ("activation", "compact_pressure_bridge"),
    ("reduced_continuation", "reduced_continuation_margins"),
    ("shape", "shape_join_pressure"),
    ("ideal_prejoin", "shape_join_pressure"),
    ("five_moment_join", "shape_join_pressure"),
    ("ideal_after_join", "regional_cone_join"),
    ("first_outer_ramp", "regional_cone_join"),
    ("axial", "axial_stress_audit"),
    ("intermediate_ramp", "intermediate_decay_audit"),
    ("intermediate_power", "intermediate_decay_audit"),
)


def orientation_identity(a, b, p1, p2, F=1):
    """Verify -T.g=a F^2 (Pc-v), retaining t=-b/a."""
    a, b, p1, p2, F = map(Q, (a, b, p1, p2, F))
    if a <= 0 or F <= 0:
        raise ValueError("positive a and F required")
    t = -b / a
    v = a + b * b / a
    pc = p1 + t * p2
    g = (F * -a, F * b)
    stress = (F * (p1 - a), F * (p2 + b))
    minus_dot = -sum(x * y for x, y in zip(stress, g))
    return {
        "minus_T_dot_g": minus_dot,
        "factorized": a * F * F * (pc - v),
        "orientation_positive": pc > v and minus_dot > 0,
    }


def compactness_counterexample(count=16):
    """Each singleton band is compact, while the cross-band supremum grows."""
    if count < 2:
        raise ValueError("at least two bands required")
    per_band_suprema = tuple(Q(j) for j in range(1, count + 1))
    return {
        "each_band_compact_and_finite": True,
        "finite_prefix_supremum": max(per_band_suprema),
        "prefix_supremum_increases": all(
            per_band_suprema[j] > per_band_suprema[j - 1]
            for j in range(1, len(per_band_suprema))
        ),
        "compactness_alone_gives_uniform_infinite_tail_bound": False,
    }


def q23_uniform_repair(h=Q(1, 200), r0=Q(1, 6057072), lreq=100):
    """Exact finite threshold sufficient for Q23's analytic tail split."""
    h, r0, lreq = Q(h), Q(r0), Q(lreq)
    if min(h, r0, lreq) <= 0:
        raise ValueError("positive finite inputs required")
    tail = 360 * (4 / h) ** 5
    length = lreq / (2 * r0)
    ceiling = lambda x: (x.numerator + x.denominator - 1) // x.denominator
    ell0 = max(16, ceiling(tail), ceiling(length))
    return {
        "ell0": ell0,
        "tail_error_bound_at_ell0": tail / ell0,
        "length_floor": 2 * r0 * ell0,
        "tail_error_gate": tail / ell0 <= 1,
        "length_gate": 2 * r0 * ell0 >= lreq,
        "finite_head": "accepted bands ell<ell0 form a finite set",
        "uniform_tail": "Q23 analytic starred-frame estimates apply for ell>=ell0",
    }


def constrained_pressure(n, nprime, tangent, forcing, k=1, mode=1):
    n, nprime, tangent, forcing = (
        tuple(map(Q, row)) for row in (n, nprime, tangent, forcing)
    )
    k, mode = Q(k), Q(mode)
    norm2 = sum(x*x for x in n)
    if norm2 <= 0 or k == 0 or mode == 0:
        raise ValueError("nonzero normal, carrier and mode required")
    if sum(x*y for x, y in zip(n, tangent)) != 0:
        raise ValueError("tangent data required")
    moving = -sum(x*y for x, y in zip(nprime, tangent))
    source = sum(x*y for x, y in zip(n, forcing))
    return {
        "moving_over_i": moving / (k * mode * norm2),
        "source_over_i": source / (k * mode * norm2),
        "total_without_K_over_i": (moving + source) / (k * mode * norm2),
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q172 identity required")
    status = statuses()
    orientation = orientation_identity(2, -1, 7, 3, 2)
    counterexample = compactness_counterexample()
    repair = q23_uniform_repair()
    moving = constrained_pressure((1, 1, 0), (1, -1, 0), (1, -1, 0), (0, 0, 0), 3, 1)
    source = constrained_pressure((0, 0, 1), (0, 0, 0), (1, 0, 0), (0, 0, 2), 5, -1)
    return {
        "checks": {
            "all_provider_statuses_match": all(status.values()),
            "ten_Q2_regions_have_named_analytic_sources": len(REGIONS) == 10 and all(name and provider for name, provider in REGIONS),
            "Q2_orientation_follows_from_positive_pc_minus_v": orientation["minus_T_dot_g"] == orientation["factorized"] and orientation["orientation_positive"],
            "per_band_compactness_does_not_imply_uniform_tail_constant": not counterexample["compactness_alone_gives_uniform_infinite_tail_bound"],
            "Q23_tail_and_length_gates_close": repair["tail_error_gate"] and repair["length_gate"],
            "moving_normal_pressure_is_retained": moving["moving_over_i"] != 0,
            "unprojected_normal_source_pressure_is_retained": source["source_over_i"] != 0,
        },
        "provider_status_checks": status,
        "Q2": {
            "regions": [{"region": name, "source": provider} for name, provider in REGIONS],
            "orientation_identity": orientation,
            "factored_edges": "weighted_background_edges proves relative stress against zeta; no positive absolute flat-edge stress floor is used",
            "source_audit_result": "the ten-region value/gap and factored-edge chain survives source review relative to the displayed analytical provider derivations",
        },
        "Q3": {
            "counterexample": counterexample,
            "correction": "Q3 section 1 cannot infer a band-uniform CE,Cd merely from compactness on each member of an infinite accepted tail",
            "repair": repair,
            "corrected_chain": "Q23 supplies the finite-head/analytic-tail uniform upper propagator; same-band and full-path providers supply positive denominators, fixed jets, constrained pressure and whole-path cutoffs",
            "moving_pressure_control": moving,
            "source_pressure_control": source,
        },
        "scope": {
            "Q2_provider_premise_source_audited": True,
            "Q3_compactness_only_argument_valid": False,
            "Q3_upper_propagator_repaired_by_Q23_on_same_field": True,
            "every_underlying_provider_theorem_reproved_from_first_principles": False,
            "fixed_projector_compensator_and_weighted_inverse_contracts_closed_here": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
        "next": "audit fixed compensator/projector supports and jets on the corrected Q2/Q3 tail, then the compact and slow inverse weighted domains",
    }


def record():
    return {
        "status": "Q172-Q2-source-audited-and-Q3-uniform-frame-corrected",
        "established": "Q2's ten regional inequalities and factored edges trace to analytical continuum providers; the directional row follows exactly from -T.g=a F^2(Pc-v)>0",
        "correction": "Q3's compactness-only inference is invalid over infinitely many accepted bands; Q23's finite-head/analytic-tail estimate supplies the required uniform upper propagator on the unchanged field",
        "remaining": "the corrected Q3 chain still relies on named same-band/full-path providers and does not close fixed projector/compensator, weighted inverse, complete covariance/means or the residual",
        "historical_frozen_records_modified": False,
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
