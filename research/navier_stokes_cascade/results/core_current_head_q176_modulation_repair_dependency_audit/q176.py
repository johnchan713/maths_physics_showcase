"""Source-level compact modulation/repair audit and dependency correction."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q176-MODULATION-REPAIR-DEPENDENCY-AUDIT"
PROVIDERS = (
    "core_current_head_q175_canonical_source_route_reconciliation",
    "prescribed_continuation_foundation",
    "loop_modulation_bounds",
    "repair_state_bounds",
    "regional_cone_join",
    "modulation_repair_all_order",
    "annular_exterior_foundation",
    "matching_root_foundation",
    "compact_global_reconciliation",
    "core_pressure_identification",
    "foundation_bundle_route",
)
EXPECTED = {
    "core_current_head_q175_canonical_source_route_reconciliation": "Q175-canonical-typed-source-route-reconciled-and-independent-frontier-restored",
    "prescribed_continuation_foundation": "pre-modulation-prescribed-continuation-through-shape-all-fixed-jets",
    "loop_modulation_bounds": "fixed-phase-loop-and-modulation-errors-bounded",
    "repair_state_bounds": "repair-state-bound-and-conditional-frequency-acceptance",
    "regional_cone_join": "regional-input-gaps-and-leading-value-cone-joined",
    "modulation_repair_all_order": "same-finite-N-modulation-moving-phase-and-exact-repair-all-order-map-estimate",
    "annular_exterior_foundation": "annular-exterior-all-order-leading-regularity-and-global-moment-foundation-relative-to-exact-roots",
    "matching_root_foundation": "actual-annular-and-compact-matching-roots-audited-relative-to-retained-core-and-compact-input-estimates",
    "compact_global_reconciliation": "explicit-cutoff-to-same-root-compact-cone-and-global-pressure-reconciliation",
    "core_pressure_identification": "leading-core-pressure-identified-with-regular-angular-inverse-and-retained-higher-memory",
    "foundation_bundle_route": "foundation-route-closed-with-independent-continuum-validation-frontier",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q176 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


NODES = (
    "continuation_shape",
    "low_order_loop_repair",
    "regional_value_selection",
    "all_order_same_map",
    "annular_exterior_composition",
    "matching_root_reaudit",
    "compact_global_reconciliation",
    "pressure_identification",
)
EDGES = (
    ("continuation_shape", "regional_value_selection"),
    ("low_order_loop_repair", "regional_value_selection"),
    ("continuation_shape", "all_order_same_map"),
    ("regional_value_selection", "all_order_same_map"),
    ("all_order_same_map", "annular_exterior_composition"),
    ("regional_value_selection", "annular_exterior_composition"),
    ("all_order_same_map", "matching_root_reaudit"),
    ("annular_exterior_composition", "matching_root_reaudit"),
    ("regional_value_selection", "compact_global_reconciliation"),
    ("matching_root_reaudit", "compact_global_reconciliation"),
    ("annular_exterior_composition", "pressure_identification"),
    ("compact_global_reconciliation", "pressure_identification"),
)


def topological_order(nodes=NODES, edges=EDGES):
    if len(set(nodes)) != len(nodes):
        raise ValueError("duplicate node")
    known = set(nodes)
    if any(a not in known or b not in known for a, b in edges):
        raise ValueError("unknown endpoint")
    incoming = {node: set() for node in nodes}
    outgoing = {node: set() for node in nodes}
    for a, b in edges:
        incoming[b].add(a)
        outgoing[a].add(b)
    ready = [node for node in nodes if not incoming[node]]
    out = []
    while ready:
        node = ready.pop(0)
        out.append(node)
        for target in tuple(outgoing[node]):
            incoming[target].remove(node)
            if not incoming[target] and target not in out and target not in ready:
                ready.append(target)
        ready.sort(key=nodes.index)
    if len(out) != len(nodes):
        raise ValueError("cycle")
    return tuple(out)


def repair_gate(beta=1000, q=100000, H_floor=8):
    if any(type(v) is not int or v <= 0 for v in (beta, q, H_floor)):
        raise ValueError("positive integer gate data required")
    upper = Q(4 * beta * beta * q, H_floor ** 16)
    return {
        "beta": beta,
        "q": q,
        "H_floor": H_floor,
        "jacobian_perturbation_upper_bound": upper,
        "strict_half_gate": upper < Q(1, 2),
        "inverse_bound": 2 * beta,
    }


def cone_gate(G_floor=4):
    if type(G_floor) is not int or G_floor < 1:
        raise ValueError("positive integer G floor required")
    left = 512 * 4 ** 10
    right = G_floor ** 30
    return {
        "left": left,
        "right_at_floor": right,
        "strict_cone_absorption": left < right,
    }


def compact_moment_inverse(rows, lam):
    if len(rows) != 5:
        raise ValueError("five transformed rows required")
    lam = Q(lam)
    if lam <= 0:
        raise ValueError("positive lambda required")
    r1, r2, r3, r4, r5 = map(Q, rows)
    return (r1, r3, r1 - lam * r2, -r4, r5)


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q176 identity required")
    status = statuses()
    regional = evidence("regional_cone_join")
    maps = evidence("modulation_repair_all_order")
    exterior = evidence("annular_exterior_foundation")
    matching = evidence("matching_root_foundation")
    compact = evidence("compact_global_reconciliation")
    pressure = evidence("core_pressure_identification")
    foundation = evidence("foundation_bundle_route")
    rg, mf, ef, xf, pf, ff = (regional["flags"], maps["flags"], exterior["flags"],
                              matching["flags"], pressure["flags"], foundation["flags"])
    order = topological_order()
    pos = {node: i for i, node in enumerate(order)}
    repair = repair_gate()
    sharp_repair = repair_gate(beta=800, q=25000)
    cone = cone_gate()
    zero_defects = compact_moment_inverse((0, 0, 0, 0, 0), Q(1, 10 ** 6))
    coarse = foundation["reviews"]["exact_route"]["topological_order"]
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "regional_selection_precedes_all_order_same_map_theorem": pos["regional_value_selection"] < pos["all_order_same_map"],
        "all_order_map_precedes_annular_exterior_composition": pos["all_order_same_map"] < pos["annular_exterior_composition"],
        "later_matching_root_reaudit_is_not_an_input_to_all_order_map": pos["all_order_same_map"] < pos["matching_root_reaudit"],
        "refined_route_is_acyclic": all(pos[a] < pos[b] for a, b in EDGES),
        "coarse_foundation_route_conflates_two_matching_and_modulation_levels": coarse.index("matching_roots") < coarse.index("modulation_repair") and pos["matching_root_reaudit"] > pos["all_order_same_map"],
        "same_selected_low_order_map_is_accepted": rg["same_finite_frequency_leading_value_cone_joined"] and mf["low_order_cone_acceptance_identifies_same_map_conditionally"],
        "conservative_all_order_Jacobian_gate_is_strict": repair["strict_half_gate"] and repair["inverse_bound"] == 2000,
        "later_matching_reaudit_sharpens_the_same_gate": sharp_repair["strict_half_gate"] and sharp_repair["inverse_bound"] == 1600,
        "cone_error_is_absorbed_at_a_very_weak_G_floor": cone["strict_cone_absorption"],
        "all_order_derivatives_use_same_root_without_Ck_reselection": mf["same_accepted_C1_repair_has_all_fixed_derivatives_by_pointwise_Jacobian"] and not mf["simultaneous_Ck_contractions_or_frequency_reselection_required"],
        "whole_compact_input_is_supplied_in_the_composition_checkpoint": ef["whole_compact_smooth_input_to_same_finite_N_maps_supplied_relative_to_exact_core_and_finite_existence_bounds"],
        "later_matching_reaudit_retains_all_five_rows_and_pressure": xf["all_five_moments_and_complete_pressure_memory_retained"],
        "zero_transformed_rows_are_exactly_zero_original_defects": zero_defects == (Q(0),) * 5,
        "global_reconciliation_keeps_same_frequency_root_and_pressure": compact["gates"]["same_N_exact_repair_bridge_available"] and compact["gates"]["global_five_moments_and_pressure_memory_available"],
        "leading_pressure_identity_keeps_positive_order_memory": pf["leading_pressure_identification_already_supplied_by_same_root_ledger"] and pf["higher_background_Pi1_Pi2_and_complete_canonical_pressure_memory_retained"],
        "independent_continuum_certificate_remains_open": not ff["formal_or_external_independent_certificate"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "refined_route": {"nodes": NODES, "edges": EDGES, "topological_order": order},
        "repair_gate": {k: str(v) if isinstance(v, Q) else v for k, v in repair.items()},
        "sharpened_matching_reaudit_gate": {k: str(v) if isinstance(v, Q) else v for k, v in sharp_repair.items()},
        "cone_gate": cone,
        "correction": {
            "foundation_route_matching_roots_node_was_too_coarse": True,
            "low_order_repair_acceptance_and_later_matching_root_reaudit_are_distinct": True,
            "regional_cone_selection_precedes_the_all_order_same_map_theorem": True,
            "annular_exterior_checkpoint_proves_the_input_then_composes_the_conditional_map": True,
            "historical_frozen_evidence_modified": False,
        },
        "established": {
            "conditional_same_map_all_order_modulation_and_exact_repair": True,
            "exact_five_row_pressure_memory": True,
            "inner_compact_outer_leading_cone_join": True,
            "radial_field_errors_uniformly_small_at_all_derivative_orders": False,
        },
        "scope": {
            "actual_selected_continuum_input_suprema_independently_recomputed_here": False,
            "finite_order_existence_inverse_and_gap_inputs_discharged_unconditionally": False,
            "uniform_actual_pulse_or_complete_residual_proved": False,
            "unconditional_force_exclusion_or_blowup_proved": False,
        },
        "next": "independently validate the retained finite-order continuation, low-order loop/repair and regional input suprema on the selected profile, then audit the weighted first-stress edge quotient and reserved supports; do not treat the later matching-root reaudit as an upstream premise",
    }


def record():
    return {
        "status": "Q176-compact-modulation-exact-repair-audited-with-refined-dependency-order",
        "established": "the fixed-N moving-phase modulation and exact five-row repair form a valid conditional same-map all-order implication; the strict Jacobian and cone absorption inequalities close and full pressure memory joins inner compact and exterior leading cones",
        "correction": "the coarse foundation graph placed matching_root_foundation upstream of modulation, but that frozen matching record is a later reaudit; the actual upstream inputs are the low-order loop/repair bounds and regional selected-map acceptance",
        "obstruction": "actual selected-profile continuation and low-order regional input suprema remain retained continuum premises; all-order finiteness does not make radial errors N-uniform or independently certify those inputs",
        "historical_frozen_records_modified": False,
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
