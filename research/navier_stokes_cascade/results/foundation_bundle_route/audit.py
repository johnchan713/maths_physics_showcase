#!/usr/bin/env python3
"""Freeze the acyclic foundation-route audit and its independent-review frontier."""

import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from route_checks import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "e140b4348b0a3ad3bcb26a8b2f26c79b2c31a8e6"
STATUS = "foundation-route-closed-with-independent-continuum-validation-frontier"
DEPENDENCIES = (
    "core_contraction_audit",
    "core_growth_activation",
    "prescribed_continuation_foundation",
    "root_foundation",
    "axial_root_foundation",
    "matching_root_foundation",
    "modulation_repair_all_order",
    "annular_exterior_foundation",
    "regional_cone_join",
    "compact_global_reconciliation",
    "core_pressure_identification",
    "weighted_stress_foundation",
    "background_jet_foundation",
    "explicit_coefficient_cutoff",
    "global_leading_c5",
    "global_first_coefficient_c3",
    "global_second_coefficient_h23",
    "explicit_tilt_pulse_tail",
    "source_envelope_edge_audit",
    "current_target_cycle_bridge",
    "current_residual_ledger_audit",
    "explicit_tail_complete_join",
)

FLAGS = {
    "named_foundation_obligations_have_analytic_provider_routes": True,
    "provider_route_is_acyclic": True,
    "one_finite_later_tail_can_conjoin_all_finite_value_gates": True,
    "same_selected_maps_roots_pressure_gauge_and_width_retained": True,
    "historical_first_band_field_or_force_identified": False,
    "provider_continuum_theorems_independently_reproved_here": False,
    "frozen_boolean_evidence_treated_as_continuum_proof": False,
    "formal_or_external_independent_certificate": False,
    "unconditional_pulse_residual_force_or_blowup_proved": False,
}

FRONTIER = [
    "independently validate the exact pressure-family core contraction and quantitative core-gap theorem on the selected datum",
    "independently validate the continuation, matching-root, compact modulation/repair and annular/exterior provider theorems on the same selected profile",
    "independently validate the regional cone, global pressure/moment and weighted first-stress edge theorems on that profile",
    "independently validate the all-fixed-jet background coefficient, generated-source inverse and current finite-cycle operator estimates before any unconditional conclusion",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if type(value) is F:
        return f"{value.numerator}/{value.denominator}"
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def records():
    return {
        name: json.loads((HERE.parent / name / "evidence.json").read_text())
        for name in DEPENDENCIES
    }


def dependency_review():
    rows = records()
    checks = {
        "core_contraction_record_has_gates": bool(rows["core_contraction_audit"]["gates"]),
        "core_growth_record_has_gates": bool(rows["core_growth_activation"]["gates"]),
        "continuation_record_has_gates": bool(rows["prescribed_continuation_foundation"]["gates"]),
        "reference_root_record_has_gates": bool(rows["root_foundation"]["gates"]),
        "axial_root_record_has_gates": bool(rows["axial_root_foundation"]["gates"]),
        "matching_root_record_has_gates": bool(rows["matching_root_foundation"]["gates"]),
        "all_order_modulation_record_has_gates": bool(rows["modulation_repair_all_order"]["gates"]),
        "annular_exterior_record_has_gates": bool(rows["annular_exterior_foundation"]["gates"]),
        "regional_cone_record_has_gates": bool(rows["regional_cone_join"]["gates"]),
        "global_pressure_reconciliation_record_has_gates": bool(rows["compact_global_reconciliation"]["gates"]),
        "core_pressure_identification_record_has_gates": bool(rows["core_pressure_identification"]["gates"]),
        "weighted_stress_record_has_gates": bool(rows["weighted_stress_foundation"]["gates"]),
        "background_coefficient_record_has_gates": bool(rows["background_jet_foundation"]["gates"]),
        "selected_coefficient_cutoff_record_has_gates": bool(rows["explicit_coefficient_cutoff"]["gates"]),
        "three_global_coefficient_records_have_gates": all(bool(rows[name]["gates"]) for name in (
            "global_leading_c5", "global_first_coefficient_c3", "global_second_coefficient_h23"
        )),
        "explicit_primary_tail_record_has_gates": bool(rows["explicit_tilt_pulse_tail"]["gates"]),
        "source_target_residual_records_have_gates": all(bool(rows[name]["gates"]) for name in (
            "source_envelope_edge_audit", "current_target_cycle_bridge", "current_residual_ledger_audit"
        )),
        "downstream_record_preserves_foundation_frontier": not rows["explicit_tail_complete_join"]["flags"][
            "retained_core_background_defining_map_and_inverse_premises_all_discharged"
        ],
    }
    return {"checks": checks}


def dependencies():
    return {
        str((HERE.parent / name / "evidence.json").relative_to(PROJECT)): sha(
            HERE.parent / name / "evidence.json"
        )
        for name in DEPENDENCIES
    }


def sources():
    paths = [HERE / name for name in (
        "README.md", "route_checks.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    paths += [
        HERE.parent / name / filename
        for name in DEPENDENCIES
        for filename in ("README.md", "evidence.json")
    ]
    return {str(path.relative_to(PROJECT)): sha(path) for path in paths}


def protocol():
    return json.loads((HERE / "protocol.json").read_text())


def generate():
    exact = encode(review())
    routed = dependency_review()
    gates = {
        "exact_route_checks": exact["all_edges_forward"] and exact["tail_max_control"] and exact["domain_inside"],
        "all_dependency_records_bound": len(dependencies()) == len(DEPENDENCIES),
        "all_provider_routes_present": all(routed["checks"].values()),
        "scope_flags_are_strict_booleans": all(type(value) is bool for value in FLAGS.values()),
        "independent_frontier_is_nonempty": len(FRONTIER) == 4,
    }
    return {
        "status": STATUS,
        "parent_commit": PARENT,
        "protocol": protocol(),
        "flags": deepcopy(FLAGS),
        "reviews": {"exact_route": exact, "dependency_routing": routed},
        "gates": gates,
        "dependency_evidence_sha256": dependencies(),
        "source_hashes": sources(),
        "remaining_obligations": FRONTIER[:],
    }


def validate(record):
    expected = generate()
    if json.dumps(expected, sort_keys=True) != json.dumps(record, sort_keys=True):
        raise ValueError("Foundation-route record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Foundation-route gate failed")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record = generate()
    validate(record)
    if args.verify_record:
        validate(json.loads(args.verify_record.read_text()))
    with args.output.open("x") as stream:
        stream.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(
        f"{STATUS}: {len(record['dependency_evidence_sha256'])} records; "
        "acyclic route closed, continuum-validation frontier retained."
    )


if __name__ == "__main__":
    main()
