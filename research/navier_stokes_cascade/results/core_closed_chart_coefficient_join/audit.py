#!/usr/bin/env python3
"""Freeze the closed-chart coefficient compactness join."""
import argparse
import hashlib
import json
from pathlib import Path

from compactness import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_corner_tail_slack_reconciliation",
    "weighted_background_edges",
    "physical_pulse_transfer",
    "matching_root_foundation",
    "core_compact_slow_inverse_domain_join",
    "core_first_real_high_jet_domain_closure",
    "core_certified_finite_stage_domain_join",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dependency_review():
    records, sources = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent / name / "evidence.json"
        record = json.loads(path.read_text())
        records[str(path.relative_to(PROJECT))] = sha(path)
        for rel, digest in record.get("source_hashes", {}).items():
            target = PROJECT / rel
            sources[name + ":" + rel] = target.exists() and sha(target) == digest
    return records, sources


def generate():
    frozen, inherited = dependency_review()
    result = review()
    own = [HERE / name for name in (
        "README.md", "compactness.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    scope = result["scope"]
    gates = {
        "eight_exact_controls": len(result["checks"]) == 8 and all(result["checks"].values()),
        "seven_frozen_dependencies": len(frozen) == len(DEPENDENCIES),
        "all_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "conditional_scope_preserved": (
            scope["closed_chart_coefficient_join_relative_to_primitive_contracts"]
            and not scope["actual_exact_core_background_and_primitive_margins_independently_proved_here"]
            and not scope["unconditional_complete_residual_compact_force_or_blowup"]
        ),
    }
    return {
        "status": "closed-chart-coefficient-finiteness-joined-relative-to-certified-primitive-contracts",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": gates,
        "dependency_evidence_sha256": frozen,
        "inherited_source_checks": inherited,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "independently verify the exact leading/background shear, growth, cone and all-fixed-order closed-chart seminorms from their defining profiles and Borel sum",
            "independently verify the global support-separation, root-input and compact-input margins used by the primitive contracts on the unchanged common tail",
            "then replay the now-joined coefficient, logarithmic-slack, Cartesian-corner, localization, force and comparison implications on the same field",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Closed-chart coefficient join, provenance or scope changed")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record = generate()
    validate(record)
    if args.verify_record:
        validate(json.loads(args.verify_record.read_text()))
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print("closed-chart coefficient join:", len(record["review"]["checks"]),
          "exact controls;", len(record["inherited_source_checks"]),
          "dependency source checks")


if __name__ == "__main__":
    main()
