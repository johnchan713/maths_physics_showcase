#!/usr/bin/env python3
"""Freeze the closed-corner flatness reduction and exact obstruction."""
import argparse
import hashlib
import json
from pathlib import Path

from corner import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_same_field_pressure_reconciliation",
    "core_certified_complete_potential_diagonal_join",
    "core_certified_complete_potential_endpoint_join",
    "core_certified_complete_potential_force_join",
    "canonical_local_realization",
    "current_endpoint_jets",
    "force_foundation",
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
        "README.md", "corner.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    gates = {
        "eight_exact_controls": len(result["checks"]) == 8 and all(result["checks"].values()),
        "seven_frozen_dependencies": len(frozen) == len(DEPENDENCIES),
        "all_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "obstruction_scope_preserved": (
            result["scope"]["coordinate_transfer_proved"]
            and not result["scope"]["closed_chart_constants_proved_finite_for_actual_field"]
            and not result["scope"]["compact_force_and_comparison_unconditional"]
            and not result["scope"]["blowup_proved"]
        ),
    }
    return {
        "status": "closed-corner-transfer-proved-and-uniform-constant-obstruction-isolated",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": gates,
        "dependency_evidence_sha256": frozen,
        "inherited_source_checks": inherited,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "prove finite fixed-stage residual constants on the full closed chart uniformly for every compact X interval",
            "prove the same closed-chart transition-uniform constants for every complete-potential cutoff increment before scale absorption",
            "then transfer arbitrary q powers to the Cartesian corner with one extra time derivative and replay localization, force and comparison",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Corner-flatness obstruction, provenance or scope changed")


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
    print("corner-flatness obstruction:", len(record["review"]["checks"]),
          "exact controls;", len(record["inherited_source_checks"]),
          "dependency source checks")


if __name__ == "__main__":
    main()
