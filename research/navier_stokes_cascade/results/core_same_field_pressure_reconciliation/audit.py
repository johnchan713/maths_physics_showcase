#!/usr/bin/env python3
"""Freeze the same-field core/global leading-pressure reconciliation."""
import argparse
import hashlib
import json
from pathlib import Path

from identities import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_certified_complete_potential_comparison_join",
    "core_pressure_identification",
    "core_fixed_point_ledger_recheck",
    "core_modulation_repair_join_recheck",
    "core_global_pressure_moment_recheck",
    "core_leading_background_borel_join",
    "core_certified_complete_potential_endpoint_join",
    "core_certified_complete_potential_force_join",
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
    protocol = json.loads((HERE / "protocol.json").read_text())
    own = [HERE / name for name in (
        "README.md", "identities.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    gates = {
        "seven_exact_controls": len(result["checks"]) == 7 and all(result["checks"].values()),
        "eight_frozen_dependencies": len(frozen) == len(DEPENDENCIES),
        "all_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "order_and_scope_preserved": (
            result["scope"]["leading_q_zero_pressure_identified_on_same_field"]
            and not result["scope"]["full_q_positive_pressure_identified_with_core_pressure"]
            and not result["scope"]["fixed_jet_closed_corner_residual_bound_proved_here"]
            and not result["scope"]["unconditional_force_residual_or_blowup_proved"]
        ),
    }
    return {
        "status": "same-field-core-pressure-identified-with-completed-leading-backward-pressure",
        "protocol": protocol,
        "review": result,
        "gates": gates,
        "dependency_evidence_sha256": frozen,
        "inherited_source_checks": inherited,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "prove for every fixed spacetime jet the closed-chart residual bound through |eta|=1 and T=0, including the extra time derivative used by localization",
            "make that corner flatness uniform on each compact X-region before promoting the compact force and global comparison",
            "retain q-positive background, wave, mean and localization pressure in the complete residual and force",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Same-field pressure reconciliation, provenance or scope changed")


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
    print("same-field pressure reconciliation:", len(record["review"]["checks"]),
          "exact controls;", len(record["inherited_source_checks"]),
          "dependency source checks")


if __name__ == "__main__":
    main()
