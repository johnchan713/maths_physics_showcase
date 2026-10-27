#!/usr/bin/env python3
"""Freeze the primitive-margin normalization reconciliation."""
import argparse
import hashlib
import json
from pathlib import Path

from margins import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_closed_chart_coefficient_join",
    "core_leading_background_borel_join",
    "core_same_band_gap_join",
    "core_canonical_full_path_ode_join",
    "core_fixed_profile_projector_join",
    "core_compact_slow_inverse_domain_join",
    "matching_root_foundation",
    "auxiliary_geometry_audit",
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
        "README.md", "margins.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    scope = result["scope"]
    gates = {
        "nine_exact_controls": len(result["checks"]) == 9 and all(result["checks"].values()),
        "eight_frozen_dependencies": len(frozen) == len(DEPENDENCIES),
        "all_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "corrected_scope_preserved": (
            scope["heterogeneous_scalar_Delta_prim_rejected"]
            and scope["scale_compatible_primitive_certificate_joined"]
            and not scope["unconditional_complete_residual_force_comparison_or_blowup"]
        ),
    }
    return {
        "status": "scale-compatible-primitive-margin-certificate-reconciled",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": gates,
        "dependency_evidence_sha256": frozen,
        "inherited_source_checks": inherited,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "replay the closed-chart coefficient, logarithmic-slack and Cartesian-corner implications using the corrected scale-compatible certificate",
            "then re-audit localization, compact force, energy and comparison on that same field without replacing weighted edge control by a positive pointwise floor",
            "retain complete pressure memory, all viscosity, cutoff and interaction terms and do not infer blow-up before the full chain is checked",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Primitive-margin reconciliation, provenance or scope changed")


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
    print("primitive-margin reconciliation:", len(record["review"]["checks"]),
          "exact controls;", len(record["inherited_source_checks"]),
          "dependency source checks")


if __name__ == "__main__":
    main()
