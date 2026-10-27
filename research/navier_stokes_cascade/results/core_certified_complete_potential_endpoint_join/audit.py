#!/usr/bin/env python3
"""Freeze the certified-diagonal punctured endpoint join."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from controls import review
from endpoint import closure_statement


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_certified_complete_potential_diagonal_join",
    "core_finite_stage_compact_recurrence",
    "current_endpoint_jets",
    "endpoint_foundation",
    "physical_pulse_transfer",
    "core_finite_support_closure",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, list):
        return [encode(item) for item in value]
    return value


def generate():
    frozen, inherited = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent / name / "evidence.json"
        record = json.loads(path.read_text())
        frozen[str(path.relative_to(PROJECT))] = sha(path)
        for rel, digest in record.get("source_hashes", {}).items():
            target = PROJECT / rel
            inherited[name + ":" + rel] = target.exists() and sha(target) == digest
    result, scope = review(), closure_statement()
    sources = [HERE / name for name in (
        "README.md", "endpoint.py", "controls.py", "audit.py",
        "test_audit.py", "protocol.json", "requirements.txt", ".gitignore",
    )]
    sources += [
        HERE.parent / "core_certified_complete_potential_diagonal_join/diagonal_join.py",
        HERE.parent / "endpoint_foundation/terminal_bounds.py",
        HERE.parent / "current_endpoint_jets/endpoint.py",
    ]
    gates = dict(
        exact_controls=all(result["checks"].values()),
        six_frozen_records=len(frozen) == 6,
        historical_sources=all(inherited.values()),
        conditional_scope=(
            scope["same_certified_complete_diagonal_field_used"]
            and scope["certified_fixed_stage_domains_retained"]
            and not scope["stage_uniform_constant_bound_required"]
            and scope["compatible_one_sided_C_infinity_potential_direct_swirl_velocity_and_pressure_jets_on_punctured_terminal_compacts_follow_conditionally"]
            and not scope["terminal_origin_regularized"]
            and not scope["primitive_all_order_endpoint_contracts_independently_proved"]
            and not scope["localization_initialization_admissible_force_energy_global_comparison_or_blowup_proved"]
        ),
    )
    return dict(
        status="certified-domain-conditional-complete-potential-punctured-endpoint-join",
        review=result,
        closure=scope,
        gates=gates,
        dependency_evidence_sha256=frozen,
        inherited_source_checks=inherited,
        source_hashes={str(path.relative_to(PROJECT)): sha(path) for path in sources},
        remaining_obligations=[
            "join localization and zero initialization to this same endpoint field without resetting gauges or cutoffs",
            "then construct the compact smooth force and verify whole-space energy and global comparison on that field",
            "retain the singular origin and every unresolved continuum premise before any blow-up claim",
        ],
    )


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Certified complete-potential endpoint join, provenance or scope changed")


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
    args.output.write_text(json.dumps(encode(record), indent=2, sort_keys=True) + "\n")
    print("certified-diagonal endpoint join:", len(record["review"]["checks"]),
          "exact controls;", len(record["inherited_source_checks"]),
          "historical source checks")


if __name__ == "__main__":
    main()
