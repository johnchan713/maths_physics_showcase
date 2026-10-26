#!/usr/bin/env python3
"""Freeze the corrected complete-potential diagonal join and provenance."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from controls import review
from diagonal import closure_statement


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_finite_stage_compact_recurrence",
    "finite_stage_induction",
    "core_finite_support_closure",
    "canonical_local_realization",
    "realization_foundation",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, list):
        return [encode(item) for item in value]
    if isinstance(value, Q):
        return str(value)
    return value


def generate():
    frozen, inherited = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent / name / "evidence.json"
        record = json.loads(path.read_text())
        frozen[str(path.relative_to(PROJECT))] = sha(path)
        for rel, digest in record.get("source_hashes", {}).items():
            inherited[name + ":" + rel] = (PROJECT / rel).exists() and sha(PROJECT / rel) == digest
    result = review()
    sources = [HERE / name for name in (
        "README.md", "diagonal.py", "controls.py", "audit.py",
        "test_audit.py", "requirements.txt", ".gitignore",
    )]
    sources += [
        HERE.parent / "core_finite_stage_compact_recurrence/recurrence.py",
        HERE.parent / "finite_stage_induction/induction.py",
        HERE.parent / "canonical_local_realization/realization.py",
        HERE.parent / "realization_foundation/realized_checks.py",
    ]
    statement = closure_statement()
    gates = dict(
        exact_controls=all(result["checks"].values()),
        five_frozen_records=len(frozen) == 5,
        historical_sources=all(inherited.values()),
        conditional_scope=(
            statement["one_locally_finite_divergence_free_field_follows_conditionally"]
            and statement["all_fixed_jets_of_complete_local_residual_are_flat_conditionally"]
            and statement["protected_inner_velocity_and_heat_exterior_preserved_conditionally"]
            and not statement["primitive_analytic_contracts_independently_proved"]
            and not statement["closed_time_localization_forcing_energy_global_comparison_or_blowup_proved"]
        ),
    )
    return dict(
        status="conditional-corrected-complete-potential-diagonal-join",
        review=result,
        closure=statement,
        gates=gates,
        dependency_evidence_sha256=frozen,
        inherited_source_checks=inherited,
        source_hashes={str(path.relative_to(PROJECT)): sha(path) for path in sources},
    )


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Complete-potential diagonal join, provenance or scope changed")


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
    print("corrected complete-potential diagonal join:",
          len(record["review"]["checks"]), "exact controls;",
          len(record["inherited_source_checks"]), "historical source checks")


if __name__ == "__main__":
    main()
