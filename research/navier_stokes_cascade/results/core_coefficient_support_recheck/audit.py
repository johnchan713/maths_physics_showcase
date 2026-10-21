#!/usr/bin/env python3
"""Freeze the coefficient-support join and reserved-pressure obstruction."""

import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "same-profile-finite-coefficient-supports-and-reserved-pressure-obstruction-rechecked"
DEPENDENCIES = (
    "core_global_pressure_moment_recheck", "coefficient_support_foundation",
    "coefficient_support_audit", "annular_exterior_foundation", "weighted_stress_foundation",
)
FLAGS = {
    "same_leading_global_identities_joined_to_finite_coefficient_supports": True,
    "ordinary_linear_rows_and_normalization_costs_rechecked": True,
    "shifted_viscosity_and_two_exceptional_pressure_orders_retained": True,
    "sharp_reserved_pressure_error_floor_derived": True,
    "stationary_piece_and_complete_cutoff_derivatives_derived": True,
    "finite_inner_volterra_existence_lemmas_remain_dependencies": True,
    "manufactured_numerical_controls_are_continuum_proofs": False,
    "uniform_actual_pulses_complete_covariance_residual_or_blowup_proved": False,
    "formal_proof_assistant_certificate": False,
}
REMAINING = [
    "prove uniform actual-pulse existence, inverse, strict-gap, value, derivative and support bounds at every physical scale with selected cutoff and tilt",
    "join complete covariance means including radial and axial velocity, curl remainders, both exceptional pressure terms and all cutoff derivatives",
    "prove complete residual gain including viscosity, every cutoff and every interaction before a force or blow-up conclusion",
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serial(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {key: serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(item) for item in value]
    return value


def dependency_review():
    bound, checks = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent / name / "evidence.json"
        record = json.loads(path.read_text())
        bound[str(path.relative_to(PROJECT))] = digest(path)
        for rel, expected in record["source_hashes"].items():
            checks[name + ":" + rel] = digest(PROJECT / rel) == expected
    return bound, checks


def sources():
    files = [HERE / n for n in ("README.md", "audit.py", "proof_checks.py", "test_audit.py", "protocol.json", "requirements.txt", ".gitignore")]
    files += [HERE.parent / name / "README.md" for name in DEPENDENCIES]
    files += [HERE.parent / p for p in (
        "matching_root_foundation/matching_bounds.py", "root_foundation/root_bounds.py",
        "coefficient_support_foundation/conservative_checks.py",
    )]
    return {str(path.relative_to(PROJECT)): digest(path) for path in files}


def generate():
    data = review()
    bound, inherited = dependency_review()
    return {
        "status": STATUS,
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "flags": deepcopy(FLAGS),
        "review": serial(data),
        "gates": {
            "all_exact_checks_pass": all(data["checks"].values()),
            "all_declared_dependency_source_hashes_match": bool(inherited) and all(inherited.values()),
            "five_frozen_dependency_records_bound": len(bound) == 5,
        },
        "dependency_evidence_sha256": bound,
        "inherited_source_checks": inherited,
        "source_hashes": sources(),
        "remaining_obligations": REMAINING[:],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Coefficient support, pressure, provenance or scope record changed")


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
    with args.output.open("x") as output:
        output.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(f"{STATUS}: exact support/pressure controls and {len(record['inherited_source_checks'])} inherited source bindings pass")


if __name__ == "__main__":
    main()
