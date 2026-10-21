#!/usr/bin/env python3
"""Freeze the repaired-profile global pressure and moment recheck."""

import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from proof_checks import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "fbefd07f4fb1b3738a00429dd282ce6f4e5c715e"
STATUS = "same-repaired-profile-global-pressure-moments-and-exterior-stress-identities-independently-rechecked"
DEPENDENCIES = (
    "core_modulation_repair_join_recheck",
    "annular_exterior_foundation",
    "pulse_moment_audit",
    "heat_exterior_audit",
    "global_moment_stress",
)

FLAGS = {
    "renormalized_heat_moment_and_endpoint_jets_rederived": True,
    "finite_endpoint_I_to_global_R_conversion_condition_made_explicit": True,
    "same_repaired_profile_all_five_global_defects_rederived": True,
    "same_gauge_backward_pressure_identity_rederived": True,
    "zero_defect_exterior_velocity_and_stress_join_rederived": True,
    "complete_pressure_memory_and_eta_derivative_retained": True,
    "coefficient_and_reserved_mean_supports_rederived_here": False,
    "uniform_actual_pulse_covariance_residual_or_blowup_proved": False,
    "formal_proof_assistant_or_external_peer_certificate": False,
}

REMAINING = [
    "independently audit finite coefficient solves, their common supports, and untouched reserved velocity and exceptional pressure means",
    "prove uniform actual-pulse existence, inverse, gap, value, derivative and support estimates at every physical scale",
    "join complete wave covariance and mean corrections including radial and axial velocity, curl remainders and both exceptional pressure tails",
    "prove the full residual gain including viscosity, every cutoff and every interaction before any force or blow-up conclusion",
]


def protocol():
    return json.loads((HERE / "protocol.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dependency_hashes():
    return {
        str((HERE.parent / name / "evidence.json").relative_to(PROJECT)):
            sha(HERE.parent / name / "evidence.json")
        for name in DEPENDENCIES
    }


def source_hashes():
    paths = [HERE / name for name in (
        "README.md", "proof_checks.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    paths += [
        HERE.parent / name / filename
        for name in DEPENDENCIES
        for filename in ("README.md", "evidence.json")
    ]
    return {str(path.relative_to(PROJECT)): sha(path) for path in paths}


def serializable(value):
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(item) for item in value]
    if isinstance(value, Fraction):
        return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    return value


def generate():
    raw = review()
    gates = {
        "all_exact_global_pressure_moment_and_stress_checks_pass": all(raw["checks"].values()),
        "all_frozen_dependencies_bound": len(dependency_hashes()) == len(DEPENDENCIES),
        "scope_flags_are_strict_booleans": all(type(value) is bool for value in FLAGS.values()),
        "coefficient_and_physical_frontier_retained": len(REMAINING) == 4,
    }
    return {
        "status": STATUS,
        "parent_commit": PARENT,
        "protocol": protocol(),
        "flags": deepcopy(FLAGS),
        "review": serializable(raw),
        "gates": gates,
        "dependency_evidence_sha256": dependency_hashes(),
        "source_hashes": source_hashes(),
        "remaining_obligations": REMAINING[:],
    }


def validate(record):
    expected = generate()
    if json.dumps(record, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise ValueError("Global pressure/moment recheck record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Global pressure/moment recheck gate failed")


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
    print(f"{STATUS}: five global defects and same-gauge exterior stress identities rechecked.")


if __name__ == "__main__":
    main()
