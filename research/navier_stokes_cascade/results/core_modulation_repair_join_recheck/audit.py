#!/usr/bin/env python3
"""Freeze the same-profile modulation, repair, and cone-join recheck."""

import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from proof_checks import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "f9c45658e75a6f0a8122875e383185679f4bec55"
STATUS = "same-profile-compact-modulation-repair-and-regional-cone-join-independently-rechecked"
DEPENDENCIES = (
    "core_continuation_matching_recheck",
    "explicit_variance_cutoff",
    "loop_modulation_bounds",
    "repair_state_bounds",
    "regional_cone_join",
    "annular_exterior_foundation",
)

FLAGS = {
    "signed_variance_root_through_zero_rederived": True,
    "moving_phase_inverse_and_slow_fast_chain_rule_rederived": True,
    "same_frequency_modulation_and_extra_moment_gain_rederived": True,
    "same_compact_five_row_root_applied_with_complete_pressure_memory": True,
    "inner_compact_annular_exterior_value_cone_join_rederived": True,
    "original_width_parameters_roots_frequency_and_gauge_unchanged": True,
    "complete_global_pressure_moment_tail_identity_rederived_here": False,
    "uniform_actual_pulse_covariance_residual_or_blowup_proved": False,
    "formal_proof_assistant_or_external_peer_certificate": False,
}

REMAINING = [
    "independently retrace the global renormalized heat moment and complete pressure/moment identities on this repaired profile",
    "audit finite coefficients, common supports, and untouched reserved velocity and exceptional pressure means",
    "prove uniform actual-pulse value, derivative, inverse, gap and support estimates at every physical scale",
    "join complete covariance means including axial velocity and pressure, then bound viscosity, all cutoffs and all interactions before any blow-up claim",
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
        "all_exact_modulation_repair_and_join_checks_pass": all(raw["checks"].values()),
        "all_frozen_dependencies_bound": len(dependency_hashes()) == len(DEPENDENCIES),
        "scope_flags_are_strict_booleans": all(type(value) is bool for value in FLAGS.values()),
        "global_and_physical_frontier_retained": len(REMAINING) == 4,
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
        raise ValueError("Modulation/repair/cone-join recheck record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Modulation/repair/cone-join recheck gate failed")


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
    print(f"{STATUS}: same N/root applied and leading regional value cone joined.")


if __name__ == "__main__":
    main()
