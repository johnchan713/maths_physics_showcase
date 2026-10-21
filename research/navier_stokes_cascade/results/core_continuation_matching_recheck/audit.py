#!/usr/bin/env python3
"""Freeze the continuation and same-profile matching-root recheck."""

import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from proof_checks import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "96865d42b7977ee08cc27ae886d6f07aff0acb8e"
STATUS = "core-continuation-and-both-same-profile-matching-roots-independently-rechecked"
DEPENDENCIES = (
    "core_activation_volterra_recheck",
    "prescribed_continuation_foundation",
    "shape_join_pressure",
    "matching_root_foundation",
    "annular_exterior_foundation",
    "global_pressure_moment_foundation",
)

FLAGS = {
    "post_taper_continuation_slopes_corrected_and_rederived": True,
    "all_five_forward_moments_and_complete_pressure_memory_retained": True,
    "annular_root_contraction_and_same_root_inclusion_rederived": True,
    "outer_prefix_supplies_same_profile_to_compact_repair_patch": True,
    "compact_root_contraction_and_same_root_inclusion_rederived": True,
    "both_transformed_zero_sets_keep_the_fifth_pressure_row": True,
    "compact_modulation_and_exact_repair_application_rederived_here": False,
    "inner_compact_outer_cone_join_rederived_here": False,
    "unconditional_pulse_residual_force_or_blowup_proved": False,
    "formal_proof_assistant_or_external_peer_certificate": False,
}

REMAINING = [
    "independently audit compact modulation, moving phase and application of the same exact repair on this fixed profile",
    "join the repaired compact cone bounds to the independently checked inner and annular/exterior bounds",
    "retrace global pressure/moment identities, coefficient supports and untouched reserved means",
    "prove uniform actual-pulse and covariance-mean estimates and the complete viscosity/cutoff/interaction residual before any blow-up claim",
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
    exact = serializable(raw)
    gates = {
        "all_exact_continuation_matching_and_pressure_checks_pass": all(raw["checks"].values()),
        "all_frozen_dependencies_bound": len(dependency_hashes()) == len(DEPENDENCIES),
        "scope_flags_are_strict_booleans": all(type(value) is bool for value in FLAGS.values()),
        "compact_and_later_frontier_retained": len(REMAINING) == 4,
    }
    return {
        "status": STATUS,
        "parent_commit": PARENT,
        "protocol": protocol(),
        "flags": deepcopy(FLAGS),
        "review": exact,
        "gates": gates,
        "dependency_evidence_sha256": dependency_hashes(),
        "source_hashes": source_hashes(),
        "remaining_obligations": REMAINING[:],
    }


def validate(record):
    expected = generate()
    if json.dumps(record, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise ValueError("Continuation/matching-root recheck record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Continuation/matching-root recheck gate failed")


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
    print(f"{STATUS}: annular and compact roots are the same roots on one retained profile.")


if __name__ == "__main__":
    main()
