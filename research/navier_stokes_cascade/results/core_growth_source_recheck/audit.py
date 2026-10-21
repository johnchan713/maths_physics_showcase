#!/usr/bin/env python3
"""Freeze the independent angular-source and two-sector exit recheck."""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "3476acc9ccdb5921e20964f86180128fac094ad6"
STATUS = "core-angular-source-and-two-sector-exit-independently-rederived"
DEPENDENCIES = (
    "core_growth_activation",
    "core_fixed_point_ledger_recheck",
    "core_pressure_gap_recheck",
)

FLAGS = {
    "alternating_scalar_bounds_independently_rederived": True,
    "weighted_eta_reciprocal_bound_keeps_Hstar_factor": True,
    "positive_source_without_division_by_chi_rederived": True,
    "forward_shear_and_angular_exit_growth_rederived": True,
    "axial_exit_pressure_growth_with_physical_scaling_rederived": True,
    "flat_activation_volterra_and_moment_bounds_rederived_here": False,
    "continuation_matching_modulation_or_exterior_rederived_here": False,
    "unconditional_pulse_residual_force_or_blowup_proved": False,
    "formal_proof_assistant_or_external_peer_certificate": False,
}

REMAINING = [
    "independently rederive the flat activation Volterra bounds and complete pressure-moment factors at the unchanged width",
    "continue through prescribed continuation and matching roots on the same selected profile",
    "audit compact modulation and exact repair with full pressure memory and join the regional cone estimates",
    "retain all later coefficient-support, actual-pulse, covariance-mean and complete residual obligations",
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


def generate():
    exact = review()
    gates = {
        "all_exact_scalar_source_and_exit_checks_pass": all(exact["checks"].values()),
        "all_frozen_dependencies_bound": len(dependency_hashes()) == len(DEPENDENCIES),
        "scope_flags_are_strict_booleans": all(type(value) is bool for value in FLAGS.values()),
        "flat_activation_and_later_frontier_retained": len(REMAINING) == 4,
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
        raise ValueError("Core growth-source recheck record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Core growth-source recheck gate failed")


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
    print(f"{STATUS}: source >= 2.5 + .95 Lambda L chi; native exit v > 2.3.")


if __name__ == "__main__":
    main()
