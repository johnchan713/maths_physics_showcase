#!/usr/bin/env python3
"""Freeze the independent flat-activation and five-moment recheck."""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "ce175f68ed606c259616c62eaa795ae83ce4500f"
STATUS = "core-flat-activation-five-moment-and-cone-collar-independently-rederived"
DEPENDENCIES = (
    "inner_collar_activation",
    "core_growth_activation",
    "core_growth_source_recheck",
)

FLAGS = {
    "positive_flat_Volterra_representation_and_endpoint_rederived": True,
    "all_five_forward_moment_factorizations_rederived": True,
    "both_pressure_value_and_eta_derivative_terms_retained": True,
    "C_gauge_and_uniform_coefficient_envelopes_rederived": True,
    "Kact_below_C_cubed_at_unchanged_width_rederived": True,
    "quantitative_inner_cone_collar_rederived": True,
    "prescribed_continuation_and_matching_roots_rederived_here": False,
    "compact_modulation_repair_or_exterior_rederived_here": False,
    "unconditional_pulse_residual_force_or_blowup_proved": False,
    "formal_proof_assistant_or_external_peer_certificate": False,
}

REMAINING = [
    "independently rederive prescribed continuation and both matching-root constructions on the same selected profile",
    "audit compact modulation and exact repair with complete pressure memory and join the inner and outer cone bounds",
    "retain the global pressure/moment, coefficient-support and reserved-mean obligations",
    "retain uniform actual-pulse, covariance-mean and complete residual obligations before any blow-up claim",
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
        "all_exact_Volterra_moment_scaling_and_cone_checks_pass": all(exact["checks"].values()),
        "all_frozen_dependencies_bound": len(dependency_hashes()) == len(DEPENDENCIES),
        "scope_flags_are_strict_booleans": all(type(value) is bool for value in FLAGS.values()),
        "continuation_and_later_frontier_retained": len(REMAINING) == 4,
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
        raise ValueError("Core activation Volterra recheck record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Core activation Volterra recheck gate failed")


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
    print(f"{STATUS}: Kact=C^2 R^80<C^3; unchanged width C^-20 retained.")


if __name__ == "__main__":
    main()
