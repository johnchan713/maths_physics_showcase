#!/usr/bin/env python3
"""Freeze the independent shared-core pressure and separation recheck."""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "0e605bb3b289ed662ea5333d0335baac8302a1ff"
STATUS = "shared-core-pressure-calculus-rederived-with-sharper-gap"
DEPENDENCIES = (
    "core_contraction_audit",
    "core_growth_activation",
    "core_pressure_identification",
    "foundation_bundle_route",
)

FLAGS = {
    "same_reference_measure_and_pressure_gauge_traced": True,
    "complex_pressure_input_bound_rederived": True,
    "real_pressure_derivative_and_atom_bounds_rederived": True,
    "velocity_remainder_bound_sharpened_to_21_eta_plus_5j": True,
    "small_source_partition_sharpened_to_eta_over_j_below_1_over_39": True,
    "angular_weight_floor_sharpened_to_7744_over_7745": True,
    "full_nonlinear_banach_map_independently_reproved_here": False,
    "flat_activation_volterra_and_moment_bounds_independently_reproved_here": False,
    "later_foundation_providers_independently_reproved_here": False,
    "unconditional_pulse_residual_force_or_blowup_proved": False,
}

REMAINING = [
    "independently rederive every nonlinear monomial and Lipschitz constant in the all-index core fixed-point map",
    "independently rederive the angular-source comparison, flat Volterra activation and complete collar moment estimates",
    "then continue the foundation frontier through continuation, matching roots, compact modulation/repair and the exterior join",
    "retain the later global pressure, weighted-edge, coefficient, source-cycle and residual obligations until independently validated",
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
        "all_exact_pressure_and_separation_checks_pass": all(exact["checks"].values()),
        "all_frozen_dependencies_bound": len(dependency_hashes()) == len(DEPENDENCIES),
        "scope_flags_are_strict_booleans": all(type(value) is bool for value in FLAGS.values()),
        "unresolved_continuum_frontier_retained": len(REMAINING) == 4,
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
        raise ValueError("Shared-core pressure recheck record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Shared-core pressure recheck gate failed")


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
    print(
        f"{STATUS}: exact shared datum retained; "
        "|eta|/j < 1/39 and chi >= 7744/7745 on the small-source region."
    )


if __name__ == "__main__":
    main()
