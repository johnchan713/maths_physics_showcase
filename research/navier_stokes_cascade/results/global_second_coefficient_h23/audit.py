#!/usr/bin/env python3
"""Freeze the shifted-viscosity second repair and explicit global H2,3."""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from global_second import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "5b875cd821b89c02f23939346c7fa0440c2cf33a"
STATUS = "shifted-viscosity-second-repair-and-global-H23-assembled"
DEPENDENCIES = (
    "global_first_coefficient_c3",
    "global_leading_c5",
    "higher_pressure_source_moments",
    "axis_second_coefficient",
    "explicit_coefficient_cutoff",
    "coefficient_support_foundation",
    "cutoff_higher_jet_interface",
)
FLAGS = {
    "whole_support_Q10_Q11_with_both_shifted_viscosity_moments_assembled": True,
    "finalized_F1_and_U1_bounds_kept_distinct": True,
    "second_cutoff_and_all_repair_radial_jets_through_three_charged": True,
    "exact_inverse_800_and_fixed_lower_first_squares_retained": True,
    "explicit_global_second_H23_closed": True,
    "historical_first_physical_band_raw_deficit_gate_verified": False,
    "uniform_actual_physical_pulse_family_verified": False,
    "complete_covariance_means_residual_or_blowup_verified": False,
}
REMAINING = [
    "insert K0,5, K1,3 and H2,3 into the actual physical pulse maps and prove uniform value derivative and support estimates at every retained scale",
    "verify the historical first-band raw-deficit inequality or construct a separately certified later tail without changing the retained proof",
    "complete covariance mean corrections including axial velocity and both exceptional pressure tails, then bound viscosity, every cutoff and all residual interactions",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def protocol():
    return json.loads((HERE / "protocol.json").read_text())


def dependencies():
    return {
        str((HERE.parent / name / "evidence.json").relative_to(PROJECT)): sha(
            HERE.parent / name / "evidence.json"
        )
        for name in DEPENDENCIES
    }


def sources():
    files = [
        HERE / name
        for name in (
            "README.md", "global_second.py", "audit.py", "test_audit.py",
            "protocol.json", "requirements.txt", ".gitignore",
        )
    ]
    files += [HERE.parent / name / "README.md" for name in DEPENDENCIES]
    files += [HERE.parent / name / "evidence.json" for name in DEPENDENCIES]
    return {str(path.relative_to(PROJECT)): sha(path) for path in files}


def generate():
    checked = review()
    gates = deepcopy(checked["checks"])
    gates["all_scope_flags_are_strict_booleans"] = all(
        type(value) is bool for value in FLAGS.values()
    )
    gates["seven_frozen_dependency_records_bound"] = len(dependencies()) == 7
    return {
        "status": STATUS,
        "protocol": protocol(),
        "flags": deepcopy(FLAGS),
        "reviews": checked,
        "gates": gates,
        "dependency_evidence_sha256": dependencies(),
        "source_hashes": sources(),
        "remaining_obligations": REMAINING[:],
    }


def validate(record):
    expected = generate()
    if json.dumps(expected, sort_keys=True) != json.dumps(record, sort_keys=True):
        raise ValueError("Global second H2,3 record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Global second H2,3 check failed")


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
        f"{STATUS}: {len(record['gates'])} checks; explicit global H2,3 closed; "
        "all-scale pulse and complete covariance/residual obligations remain open."
    )


if __name__ == "__main__":
    main()
