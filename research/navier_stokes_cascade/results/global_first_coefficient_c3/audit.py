#!/usr/bin/env python3
"""Freeze the full-pressure-memory first repair and explicit global K1,3."""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from global_first import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "6ef0567c8ef3373735d956dc0fec59e864878205"
STATUS = "full-pressure-memory-first-repair-and-global-K13-assembled"
DEPENDENCIES = (
    "global_leading_c5",
    "first_coefficient_c3_boundary",
    "axis_first_coefficient",
    "collar_operator_envelope",
    "coefficient_support_foundation",
    "cutoff_higher_jet_interface",
    "global_swirl_envelope",
)
FLAGS = {
    "selected_whole_support_K05_inserted_in_signed_Ceta3_sources": True,
    "axis_viscosity_and_complete_pressure_memory_retained": True,
    "local_core_collar_and_selected_cutoff_value_constant_assembled": True,
    "exact_inverse_800_and_actual_disjoint_U_bump_value_used": True,
    "finalized_global_first_K13_closed": True,
    "second_signed_moments_repair_and_H23_closed": False,
    "uniform_actual_physical_pulse_family_verified": False,
    "complete_covariance_means_residual_or_blowup_verified": False,
}
REMAINING = [
    "use the finalized K1,3 formula and leading Ceta4 row in the full-support shifted-viscosity Q10,Q11 identities, then assemble the second repair and H2,3",
    "audit uniform actual-pulse value/derivative/support estimates on every physical scale",
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
            "README.md",
            "global_first.py",
            "audit.py",
            "test_audit.py",
            "protocol.json",
            "requirements.txt",
            ".gitignore",
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
        raise ValueError("Global first C3 record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Global first C3 check failed")


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
        f"{STATUS}: {len(record['gates'])} checks; explicit global K1,3 closed; "
        "second repair and H2,3 remain open."
    )


if __name__ == "__main__":
    main()
