#!/usr/bin/env python3
"""Freeze the later reference C5 lift and selected whole-support K0,5."""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from leading5 import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "a0c7b2a48a0f36a3adf73a89624de29c05666b6c"
STATUS = "same-reference-fifth-angular-and-selected-whole-support-leading-K05"
DEPENDENCIES = (
    "native_sixth_angular",
    "final_fifth_moment_interface",
    "compact_global_reconciliation",
    "fourth_angular_axial_interface",
    "global_axial_source_envelope",
    "regional_cone_join",
    "first_coefficient_c3_boundary",
)
FLAGS = {
    "same_later_reference_angular_root_lifted_to_Ceta5": True,
    "complete_later_reference_energy_and_strength_Ceta5_bounded": True,
    "later_reference_axial_pulse_Ceta5_bounded": True,
    "native_FINAL_compact_and_later_regions_joined_in_whole_support_K05": True,
    "selected_parameter_whole_support_leading_K05_quantified": True,
    "finalized_global_first_K13_closed": False,
    "second_signed_moments_repair_and_H23_closed": False,
    "uniform_actual_physical_pulse_family_verified": False,
    "complete_covariance_means_residual_or_blowup_verified": False,
}
REMAINING = [
    "insert K0,5=H^262 into the full-support Q03,Q13 formulas and assemble the explicit local/cutoff and inverse-800 constants into finalized K1,3",
    "use K1,3 with the leading Ceta4 row to evaluate Q10,Q11, complete the second repair and quantify H2,3",
    "audit all physical-scale pulse value/derivative/support bounds, complete-wave covariance means including axial velocity and exceptional pressure tails, and every viscosity/cutoff/interaction residual term",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def protocol():
    return {
        "parent_commit": PARENT,
        "later_root": "The same accepted angular correction root has Ceta5 norm at most 2^32 lambda^30; no root, frequency, width or schedule is reselected.",
        "energy": "The exact fifth interpolation derivative has coefficients 120,160,32 and raw cap 3936, giving |g_theta|5<72; complete reference energy remains |R|5<=8192T.",
        "strength": "The unchanged positive quadratic strength recurrence has all five normalized jets <=12 epsilon, hence |Amp|5<2 and |ci|5<=4lambda^20.",
        "later_cover": "The complete later reference axial cover 4eta, taper, zero pieces and pulse has |U0|Ceta5<=132.",
        "global_join": "Native PREMOD gives U Ceta6<A; the selected cutoff discharges F5 so the FINAL compact difference is <H^261; adding the native field and joining the later reference gives K0,5<H^262 on [0,Xv].",
        "scope": "Analytical selected-parameter leading fifth-angular join relative to retained foundations; global K1,3, second repair H2,3, all-scale pulses, covariance means, complete residual and blowup remain open.",
    }


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
            "leading5.py",
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


def serial(value):
    from fractions import Fraction

    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {key: serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(item) for item in value]
    return value


def generate():
    if json.loads((HERE / "protocol.json").read_text()) != protocol():
        raise ValueError("Protocol changed")
    checked = serial(review())
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
        raise ValueError("Global leading C5 record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Global leading C5 check failed")


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
        f"{STATUS}: {len(record['gates'])} checks; selected K0,5<H^262 closed; "
        "global K1,3 and H2,3 remain open."
    )


if __name__ == "__main__":
    main()
