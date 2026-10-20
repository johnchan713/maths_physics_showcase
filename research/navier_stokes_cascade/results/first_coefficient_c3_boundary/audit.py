#!/usr/bin/env python3
"""Freeze the local first C3 estimate and exact global leading-C5 boundary."""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from boundary import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "df9bd482678a6c9cf502233bbc080c455c902630"
STATUS = "local-first-C3-and-exact-global-leading-C5-repair-boundary"
DEPENDENCIES = (
    "axis_first_coefficient",
    "signed_pressure_source_moments",
    "coefficient_support_foundation",
    "explicit_coefficient_cutoff",
    "native_sixth_angular",
    "final_fifth_moment_interface",
    "compact_global_reconciliation",
    "fourth_angular_axial_interface",
)
FLAGS = {
    "whole_inner_first_coefficient_X3_Ceta3_bound_derived": True,
    "selected_eta_independent_coefficient_cutoff_propagated_through_Ceta3": True,
    "full_support_signed_leading_Ceta3_pressure_moment_bounds_reduced_to_K05": True,
    "unchanged_linear_inverse800_first_repair_Ceta3_formula_retained": True,
    "prior_roadmap_Ceta4_sufficiency_for_constructing_K13_corrected": True,
    "selected_parameter_whole_support_leading_K05_quantified": False,
    "finalized_global_first_K13_second_repair_and_H23_closed": False,
    "raw_deficit_first_band_and_uniform_pulse_family_verified": False,
    "complete_covariance_axial_pressure_means_residual_or_blowup_verified": False,
}
REMAINING = [
    "extend the finalized later axial/reference leading stage from C_eta^4 to C_eta^5 and join it to native and FINAL compact fifth rows into one whole-support K0,5 majorant",
    "insert K0,5 into Q03,Q13 and the inverse-800 repair to close global K1,3, then evaluate full-support Q10,Q11 and the second repair H2,3",
    "verify the separate raw-deficit first-band inequality or a certified later tail, all pulse defining/inverse/gap/value/derivative/support estimates, complete covariance means and the full residual",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def protocol():
    return {
        "parent_commit": PARENT,
        "correction": "The finalized leading C_eta^4 bound is sufficient for the k=1 shifted-viscosity source after K1,3 is known, but not for constructing the first repair in C_eta^3; the signed leading source moments require whole-support U0 C_eta^5.",
        "local_result": "The existing 3delta/2 complex first-coefficient strips give Ceta3 factor 1+2/delta+4/delta^2+8/delta^3 through radial order three; the eta-independent selected cutoff propagates it by the exact binomial product and retains kappa_3 J0.",
        "derivative_ladder": "Signed leading I0,I1 in C_eta^k need a,b in C_eta^(k+1), hence v in C_eta^(k+1), hence U0 in C_eta^(k+2) because v contains the eta derivative of the regular radial average.",
        "pressure_memory": "Q03,Q13 are full-support signed totals through Xv; axis viscosity remains in Q03, weighted radial viscosity cancels only in Q13, and neither is an absolute-source bound.",
        "global_boundary": "Supply an explicit selected-parameter K0,5=sup_[0,Xv]|U0|Ceta5; then V4<29K0,5, the displayed Q03,Q13 bounds, five ordinary discrepancies and the unchanged inverse-800 formula quantify the full first repair and K1,3.",
        "scope": "Analytical local C3 extension and exact conditional global-repair boundary; no selected whole-support K0,5, finalized K1,3, H2,3, all-family pulse, covariance/mean, residual or blowup conclusion.",
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
            "boundary.py",
            "audit.py",
            "test_audit.py",
            "protocol.json",
            "requirements.txt",
            ".gitignore",
        )
    ]
    files += [HERE.parent / name / "README.md" for name in DEPENDENCIES]
    files += [HERE.parent / name / "evidence.json" for name in DEPENDENCIES]
    files += [
        HERE.parent / "axis_first_coefficient" / "axis_bounds.py",
        HERE.parent / "signed_pressure_source_moments" / "moments.py",
        HERE.parent / "coefficient_support_foundation" / "support_bounds.py",
        HERE.parent / "explicit_coefficient_cutoff" / "cutoff.py",
    ]
    return {str(path.relative_to(PROJECT)): sha(path) for path in files}


def generate():
    if json.loads((HERE / "protocol.json").read_text()) != protocol():
        raise ValueError("Protocol changed")
    checked = review()
    gates = deepcopy(checked["checks"])
    gates["all_scope_flags_are_strict_booleans"] = all(
        type(value) is bool for value in FLAGS.values()
    )
    gates["eight_frozen_dependency_records_bound"] = len(dependencies()) == 8
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
        raise ValueError("First-coefficient C3 boundary record changed")
    if not all(expected["gates"].values()):
        raise ValueError("First-coefficient C3 boundary check failed")


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
        f"{STATUS}: {len(record['gates'])} checks; local C3 closed; "
        "whole-support leading K0,5 and global K1,3 remain open."
    )


if __name__ == "__main__":
    main()
