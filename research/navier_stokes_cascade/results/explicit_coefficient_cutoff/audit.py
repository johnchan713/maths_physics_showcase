#!/usr/bin/env python3
"""Freeze the missing coefficient cutoff shape and bind its consequences."""

import argparse
import hashlib
import json
from copy import deepcopy
from pathlib import Path

from cutoff import consequences, review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "61ee5de6493e359c9a4e6056d5923ef9dc794d89"
STATUS = "explicit-coefficient-cutoff-all-fixed-jets-and-pulse-gate-boundary"
DEPENDENCIES = (
    "compact_global_reconciliation",
    "cutoff_repair_derivative_budget",
    "cutoff_higher_jet_interface",
    "sharp_current_scale_floor",
    "second_coefficient_source",
    "axis_second_coefficient",
    "higher_pressure_source_moments",
    "finite_head_frame_jets",
    "physical_pulse_foundation",
)
FLAGS = {
    "qualitative_coefficient_cutoff_is_now_reproducibly_specialized": True,
    "same_Cminus20_width_keep_cut_support_roots_gauges_and_pressure_memory_retained": True,
    "coefficient_and_variance_cutoffs_are_not_conflated": True,
    "all_fixed_normalized_and_physical_coefficient_cutoff_jets_enclosed": True,
    "cutoff_only_input_to_H23_and_pulse_jet_dictionary_discharged": True,
    "actual_first_scale_raw_deficit_and_mixed_gate_verified": False,
    "global_first_U1_Ceta3_and_finalized_second_H23_verified": False,
    "uniform_complete_pulse_family_covariance_means_residual_or_blowup_verified": False,
}
REMAINING = [
    "supply the actual raw same-tilt deficit d0 and endpoint constants and verify the displayed mixed first-scale inequality with Cstar(729,2916)",
    "quantify the finalized global first U1 C_eta^3 norm, feed the complete full-support Q10,Q11 moments to the second repair, and close H2,3",
    "then audit every pulse-family defining inverse gap value derivative and support estimate, complete covariance means including axial velocity and both pressure tails, and the full viscous cutoff interaction residual",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def protocol():
    return {
        "parent_commit": PARENT,
        "manuscript_sha256": "0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f",
        "correction": "The coefficient cutoff in Lemma 5.2 Step 1 was qualitatively chosen but not assigned a formula. It is not the variance cutoff fixed in the parent stream. The repeated wording unchanged cutoff is specialized here to one reproducible admissible member.",
        "choice": "phi=exp(-1/s-1/(1-s)); K=1-integral(phi)/total; kappa(X)=K(log(X/Xkeep)/Delta), Delta=tc/16, with the same Xkeep,Xcut and tc=C^-20",
        "normalized_bounds": "sup|K'|<729, sup|K''|<2916, sup|K'''|<58320; every fixed higher jet follows from the same exact polynomial recurrence",
        "physical_bounds": "partial_X^m kappa=X^-m sum_j s(m,j) Delta^-j K^(j); replacing X by Xkeep gives an explicit finite bound at every fixed m",
        "first_blocker": "H=gammaB/8-Cm ell^-6 must be positive and Cstar(729,2916)<=2^(2h(ell-2))*min(rb/8,omega^2 d0/(16Lpsi),H/(8Cp)); neither d0 nor the old first-band endpoint constants are supplied by cutoff smoothness",
        "H23_boundary": "the cutoff product part of H2,3 is now explicit, but finalized H2,3 still requires global first U1 C_eta^3, full-support Q10,Q11, the second five-row repair and its bump jets",
        "scope": "analytical construction specialization and exact dependency boundary; not a numerical actual-profile evaluation, all-family pulse proof, complete residual theorem or blowup proof",
    }


def source_hashes():
    files = [
        HERE / name
        for name in (
            "README.md",
            "cutoff.py",
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


def dependency_bindings():
    return {
        str((HERE.parent / name / "evidence.json").relative_to(PROJECT)): sha(
            HERE.parent / name / "evidence.json"
        )
        for name in DEPENDENCIES
    }


def generate():
    if json.loads((HERE / "protocol.json").read_text()) != protocol():
        raise ValueError("Protocol changed")
    selected = review()
    gates = deepcopy(selected["checks"])
    gates["all_flags_are_strict_booleans"] = all(
        type(value) is bool for value in FLAGS.values()
    )
    gates["all_consequences_are_strict_booleans"] = all(
        type(value) is bool for value in consequences().values()
    )
    gates["nine_frozen_dependency_records_bound"] = len(dependency_bindings()) == 9
    return {
        "status": STATUS,
        "protocol": protocol(),
        "flags": deepcopy(FLAGS),
        "reviews": {
            "selected_coefficient_cutoff": selected,
            "consequences": consequences(),
            "first_unclosed_mixed_gate": {
                "M": 729,
                "N": 2916,
                "H": "gammaB/8-Cm*ell^-6",
                "budget": "2^(2h(ell-2))*min(rb/8,omega^2*d0/(16Lpsi),H/(8Cp))",
                "required": "H>0 and Cstar(729,2916)<=budget",
            },
        },
        "gates": gates,
        "dependency_evidence_sha256": dependency_bindings(),
        "source_hashes": source_hashes(),
        "remaining_obligations": REMAINING[:],
    }


def validate(record):
    expected = generate()
    if json.dumps(expected, sort_keys=True) != json.dumps(record, sort_keys=True):
        raise ValueError("Explicit coefficient-cutoff record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Explicit coefficient-cutoff check failed")


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
        f"{STATUS}: {len(record['gates'])} checks; coefficient cutoff fixed "
        "with all fixed jets; first-scale and H2,3 gates remain open."
    )


if __name__ == "__main__":
    main()
