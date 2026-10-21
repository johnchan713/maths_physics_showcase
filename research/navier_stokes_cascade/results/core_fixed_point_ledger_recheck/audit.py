#!/usr/bin/env python3
"""Freeze the independent all-index core fixed-point ledger recheck."""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "6582f51d627b46fdcf244335b880223fcb6b59e3"
STATUS = "core-fixed-point-all-index-ledger-independently-rederived"
DEPENDENCIES = (
    "core_contraction_audit",
    "core_pressure_gap_recheck",
    "foundation_bundle_route",
)

FLAGS = {
    "all_index_product_and_integrated_derivative_bounds_rederived": True,
    "order_one_angular_inverse_bound_rederived": True,
    "all_normalized_R1_R2_monomials_enumerated": True,
    "complete_pressure_primitive_and_three_pressure_terms_retained": True,
    "map_invariance_and_lipschitz_majorants_rederived": True,
    "selected_Lambda_formula_gives_strict_contraction": True,
    "full_core_growth_and_activation_independently_reproved_here": False,
    "later_continuation_roots_modulation_or_exterior_reproved_here": False,
    "unconditional_pulse_residual_force_or_blowup_proved": False,
    "formal_proof_assistant_or_external_peer_certificate": False,
}

REMAINING = [
    "independently rederive the angular scalar comparison and positive source bound without division by chi",
    "independently rederive the axial-sector exit growth and flat Volterra activation with complete pressure moments",
    "continue through prescribed continuation, matching roots, compact modulation and the annular/exterior join on the same profile",
    "retain the later regional cone, global pressure, weighted-edge, coefficient, source-cycle and residual obligations",
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
        "all_exact_operator_source_and_contraction_checks_pass": all(exact["checks"].values()),
        "all_frozen_dependencies_bound": len(dependency_hashes()) == len(DEPENDENCIES),
        "scope_flags_are_strict_booleans": all(type(value) is bool for value in FLAGS.values()),
        "growth_and_later_frontier_retained": len(REMAINING) == 4,
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
        raise ValueError("Core fixed-point ledger record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Core fixed-point ledger gate failed")


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
        f"{STATUS}: 21 source monomials; pressure coefficient "
        "2720*A0^4 and Lipschitz coefficient 5440*A0^4 retained."
    )


if __name__ == "__main__":
    main()
