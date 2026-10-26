#!/usr/bin/env python3
"""Freeze the compact/slow weighted-domain primitive join."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from inverse_join import review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "compact-slow-weighted-domain-and-support-contract-joined"
DEPENDENCIES = (
    "core_fixed_profile_projector_join",
    "core_rho_projector_signed_recompute",
    "core_initial_haar_row_budget",
    "core_native_haar_primitive_join",
    "core_symbolic_weighted_maps",
    "core_fast_inverse_budget",
    "core_slow_inverse_budget",
    "core_real_refresh_obligation_join",
    "core_initial_final_row_obstruction",
    "core_finite_support_closure",
)
FLAGS = {
    "same_choice_projector_and_current_weighted_rows_used": True,
    "shifted_compact_maps_joined_on_complete_weighted_domain": True,
    "zero_haar_temporal_inverse_and_actual_A1_velocity_retained": True,
    "fixed_slow_response_joined_on_complete_reserved_potential_hull": True,
    "whole_moving_pressure_product_and_both_pressure_tails_retained": True,
    "all_fixed_output_request_graphs_are_finite": True,
    "whole_path_and_zero_extension_support_contract_joined": True,
    "strict_final_row_gain_from_positive_majorants_proved": False,
    "uniformity_over_derivative_order_stage_width_or_tilt_proved": False,
    "later_generated_harmonics_or_complete_residual_certified": False,
    "unconditional_pulse_force_residual_improvement_or_blowup_proved": False,
}
REMAINING = [
    "prove delta_source+delta_viscosity>0 for the first-real final row, or derive an exact cancellation replacing the zero-margin estimate",
    "join every later generated harmonic source and complete covariance/mean correction on these same weighted domains",
    "close the full pointwise residual with viscosity, pressure, cutoffs and every old/new interaction before any force or blow-up claim",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dependency_review():
    records, sources = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent / name / "evidence.json"
        record = json.loads(path.read_text())
        records[str(path.relative_to(PROJECT))] = sha(path)
        for relative, digest in record.get("source_hashes", {}).items():
            source = PROJECT / relative
            sources[name + ":" + relative] = source.exists() and sha(source) == digest
    return records, sources


def generate():
    records, inherited = dependency_review()
    result = review()
    own = [HERE / name for name in (
        "README.md", "inverse_join.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_symbolic_weighted_maps/README.md",
        HERE.parent / "core_fast_inverse_budget/README.md",
        HERE.parent / "core_slow_inverse_budget/README.md",
        HERE.parent / "core_real_refresh_obligation_join/README.md",
        HERE.parent / "core_initial_final_row_obstruction/README.md",
    ]
    gates = {
        "all_exact_domain_request_and_obstruction_controls": bool(result["checks"]) and all(result["checks"].values()),
        "ten_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scientific_scope_retained": (
            FLAGS["shifted_compact_maps_joined_on_complete_weighted_domain"]
            and FLAGS["all_fixed_output_request_graphs_are_finite"]
            and not FLAGS["strict_final_row_gain_from_positive_majorants_proved"]
            and not FLAGS["unconditional_pulse_force_residual_improvement_or_blowup_proved"]
        ),
    }
    return {
        "status": STATUS,
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "flags": deepcopy(FLAGS),
        "review": result,
        "gates": gates,
        "dependency_evidence_sha256": records,
        "inherited_source_checks": inherited,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own + providers},
        "remaining_obligations": list(REMAINING),
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Compact/slow inverse domain, provenance, obstruction or scope changed")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record = generate()
    validate(record)
    if args.verify_record:
        validate(json.loads(args.verify_record.read_text()))
    with args.output.open("x") as output:
        output.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print("compact/slow domain join:", len(record["review"]["checks"]), "exact controls;", len(record["inherited_source_checks"]), "dependency source checks")


if __name__ == "__main__":
    main()
