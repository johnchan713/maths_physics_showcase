#!/usr/bin/env python3
"""Freeze the same-tail canonical full-path ODE, pressure and cutoff join."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from ode_join import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_same_band_gap_join",
    "core_phase_coefficient_budget",
    "core_tangent_propagator_closure",
    "canonical_primary_foundation",
    "core_canonical_entrance_jet_closure",
    "core_canonical_haar_column_jet_closure",
    "sourced_inverse_foundation",
    "core_finite_support_closure",
)
STATUS = "same-tail-canonical-full-path-ode-pressure-and-cutoff-contract-joined"
FLAGS = {
    "same_rule_defined_background_tail_and_original_physical_representatives_used": True,
    "canonical_reference_mixing_entrance_retained_and_not_relabelled": True,
    "finite_interval_homogeneous_primary_existence_uniqueness_joined": True,
    "full_moving_normal_pressure_and_potential_retained": True,
    "adjacent_sourced_inverse_keeps_unprojected_normal_source_pressure": True,
    "varying_reference_growth_retained_at_every_fixed_derivative_order": True,
    "whole_path_raw_solution_bounds_and_gaussian_ratio_cancellation_joined": True,
    "temporal_and_transverse_cutoffs_apply_to_potential_and_full_pressure": True,
    "raw_outgoing_response_declared_zero_or_compactly_supported": False,
    "canonical_and_inherited_primary_families_identified": False,
    "unspecified_historical_first_band_or_old_localized_force_certified": False,
    "fixed_compensator_projector_or_weighted_inverse_contracts_closed_here": False,
    "unconditional_pulse_residual_force_or_blowup_proved": False,
}
REMAINING = [
    "independently validate the fixed compensator/projector profiles, their actual supports and fixed jets on this same accepted domain",
    "join the compact and slow inverse maps with their weighted-domain and whole-path support estimates",
    "join generated harmonic sources, complete covariance and fresh axial/full-pressure means, then the full viscosity/cutoff/interaction residual",
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
    protocol = json.loads((HERE / "protocol.json").read_text())
    own = [HERE / name for name in (
        "README.md", "ode_join.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    gates = {
        "all_exact_ode_pressure_cutoff_controls": bool(result["checks"]) and all(result["checks"].values()),
        "eight_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scientific_scope_retained": (
            FLAGS["finite_interval_homogeneous_primary_existence_uniqueness_joined"]
            and FLAGS["adjacent_sourced_inverse_keeps_unprojected_normal_source_pressure"]
            and not FLAGS["fixed_compensator_projector_or_weighted_inverse_contracts_closed_here"]
            and not FLAGS["unconditional_pulse_residual_force_or_blowup_proved"]
        ),
    }
    return {
        "status": STATUS,
        "protocol": protocol,
        "flags": deepcopy(FLAGS),
        "review": result,
        "gates": gates,
        "dependency_evidence_sha256": records,
        "inherited_source_checks": inherited,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": list(REMAINING),
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Canonical full-path ODE, pressure, cutoff join, provenance or scope changed")


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
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(
        "canonical full-path ODE join:", len(record["review"]["checks"]),
        "exact controls;", len(record["inherited_source_checks"]),
        "dependency source checks",
    )


if __name__ == "__main__":
    main()
