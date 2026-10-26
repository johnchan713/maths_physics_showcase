#!/usr/bin/env python3
"""Freeze the certified complete-potential diagonal join and correction."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from diagonal_join import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "certified-domain-conditional-complete-potential-diagonal-join"
DEPENDENCIES = (
    "core_certified_finite_stage_domain_join",
    "core_complete_potential_diagonal_join",
    "core_finite_stage_compact_recurrence",
    "finite_stage_induction",
    "core_finite_support_closure",
    "canonical_local_realization",
    "realization_foundation",
    "core_first_sourced_high_jet_domain_closure",
)
FLAGS = {
    "certified_fixed_stage_domains_joined_to_complete_potential_diagonal": True,
    "previous_stage_growth_obstruction_corrected": True,
    "one_locally_finite_divergence_free_field_constructed_conditionally": True,
    "all_fixed_Cartesian_jets_of_local_complete_residual_flat_conditionally": True,
    "complete_pressure_viscosity_cutoff_crosses_and_square_retained": True,
    "stage_uniform_constants_or_raw_series_convergence_required": False,
    "actual_infinite_coefficients_or_cutoff_scales_numerically_evaluated": False,
    "closed_time_endpoint_and_localization_joined": False,
    "global_force_energy_comparison_or_blowup_proved": False,
}
REMAINING = [
    "join this same certified-domain diagonal to the closed-time endpoint atlas and complete one-sided potential, velocity and pressure jets",
    "then join localization, zero initialization, compact smooth force, whole-space energy and global comparison on that same field and gauge",
    "retain all unresolved continuum premises and do not promote the conditional local residual statement to unconditional blow-up",
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
        "README.md", "diagonal_join.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_certified_finite_stage_domain_join/certified_recurrence.py",
        HERE.parent / "core_complete_potential_diagonal_join/diagonal.py",
        HERE.parent / "canonical_local_realization/realization.py",
        HERE.parent / "finite_stage_induction/induction.py",
        HERE.parent / "realization_foundation/realized_checks.py",
    ]
    gates = {
        "all_exact_join_controls": bool(result["checks"]) and all(result["checks"].values()),
        "eight_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scientific_scope_retained": (
            FLAGS["certified_fixed_stage_domains_joined_to_complete_potential_diagonal"]
            and FLAGS["previous_stage_growth_obstruction_corrected"]
            and FLAGS["one_locally_finite_divergence_free_field_constructed_conditionally"]
            and not FLAGS["actual_infinite_coefficients_or_cutoff_scales_numerically_evaluated"]
            and not FLAGS["closed_time_endpoint_and_localization_joined"]
            and not FLAGS["global_force_energy_comparison_or_blowup_proved"]
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
        raise ValueError("Certified complete-potential diagonal join, provenance or scope changed")


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
    print(
        "certified complete-potential diagonal join:",
        len(record["review"]["checks"]), "exact controls;",
        len(record["inherited_source_checks"]), "dependency source checks",
    )


if __name__ == "__main__":
    main()
