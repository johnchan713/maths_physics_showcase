#!/usr/bin/env python3
"""Freeze the certified finite-stage domain join and diagonal obstruction."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from certified_recurrence import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "certified-fixed-stage-high-jet-domain-join-with-diagonal-obstruction"
DEPENDENCIES = (
    "core_first_sourced_high_jet_domain_closure",
    "core_first_real_high_jet_domain_closure",
    "core_finite_stage_compact_recurrence",
    "core_post_signed_seed_sync",
    "core_next_sourced_wave_join",
    "core_next_sourced_mean_tail_join",
    "core_post_signed_mean_closure",
    "core_post_signed_residual_join",
    "finite_stage_induction",
    "core_finite_support_closure",
)
FLAGS = {
    "fixed_stage_compact_request_graphs_joined_to_certified_primitive_domains": True,
    "historical_first_real_interfaces_distinguished_from_primitive_ceilings": True,
    "common_poles_recomputed_from_fully_expanded_primitive_requests": True,
    "complete_covariance_pressure_axial_velocity_means_and_supports_retained": True,
    "every_fixed_stage_jet_and_power_has_a_finite_analytical_compact_graph": True,
    "giant_jets_or_operator_constants_numerically_materialized": False,
    "stage_uniform_constant_growth_bound_proved": False,
    "q_dependent_one_field_diagonal_proved": False,
    "complete_forcing_energy_or_blowup_proved": False,
}
REMAINING = [
    "bound the growth of C_{j,m} strongly enough to verify log C_{j(q),m} <= (h B_{j(q)}-loss_m-N) log(1/q)",
    "construct and converge one q-dependent complete-potential field rather than a separate finite state for each target",
    "then re-audit every viscosity, cutoff, interaction, localization, force and energy term on that same diagonal before any blow-up claim",
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
        "README.md", "certified_recurrence.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_finite_stage_compact_recurrence/recurrence.py",
        HERE.parent / "core_first_sourced_high_jet_domain_closure/sourced_join.py",
        HERE.parent / "core_first_real_high_jet_domain_closure/high_jet_join.py",
        HERE.parent / "finite_stage_induction/induction.py",
        HERE.parent / "core_finite_support_closure/README.md",
    ]
    gates = {
        "all_exact_join_controls": bool(result["checks"]) and all(result["checks"].values()),
        "ten_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scientific_scope_retained": (
            FLAGS["fixed_stage_compact_request_graphs_joined_to_certified_primitive_domains"]
            and not FLAGS["giant_jets_or_operator_constants_numerically_materialized"]
            and not FLAGS["stage_uniform_constant_growth_bound_proved"]
            and not FLAGS["q_dependent_one_field_diagonal_proved"]
            and not FLAGS["complete_forcing_energy_or_blowup_proved"]
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
        raise ValueError("Certified finite-stage domain join, provenance or scope changed")


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
        "certified finite-stage domain join:",
        len(record["review"]["checks"]), "exact controls;",
        len(record["inherited_source_checks"]), "dependency source checks",
    )


if __name__ == "__main__":
    main()
