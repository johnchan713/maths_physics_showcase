#!/usr/bin/env python3
"""Freeze the first-sourced high-jet/domain closure and residual obstruction."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from sourced_join import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "certified-first-sourced-high-jet-and-compact-tail-domain-closure"
DEPENDENCIES = (
    "core_first_real_high_jet_domain_closure",
    "core_first_sourced_physical_tail_closure",
    "first_sourced_foundation",
    "core_first_sourced_inverse_budget",
    "core_sourced_normal_transfer",
    "core_sourced_covariance_moment_join",
    "core_reweighted_fresh_mean_closure",
    "core_finite_support_closure",
    "core_post_signed_seed_sync",
    "core_compact_slow_inverse_domain_join",
)
FLAGS = {
    "six_first_sourced_primitive_contracts_joined": True,
    "historical_first_real_interface_ceiling_distinguished_from_primitive_ceiling": True,
    "actual_first_sourced_A0_and_A1_fixed_order_inequalities_proved_on_certified_tail": True,
    "every_fixed_physical_jet_and_power_has_first_sourced_compact_tail_bound": True,
    "complete_complex_covariance_pressure_axial_velocity_and_mean_refresh_retained": True,
    "uniformity_over_all_orders_stages_or_physical_scales_proved": False,
    "later_stage_compact_graph_recurrence_proved": False,
    "complete_regular_or_physical_residual_decay_proved": False,
    "force_or_blowup_certified": False,
}
REMAINING = [
    "construct the finite stage-indexed backward high-jet and common-pole recurrence for every retained later sourced state",
    "synchronize each later compact A0/A1 graph with the matching regular state and prove orders satisfying h*a > N+Lambda+2*m+eta",
    "then join the controlled physical diagonal, all viscosity/cutoff/interactions and forcing estimates before any blow-up claim",
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
        "README.md", "sourced_join.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_first_sourced_physical_tail_closure/closure.py",
        HERE.parent / "core_first_real_high_jet_domain_closure/high_jet_join.py",
        HERE.parent / "first_sourced_foundation/README.md",
        HERE.parent / "core_sourced_covariance_moment_join/README.md",
        HERE.parent / "core_reweighted_fresh_mean_closure/README.md",
    ]
    gates = {
        "all_exact_join_controls": bool(result["checks"]) and all(result["checks"].values()),
        "ten_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scientific_scope_retained": (
            FLAGS["six_first_sourced_primitive_contracts_joined"]
            and FLAGS["actual_first_sourced_A0_and_A1_fixed_order_inequalities_proved_on_certified_tail"]
            and not FLAGS["uniformity_over_all_orders_stages_or_physical_scales_proved"]
            and not FLAGS["later_stage_compact_graph_recurrence_proved"]
            and not FLAGS["complete_regular_or_physical_residual_decay_proved"]
            and not FLAGS["force_or_blowup_certified"]
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
        raise ValueError("First-sourced high-jet/domain closure, provenance or scope changed")


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
        "first-sourced high-jet/domain closure:",
        len(record["review"]["checks"]),
        "exact controls;",
        len(record["inherited_source_checks"]),
        "dependency source checks",
    )


if __name__ == "__main__":
    main()
