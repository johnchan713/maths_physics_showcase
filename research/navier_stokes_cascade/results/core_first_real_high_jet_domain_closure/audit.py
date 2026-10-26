#!/usr/bin/env python3
"""Freeze the certified first-real high-jet/domain closure."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from high_jet_join import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "certified-first-real-high-jet-and-compact-tail-domain-closure"
DEPENDENCIES = (
    "core_sharp_initial_target_domain_join",
    "core_leading_background_borel_join",
    "core_same_band_gap_join",
    "core_canonical_full_path_ode_join",
    "core_fixed_profile_projector_join",
    "core_compact_slow_inverse_domain_join",
    "core_initial_phi_high_jet_closure",
    "core_initial_a0_physical_tail_closure",
    "core_initial_physical_stage_sync",
)
FLAGS = {
    "five_primitive_accepted_tail_contracts_joined": True,
    "actual_current_Phi_and_DzPhi_fixed_order_inequalities_proved_on_certified_tail": True,
    "actual_current_radial_source_fixed_order_inequality_proved_on_certified_tail": True,
    "current_A1_and_A0_compact_tails_have_strict_kappa_margin": True,
    "every_fixed_physical_jet_and_power_has_finite_compact_tail_bound": True,
    "whole_support_seams_and_complete_pressure_memory_retained": True,
    "uniformity_over_all_orders_stages_or_physical_scales_proved": False,
    "post_sourced_changed_state_compact_DAG_proved": False,
    "complete_regular_or_physical_residual_improvement_proved": False,
    "force_or_blowup_certified": False,
}
REMAINING = [
    "carry the high-jet graph through the first sourced complex wave, its complete six-entry covariance, changed pressure moments and refreshed fast/slow means",
    "join the changed-state A0/A1 compact tails to the matching post-sourced regular residual without mixing states",
    "then establish uniform actual-pulse scale/support bounds and the complete viscosity, cutoff and interaction residual before any force or blow-up claim",
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
        "README.md", "high_jet_join.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_initial_phi_high_jet_closure/closure.py",
        HERE.parent / "core_initial_a0_physical_tail_closure/closure.py",
        HERE.parent / "core_sharp_initial_target_domain_join/rate_join.py",
        HERE.parent / "core_fixed_profile_projector_join/README.md",
        HERE.parent / "core_compact_slow_inverse_domain_join/README.md",
    ]
    gates = {
        "all_exact_join_controls": bool(result["checks"]) and all(result["checks"].values()),
        "nine_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scientific_scope_retained": (
            FLAGS["five_primitive_accepted_tail_contracts_joined"]
            and FLAGS["actual_current_Phi_and_DzPhi_fixed_order_inequalities_proved_on_certified_tail"]
            and not FLAGS["uniformity_over_all_orders_stages_or_physical_scales_proved"]
            and not FLAGS["post_sourced_changed_state_compact_DAG_proved"]
            and not FLAGS["complete_regular_or_physical_residual_improvement_proved"]
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
        raise ValueError("First-real high-jet/domain closure, provenance or scope changed")


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
        "first-real high-jet/domain closure:",
        len(record["review"]["checks"]),
        "exact controls;",
        len(record["inherited_source_checks"]),
        "dependency source checks",
    )


if __name__ == "__main__":
    main()
