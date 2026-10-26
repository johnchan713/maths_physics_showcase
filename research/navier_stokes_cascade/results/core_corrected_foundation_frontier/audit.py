#!/usr/bin/env python3
"""Freeze the corrected-field foundation-frontier reconciliation."""
import argparse
import hashlib
import json
from pathlib import Path

from frontier import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_corrected_comparison_join", "core_pressure_identification",
    "canonical_value_domain", "core_fixed_point_ledger_recheck",
    "core_growth_source_recheck", "core_activation_volterra_recheck",
    "core_continuation_matching_recheck", "core_modulation_repair_join_recheck",
    "core_global_pressure_moment_recheck", "core_coefficient_support_recheck",
    "core_primary_tail_recheck", "core_covariance_mean_recheck",
    "core_dynamic_mean_recheck", "core_slow_jet_loss_recheck",
    "core_pressure_shape_recheck", "core_actual_compensator_gate",
    "core_fixed_profile_tame", "core_initial_phi_high_jet_closure",
    "core_finite_stage_compact_recurrence", "core_complete_potential_diagonal_join",
)

FLAGS = {
    "leading_core_pressure_is_identified_with_same_completed_backward_pressure": True,
    "old_first_band_acceptance_is_required_for_current_later_tail": False,
    "compact_modulation_repair_and_global_pressure_rechecks_already_exist": True,
    "coefficient_support_primary_tail_and_complete_mean_rechecks_already_exist": True,
    "corrected_head_next_obligation_text_is_current_frontier": False,
    "five_primitive_accepted_tail_contracts_independently_certified_here": False,
    "passing_dependency_tests_is_a_continuum_proof": False,
    "unconditional_pulse_residual_force_or_blowup_proved": False,
}

REMAINING = [
    "independently validate the actual analytic leading/background and Borel identification",
    "independently validate the same-band phase/normal/value/cone gaps and canonical full-path ODE",
    "independently validate the fixed compensator/projector supports and exact compact/slow inverse domains",
    "then replay the complete current residual with viscosity, cutoffs and all interactions on those certified inputs",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dependency_review():
    records, sources = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent / name / "evidence.json"
        record = json.loads(path.read_text())
        records[str(path.relative_to(PROJECT))] = sha(path)
        for rel, digest in record.get("source_hashes", {}).items():
            target = PROJECT / rel
            sources[name + ":" + rel] = target.exists() and sha(target) == digest
    return records, sources


def generate():
    records, inherited = dependency_review()
    result = review()
    protocol = json.loads((HERE / "protocol.json").read_text())
    own = [HERE / name for name in (
        "README.md", "frontier.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    gates = {
        "exact_frontier_controls": all(result["checks"].values()),
        "twenty_dependency_records_bound": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scope_remains_conditional": not FLAGS["five_primitive_accepted_tail_contracts_independently_certified_here"] and not FLAGS["unconditional_pulse_residual_force_or_blowup_proved"],
    }
    return {
        "status": "corrected-foundation-frontier-reconciled-with-pressure-and-later-tail-not-repeated",
        "protocol": protocol, "flags": dict(FLAGS), "review": result,
        "gates": gates, "dependency_evidence_sha256": records,
        "inherited_source_checks": inherited,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": list(REMAINING),
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Corrected foundation frontier, provenance or scope changed")


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
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print("corrected foundation frontier:", len(record["review"]["checks"]),
          "exact controls;", len(record["inherited_source_checks"]),
          "dependency source checks")


if __name__ == "__main__":
    main()
