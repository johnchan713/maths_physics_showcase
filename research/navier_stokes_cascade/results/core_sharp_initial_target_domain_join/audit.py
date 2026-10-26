#!/usr/bin/env python3
"""Freeze the sharper initial target rate on the joined inverse domain."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from rate_join import review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "sharper-initial-target-rate-and-strict-first-real-domain-join"
DEPENDENCIES = (
    "core_compact_slow_inverse_domain_join",
    "core_initial_target_rate_gain",
    "core_initial_actual_profile_gate",
    "core_initial_viscosity_sharpness_witness",
    "core_initial_haar_row_budget",
    "core_symbolic_initial_haar_join",
    "core_rho_projector_signed_recompute",
    "core_initial_final_row_obstruction",
)
FLAGS = {
    "stale_alpha_target_downgrade_corrected": True,
    "initial_Haar_rate_beta_retained_through_linear_signed_maps": True,
    "complete_first_real_fast_source_order_alpha_joined": True,
    "radial_viscosity_order_three_minus_five_kappa_joined": True,
    "conditional_strict_first_real_margin_kappa_proved": True,
    "same_choice_projector_weighted_domain_and_pressure_memory_retained": True,
    "conservative_finite_request_graph_remains_sufficient": True,
    "actual_same_state_high_jet_enclosures_proved": False,
    "later_complex_wave_rate_improved": False,
    "complete_residual_force_or_blowup_certified": False,
}
REMAINING = [
    "supply the actual same-state high-jet and whole-support inequalities for the corrected beta-order target and retained A1 recipe",
    "re-evaluate the radial A0 and global pressure/moment support interfaces on the corrected first-real state",
    "join later complex harmonics, all covariance/mean corrections and the complete viscosity/cutoff/interaction residual",
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
        "README.md", "rate_join.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_initial_target_rate_gain/rate_gain.py",
        HERE.parent / "core_initial_haar_row_budget/README.md",
        HERE.parent / "core_rho_projector_signed_recompute/README.md",
        HERE.parent / "core_compact_slow_inverse_domain_join/README.md",
    ]
    gates = {
        "all_exact_rate_correction_controls": bool(result["checks"]) and all(result["checks"].values()),
        "eight_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scientific_scope_retained": (
            FLAGS["conditional_strict_first_real_margin_kappa_proved"]
            and not FLAGS["actual_same_state_high_jet_enclosures_proved"]
            and not FLAGS["later_complex_wave_rate_improved"]
            and not FLAGS["complete_residual_force_or_blowup_certified"]
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
        raise ValueError("Sharper initial target rate, provenance or scope changed")


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
    print("sharp initial target/domain join:", len(record["review"]["checks"]), "exact controls;", len(record["inherited_source_checks"]), "dependency source checks")


if __name__ == "__main__":
    main()
