#!/usr/bin/env python3
"""Freeze the rule-defined-background same-band phase/value/cone gap join."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from gap_join import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_leading_background_borel_join",
    "regional_cone_join",
    "core_background_finite_head_budget",
    "weighted_background_edges",
    "current_frame_freezing",
    "core_primary_tail_recheck",
    "core_physical_representative_closure",
    "core_reference_scaled_normal_floor",
)
STATUS = "same-rule-defined-background-tail-phase-normal-value-cone-contract-joined"
FLAGS = {
    "same_rule_defined_background_and_completed_leading_profile_used": True,
    "leading_regional_value_orientation_and_cone_gaps_retained": True,
    "finite_head_and_weighted_stress_errors_are_order_q_to_2h": True,
    "one_derivative_independent_integer_selects_the_certified_later_tail": True,
    "same_tail_has_positive_background_spectral_and_same_tilt_cone_gaps": True,
    "phase_uses_full_F_G_axial_transport_and_original_carriers": True,
    "reference_relative_full_phase_normal_floor_joined_on_that_tail": True,
    "optional_native_Ctan_or_Cphi_margin_required_for_original_full_frame_route": False,
    "full_pressure_memory_and_radial_cutoff_curl_retained": True,
    "unspecified_historical_first_band_or_old_localized_force_certified": False,
    "canonical_homogeneous_full_path_ODE_independently_closed_here": False,
    "unconditional_pulse_residual_force_or_blowup_proved": False,
}
REMAINING = [
    "independently validate the canonical homogeneous full-path ODE with its moving normal, unprojected pressure source and every cutoff on this accepted tail",
    "audit fixed compensator/projector supports and jets and the compact/slow inverse weighted domains",
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
        "README.md", "gap_join.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    gates = {
        "all_exact_join_controls": bool(result["checks"]) and all(result["checks"].values()),
        "eight_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scientific_scope_retained": (
            FLAGS["reference_relative_full_phase_normal_floor_joined_on_that_tail"]
            and not FLAGS["canonical_homogeneous_full_path_ODE_independently_closed_here"]
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
        raise ValueError("Same-band phase/value/cone join, provenance or scope changed")


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
        "same-band phase/value/cone join:", len(record["review"]["checks"]),
        "exact controls;", len(record["inherited_source_checks"]),
        "dependency source checks",
    )


if __name__ == "__main__":
    main()
