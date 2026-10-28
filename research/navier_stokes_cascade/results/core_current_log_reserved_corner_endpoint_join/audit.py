#!/usr/bin/env python3
"""Freeze the current corner and punctured-endpoint compatibility join."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from transfer import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "current-log-reserved-corner-and-punctured-endpoint-joined-conditionally"
DEPENDENCIES = (
    "core_current_log_reserved_residual_reconciliation",
    "current_endpoint_jets",
    "core_log_reserved_cartesian_corner_replay",
    "core_log_reserved_endpoint_compatibility_join",
    "core_corner_flatness_obstruction",
    "core_primitive_margin_normalization_reconciliation",
)
FLAGS = {
    "current_request_dependent_residual_stage_transferred_to_corner": True,
    "extra_A0_A1_source_derivatives_retained_at_corner": True,
    "current_punctured_endpoint_derivative_ledger_joined": True,
    "complete_pressure_covariance_axial_and_exterior_memory_retained": True,
    "endpoint_uniform_all_order_physical_bound_proved_here": False,
    "terminal_velocity_or_pressure_regular_at_origin_proved": False,
    "localization_force_energy_or_comparison_replayed_here": False,
    "unconditional_pulse_residual_force_exclusion_or_blowup_proved": False,
}
REMAINING = [
    "replay all twelve localization terms on this exact current corner/endpoint package",
    "construct the compatible compact future force and revalidate preterminal energy",
    "then replay fixed-gauge pressure recovery and expanding-ball comparison before any exclusion or blow-up claim",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if isinstance(value, Fraction):
        return str(value)
    return value


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
        "README.md", "transfer.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_current_log_reserved_residual_reconciliation/residual.py",
        HERE.parent / "current_endpoint_jets/endpoint.py",
        HERE.parent / "core_corner_flatness_obstruction/corner.py",
        HERE.parent / "core_primitive_margin_normalization_reconciliation/margins.py",
    ]
    gates = {
        "ten_exact_transfer_controls":
            len(result["checks"]) == 10 and all(result["checks"].values()),
        "six_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scientific_scope_retained": (
            result["scope"]["current_same_field_residual_cartesian_corner_flat_conditionally"]
            and result["scope"]["current_closed_time_localization_inputs_joined_conditionally"]
            and not result["scope"]["endpoint_uniform_all_order_physical_bound_proved_here"]
            and not result["scope"]["terminal_velocity_or_pressure_regular_at_origin_proved"]
            and not result["scope"]["localization_force_energy_comparison_replayed_here"]
            and not result["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"]
        ),
    }
    return {
        "status": STATUS,
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "flags": deepcopy(FLAGS),
        "review": encode(result),
        "gates": gates,
        "dependency_evidence_sha256": records,
        "inherited_source_checks": inherited,
        "source_hashes": {
            str(path.relative_to(PROJECT)): sha(path) for path in own + providers
        },
        "remaining_obligations": list(REMAINING),
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Current corner/endpoint transfer, provenance or scope changed")


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
        "current corner/endpoint join:", len(record["review"]["checks"]),
        "exact controls;", len(record["inherited_source_checks"]),
        "dependency source checks",
    )


if __name__ == "__main__":
    main()
