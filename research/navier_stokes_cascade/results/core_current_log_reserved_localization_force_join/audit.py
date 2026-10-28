#!/usr/bin/env python3
"""Freeze current log-reserved localization and compact-force gluing."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from localization_force import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "current-log-reserved-localization-and-compact-force-joined-conditionally"
DEPENDENCIES = (
    "core_current_log_reserved_corner_endpoint_join",
    "core_current_log_reserved_residual_reconciliation",
    "current_endpoint_jets",
    "core_log_reserved_localization_force_replay",
    "core_corrected_force_join",
    "force_foundation",
    "current_localized_force",
)
FLAGS = {
    "current_same_field_axisymmetric_localization_replayed": True,
    "all_twelve_time_viscous_pressure_and_nonlinear_groups_retained": True,
    "force_divergence_in_pressure_equation_retained": True,
    "current_terminal_force_jets_joined_conditionally": True,
    "one_smooth_compact_future_force_constructed_conditionally": True,
    "actual_terminal_norms_or_extension_scales_numerically_evaluated": False,
    "terminal_velocity_or_pressure_regular_at_origin_proved": False,
    "future_velocity_or_pressure_constructed": False,
    "preterminal_energy_or_comparison_replayed_here": False,
    "unconditional_pulse_residual_force_exclusion_or_blowup_proved": False,
}
REMAINING = [
    "replay the fixed-time preterminal energy identity for this exact localized current field and force",
    "retain the factor-two dissipation bound and take energy and dissipation terminal limits separately",
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
        "README.md", "localization_force.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_current_log_reserved_corner_endpoint_join/transfer.py",
        HERE.parent / "core_current_log_reserved_residual_reconciliation/residual.py",
        HERE.parent / "current_endpoint_jets/endpoint.py",
        HERE.parent / "core_corrected_force_join/localization.py",
    ]
    scope = result["scope"]
    gates = {
        "ten_exact_localization_force_controls":
            len(result["checks"]) == 10 and all(result["checks"].values()),
        "seven_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "localization_and_force_joined": (
            scope["current_same_field_axisymmetric_localization_replayed_conditionally"]
            and scope["one_smooth_compact_force_extension_constructed_conditionally"]
        ),
        "scientific_scope_retained": (
            not scope["terminal_velocity_or_pressure_regular_at_origin_proved"]
            and not scope["future_velocity_or_pressure_constructed"]
            and not scope["preterminal_energy_replayed_here"]
            and not scope["pressure_recovery_comparison_or_blowup_proved"]
            and not scope["unconditional_pulse_residual_or_force_claimed"]
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
        raise ValueError("Current localization/force join, provenance or scope changed")


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
        "current localization/force join:", len(record["review"]["checks"]),
        "exact controls;", len(record["inherited_source_checks"]),
        "dependency source checks",
    )


if __name__ == "__main__":
    main()
