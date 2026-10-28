#!/usr/bin/env python3
"""Freeze the exact-current log-reserved preterminal energy join."""

import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from preterminal_energy import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "current-log-reserved-preterminal-energy-joined-conditionally"
DEPENDENCIES = (
    "core_current_log_reserved_localization_force_join",
    "core_current_log_reserved_corner_endpoint_join",
    "core_current_log_reserved_residual_reconciliation",
    "core_log_reserved_preterminal_energy_replay",
    "core_corrected_force_join",
    "force_foundation",
)
FLAGS = {
    "exact_current_localized_field_and_force_used": True,
    "complete_pressure_memory_and_unit_viscosity_retained": True,
    "zero_norm_regularization_retained": True,
    "factor_two_dissipation_retained": True,
    "energy_and_dissipation_terminal_limits_taken_separately": True,
    "terminal_velocity_or_pressure_regular_at_origin_proved": False,
    "future_velocity_or_pressure_constructed": False,
    "pressure_recovery_or_comparison_replayed_here": False,
    "unconditional_pulse_residual_force_exclusion_or_blowup_proved": False,
}
REMAINING = [
    "replay fixed-gauge pressure recovery for the exact current field and compact force with -div(f) retained",
    "prove expanding-ball comparison for every admissible smooth finite-energy competitor with the same datum and force",
    "only then revisit competitor exclusion; no unconditional blow-up claim is available",
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
        "README.md", "preterminal_energy.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_current_log_reserved_localization_force_join/localization_force.py",
        HERE.parent / "core_log_reserved_preterminal_energy_replay/preterminal_energy.py",
        HERE.parent / "core_corrected_force_join/localization.py",
    ]
    scope = result["scope"]
    gates = {
        "eleven_exact_energy_controls":
            len(result["checks"]) == 11 and all(result["checks"].values()),
        "six_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match":
            bool(inherited) and all(inherited.values()),
        "exact_current_preterminal_energy_joined": (
            scope["exact_current_field_preterminal_energy_identity_replayed_conditionally"]
            and scope["separate_finite_energy_and_total_dissipation_bounds_obtained_conditionally"]
        ),
        "scientific_scope_retained": (
            not scope["terminal_velocity_or_pressure_regular_at_origin_proved"]
            and not scope["future_velocity_or_pressure_constructed"]
            and not scope["fixed_gauge_pressure_recovery_replayed_here"]
            and not scope["expanding_ball_comparison_or_competitor_exclusion_proved"]
            and not scope["unconditional_pulse_residual_force_or_blowup_claimed"]
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
        raise ValueError("Current energy join, provenance or scope changed")


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
        "current preterminal energy join:", len(record["review"]["checks"]),
        "exact controls;", len(record["inherited_source_checks"]),
        "dependency source checks",
    )


if __name__ == "__main__":
    main()
