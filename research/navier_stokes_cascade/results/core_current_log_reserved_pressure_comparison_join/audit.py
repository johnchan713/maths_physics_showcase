#!/usr/bin/env python3
"""Freeze current fixed-gauge pressure recovery and comparison."""

import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from pressure_comparison import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "current-log-reserved-pressure-and-expanding-ball-comparison-joined-conditionally"
DEPENDENCIES = (
    "core_current_log_reserved_preterminal_energy_join",
    "core_current_log_reserved_localization_force_join",
    "core_log_reserved_pressure_comparison_replay",
    "core_corrected_comparison_join",
    "comparison_foundation",
    "core_same_field_positive_core_reconciliation",
    "core_leading_profile_certificate",
    "core_current_log_source_reconciliation",
)
FLAGS = {
    "exact_current_field_force_and_zero_datum_used": True,
    "individual_force_divergence_retained_and_difference_force_cancelled_once": True,
    "harmonic_pressure_gradient_removed_by_finite_energy_and_momentum": True,
    "all_expanding_ball_remainders_decay": True,
    "current_protected_core_identity_joined_conditionally": True,
    "finite_energy_hypothesis_removed": False,
    "future_velocity_or_pressure_constructed": False,
    "all_continuum_contracts_independently_reproved_here": False,
    "unconditional_admissible_force_exclusion_or_blowup_proved": False,
}
REMAINING = [
    "audit the new current-head dependency graph and enumerate every still-conditional primitive or operator contract",
    "independently discharge each remaining continuum contract before promoting admissible-force or exclusion statements",
    "do not infer an unconditional Navier-Stokes blow-up theorem from hashes, tests or conditional comparison",
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
            source_path = PROJECT / relative
            sources[name + ":" + relative] = (
                source_path.exists() and sha(source_path) == digest
            )
    return records, sources


def generate():
    records, inherited = dependency_review()
    result = review()
    own = [HERE / name for name in (
        "README.md", "pressure_comparison.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_current_log_reserved_preterminal_energy_join/preterminal_energy.py",
        HERE.parent / "core_corrected_comparison_join/comparison.py",
        HERE.parent / "core_same_field_positive_core_reconciliation/core_path.py",
        HERE.parent / "core_leading_profile_certificate/certificate.py",
        HERE.parent / "core_current_log_source_reconciliation/reconcile.py",
    ]
    scope = result["scope"]
    gates = {
        "eleven_exact_comparison_controls":
            len(result["checks"]) == 11 and all(result["checks"].values()),
        "eight_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match":
            bool(inherited) and all(inherited.values()),
        "conditional_current_comparison_joined": (
            scope["fixed_gauge_pressure_gradient_recovered_conditionally"]
            and scope["exact_current_field_expanding_ball_comparison_replayed_conditionally"]
            and scope["global_smooth_finite_energy_competitor_exclusion_follows_conditionally"]
        ),
        "scientific_scope_retained": (
            not scope["finite_energy_hypothesis_removed"]
            and not scope["mixed_pressure_moment_isotropy_assumed"]
            and not scope["future_velocity_or_pressure_constructed"]
            and not scope["all_frozen_continuum_contracts_independently_reproved_here"]
            and not scope["unconditional_admissible_force_or_blowup_proved"]
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
        raise ValueError("Current comparison join, provenance or scope changed")


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
        "current pressure/comparison join:", len(record["review"]["checks"]),
        "exact controls;", len(record["inherited_source_checks"]),
        "dependency source checks",
    )


if __name__ == "__main__":
    main()
