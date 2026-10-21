#!/usr/bin/env python3
"""Freeze the common primary tail and all-band edge/flatness recheck."""

import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "same-profile-primary-tail-value-gates-and-dyadic-edge-flatness-rechecked"
DEPENDENCIES = (
    "core_coefficient_support_recheck", "explicit_tilt_pulse_tail",
    "physical_pulse_foundation", "auxiliary_geometry_audit",
    "canonical_value_domain", "sourced_inverse_foundation",
)
FLAGS = {
    "one_integer_joins_all_finite_primary_value_gates": True,
    "binary_scale_exponent_certified_without_materializing_its_power": True,
    "all_fixed_primary_physical_jets_joined_on_same_certified_tail": True,
    "edge_width_cost_and_retained_half_edge_factor_explicit": True,
    "additive_envelope_and_support_controlled_cutoff_errors_have_all_band_flatness_bound": True,
    "complete_radial_axial_and_two_exceptional_pressure_terms_retained": True,
    "retained_core_frame_inverse_and_source_envelope_lemmas_are_dependencies": True,
    "all_generated_sources_satisfy_inverse_input_hypotheses": False,
    "complete_covariance_means_full_residual_or_blowup_proved": False,
    "old_first_band_or_old_localized_force_reidentified": False,
    "manufactured_controls_are_continuum_suprema": False,
}
REMAINING = [
    "join complete-wave covariance and current mean corrections on the common tail, including axial velocity, curl remainders, signed squares and Pi1/Pi2 cutoff derivatives",
    "prove full lifted-path envelope, endpoint-extension, edge/jet and harmonic-moment bounds for every generated nonlinear source",
    "prove complete residual with viscosity, every cutoff and interaction, and propagate any new field/force without identifying an unverified historical first band",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serial(x):
    if isinstance(x, Fraction):
        return str(x)
    if isinstance(x, dict):
        return {key: serial(value) for key, value in x.items()}
    if isinstance(x, (tuple, list)):
        return [serial(value) for value in x]
    return x


def provenance():
    bound, checks = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent / name / "evidence.json"
        frozen = json.loads(path.read_text())
        bound[str(path.relative_to(PROJECT))] = sha(path)
        for rel, expected in frozen["source_hashes"].items():
            checks[name + ":" + rel] = sha(PROJECT / rel) == expected
    return bound, checks


def sources():
    files = [HERE / n for n in ("README.md", "proof_checks.py", "audit.py", "test_audit.py", "protocol.json", "requirements.txt", ".gitignore")]
    files += [HERE.parent / n / "README.md" for n in DEPENDENCIES]
    files += [HERE.parent / "explicit_tilt_pulse_tail/tail_bounds.py"]
    return {str(path.relative_to(PROJECT)): sha(path) for path in files}


def generate():
    result = review()
    bound, checks = provenance()
    return {
        "status": STATUS,
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "flags": deepcopy(FLAGS),
        "review": serial(result),
        "gates": {
            "all_exact_primary_tail_controls_pass": all(result["checks"].values()),
            "declared_dependency_source_hashes_match": bool(checks) and all(checks.values()),
            "six_frozen_records_bound": len(bound) == 6,
        },
        "dependency_evidence_sha256": bound,
        "inherited_source_checks": checks,
        "source_hashes": sources(),
        "remaining_obligations": REMAINING[:],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Primary tail, flatness, source binding or scientific scope changed")


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
    with args.output.open("x") as output:
        output.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(f"{STATUS}: one primary tail, explicit dyadic/edge constants, {len(record['inherited_source_checks'])} inherited source checks")


if __name__ == "__main__":
    main()
