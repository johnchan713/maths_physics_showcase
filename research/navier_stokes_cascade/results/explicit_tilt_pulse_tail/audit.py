#!/usr/bin/env python3
"""Freeze explicit tilt/raw-deficit and separately certified pulse-tail data."""

import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tail_bounds import physical_exponents, tail_certificate, tilt_certificate


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "9ff05e276d247ec23c15f08f22906234d7037c5a"
STATUS = "explicit-tilt-raw-deficit-and-primary-pulse-tail-certificate"
DEPENDENCIES = (
    "global_second_coefficient_h23",
    "explicit_coefficient_cutoff",
    "leading_cone_margin",
    "same_tilt_deficit",
    "leading_frame_freezing",
    "current_frame_freezing",
    "physical_pulse_foundation",
)
FLAGS = {
    "previously_existential_tilt_now_has_explicit_specialization": True,
    "explicit_specialization_has_quantitative_raw_deficit": True,
    "selected_cutoff_mixed_gate_has_finite_later_tail": True,
    "primary_physical_pulse_value_derivative_and_support_transfer_on_new_tail": True,
    "unspecified_historical_first_integer_verified": False,
    "old_localized_force_automatically_reidentified_with_new_tail": False,
    "all_generated_sources_complete_covariance_residual_or_blowup_verified": False,
}
REMAINING = [
    "propagate the separately certified tail through localization and force without identifying it with an absent historical integer",
    "verify complete-wave covariance and mean corrections including axial velocity, curl remainders and both exceptional pressure tails",
    "prove every generated source satisfies the sourced inverse hypotheses and bound viscosity, cutoffs and all residual interactions",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if type(value) is F:
        return f"{value.numerator}/{value.denominator}"
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def protocol():
    return json.loads((HERE / "protocol.json").read_text())


def dependencies():
    return {
        str((HERE.parent / name / "evidence.json").relative_to(PROJECT)): sha(
            HERE.parent / name / "evidence.json"
        )
        for name in DEPENDENCIES
    }


def sources():
    files = [HERE / name for name in (
        "README.md", "tail_bounds.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    files += [HERE.parent / name / "README.md" for name in DEPENDENCIES]
    files += [HERE.parent / name / "evidence.json" for name in DEPENDENCIES]
    return {str(path.relative_to(PROJECT)): sha(path) for path in files}


def review():
    tilt_rows = []
    for power in range(1, 13):
        item = tilt_certificate(F(1, 2**power))
        tilt_rows.append({
            "power": power,
            "valid": item["chosen_d0"] <= item["actual_deficit_floor"],
            "strict": item["unused_slack"] > 0,
            "ratio": item["coarse_ratio_upper"],
        })
    fixture = tail_certificate(
        F(1, 200), F(1000), F(1, 9), F(1, 8), F(13), F(17),
        F(19), F(6), F(1, 23), F(1, 64),
    )
    exponent_rows = [
        {"a": a, "b": b, **physical_exponents(F(1, 200), a, b)}
        for a in range(7) for b in range(7 - a)
    ]
    checks = {
        "twelve_exact_tilt_rows_pass": all(
            row["valid"] and row["strict"] and row["ratio"] < 1
            for row in tilt_rows
        ),
        "manufactured_tail_fixture_passes_all_gates": all(fixture["checks"].values()),
        "twenty_eight_physical_exponent_rows_recorded": len(exponent_rows) == 28,
        "pressure_has_one_more_spatial_loss_than_velocity": all(
            row["pressure"] - row["velocity"] == row["sx"]
            for row in exponent_rows
        ),
        "potential_has_one_less_spatial_loss_than_velocity": all(
            row["velocity"] - row["potential"] == row["sx"]
            for row in exponent_rows
        ),
    }
    return {
        "tilt_rows": tilt_rows,
        "tail_fixture": fixture,
        "physical_exponent_rows": exponent_rows,
        "checks": checks,
    }


def generate():
    reviewed = review()
    gates = deepcopy(reviewed["checks"])
    gates["all_scope_flags_are_strict_booleans"] = all(
        type(value) is bool for value in FLAGS.values()
    )
    gates["seven_frozen_dependency_records_bound"] = len(dependencies()) == 7
    return encode({
        "status": STATUS,
        "protocol": protocol(),
        "flags": deepcopy(FLAGS),
        "reviews": reviewed,
        "gates": gates,
        "dependency_evidence_sha256": dependencies(),
        "source_hashes": sources(),
        "remaining_obligations": REMAINING[:],
    })


def validate(record):
    expected = generate()
    if json.dumps(expected, sort_keys=True) != json.dumps(record, sort_keys=True):
        raise ValueError("Explicit tilt/pulse-tail record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Explicit tilt/pulse-tail check failed")


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
    with args.output.open("x") as stream:
        stream.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(
        f"{STATUS}: {len(record['gates'])} checks; explicit raw deficit and "
        "separate primary-pulse tail certified; complete covariance/residual remain open."
    )


if __name__ == "__main__":
    main()
