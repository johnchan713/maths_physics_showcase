#!/usr/bin/env python3
"""Freeze Gaussian-adapted ambient forced-inverse norm evidence."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=("core_initial_haar_row_budget","core_primary_root_budget","first_sourced_foundation",
              "regular_source_foundation","actual_band_pulse_bounds")
STATUS="first-forced-inverse-Gaussian-adapted-ambient-mixed-jet-norm-budget"
FLAGS={
    "moving_constraint_requires_more_than_direct_tangent_derivative_bound":True,
    "Gaussian_adapted_extension_preserves_actual_tangent_solution":True,
    "full_profile_weight_and_polynomial_ambient_propagator_retained":True,
    "explicit_finite_mixed_velocity_pressure_potential_polynomials":True,
    "moving_normal_and_unprojected_source_pressure_terms_retained":True,
    "actual_current_full_path_source_norm_inequality_independently_verified":False,
    "actual_core_coefficient_and_tangent_propagator_inputs_independently_verified":False,
    "complete_mean_pressure_memory_replaced_by_wave_pressure":False,
    "historical_evidence_or_fixed_physical_parameters_changed":False,
    "unconditional_pulse_bounds_complete_residual_or_blowup_proved":False,
}
REMAINING=[
    "derive actual current signed-prefix regular source polynomials in inequality (9) with fresh mean maps full viscosity interactions and whole-path supports",
    "substitute actual normal matrix logarithmic profile jets and full-core tangent propagator constants into (5) and (8)",
    "join every physical scale complete covariance mean corrections and complete residual including additive cutoff pressure tails before unconditional claims",
]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        path=HERE.parent/name/"evidence.json";record=json.loads(path.read_text())
        bound[str(path.relative_to(PROJECT))]=sha(path)
        for rel,expected in record["source_hashes"].items():checks[name+":"+rel]=sha(PROJECT/rel)==expected
    return bound,checks


def sources():
    paths=[HERE/n for n in ("README.md","budgets.py","controls.py","audit.py","test_audit.py","protocol.json","requirements.txt",".gitignore")]
    paths += [HERE.parent/n/"README.md" for n in DEPENDENCIES]
    paths += [HERE.parent/"core_primary_root_budget/budgets.py",HERE.parent/"first_sourced_foundation/sourced_checks.py",
              HERE.parent/"modulation_repair_all_order/map_checks.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_ambient_constraint_and_norm_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "five_historical_records_bound":len(bound)==5,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Forced inverse Gaussian adaptation, current source scope, norm recurrence or historical bindings changed")


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True,type=Path);parser.add_argument("--verify-record",type=Path)
    args=parser.parse_args(argv)
    if args.output.exists():raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record=generate();validate(record)
    if args.verify_record:validate(json.loads(args.verify_record.read_text()))
    with args.output.open("x") as output:output.write(json.dumps(record,indent=2,sort_keys=True)+"\n")
    print(f"{STATUS}: {len(record['review']['checks'])} exact controls; {len(record['inherited_source_checks'])} inherited source checks")


if __name__=="__main__":main()
