#!/usr/bin/env python3
"""Freeze exact initial mean-defect budgets without promoting later closure."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "complete-initial-mean-defects-and-projected-quadratic-conditioning-rechecked"
DEPENDENCIES = (
    "core_primary_tail_recheck", "core_coefficient_support_recheck",
    "compact_mean_foundation", "complete_covariance_initial_means",
    "canonical_initial_mean_join", "current_target_moment_audit",
)
FLAGS = {
    "complete_covariance_axial_velocity_and_curl_squares_retained": True,
    "both_pressure_tails_and_dynamic_mean_terms_retained": True,
    "initial_nonlinear_defects_recomputed_after_linear_five_map": True,
    "projected_angular_and_axial_energy_lower_bounds_derived": True,
    "fixed_hull_axial_quadratic_inverse_tilt_cost_derived": True,
    "joint_pressure_axial_dynamic_cancellation_gate_derived": True,
    "same_common_primary_tail_preserved": True,
    "actual_dynamic_terms_cannot_cancel_energy": False,
    "bounds_uniform_across_all_actual_physical_scales": False,
    "later_nonzero_means_or_changing_wave_covered_by_initial_formula": False,
    "all_generated_sources_admissible": False,
    "complete_residual_improvement_or_blowup_proved": False,
    "manufactured_tests_are_continuum_proof": False,
}
REMAINING = [
    "bound actual dynamic C_P,C_Q and g*u plus every later old-mean and changing-wave cross term in the same weighted source norms on the common tail",
    "close the current signed-target cancellation gate with complete pressure and axial velocity, rather than inferring a new exponent from the linear five-map",
    "prove uniform generated-source supports/jets and the complete residual including viscosity, all cutoffs and interactions",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provenance():
    bound, checks = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent / name / "evidence.json"
        frozen = json.loads(path.read_text())
        bound[str(path.relative_to(PROJECT))] = sha(path)
        for rel, expected in frozen["source_hashes"].items():
            checks[name+":"+rel] = sha(PROJECT/rel)==expected
    return bound, checks


def sources():
    files = [HERE/n for n in ("README.md","proof_checks.py","audit.py","test_audit.py","protocol.json","requirements.txt",".gitignore")]
    files += [HERE.parent/n/"README.md" for n in DEPENDENCIES]
    files += [HERE.parent/"complete_covariance_initial_means/means.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in files}


def generate():
    result = review()
    bound, checks = provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),
                flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_mean_controls_pass":all(result["checks"].values()),
                       "declared_dependency_source_hashes_match":bool(checks) and all(checks.values()),
                       "six_frozen_records_bound":len(bound)==6},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected = generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Mean budgets, source binding or scientific scope changed")


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--verify-record",type=Path)
    args=parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record=generate()
    validate(record)
    if args.verify_record:
        validate(json.loads(args.verify_record.read_text()))
    with args.output.open("x") as output:
        output.write(json.dumps(record,indent=2,sort_keys=True)+"\n")
    print(f"{STATUS}: {len(record['review']['checks'])} exact controls, {len(record['inherited_source_checks'])} inherited source checks")


if __name__=="__main__":
    main()
