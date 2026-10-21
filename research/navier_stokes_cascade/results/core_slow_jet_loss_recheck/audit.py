#!/usr/bin/env python3
"""Freeze exact axial inverse moments and the isolated target-jet obstruction."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS="exact-slow-potential-moments-and-target-jet-pressure-response-rechecked"
DEPENDENCIES=("core_reserved_geometry_recheck","core_dynamic_mean_recheck",
              "compact_mean_foundation","current_target_cycle_bridge",
              "complex_mean_relaxation","complete_covariance_initial_means")
FLAGS={
    "prescribed_two_bump_inverse_logarithmic_moment_derived":True,
    "physical_Q_chart_normalization_retained_before_differentiation":True,
    "nonzero_cubic_target_jet_pressure_coefficient_identified":True,
    "whole_axial_quartic_coefficient_and_cancellation_condition_derived":True,
    "moving_pressure_commutator_retained_after_principal_cancellation":True,
    "same_finite_order_target_norm_obstruction_is_scoped_to_isolated_remainder":True,
    "complete_fixed_wave_and_both_base_pressure_tails_retained":True,
    "all_smooth_finite_cycles_or_actual_canonical_sequence_refuted":False,
    "actual_compensator_special_cancellation_identity_verified":False,
    "actual_uniform_high_target_jet_budgets_proved":False,
    "new_tail_profile_or_tilt_chosen":False,
    "unconditional_residual_improvement_or_blowup_proved":False,
    "finite_controls_are_actual_continuum_suprema":False,
}
REMAINING=[
    "prove actual current-target higher-jet inequalities for the exact cubic pressure and whole moving axial-moment response, or show the full linked cycle cancellation that replaces them",
    "evaluate or bound the actual compensator/profile moment combination before claiming special quartic cancellation",
    "join all lower transport, old means, generated-wave supports and complete viscosity/cutoff interactions into the full pointwise residual estimate",
]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        path=HERE.parent/name/"evidence.json";frozen=json.loads(path.read_text())
        bound[str(path.relative_to(PROJECT))]=sha(path)
        for rel,expected in frozen["source_hashes"].items():checks[name+":"+rel]=sha(PROJECT/rel)==expected
    return bound,checks


def sources():
    paths=[HERE/n for n in ("README.md","proof_checks.py","audit.py","test_audit.py","protocol.json","requirements.txt",".gitignore")]
    paths += [HERE.parent/n/"README.md" for n in DEPENDENCIES]
    paths += [HERE.parent/"complete_covariance_initial_means/means.py",HERE.parent/"compact_mean_corrections/README.md"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_moment_and_operator_controls_pass":all(result["checks"].values()),
                       "six_historical_records_bound":len(bound)==6,
                       "all_declared_historical_source_hashes_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Slow moments, high-jet response, source binding or scientific scope changed")


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
