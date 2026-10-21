#!/usr/bin/env python3
"""Freeze the conditional current scalar recurrence and retained scale losses."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS="finite-current-scalar-cycle-and-polynomial-pressure-input-budget"
DEPENDENCIES=("core_scalar_target_pressure_join","core_fixed_profile_tame","core_dynamic_mean_recheck",
              "current_target_cycle_bridge","complex_mean_relaxation","compact_mean_foundation")
FLAGS={
    "existing_postslow_order_not_presented_as_new":True,
    "current_wave_fast_and_fresh_slow_scalar_norm_recurrence_derived":True,
    "sharper_accumulated_Haar_rate_used_with_finite_sums":True,
    "polynomial_scale_losses_retained_in_pressure_input_bounds":True,
    "positive_margin_given_tail_test_derived_without_new_tail_selection":True,
    "actual_seed_wave_mean_inverse_and_other_defect_inputs_independently_verified":False,
    "uniform_infinite_stage_or_all_order_constant_proved":False,
    "actual_first_band_or_full_residual_certified":False,
    "actual_profile_or_accepted_tail_changed":False,
    "unconditional_pulse_or_blowup_proved":False,
}
REMAINING=[
    "independently discharge the retained canonical seed and actual current complete-wave, fast/slow inverse and other scalar-defect norm inputs",
    "check actual polynomial scale-loss and pressure budgets against the required tolerances on the existing common domain",
    "join every generated support, full pressure memory, complete covariance, viscosity, cutoffs and nonlinear interaction in the complete residual",
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
    paths=[HERE/n for n in ("README.md","proof_checks.py","audit.py","test_audit.py","protocol.json","requirements.txt",".gitignore")]
    paths += [HERE.parent/n/"README.md" for n in DEPENDENCIES]
    paths += [HERE.parent/"core_dynamic_mean_recheck/proof_checks.py",HERE.parent/"complete_covariance_initial_means/means.py",
              HERE.parent/"canonical_finite_cycle_join/README.md",HERE.parent/"physical_pulse_transfer/README.md"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_current_scalar_and_scale_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "six_historical_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Current scalar recurrence, retained scale losses, source bindings or scope changed")


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
