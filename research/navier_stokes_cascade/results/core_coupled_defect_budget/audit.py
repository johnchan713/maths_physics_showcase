#!/usr/bin/env python3
"""Freeze conditional coupled scalar budgets and their derivative costs."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS="conditional-coupled-pressure-axial-scalar-and-whole-moment-budget"
DEPENDENCIES=("core_scalar_cycle_budget","core_dynamic_mean_recheck","core_scalar_target_pressure_join",
              "core_reserved_geometry_recheck","current_target_cycle_bridge","compact_mean_foundation")
FLAGS={
    "coupled_pressure_and_axial_scalar_polynomial_budgets_derived":True,
    "complete_covariance_and_fast_radial_correlations_retained":True,
    "whole_moving_pressure_product_differentiated":True,
    "axial_moment_requires_up_to_four_extra_incoming_target_jets":True,
    "existing_postslow_order_not_presented_as_new":True,
    "actual_seed_and_inverse_inputs_independently_verified":False,
    "uniform_infinite_stage_or_all_order_constant_proved":False,
    "complete_pointwise_residual_certified":False,
    "actual_profile_or_accepted_tail_changed":False,
    "unconditional_pulse_or_blowup_proved":False,
}
REMAINING=[
    "independently verify the canonical seed and actual current wave/inverse/source/edge hypotheses and their fixed-profile polynomial constants",
    "check the coupled derivative and scale-loss budgets against actual tolerances on the existing common domain",
    "join complete physical support, pressure memory, covariance, viscosity, cutoffs and all nonlinear interactions in the pointwise residual",
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
    paths += [HERE.parent/"core_scalar_cycle_budget/proof_checks.py",HERE.parent/"core_dynamic_mean_recheck/proof_checks.py",
              HERE.parent/"complete_covariance_initial_means/means.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_coupled_moment_and_budget_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "six_historical_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Coupled current moments, polynomial budgets, derivative ledger, bindings or scope changed")


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
