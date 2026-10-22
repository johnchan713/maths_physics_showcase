#!/usr/bin/env python3
"""Freeze the complete-primary to initial scalar-prefix quantitative budget."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS="quantitative-complete-primary-initial-scalar-prefix-with-pressure-channel-split"
DEPENDENCIES=("core_fast_inverse_budget","core_slow_inverse_budget","canonical_initial_mean_join",
              "canonical_primary_foundation","core_coupled_defect_budget","compact_mean_corrections")
FLAGS={
    "initial_fast_residual_inputs_bounded_from_complete_primary_norms":True,
    "slow_compensator_removed_only_from_zero_Haar_fast_channel":True,
    "whole_moving_Haar_pressure_moment_and_actual_A1_retained":True,
    "explicit_initial_S2_and_sharper_Haar_velocity_budgets_derived":True,
    "C1_initial_case_checked_separately_from_later_cycle_halfline":True,
    "known_qualitative_initial_S2_order_not_presented_as_new":True,
    "primary_aggregate_root_and_edge_hypotheses_independently_verified":False,
    "full_signed_sourced_canonical_finite_cycle_seed_verified":False,
    "complete_pointwise_mean_or_full_residual_verified":False,
    "uniform_infinite_stage_or_all_order_constant_proved":False,
    "actual_profile_tilt_width_or_accepted_tail_changed":False,
    "unconditional_pulse_or_blowup_proved":False,
}
REMAINING=[
    "verify actual primary aggregate/root/edge bounds and retained core/background floors at the displayed coefficient orders",
    "extend the initial scalar-prefix chain through the actual signed/sourced canonical prefix with current source and support bounds",
    "join complete pressure memory, actual A0/A1 fields, viscosity, cutoffs and old/new interactions in the pointwise residual before infinite realization or blow-up claims",
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
    paths += [HERE.parent/"core_fast_inverse_budget/proof_checks.py",HERE.parent/"core_slow_inverse_budget/proof_checks.py",
              HERE.parent/"core_coupled_defect_budget/proof_checks.py",HERE.parent/"complete_covariance_initial_means/means.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_initial_pressure_prefix_and_norm_chain_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "six_historical_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Initial pressure channels, current prefix, primary norm chain, bindings or scope changed")


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
