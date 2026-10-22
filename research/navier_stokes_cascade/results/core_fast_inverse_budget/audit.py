#!/usr/bin/env python3
"""Freeze explicit current-residual fast-map and actual axial-error budgets."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS="explicit-fast-inverse-actual-axial-repair-and-current-residual-jet-budget"
DEPENDENCIES=("core_slow_inverse_budget","core_coupled_defect_budget","core_scalar_cycle_budget",
              "compact_mean_foundation","compact_mean_corrections","current_target_cycle_bridge")
FLAGS={
    "explicit_current_residual_temporal_and_compact_norm_bounds_derived":True,
    "actual_A1_velocity_retained_without_spurious_graph_loss":True,
    "fast_time_error_is_A1_of_current_residual_with_correct_sign":True,
    "direct_flat_error_bound_avoids_temporal_inverse_loss":True,
    "channel_specific_fast_and_joined_slow_input_jet_counts_retained":True,
    "qualitative_inverse_existence_and_flatness_not_presented_as_new":True,
    "actual_canonical_seed_and_current_residual_input_norms_verified":False,
    "actual_sourced_inverse_inputs_independently_verified":False,
    "uniform_infinite_stage_or_all_order_constant_proved":False,
    "actual_profile_tilt_width_or_accepted_tail_changed":False,
    "complete_residual_unconditional_pulse_or_blowup_proved":False,
}
REMAINING=[
    "verify actual canonical seed and sourced/current prefast residual seminorms and edge supports at the displayed derivative orders",
    "join those inputs to the explicit coupled fast/slow scalar budgets on the existing common domain with retained core/background floors",
    "verify complete pressure memory, all actual A0/A1 fields, viscosity, cutoffs and old/new interactions in the pointwise residual before any infinite realization",
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
    paths += [HERE.parent/"core_coupled_defect_budget/proof_checks.py",HERE.parent/"core_scalar_cycle_budget/proof_checks.py",
              HERE.parent/"compact_mean_corrections/operators.py",HERE.parent/"complete_covariance_initial_means/means.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_fast_map_current_error_and_budget_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "six_historical_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Actual fast velocity, current residual error, derivative budget, bindings or scope changed")


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
