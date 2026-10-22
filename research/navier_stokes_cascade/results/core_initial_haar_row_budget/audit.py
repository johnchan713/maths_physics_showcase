#!/usr/bin/env python3
"""Freeze conditional complete initial Haar-row and weighted target budgets."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=("core_signed_target_edge_budget","core_primary_root_budget","core_fast_inverse_budget",
              "core_initial_scalar_seed_budget","compact_mean_foundation","canonical_initial_mean_join")
STATUS="initial-complete-weighted-Haar-row-budget-and-sharper-average-rate"
FLAGS={
    "current_initial_Haar_rows_derived_from_complete_fast_slow_covariance_pressure_state":True,
    "sharper_Haar_order_two_minus_two_kappa_retains_full_row_two_minus_three_kappa":True,
    "weighted_operator_route_uses_edge_directed_primitives_and_interior_A1_support":True,
    "actual_cutoff_plateau_gap_retained_in_weighted_primitive_constant":True,
    "independent_initial_current_row_norm_replaced_by_lower_level_budgets":True,
    "whole_current_pressure_and_nonzero_fast_correlations_retained":True,
    "canonical_ODE_core_background_and_all_fixed_map_premises_independently_closed":False,
    "sharper_Haar_rate_promoted_to_full_unaveraged_row":False,
    "initial_real_and_zero_Haar_reductions_extended_to_later_complex_states":False,
    "actual_nonzero_angular_sources_or_complete_pointwise_residual_verified":False,
    "historical_evidence_or_fixed_physical_parameters_changed":False,
    "unconditional_pulse_bounds_infinite_realization_or_blowup_proved":False,
}
REMAINING=[
    "verify the actual nonzero-angular regular source through the signed prefix and first sourced inverse, including full-path supports and all viscosity/cutoff/interaction terms",
    "retain and independently discharge canonical root/ODE, exact core/background and fixed profile/map construction premises before an unconditional all-pulse statement",
    "join the complete pointwise residual and infinite realization without promoting the initial Haar rate to a full or later-complex-state improvement",
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
    paths += [HERE.parent/"core_primary_root_budget/budgets.py",HERE.parent/"core_signed_target_edge_budget/budgets.py",
              HERE.parent/"core_initial_scalar_seed_budget/proof_checks.py",HERE.parent/"complete_covariance_initial_means/means.py",
              HERE.parent/"core_scalar_target_pressure_join/proof_checks.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_current_Haar_row_and_weighted_norm_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "six_historical_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Initial Haar-row/pressure equations, weighted route, historical bindings or scoped rate changed")


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
