#!/usr/bin/env python3
"""Freeze explicit conditional canonical-root and signed coefficient budgets."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=("core_primary_aggregate_atlas","core_initial_scalar_seed_budget",
              "canonical_primary_foundation","primary_covariance_foundation","weighted_stress_foundation")
STATUS="quantitative-canonical-root-edge-curl-and-original-root-signed-budget"
FLAGS={
    "root_reciprocal_and_flat_constants_given_by_positive_recurrences":True,
    "normalized_outer_target_allowed_polynomial_edge_pole":True,
    "complete_curl_and_raw_factorial_conversion_retained":True,
    "root_budgets_substituted_into_aggregate_and_whole_pressure_moment":True,
    "whole_moment_requires_underlying_inputs_through_n_plus_12":True,
    "signed_mapping_uses_original_canonical_root_and_keeps_signed_square":True,
    "canonical_ODE_column_normal_and_exact_core_inputs_independently_discharged":False,
    "actual_current_signed_target_or_all_sourced_inputs_verified":False,
    "complete_pointwise_residual_or_infinite_realization_verified":False,
    "historical_evidence_or_frozen_physical_parameters_changed":False,
    "unconditional_pulse_bounds_or_blowup_proved":False,
}
REMAINING=[
    "discharge the actual canonical column H_j, target T_j, homogeneous/cutoff W_j, full normal/curl D_j and strict value-floor inputs against retained exact core/background premises",
    "verify the current signed target and then the actual sourced residual weighted norms and complete path/support hypotheses",
    "close the complete pointwise residual with moving pressure memory, A0/A1, axial means, viscosity, cutoffs and nonlinear interactions before infinite realization or blow-up",
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
    paths += [HERE.parent/"core_primary_aggregate_atlas/proof_checks.py",HERE.parent/"core_initial_scalar_seed_budget/proof_checks.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_root_edge_signed_and_chain_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "five_historical_records_bound":len(bound)==5,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Root/edge/curl budgets, original-root map, provenance or conditional scope changed")


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
