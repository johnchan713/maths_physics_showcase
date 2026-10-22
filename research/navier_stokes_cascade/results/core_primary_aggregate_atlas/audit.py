#!/usr/bin/env python3
"""Freeze common-atlas aggregate bounds and matched physical normalizations."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS="common-atlas-primary-aggregate-and-matched-band-normalization-budget"
DEPENDENCIES=("core_initial_scalar_seed_budget","auxiliary_geometry_audit","canonical_primary_foundation",
              "physical_pulse_foundation","core_fast_inverse_budget","core_slow_inverse_budget")
FLAGS={
    "aggregate_norms_derived_from_actual_closed_separation_and_native_bounds":True,
    "fixed_width_edge_unweighting_constants_retained":True,
    "common_covering_raw_derivative_cost_quantified":True,
    "coarsest_normalization_matches_carrier_and_fast_time_floor":True,
    "arbitrary_reference_requires_explicit_carrier_and_time_constants":True,
    "actual_atlas_geometry_not_presented_as_new":True,
    "per_label_canonical_ODE_root_and_edge_estimates_independently_reproved":False,
    "full_signed_sourced_seed_or_pointwise_residual_verified":False,
    "uniform_infinite_stage_or_all_order_constant_proved":False,
    "actual_profile_radius_tilt_width_or_accepted_tail_changed":False,
    "unconditional_pulse_or_blowup_proved":False,
}
REMAINING=[
    "retain and independently discharge the per-label canonical ODE/root/edge bounds relative to the exact core/background premises",
    "extend the now aggregated initial scalar estimates through the actual signed/sourced seed and source support obligations",
    "verify full pressure memory, all A0/A1 repairs, viscosity, cutoffs and interactions in the complete pointwise residual before infinite realization or blow-up claims",
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
    paths += [HERE.parent/"auxiliary_geometry_audit/geometry.py",HERE.parent/"core_initial_scalar_seed_budget/proof_checks.py",
              HERE.parent/"core_fast_inverse_budget/proof_checks.py",HERE.parent/"core_slow_inverse_budget/proof_checks.py",
              HERE.parent/"complete_covariance_initial_means/means.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_atlas_edge_normalization_and_aggregate_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "six_historical_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Aggregate atlas bounds, width/reference constants, bindings or scope changed")


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
