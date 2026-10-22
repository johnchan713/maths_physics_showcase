#!/usr/bin/env python3
"""Freeze explicit fixed-profile slow-map constants and coupled substitution."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS="explicit-fixed-profile-slow-inverse-and-coupled-scalar-substitution"
DEPENDENCIES=("core_coupled_defect_budget","core_fixed_profile_tame","core_reserved_geometry_recheck",
              "compact_mean_foundation","compact_mean_corrections","canonical_initial_mean_join")
FLAGS={
    "three_distinct_physical_response_scale_factors_retained":True,
    "five_current_signed_rows_and_radial_potential_derived":True,
    "explicit_finite_jet_profile_and_reciprocal_bounds_derived":True,
    "abstract_slow_map_constants_replaced_in_coupled_scalar_budget":True,
    "qualitative_inverse_existence_not_presented_as_new":True,
    "one_extra_radial_input_jet_and_full_pressure_moment_cost_retained":True,
    "actual_core_background_floor_or_canonical_seed_independently_verified":False,
    "actual_sourced_and_fast_inverse_inputs_independently_verified":False,
    "uniform_infinite_stage_or_all_order_constant_proved":False,
    "actual_profile_tilt_width_or_accepted_tail_changed":False,
    "complete_residual_unconditional_pulse_or_blowup_proved":False,
}
REMAINING=[
    "verify actual canonical seed and sourced/fast input norms against the explicit coupled scalar and slow-map budgets",
    "retain the actual core/background floor, fixed profile and support hypotheses and check required tolerances on the unchanged common tail",
    "join complete pressure memory, covariance, support edges, viscosity, cutoffs and every interaction in the pointwise residual before any infinite realization",
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
    paths += [HERE.parent/"core_coupled_defect_budget/proof_checks.py",HERE.parent/"core_fixed_profile_tame/proof_checks.py",
              HERE.parent/"core_reserved_geometry_recheck/proof_checks.py",HERE.parent/"compact_mean_foundation/mean_maps.py",
              HERE.parent/"complete_covariance_initial_means/means.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_response_profile_scaling_jet_and_budget_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "six_historical_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Explicit inverse profiles, physical normalization, jet costs, bindings or scope changed")


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
