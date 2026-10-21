#!/usr/bin/env python3
"""Freeze the actual-compensator data gate and its exact conditional bounds."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS="actual-compensator-profile-data-gate-and-conditional-shape-bounds"
DEPENDENCIES=("core_pressure_shape_recheck","core_reserved_geometry_recheck","core_slow_jet_loss_recheck")
FLAGS={
    "actual_profile_contract_traced":True,
    "actual_profile_formula_or_numeric_seminorms_found":False,
    "scalar_mass_moment_and_support_suffice_for_shape_bound":False,
    "positivity_alone_suffices_for_shape_derivative_bounds":False,
    "conditional_first_and_second_slow_shape_bounds_derived":True,
    "parent_pointwise_pressure_and_radial_coefficients_retained":True,
    "diagnostic_pressure_bump_promoted_to_canonical_data":False,
    "actual_high_target_jet_budget_closed":False,
    "complete_residual_improvement_or_blowup_proved":False,
}
REMAINING=[
    "freeze the actual physical rhohat formula or supply its Euler-shape seminorms B0,B1,... at every required fixed order",
    "insert those constants and the actual current-target derivatives into the pointwise pressure, radial and moving-commutator budgets",
    "join generated-source supports, viscosity, cutoffs and all interactions in the complete physical residual",
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
    paths += [HERE.parent/"compact_mean_corrections/README.md",HERE.parent/"compact_mean_corrections/reviews.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_profile_obstruction_and_budget_controls_pass":all(result["checks"].values()),
                       "three_historical_records_bound":len(bound)==3,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Actual compensator data gate, conditional bounds, source binding or scope changed")


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
