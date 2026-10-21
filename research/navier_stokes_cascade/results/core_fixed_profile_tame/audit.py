#!/usr/bin/env python3
"""Freeze the fixed-profile quantifier correction and tame estimates."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS="fixed-profile-all-order-tame-bounds-and-quantifier-correction"
DEPENDENCIES=("core_actual_compensator_gate","core_pressure_shape_recheck","core_reserved_geometry_recheck","core_slow_jet_loss_recheck")
FLAGS={
    "prior_profile_independent_counterfamilies_remain_valid":True,
    "numerical_profile_required_for_asymptotic_existence_corrected":True,
    "fixed_smooth_profile_gives_finite_constants_at_every_fixed_order":True,
    "all_fixed_radial_and_slow_composition_majorants_derived":True,
    "all_six_moving_density_pressure_terms_bounded_conditionally":True,
    "slow_time_mismatch_and_axial_viscosity_both_retained":True,
    "actual_profile_or_accepted_tail_changed":False,
    "uniformity_over_all_profiles_or_infinite_orders_proved":False,
    "actual_current_target_higher_jets_verified":False,
    "complete_residual_improvement_or_blowup_proved":False,
}
REMAINING=[
    "prove actual linked current-target and incoming-field higher derivative bounds with band and stage dependence retained",
    "evaluate their strict epsilon margins or complete linked cancellations on the existing common domain",
    "join all generated supports, complete covariance, full pressure memory, viscosity, cutoffs and nonlinear interactions in the full residual",
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
    paths += [HERE.parent/"core_reserved_geometry_recheck/proof_checks.py",HERE.parent/"core_pressure_shape_recheck/proof_checks.py",HERE.parent/"complete_covariance_initial_means/means.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_tame_and_pressure_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "four_historical_records_bound":len(bound)==4,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Fixed-profile estimates, quantifiers, source bindings or scope changed")


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
