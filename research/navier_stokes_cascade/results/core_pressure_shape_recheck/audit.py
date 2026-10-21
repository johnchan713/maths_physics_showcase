#!/usr/bin/env python3
"""Freeze the pointwise pressure-shape gate without changing the compensator."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS="pointwise-pressure-shape-gate-and-persistent-radial-target-jet-rechecked"
DEPENDENCIES=("core_slow_jet_loss_recheck","core_reserved_geometry_recheck",
              "core_dynamic_mean_recheck","compact_mean_foundation",
              "current_target_moment_audit","complete_covariance_initial_means")
FLAGS={
    "target_independent_positive_reference_density_derived":True,
    "pointwise_pressure_high_jet_requires_full_density_shape_match":True,
    "same_mass_second_moment_and_positivity_do_not_imply_shape_match":True,
    "nonzero_signed_target_response_can_survive_scalar_cancellation":True,
    "whole_moving_density_commutator_retained":True,
    "unit_compensator_cannot_remove_radial_cubic_response_coefficient":True,
    "sharp_half_L1_shape_primitive_estimate_derived":True,
    "complete_covariance_and_both_base_pressure_tails_retained":True,
    "actual_compensator_or_any_frozen_profile_changed":False,
    "actual_compensator_shape_match_verified":False,
    "actual_linked_canonical_sequence_refuted":False,
    "actual_high_target_jet_bounds_or_full_residual_closed":False,
    "unconditional_blowup_proved":False,
    "finite_diagnostics_are_actual_continuum_suprema":False,
}
REMAINING=[
    "identify or bound the actual compensator's full density-shape primitive, beyond its scalar second moment",
    "prove actual current-target high-jet budgets for the pointwise mismatch, moving commutators and unavoidable radial coefficient, or their complete linked-cycle cancellations",
    "join all generated-source support/edge inputs and complete viscosity/cutoff interactions to the full pointwise residual",
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
    paths += [HERE.parent/"core_slow_jet_loss_recheck/proof_checks.py",HERE.parent/"complete_covariance_initial_means/means.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_shape_and_operator_controls_pass":all(result["checks"].values()),
                       "six_historical_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Pressure shape, radial response, source binding or scientific scope changed")


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
