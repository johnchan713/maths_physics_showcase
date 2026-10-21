#!/usr/bin/env python3
"""Freeze scalar/current-target routing into the higher-jet pressure budget."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS="current-scalar-projected-target-and-higher-jet-pressure-join"
DEPENDENCIES=("core_fixed_profile_tame","core_dynamic_mean_recheck","current_target_moment_audit",
              "current_target_cycle_bridge","core_pressure_shape_recheck","core_slow_jet_loss_recheck")
FLAGS={
    "current_scalar_and_projected_signed_target_distinguished":True,
    "two_channels_reconstruct_full_Haar_angular_equation":True,
    "current_scalar_higher_jets_routed_from_pure_T_flux_and_integrated_residual":True,
    "epsilon_inverse_cost_of_scalar_jet_recovery_retained":True,
    "slow_time_and_axial_viscous_pressure_budgets_joined_conditionally":True,
    "postfast_scalar_correlations_and_complete_graph_curl_retained":True,
    "projected_target_bound_alone_controls_discarded_scalar_channel":False,
    "operator_projection_control_is_actual_canonical_counterexample":False,
    "actual_incoming_norms_or_complete_cycle_inputs_verified":False,
    "actual_profile_or_accepted_tail_changed":False,
    "unconditional_pulse_residual_improvement_or_blowup_proved":False,
}
REMAINING=[
    "prove actual current pure-T scalar flux/cancellation and integrated angular residual jets in the displayed higher-jet pressure inequalities",
    "retain the separate projected-target recurrence and full incoming current-field, edge and inverse estimates at each canonical substep",
    "join all generated supports, full covariance and pressure memory, viscosity, cutoffs and nonlinear interactions in the complete residual",
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
    paths += [HERE.parent/"core_dynamic_mean_recheck/proof_checks.py",HERE.parent/"complete_covariance_initial_means/means.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_scalar_and_pressure_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "six_historical_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Scalar/current-target routing, source bindings or scope changed")


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
