#!/usr/bin/env python3
"""Freeze current signed-target two-edge and scalar-channel estimates."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=("core_primary_root_budget","core_initial_scalar_seed_budget",
              "core_scalar_target_pressure_join","compact_mean_foundation","canonical_initial_mean_join")
STATUS="current-signed-target-two-edge-primitive-and-weighted-row-input-budget"
FLAGS={
    "current_signed_map_identified_by_exact_Haar_zero_moment":True,
    "explicit_two_edge_primitive_constant_preserves_full_flat_weight":True,
    "current_angular_and_whole_axial_pressure_scalar_channels_retained":True,
    "nonlinear_chart_jets_and_reciprocal_weight_cost_retained":True,
    "unweighted_bounds_even_with_exact_moments_do_not_imply_weighted_target":True,
    "counterexample_is_actual_canonical_field_or_disproof":False,
    "actual_weighted_current_Haar_row_norms_independently_verified":False,
    "actual_all_sourced_inputs_and_supports_verified":False,
    "complete_residual_or_infinite_realization_verified":False,
    "historical_evidence_or_fixed_physical_parameters_changed":False,
    "unconditional_pulse_bounds_or_blowup_proved":False,
}
REMAINING=[
    "prove inequality (8): weighted current complete regular Haar-row seminorms from the actual pressure/fast/slow equations on every accepted band",
    "retain exact projected moments, reserved-hull and whole-path support hypotheses while verifying the actual sourced nonzero rows",
    "discharge canonical ODE/core/background inputs and the full moving-pressure, axial, viscosity, A0/A1, cutoff and nonlinear residual before realization or blow-up",
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
    paths += [HERE.parent/"core_primary_root_budget/budgets.py",HERE.parent/"core_scalar_target_pressure_join/proof_checks.py",
              HERE.parent/"core_initial_scalar_seed_budget/proof_checks.py",HERE.parent/"complete_covariance_initial_means/means.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_current_projection_edge_and_budget_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "five_historical_records_bound":len(bound)==5,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Current signed target edge/scalar/chart budget, provenance or weighted-input scope changed")


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
