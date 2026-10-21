#!/usr/bin/env python3
"""Freeze band-uniform reserved geometry without certifying all field norms."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS="uniform-reserved-coordinate-compensator-jets-and-moving-support-costs-rechecked"
DEPENDENCIES=("core_dynamic_mean_recheck","compact_mean_foundation",
              "physical_pulse_foundation","core_primary_tail_recheck",
              "coefficient_support_audit","source_envelope_edge_audit")
FLAGS={
    "explicit_band_uniform_coordinate_jet_constants_derived":True,
    "complete_moving_pressure_compensator_jet_cost_supplied":True,
    "reserved_support_shell_and_log_weight_constants_supplied":True,
    "dynamic_budget_geometry_and_compensator_inputs_joined":True,
    "whole_pressure_product_and_both_background_tails_retained":True,
    "common_primary_tail_fixed_tilt_and_width_unchanged":True,
    "rectangular_union_hull_has_positive_weight_floor":False,
    "coordinate_analytic_extension_is_heat_profile_analytic_extension":False,
    "all_actual_generated_field_norms_certified":False,
    "reserved_hull_gate_covers_wider_changing_wave_supports":False,
    "complete_residual_improvement_or_blowup_proved":False,
    "order_six_diagnostics_prove_all_orders":False,
}
REMAINING=[
    "supply band-uniform actual current-target and generated-field input norms to the dynamic budget with its now explicit reserved geometry and compensator constants",
    "handle wider wave supports with full-path edge/envelope bounds and retain all fast correlations and complete signed-target cancellation",
    "prove complete pointwise residual with viscosity, every cutoff and interaction on the same tail",
]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        path=HERE.parent/name/"evidence.json";record=json.loads(path.read_text())
        bound[str(path.relative_to(PROJECT))]=sha(path)
        for rel,expected in record["source_hashes"].items():
            checks[name+":"+rel]=sha(PROJECT/rel)==expected
    return bound,checks


def sources():
    paths=[HERE/n for n in ("README.md","proof_checks.py","audit.py","test_audit.py","protocol.json","requirements.txt",".gitignore")]
    paths += [HERE.parent/n/"README.md" for n in DEPENDENCIES]
    paths += [HERE.parent/"core_dynamic_mean_recheck/proof_checks.py",HERE.parent/"complete_covariance_initial_means/means.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),
                flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_coordinate_and_support_controls_pass":all(result["checks"].values()),
                       "six_frozen_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Coordinate jets, support costs, source binding or scientific scope changed")


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--verify-record",type=Path)
    args=parser.parse_args(argv)
    if args.output.exists():raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record=generate();validate(record)
    if args.verify_record:validate(json.loads(args.verify_record.read_text()))
    with args.output.open("x") as output:output.write(json.dumps(record,indent=2,sort_keys=True)+"\n")
    print(f"{STATUS}: {len(record['review']['checks'])} exact controls; {len(record['inherited_source_checks'])} inherited source checks")


if __name__=="__main__":main()
