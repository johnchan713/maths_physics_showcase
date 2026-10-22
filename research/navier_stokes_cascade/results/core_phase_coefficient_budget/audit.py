#!/usr/bin/env python3
"""Freeze full-background to frozen-phase coefficient and value-gate recipes."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=("core_signed_prefix_source_budget","core_first_sourced_inverse_budget","actual_band_pulse_bounds",
              "core_primary_tail_recheck","core_primary_aggregate_atlas")
STATUS="full-background-frozen-phase-coefficient-recipes-varying-growth-and-normal-floor-gate"
FLAGS={
    "normal_moving_normal_K_and_eikonal_jets_derived_from_full_background_inputs":True,
    "logarithmic_profile_rate_retains_varying_growth_at_all_fixed_orders":True,
    "full_normal_floor_reduced_to_explicit_whole_path_value_inequalities":True,
    "five_variable_raw_to_factorial_conversion_and_inverse_source_interfaces_joined":True,
    "frozen_nonzero_rounding_radial_flow_and_axial_transport_retained":True,
    "actual_core_background_and_phase_value_constants_independently_verified":False,
    "actual_unchanged_accepted_domain_value_gate_substitution_verified":False,
    "homogeneous_tangent_propagator_or_all_support_premises_independently_closed":False,
    "historical_evidence_or_actual_profiles_carriers_or_accepted_bands_changed":False,
    "complete_residual_unconditional_pulses_or_blowup_proved":False,
}
REMAINING=[
    "derive actual full-background jet and phase-comparison constants and verify both value inequalities (10) on the unchanged accepted domain",
    "discharge the homogeneous tangent propagator original-root wave and remaining full-path support and fixed-map premises",
    "join quantitative first sourced complete covariance fresh full pressure means and the complete physical residual including additive tails before unconditional claims",
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
    paths += [HERE.parent/"core_signed_prefix_source_budget/budgets.py",HERE.parent/"core_first_sourced_inverse_budget/budgets.py",
              HERE.parent/"modulation_repair_all_order/map_checks.py",HERE.parent/"core_all_order_foundation/coordinate_checks.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_phase_growth_floor_and_box_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "five_historical_records_bound":len(bound)==5,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Actual phase coefficients varying growth normal value gate lower premises or historical bindings changed")


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
