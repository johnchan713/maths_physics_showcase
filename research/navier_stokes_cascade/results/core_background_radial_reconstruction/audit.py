#!/usr/bin/env python3
"""Freeze full radial-background reconstruction and low-order phase input joins."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=("core_background_finite_head_budget","global_first_coefficient_c3","global_second_coefficient_h23",
              "global_axial_source_envelope","global_swirl_envelope","background_jet_foundation","first_coefficient_c3_boundary")
STATUS="full-radial-primitive-reconstruction-and-finalized-first-coefficient-phase-input-join"
FLAGS={
    "full_E_U_radial_velocity_and_primitive_low_jets_derived":True,
    "first_angular_radial_jet_retains_the_second_angular_U_input":True,
    "factorial_to_raw_cost_and_all_inverse_radius_terms_retained":True,
    "first_tuple_norm_uses_same_global_pressure_viscosity_and_exact_linear_repair":True,
    "entire_native_band_chart_coverage_is_checked_before_phase_substitution":True,
    "actual_Borel_Euler_norm_and_normalized_schedule_certified_here":False,
    "actual_unchanged_domain_normal_floor_substitution_verified":False,
    "tangent_propagator_all_supports_or_complete_covariance_newly_closed":False,
    "historical_evidence_profiles_carriers_cutoffs_or_accepted_bands_changed":False,
    "complete_residual_unconditional_pulses_or_blowup_proved":False,
}
REMAINING=[
    "verify actual full native-to-band domain coverage for the existing leading and finalized coefficient rows",
    "certify the same Borel cutoff Euler K2 and normalized retained diagonal without identifying it with the radial coefficient cutoff",
    "substitute explicit Cphi and actual growth/shear/geometry data into beta*Lreq^2>=144*r0^2*Cphi^2 on the unchanged selector with every original check",
    "discharge tangent propagator actual-root waves supports and fixed maps before the first sourced complete covariance fresh pressure/axial means and complete residual join",
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
    paths += [HERE.parent/"core_background_finite_head_budget/budgets.py",
              HERE.parent/"global_first_coefficient_c3/global_first.py",
              HERE.parent/"global_second_coefficient_h23/global_second.py",
              HERE.parent/"modulation_repair_all_order/map_checks.py",
              HERE.parent/"core_all_order_foundation/coordinate_checks.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_radial_reconstruction_memory_and_join_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "seven_historical_records_bound":len(bound)==7,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Radial reconstruction primitive memory domain or lower-premise bindings changed")


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
