#!/usr/bin/env python3
"""Freeze existing native chart coverage and full-pressure support-aware norms."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=("core_borel_cutoff_certificate","core_background_radial_reconstruction","global_swirl_envelope",
              "global_axial_source_envelope","background_jet_foundation","auxiliary_geometry_audit","core_primary_tail_recheck")
STATUS="existing-native-chart-cover-and-support-aware-full-velocity-phase-inputs"
FLAGS={
    "existing_doubled_box_s_and_native_log_distance_bounds_explicit":True,
    "global_leading_velocity_rows_cover_the_entire_native_norm_hull":True,
    "support_aware_primitive_volume_retains_all_inverse_radius_and_raw_factorial_costs":True,
    "first_pressure_moment_repair_keeps_original_full_future_Xv":True,
    "derived_native_band_enclosures_join_complete_Cphi_without_reselection":True,
    "actual_existing_native_neighborhood_and_negative_T_extension_verified":False,
    "growth_stress_or_cone_floors_extended_to_cartesian_norm_hull":False,
    "actual_infinite_Borel_schedule_or_unchanged_selector_verified":False,
    "homogeneous_propagator_all_supports_or_full_covariance_newly_closed":False,
    "historical_evidence_profiles_cutoffs_scales_or_accepted_bands_changed":False,
    "complete_residual_unconditional_pulses_or_blowup_proved":False,
}
REMAINING=[
    "verify actual unchanged band radial-clearance/native-log-neighborhood inequalities and frozen physical representative coverage or original negative-T extension",
    "retain actual leading value and edge floors on the covered coupled domain, not the larger Cartesian norm hull",
    "identify the original infinite Borel prescription with its certified cutoff/normalized bounds and retain all original physical weighted-stress and residual schedule gates",
    "substitute complete Cphi and actual gaps into the unchanged normal-floor selector inequality, then discharge propagator supports fixed maps and complete sourced covariance fresh pressure/axial means before residual claims",
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
    paths += [HERE.parent/"core_borel_cutoff_certificate/budgets.py",
              HERE.parent/"core_background_radial_reconstruction/budgets.py",
              HERE.parent/"global_swirl_envelope/swirl_bounds.py",
              HERE.parent/"global_axial_source_envelope/envelopes.py",
              HERE.parent/"core_primary_tail_recheck/proof_checks.py",
              HERE.parent/"annular_exterior_foundation/README.md"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_native_geometry_support_and_full_pressure_join_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "seven_historical_records_bound":len(bound)==7,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Native cover support-aware norm pressure endpoint or historical scope changed")


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
