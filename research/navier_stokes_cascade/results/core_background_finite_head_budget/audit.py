#!/usr/bin/env python3
"""Freeze explicit finite-head background and unchanged selector comparison budgets."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=("core_phase_coefficient_budget","background_jet_foundation","background_summation",
              "core_fixed_profile_tame","core_primary_tail_recheck","core_signed_prefix_source_budget")
STATUS="finite-head-full-background-jet-constants-and-unchanged-normal-selector-compatibility"
FLAGS={
    "whole_q_interval_finite_head_and_geometric_tail_constants_explicit":True,
    "radial_cutoff_curl_commutator_and_Euler_conversion_cost_retained":True,
    "finite_real_chart_composition_and_full_background_coefficient_substitution_explicit":True,
    "Cphi_contains_all_seven_full_normal_error_groups":True,
    "unchanged_selector_join_reduced_to_explicit_original_length_compatibility":True,
    "actual_leading_coefficient_cutoff_and_schedule_inputs_certified_here":False,
    "actual_unchanged_domain_normal_floor_substitution_verified":False,
    "global_pressure_moment_identities_or_tangent_propagator_newly_closed":False,
    "historical_evidence_profiles_carriers_cutoffs_or_accepted_bands_changed":False,
    "complete_residual_unconditional_pulses_or_blowup_proved":False,
}
REMAINING=[
    "certify actual leading through-two and finalized coefficient-one through-one tuple jets, cutoff Euler norms, normalized retained schedule, growth/shear and geometry data",
    "evaluate all seven Cphi groups and beta*Lreq^2>=144*r0^2*Cphi^2 on the unchanged existing selector, retaining all original checks",
    "discharge homogeneous tangent propagator actual-root wave and remaining whole-path support and fixed-map premises",
    "join first sourced complete covariance fresh full-pressure moments and complete physical residual including all additive tails before unconditional claims",
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
    paths += [HERE.parent/"core_phase_coefficient_budget/budgets.py",
              HERE.parent/"core_reserved_geometry_recheck/proof_checks.py",
              HERE.parent/"core_primary_tail_recheck/proof_checks.py",
              HERE.parent/"modulation_repair_all_order/map_checks.py",
              HERE.parent/"core_all_order_foundation/coordinate_checks.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_finite_head_cutoff_chart_and_selector_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "six_historical_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Finite-head coefficient comparison selector lower premises or historical bindings changed")


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
