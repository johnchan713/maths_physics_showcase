#!/usr/bin/env python3
"""Freeze implemented Borel cutoff Euler jets and unchanged-scale certificates."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=("core_background_radial_reconstruction","core_background_finite_head_budget","background_summation",
              "background_jet_foundation","cutoff_repair_derivative_budget","core_primary_tail_recheck")
STATUS="implemented-borel-cutoff-Euler-certificate-and-unchanged-schedule-phase-gates"
FLAGS={
    "implemented_squared_flat_Borel_cutoff_Euler_K2_at_most_4928":True,
    "all_fixed_order_real_cutoff_constants_have_explicit_finite_recurrences":True,
    "ordinary_second_derivative_not_substituted_for_second_Euler_jet":True,
    "supplied_schedule_prefix_verifier_keeps_every_row_and_never_selects_scales":True,
    "full_Cphi_cutoff_dependence_and_necessary_sufficient_recipe_tests_explicit":True,
    "every_qualitative_historical_chi_identified_with_this_implementation":False,
    "actual_infinite_cutoff_schedule_or_all_original_schedule_gates_verified":False,
    "actual_native_coverage_gaps_and_unchanged_selector_substitution_verified":False,
    "tangent_propagator_supports_or_full_covariance_newly_closed":False,
    "historical_evidence_profiles_cutoffs_scales_or_accepted_bands_changed":False,
    "complete_residual_unconditional_pulses_or_blowup_proved":False,
}
REMAINING=[
    "identify the actual summed background chi and scale prescription with the implemented analytic certificate without changing either",
    "verify the retained infinite normalized schedule with all original physical potential weighted-stress and residual requirements, not a finite diagnostic prefix",
    "verify actual whole-chart native coverage and leading gaps, then substitute Cphi(4928) or a sharper certified norm in the unchanged selector inequality",
    "discharge tangent propagator supports fixed maps and original-root wave inputs before complete sourced covariance fresh pressure/axial means and the full residual join",
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
    paths += [HERE.parent/"core_background_radial_reconstruction/budgets.py",
              HERE.parent/"core_background_finite_head_budget/budgets.py",
              HERE.parent/"background_summation/checks.py",
              HERE.parent/"outer_pressure_pilot/schedule.py",
              HERE.parent/"modulation_repair_all_order/map_checks.py",
              HERE.parent/"core_all_order_foundation/coordinate_checks.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_cutoff_Euler_source_schedule_and_phase_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "six_historical_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Implemented Borel cutoff Euler schedule scope or historical bindings changed")


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
