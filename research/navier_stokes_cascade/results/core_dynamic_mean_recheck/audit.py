#!/usr/bin/env python3
"""Freeze dynamic mean upper budgets with current correlation limitations."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from proof_checks import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS="current-dynamic-mean-norm-budgets-and-correlation-rate-gates-rechecked"
DEPENDENCIES=("core_covariance_mean_recheck","current_target_moment_audit",
              "current_target_cycle_bridge","complex_mean_relaxation",
              "complete_covariance_initial_means","core_primary_tail_recheck")
FLAGS={
    "nonzero_old_Haar_means_and_complete_covariance_changes_retained":True,
    "explicit_finite_jet_scalar_upper_budgets_derived":True,
    "whole_moving_pressure_product_and_extra_derivative_cost_retained":True,
    "old_mean_and_same_phase_wave_cross_rate_controls_derived":True,
    "all_six_covariance_entries_kept_in_moment_or_target_flux":True,
    "fresh_fixed_wave_slow_map_has_zero_covariance_difference":True,
    "both_complete_base_pressure_tails_retained":True,
    "common_primary_tail_unchanged":True,
    "common_radial_hull_required_for_all_covariance_changes":True,
    "generic_correlations_are_new_square_order":False,
    "norm_gate_verified_uniformly_on_all_actual_scales":False,
    "actual_five_map_or_conditional_cycle_refuted_by_control":False,
    "pointwise_full_residual_improvement_or_blowup_proved":False,
    "polynomial_diagnostics_are_continuum_suprema":False,
}
REMAINING=[
    "establish uniform actual input-norm, shell and zeta/common-weight comparison bounds satisfying the explicit dynamic budgets on the common tail",
    "apply retained current signed cancellation to the complete target including Drt, full pointwise pressure and all fast correlations",
    "prove generated-source supports/jets and the complete residual including viscosity, every cutoff and interaction",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        path=HERE.parent/name/"evidence.json"
        record=json.loads(path.read_text())
        bound[str(path.relative_to(PROJECT))]=sha(path)
        for rel,expected in record["source_hashes"].items():
            checks[name+":"+rel]=sha(PROJECT/rel)==expected
    return bound,checks


def sources():
    files=[HERE/n for n in ("README.md","proof_checks.py","audit.py","test_audit.py","protocol.json","requirements.txt",".gitignore")]
    files += [HERE.parent/n/"README.md" for n in DEPENDENCIES]
    files += [HERE.parent/"complete_covariance_initial_means/means.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in files}


def generate():
    result=review()
    bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),
                flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_dynamic_controls_pass":all(result["checks"].values()),
                       "six_frozen_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Dynamic budgets, inherited sources or scientific scope changed")


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--verify-record",type=Path)
    args=parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record=generate()
    validate(record)
    if args.verify_record:
        validate(json.loads(args.verify_record.read_text()))
    with args.output.open("x") as output:
        output.write(json.dumps(record,indent=2,sort_keys=True)+"\n")
    print(f"{STATUS}: {len(record['review']['checks'])} exact controls; {len(record['inherited_source_checks'])} inherited source checks")


if __name__=="__main__":
    main()
