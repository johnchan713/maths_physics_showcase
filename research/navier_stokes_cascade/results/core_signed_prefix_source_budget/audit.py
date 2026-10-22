#!/usr/bin/env python3
"""Freeze conservative actual signed-prefix source-budget evidence."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=("core_first_sourced_inverse_budget","core_initial_haar_row_budget","core_fast_inverse_budget",
              "core_slow_inverse_budget","regular_source_foundation","core_initial_scalar_seed_budget")
STATUS="first-real-signed-prefix-conservative-fresh-mean-and-complete-regular-source-norm-budget"
FLAGS={
    "full_pointwise_mean_size_not_improved_Haar_rate_used":True,
    "fresh_fast_then_recomputed_pressure_scalars_then_slow_order_retained":True,
    "actual_A1_radial_mean_complete_covariance_and_pressure_memory_retained":True,
    "conditional_actual_first_source_polynomial_derived_from_lower_inputs":True,
    "raw_source_norm_interface_to_Gaussian_adapted_inverse_joined":True,
    "full_viscosity_complete_longitudinal_transport_and_mean_phase_loss_retained":True,
    "underlying_actual_core_ODE_coefficient_support_and_fixed_operator_inputs_independently_closed":False,
    "initial_Haar_improvement_extended_to_full_or_later_complex_rows":False,
    "all_additive_cutoff_tails_or_complete_pointwise_residual_verified":False,
    "historical_evidence_or_fixed_parameters_changed":False,
    "unconditional_pulse_bounds_infinite_realization_or_blowup_proved":False,
}
REMAINING=[
    "substitute and independently discharge actual core phase normal root homogeneous pulse and weighted fixed-map coefficient and whole-path support inputs",
    "carry the first sourced correction complete complex covariance through fresh full pressure and mean updates with quantitative norms",
    "join every physical scale and the complete residual including additive flat pressure cutoff tails viscosity and all interactions before unconditional claims",
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
    paths += [HERE.parent/"core_initial_haar_row_budget/budgets.py",HERE.parent/"core_fast_inverse_budget/proof_checks.py",
              HERE.parent/"regular_source_foundation/source_checks.py",HERE.parent/"core_first_sourced_inverse_budget/budgets.py",
              HERE.parent/"complete_covariance_initial_means/means.py"]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    result=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/"protocol.json").read_text()),flags=deepcopy(FLAGS),review=result,
                gates={"all_exact_full_mean_source_and_fresh_prefix_controls_pass":bool(result["checks"]) and all(result["checks"].values()),
                       "six_historical_records_bound":len(bound)==6,
                       "all_declared_historical_sources_match":bool(checks) and all(checks.values())},
                dependency_evidence_sha256=bound,inherited_source_checks=checks,
                source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected["gates"].values()):
        raise ValueError("Signed-prefix complete source recipe fresh pressure order lower premises or historical bindings changed")


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
