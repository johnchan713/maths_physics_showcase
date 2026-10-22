#!/usr/bin/env python3
"""Freeze original Riccati-length normal-floor closure and quantitative limits."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('frame_column_foundation','core_primary_tail_recheck','core_physical_representative_closure','core_tangential_normal_closure','leading_frame_foundation')
STATUS='original-rule-defined-primary-tail-normal-nonvanishing-from-retained-Riccati-length'
FLAGS={
    'full_tangential_error_bounded_by_original_Dt_over_actual_L':True,
    'physical_reference_growth_formula_implies_primitive_consistency':True,
    'original_Riccati_length_dominates_original_normal_threshold_before_ODE_use':True,
    'independent_Ctan_Cphi_compatibility_gate_removed_for_original_normal_nonvanishing':True,
    'original_reference_Bs_inverse_bound_kept_distinct_from_phase_beta_bound':True,
    'full_axial_transport_pressure_and_old_frame_budget_premises_retained':True,
    'new_analytical_assumptions_added':False,
    'actual_primitive_enclosures_or_arbitrary_historical_band_certified':False,
    'stronger_phase_beta_inverse_bound_inferred_from_coarse_Bs_floor':False,
    'full_propagator_covariance_mean_or_physical_residual_closed':False,
    'historical_fields_bands_lengths_supports_or_evidence_changed':False,
    'unconditional_pulses_or_blowup_proved':False,
}
REMAINING=[
    'retain actual validity of original leading background and frame primitive bounds including Mz and Db, physical coverage q-domain and original selector inputs',
    'use original inverse-normal constant 4/bmin^2 in sourced inverse budgets; any stronger explicitly requested phase-beta constant needs its separate justification',
    'keep unknown historical numeric-band membership separate from the rule-defined accepted tail and reject unlinked or newly enlarged length-budget inputs',
    'close full sourced propagator supports generated complete covariance fresh pressure/axial means and all viscous cutoff interaction residual estimates',
]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,expected in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==expected
    return bound,checks


def sources():
    paths=[HERE/n for n in ('README.md','budgets.py','controls.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
    paths += [HERE.parent/'frame_column_foundation/frame_budget.py',HERE.parent/'core_primary_tail_recheck/proof_checks.py',HERE.parent/'leading_frame_foundation/frame_floors.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
                gates=dict(all_exact_original_length_normal_and_scope_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                           five_historical_records_bound=len(bound)==5,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
                dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Original length normal-floor closure or original-denominator scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} inherited source checks")


if __name__=='__main__':main()
