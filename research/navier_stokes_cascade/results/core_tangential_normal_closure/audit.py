#!/usr/bin/env python3
"""Freeze low-order diagonal closure and the tangential normal margin."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_native_chart_cover','core_background_finite_head_budget','core_borel_cutoff_certificate','core_phase_coefficient_budget')
STATUS='low-order-diagonal-input-closure-and-complete-tangential-normal-margin'
FLAGS={
    'same_complete_diagonal_implies_qh_correction_through_first_order':True,
    'independent_first_coefficient_K2_and_leading_second_inputs_removed_for_floor':True,
    'full_axial_transport_and_both_tangential_tilt_errors_retained':True,
    'original_slow_gate_implies_Ctan_over_S_without_extra_factor':True,
    'normal_floor_recipe_reduced_to_4_Ctan_squared_and_length_factor_16':True,
    'new_analytical_assumptions_added':False,
    'actual_infinite_diagonal_cutoff_identity_or_first_band_certified':False,
    'higher_order_or_radial_phase_obligations_removed':False,
    'pressure_source_supports_or_complete_covariance_newly_closed':False,
    'historical_evidence_parameters_schedule_or_accepted_bands_changed':False,
    'complete_residual_unconditional_pulses_or_blowup_proved':False,
}
REMAINING=[
    'justify the SAME complete infinite Euler diagonal and actual cutoff by the original prescription, not finite samples',
    'verify actual physical representatives native value neighborhoods leading first tuple bounds and positive growth/shear gaps',
    'substitute those inputs into 4*Ctan^2<=beta*ell_existing^4 and the original slow gate, or unchanged length compatibility with factor 16',
    'retain independent higher phase radial propagator support full pressure covariance mean and complete residual obligations',
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
    paths += [HERE.parent/n/'budgets.py' for n in DEPENDENCIES]
    paths += [HERE.parent/'background_jet_foundation/README.md',HERE.parent/'background_summation/README.md']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
                gates=dict(all_exact_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                           four_historical_records_bound=len(bound)==4,
                           all_declared_historical_sources_match=bool(checks) and all(checks.values())),
                dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(record):
    expected=generate()
    if record!=expected or not all(expected['gates'].values()):raise ValueError('Low-order diagonal normal-margin source or scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} inherited source checks")


if __name__=='__main__':main()
