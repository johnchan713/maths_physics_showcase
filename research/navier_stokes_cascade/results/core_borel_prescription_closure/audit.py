#!/usr/bin/env python3
"""Freeze the original rule's infinite diagonal implication and its limits."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('background_jet_foundation','background_summation','core_tangential_normal_closure','core_borel_cutoff_certificate')
STATUS='original-real-borel-prescription-infinite-diagonal-closure-relative-to-full-coefficient-regularity'
FLAGS={
    'all_sequences_admitted_by_original_real_rule_satisfy_all_order_diagonal':True,
    'strict_normalized_gain_margin_and_finite_stage_induction_supplied':True,
    'inherited_finite_physical_log_power_requirements_retained':True,
    'independent_diagonal_and_diagnostic_cutoff_identity_removed_for_rule_defined_floor':True,
    'full_pressure_radial_curl_and_Euler_conversion_costs_retained':True,
    'new_analytical_assumptions_added':False,
    'actual_separately_fixed_numeric_sequence_or_mpmath_rounding_certified':False,
    'actual_leading_coefficient_certificate_chart_gaps_or_band_newly_verified':False,
    'historical_scales_cutoffs_supports_or_evidence_changed':False,
    'complete_residual_unconditional_pulses_or_blowup_proved':False,
}
REMAINING=[
    'verify retained leading and full coefficient regularity premises on actual required domains',
    'verify actual physical representative and native neighborhood coverage and positive leading growth/shear gaps',
    'substitute actual leading first tuple bounds into the same Ctan slow and normal inequalities at the unchanged band',
    'if an independently fixed numeric sequence is intended instead of the rule-defined family identify its same-rule adherence or prove direct diagonal inequalities without reselection',
    'retain remaining radial phase propagator supports full covariance fresh pressure/axial means and complete physical residual obligations',
]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,expected in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==expected
    return bound,checks


def sources():
    paths=[HERE/n for n in ('README.md','proof.py','controls.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
    paths += [HERE.parent/'background_jet_foundation/foundation.py',HERE.parent/'background_summation/summation.py',
              HERE.parent/'core_tangential_normal_closure/budgets.py',HERE.parent/'core_background_finite_head_budget/budgets.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
                gates=dict(all_exact_universal_source_and_transfer_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                           four_frozen_records_bound=len(bound)==4,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
                dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Original prescription universal implication or retained scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} inherited source checks")


if __name__=='__main__':main()
