#!/usr/bin/env python3
"""Freeze the scoped axial potential and same-postwave pressure-complement gain."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_angular_pressure_normal_form','core_compact_source_quotient',
    'core_next_sourced_wave_join','core_next_sourced_mean_tail_join',
    'core_preserved_haar_scalar_join','compact_mean_foundation')
STATUS='axial-potential-pressure-complement-and-same-postwave-A1-order-two'
FLAGS={
    'exact_axial_potential_normal_form_uses_only_old_coarse_RHS_inputs':True,
    'same_postwave_pressure_complement_proved_at_6q_from_angular_row_E':True,
    'projected_A1_error_order_two_needs_no_additional_fine_Haar_callback_here':True,
    'upstream_fine_row_Haar_dependencies_remain_in_older_graph':True,
    'full_potential_one_and_beta_two_require_preserved_Haar_gamma_history':True,
    'full_pressure_viscous_memory_factor_two_and_moving_cutoffs_retained':True,
    'variable_G_R_and_whole_potential_hull_retained':True,
    'homogeneous_Haar_comparison_exposes_independent_history_gate':True,
    'actual_all_scale_high_jet_Haar_profile_or_support_inputs_validated':False,
    'final_mean_fields_substituted_for_postwave_fast_source':False,
    'full_gamma_total_pressure_or_whole_regular_E_order_promoted':False,
    'complete_residual_absolute_decay_infinite_field_or_blowup_proved':False,
    'historical_frozen_evidence_fields_width_band_roots_or_carriers_changed':False
}
REMAINING=[
    'validate actual same-postwave Ez/Etheta/coarse pressure/Psi inputs and Q Phi norms at all requested jets, preserving every upstream original Haar/root/source dependency',
    'validate whole compact potential and stress-flux support traces, moving T1/A1 profiles, fixed G_R jets and original 6q/12q weighted/raw conversion constants',
    'independently validate the preserved finer Haar gamma history for the full Psi1/beta2 conclusion; the Q-only tail route does not supply it',
    'then advance the next generic signed cycle with complete pressure/mean memory; sufficient full regular orders, infinite diagonal, localization/force and growth remain before blow-up claims'
]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,expected in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==expected
    return bound,checks
def sources():
    paths=[HERE/n for n in ('README.md','bounds.py','controls.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
    paths += [HERE.parent/'core_angular_pressure_normal_form/bounds.py',HERE.parent/'core_compact_source_quotient/bounds.py',
              HERE.parent/'core_next_sourced_wave_join/budgets.py',HERE.parent/'core_next_sourced_mean_tail_join/budgets.py',
              HERE.parent/'compact_mean_corrections/mean_system.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}
def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {k:canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v
def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_axial_pressure_cutoff_and_budget_controls_pass=bool(r['checks']) and all(r['checks'].values()),
            six_historical_records_bound=len(bound)==6,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Axial potential, same-stage pressure complement, scope or provenance changed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")
if __name__=='__main__':main()
