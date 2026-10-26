import argparse,hashlib,json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_initial_actual_profile_gate','core_initial_haar_row_budget',
              'core_symbolic_initial_haar_join','core_log_pullback_projector_join',
              'core_symbolic_signed_target_join','core_rho_projector_signed_recompute',
              'core_initial_final_row_obstruction')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
    bound={};checks={}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,h in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==h
    r=json.loads(json.dumps(review(),default=str));paths=[HERE/n for n in ('README.md','rate_gain.py','controls.py','audit.py','test_audit.py','requirements.txt','.gitignore')]
    paths += [HERE.parent/'core_initial_actual_profile_gate/profile_gate.py',HERE.parent/'real_signed_covariance_update/signed.py',HERE.parent/'core_initial_final_row_obstruction/bounds.py']
    return dict(status='conditional-sharper-initial-target-and-strict-first-real-final-row-rate',review=r,
      gates=dict(exact_controls=all(r['checks'].values()),seven_frozen_records=len(bound)==7,historical_sources=all(checks.values())),
      dependency_evidence_sha256=bound,inherited_source_checks=checks,
      source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in paths},
      flags=dict(sharper_initial_target_rate_beta_retained=True,
        conditional_strict_first_real_final_row_margin_proved=True,
        actual_same_state_high_jets_proved=False,complete_residual_and_blowup_certified=False))
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Initial target-rate gain, provenance or scope changed')
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args()
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print('initial target-rate gain:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'historical source checks')
if __name__=='__main__':main()
