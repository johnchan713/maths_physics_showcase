import argparse,hashlib,json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_initial_viscosity_sharpness_witness','real_signed_covariance_update',
              'core_localized_stress_assembly','compact_mean_foundation','physical_pulse_foundation')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
    bound={};checks={}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,h in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==h
    r=json.loads(json.dumps(review(),default=str));paths=[HERE/n for n in ('README.md','profile_gate.py','controls.py','audit.py','test_audit.py','requirements.txt','.gitignore')]
    paths += [HERE.parent/'real_signed_covariance_update/signed.py',HERE.parent/'core_initial_viscosity_sharpness_witness/witness.py']
    return dict(status='conditional-actual-profile-injectivity-and-quantitative-rate-gate',review=r,
      gates=dict(exact_controls=all(r['checks'].values()),five_frozen_records=len(bound)==5,historical_sources=all(checks.values())),
      dependency_evidence_sha256=bound,inherited_source_checks=checks,
      source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in paths},
      flags=dict(nonzero_selected_actual_profile_has_internal_source_to_viscosity_cancellation=False,
        positive_epsilon_order_gain_proved=False,actual_same_state_high_jets_proved=False,
        complete_residual_and_blowup_certified=False))
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Actual-profile gate, provenance or scope changed')
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args()
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print('initial actual-profile gate:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'historical source checks')
if __name__=='__main__':main()
