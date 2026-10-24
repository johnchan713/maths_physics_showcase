import argparse,hashlib,json
from pathlib import Path
from controls import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_native_haar_primitive_join','core_post_signed_mean_closure','core_actual_scalar_gain_join','core_fast_row_gain_join','core_symbolic_real_covariance_join')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
 bound={};checks={}
 for name in DEPENDENCIES:
  p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
  for rel,h in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==h
 r=review()
 paths=[HERE/n for n in ('README.md','bounds.py','controls.py','audit.py','test_audit.py','requirements.txt','.gitignore')]
 paths += [HERE.parent/'core_post_signed_mean_closure/budgets.py',HERE.parent/'core_actual_scalar_gain_join/budgets.py',HERE.parent/'core_fast_row_gain_join/budgets.py']
 return dict(status='conditional-first-real-refresh-obligations-and-current-scalars',review=r,
  gates=dict(exact_controls=all(r['checks'].values()),five_frozen_records=len(bound)==5,historical_sources=all(checks.values())),
  dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in paths},
  flags=dict(first_real_requests_and_current_scalar_protocol_derived=True,actual_fixed_profile_support_mass_and_jets_required=True,
   actual_sigma_family_still_required=True,actual_pulse_full_residual_or_blowup_certified=False))
def validate(r):
 expected=generate()
 if r!=expected or not all(expected['gates'].values()):raise ValueError('Compensator composition, moment hull, pressure memory or frozen source changed')
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args()
 if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
 r=generate();validate(r)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print('real refresh obligations:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'historical source checks')
if __name__=='__main__':main()
