import argparse,hashlib,json
from pathlib import Path
from controls import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_symbolic_signed_target_join','core_weighted_primary_covariance_join','core_generic_signed_covariance_join','real_signed_covariance_update','core_weighted_reverse_chart')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
 bound={};checks={}
 for name in DEPENDENCIES:
  p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
  for rel,h in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==h
 r=review()
 paths=[HERE/n for n in ('README.md','bounds.py','controls.py','audit.py','test_audit.py','requirements.txt','.gitignore')]
 paths += [HERE.parent/'core_symbolic_signed_target_join/bounds.py',HERE.parent/'core_generic_signed_covariance_join/budgets.py',HERE.parent/'core_weighted_reverse_chart/bounds.py']
 return dict(status='conditional-symbolic-first-real-signed-covariance-and-full-pressure-join',review=r,
  gates=dict(exact_controls=all(r['checks'].values()),five_frozen_records=len(bound)==5,historical_sources=all(checks.values())),
  dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in paths},
  flags=dict(all_six_covariances_signed_squares_and_whole_pressure_retained=True,new_operator_domain_and_actual_projector_inputs_still_required=True,actual_source_slow_mean_full_residual_or_blowup_certified=False))

def validate(r):
 expected=generate()
 if r!=expected or not all(expected['gates'].values()):raise ValueError('Chart transfer, pressure input, scope or frozen source changed')
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args()
 if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
 r=generate();validate(r)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print('symbolic real covariance join:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'historical source checks')
if __name__=='__main__':main()
