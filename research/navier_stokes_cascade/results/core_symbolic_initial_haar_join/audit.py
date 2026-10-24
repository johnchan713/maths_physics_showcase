import argparse,hashlib,json
from pathlib import Path
from controls import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_sharp_initial_scalar_join','core_symbolic_fresh_mean_join','core_initial_haar_row_budget','core_signed_target_edge_budget','core_weighted_primary_covariance_join')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
 bound={};checks={}
 for name in DEPENDENCIES:
  p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
  for rel,h in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==h
 r=review()
 paths=[HERE/n for n in ('README.md','bounds.py','controls.py','audit.py','test_audit.py','requirements.txt','.gitignore')]
 paths += [HERE.parent/'core_sharp_initial_scalar_join/bounds.py',HERE.parent/'core_symbolic_fresh_mean_join/bounds.py',HERE.parent/'core_initial_haar_row_budget/budgets.py',HERE.parent/'core_signed_target_edge_budget/budgets.py']
 return dict(status='conditional-symbolic-complete-initial-Haar-and-projected-source-join',review=r,
  gates=dict(exact_controls=all(r['checks'].values()),five_frozen_records=len(bound)==5,historical_sources=all(checks.values())),
  dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in paths},
  flags=dict(whole_pressure_A1_all_covariance_and_sharp_means_retained=True,actual_log_pullback_and_projector_bump_inputs_still_required=True,actual_source_slow_mean_full_residual_or_blowup_certified=False))

def validate(r):
 expected=generate()
 if r!=expected or not all(expected['gates'].values()):raise ValueError('Chart transfer, pressure input, scope or frozen source changed')
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args()
 if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
 r=generate();validate(r)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print('symbolic initial Haar join:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'historical source checks')
if __name__=='__main__':main()
