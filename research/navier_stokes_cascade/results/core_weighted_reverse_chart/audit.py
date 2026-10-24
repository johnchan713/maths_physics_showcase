import argparse,hashlib,json
from pathlib import Path
from controls import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_weighted_primary_covariance_join','core_complete_curl_operator_budget','core_primary_target_cutoff_budget','core_initial_haar_row_budget','core_primary_aggregate_atlas')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
 bound={};checks={}
 for name in DEPENDENCIES:
  p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
  for rel,h in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==h
 r=review()
 paths=[HERE/n for n in ('README.md','bounds.py','controls.py','audit.py','test_audit.py','requirements.txt','.gitignore')]
 paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
 paths += [HERE.parent/n/'bounds.py' for n in DEPENDENCIES if n not in ('core_initial_haar_row_budget','core_primary_aggregate_atlas')]
 paths += [HERE.parent/'core_initial_haar_row_budget/budgets.py']
 return dict(status='conditional-same-weight-reverse-radial-chart-and-initial-row-interface',review=r,
  gates=dict(exact_controls=all(r['checks'].values()),five_frozen_records=len(bound)==5,historical_sources=all(checks.values())),
  dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in paths},
  flags=dict(same_weight_U_to_R_tensor_transfer_derived=True,raw_and_factorial_conventions_both_enclosed=True,
   no_new_derivative_order_edge_pole_or_epsilon_loss=True,all_six_covariances_and_whole_pressure_retained=True,
   other_same_hull_weighted_mean_map_inputs_still_required=True,actual_full_mean_residual_or_blowup_certified=False))
def validate(r):
 expected=generate()
 if r!=expected or not all(expected['gates'].values()):raise ValueError('Chart transfer, pressure input, scope or frozen source changed')
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args()
 if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
 r=generate();validate(r)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print('reverse weighted chart:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'historical source checks')
if __name__=='__main__':main()
