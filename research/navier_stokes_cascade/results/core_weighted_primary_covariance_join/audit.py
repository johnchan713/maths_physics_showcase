import argparse,hashlib,json
from pathlib import Path
from bounds import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_primary_target_cutoff_budget','core_primary_aggregate_atlas','core_initial_haar_row_budget','core_localized_stress_assembly','canonical_initial_mean_join')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
 bound={};checks={}
 for name in DEPENDENCIES:
  p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
  for rel,h in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==h
 r=review()
 paths=[HERE/n for n in ('README.md','bounds.py','audit.py','test_audit.py','requirements.txt','.gitignore')]
 paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
 paths += [HERE.parent/n for n in ('core_primary_target_cutoff_budget/bounds.py','core_primary_aggregate_atlas/proof_checks.py','core_initial_haar_row_budget/budgets.py','core_complete_curl_operator_budget/bounds.py')]
 return dict(status='conditional-same-weight-primary-covariance-and-mean-row-join',review=r,
  gates=dict(exact_controls=all(r['checks'].values()),five_frozen_records=len(bound)==5,historical_sources=all(checks.values())),
  dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in paths},
  flags=dict(all_six_entries_full_curl_and_pressure_retained=True,weight_preserved_before_aggregation=True,original_mean_inputs_still_required=True,same_weight_U_to_R_row_transfer_still_required=True,actual_full_mean_residual_or_blowup_certified=False))
def validate(r):
 expected=generate()
 if r!=expected or not all(expected['gates'].values()):raise ValueError('Weighted bridge, pressure input, scope or frozen source changed')
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args()
 if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
 r=generate();validate(r)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print('weighted primary bridge:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'historical source checks')
if __name__=='__main__':main()
