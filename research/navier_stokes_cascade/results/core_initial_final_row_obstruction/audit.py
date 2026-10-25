import argparse,hashlib,json
from pathlib import Path
from controls import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_initial_postfast_refresh','core_post_signed_mean_closure','core_axial_potential_pressure_gain','core_final_mean_row_closure','post_signed_mean_relaxation')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
 bound={};checks={}
 for name in DEPENDENCIES:
  p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
  for rel,h in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==h
 r=review();paths=[HERE/n for n in ('README.md','bounds.py','controls.py','audit.py','test_audit.py','requirements.txt','.gitignore')]
 paths += [HERE.parent/'core_initial_postfast_refresh/bounds.py',HERE.parent/'core_post_signed_mean_closure/budgets.py',HERE.parent/'core_axial_potential_pressure_gain/bounds.py',HERE.parent/'core_final_mean_row_closure/budgets.py',HERE.parent/'post_signed_mean_relaxation/relax.py']
 return dict(status='conditional-initial-final-row-and-sharp-radial-viscosity-obstruction',review=r,
  gates=dict(exact_controls=all(r['checks'].values()),five_frozen_records=len(bound)==5,historical_sources=all(checks.values())),
  dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in paths},
  flags=dict(complete_conditional_final_row_and_pressure_assembled=True,same_state_A1_raw_recipe_materialized=True,
   unique_zero_margin_radial_viscosity_group_identified=True,actual_high_jets_support_positive_gain_complete_residual_and_blowup_certified=False))
def validate(r):
 expected=generate()
 if r!=expected or not all(expected['gates'].values()):raise ValueError('Initial final-row obstruction, provenance or scope changed')
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args()
 if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
 r=generate();validate(r)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print('initial final row obstruction:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'historical source checks')
if __name__=='__main__':main()
