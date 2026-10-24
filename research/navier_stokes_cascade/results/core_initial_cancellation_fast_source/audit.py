import argparse,hashlib,json
from pathlib import Path
from controls import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_full_compact_pressure_response','core_fast_row_gain_join','core_fast_inverse_budget','core_sharp_initial_scalar_join','core_symbolic_fresh_mean_join')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
 bound={};checks={}
 for name in DEPENDENCIES:
  p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
  for rel,h in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==h
 r=review()
 paths=[HERE/n for n in ('README.md','bounds.py','controls.py','audit.py','test_audit.py','requirements.txt','.gitignore')]
 paths += [HERE.parent/'core_full_compact_pressure_response/bounds.py',HERE.parent/'core_fast_row_gain_join/budgets.py',HERE.parent/'core_sharp_initial_scalar_join/bounds.py',HERE.parent/'core_symbolic_fresh_mean_join/bounds.py']
 return dict(status='conditional-initial-cancellation-current-regular-fast-source-with-deferred-A1',review=r,
  gates=dict(exact_controls=all(r['checks'].values()),five_frozen_records=len(bound)==5,historical_sources=all(checks.values())),
  dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in paths},
  flags=dict(initial_cancellation_and_current_regular_source_program_derived=True,actual_fixed_profile_support_mass_and_jets_required=True,actual_raw_T0_A0_profiles_charts_and_plateau_required=True,
   inherited_sigma_identity_or_actual_continuum_fields_certified=False,actual_fine_axial_source_fast_slow_pulse_or_blowup_certified=False))
def validate(r):
 expected=generate()
 if r!=expected or not all(expected['gates'].values()):raise ValueError('Initial cancellation, original A1 obligation, source domain or frozen source changed')
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args()
 if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
 r=generate();validate(r)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print('initial cancellation fast source:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'historical source checks')
if __name__=='__main__':main()
