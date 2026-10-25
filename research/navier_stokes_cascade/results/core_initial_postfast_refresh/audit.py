import argparse,hashlib,json
from pathlib import Path
from controls import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_initial_quotient_fast_inverse','core_real_refresh_obligation_join','core_axial_potential_pressure_gain','core_full_compact_pressure_response')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
 bound={};checks={}
 for name in DEPENDENCIES:
  p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
  for rel,h in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==h
 r=review();paths=[HERE/n for n in ('README.md','bounds.py','controls.py','audit.py','test_audit.py','requirements.txt','.gitignore')]
 paths += [HERE.parent/'core_initial_quotient_fast_inverse/bounds.py',HERE.parent/'core_real_refresh_obligation_join/bounds.py',HERE.parent/'core_axial_potential_pressure_gain/bounds.py',HERE.parent/'core_full_compact_pressure_response/bounds.py']
 return dict(status='conditional-complete-initial-postfast-pressure-scalar-refresh-and-sharp-slow-input',review=r,
  gates=dict(exact_controls=all(r['checks'].values()),four_frozen_records=len(bound)==4,historical_sources=all(checks.values())),
  dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in paths},
  flags=dict(complete_postfast_pressure_memory_assembled=True,signed_scalar_change_not_double_counted=True,
   conditional_sharp_slow_input_materialized=True,current_A1_strict_request_derived=True,
   actual_current_Phi_high_jets_support_complete_residual_and_blowup_certified=False))
def validate(r):
 expected=generate()
 if r!=expected or not all(expected['gates'].values()):raise ValueError('Postfast refresh, provenance or scope changed')
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args()
 if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
 r=generate();validate(r)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print('initial postfast refresh:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'historical source checks')
if __name__=='__main__':main()
