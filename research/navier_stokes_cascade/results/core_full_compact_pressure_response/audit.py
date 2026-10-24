import argparse,hashlib,json
from pathlib import Path
from controls import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_rho_projector_signed_recompute','core_symbolic_weighted_maps','core_fast_inverse_budget','core_fast_row_gain_join','core_native_haar_primitive_join')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
 bound={};checks={}
 for name in DEPENDENCIES:
  p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
  for rel,h in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==h
 r=review()
 paths=[HERE/n for n in ('README.md','bounds.py','controls.py','audit.py','test_audit.py','requirements.txt','.gitignore')]
 paths += [HERE.parent/'core_symbolic_weighted_maps/bounds.py',HERE.parent/'core_fast_inverse_budget/README.md',HERE.parent/'core_fast_row_gain_join/budgets.py',HERE.parent/'core_rho_projector_signed_recompute/bounds.py']
 return dict(status='conditional-full-auxiliary-compact-pressure-and-changed-covariance-response',review=r,
  gates=dict(exact_controls=all(r['checks'].values()),five_frozen_records=len(bound)==5,historical_sources=all(checks.values())),
  dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in paths},
  flags=dict(full_auxiliary_T0_A0_and_changed_covariance_response_derived=True,actual_fixed_profile_support_mass_and_jets_required=True,actual_raw_T0_A0_profiles_charts_and_plateau_required=True,
   inherited_sigma_identity_or_actual_continuum_fields_certified=False,actual_pulse_full_residual_or_blowup_certified=False))
def validate(r):
 expected=generate()
 if r!=expected or not all(expected['gates'].values()):raise ValueError('Full shifted compact map, covariance response, scope or frozen source changed')
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args()
 if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
 r=generate();validate(r)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print('full compact pressure response:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'historical source checks')
if __name__=='__main__':main()
