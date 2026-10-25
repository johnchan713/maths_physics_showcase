import argparse,hashlib,json
from pathlib import Path
from controls import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_initial_cancellation_fast_source','core_compact_source_quotient','core_axial_potential_pressure_gain','core_symbolic_weighted_maps','core_prefast_weighted_source_join')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
 bound={};checks={}
 for name in DEPENDENCIES:
  p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
  for rel,h in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==h
 r=review();paths=[HERE/n for n in ('README.md','bounds.py','controls.py','audit.py','test_audit.py','requirements.txt','.gitignore')]
 paths += [HERE.parent/'core_initial_cancellation_fast_source/bounds.py',HERE.parent/'core_compact_source_quotient/bounds.py',HERE.parent/'core_axial_potential_pressure_gain/bounds.py',HERE.parent/'core_symbolic_weighted_maps/bounds.py']
 return dict(status='conditional-original-prefast-quotient-correction-and-whole-source-fast-inverse',review=r,
  gates=dict(exact_controls=all(r['checks'].values()),five_frozen_records=len(bound)==5,historical_sources=all(checks.values())),
  dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in paths},
  flags=dict(original_prefast_A1_direct_quotient_order_two_derived=True,normalized_200001_jet_requirement_corrected=True,
   whole_current_source_and_conditional_fast_map_coefficients_derived=True,actual_fixed_profile_support_and_field_enclosures_required=True,
   postfast_pressure_scalar_slow_complete_residual_and_blowup_certified=False))
def validate(r):
 expected=generate()
 if r!=expected or not all(expected['gates'].values()):raise ValueError('Initial quotient correction, fast inverse, provenance or scope changed')
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args()
 if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
 r=generate();validate(r)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print('initial quotient fast inverse:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'historical source checks')
if __name__=='__main__':main()
