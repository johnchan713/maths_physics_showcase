#!/usr/bin/env python3
"""Freeze an explicit admissible variance cutoff and discharge Z4/Z5."""
import argparse,hashlib,importlib.util,json,sys
from copy import deepcopy
from functools import lru_cache
from pathlib import Path
from cutoff import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
PARENT='e1fbc36677e27239c2da686d75080c3e1d74c068'
STATUS='explicit-source-admissible-variance-cutoff-and-Z4-Z5-gates'
MANIFEST_SHA='364e9fef3314e8eae42b27a593b2e6bc9156188cf547d8b36896672c94a7ed96'
DEPENDENCIES=('native_sixth_angular','fifth_loop_input_gate','fourth_loop_shape_gate',
 'actual_cutoff_pullback','cutoff_composition_obstruction','joint_fourth_loop_gamma')
FLAGS=dict(
 inherited_claim_that_a_quantitative_actual_cutoff_was_already_fixed_corrected=True,
 explicit_source_admissible_Cinfinity_variance_cutoff_now_frozen=True,
 same_deltaL_window_plateaus_widths_roots_gauges_supports_and_pressure_memory_retained=True,
 actual_raw_Z4_Z5_enclosed_for_the_now_selected_cutoff=True,
 fourth_target_J_and_reduced_fifth_cutoff_gate_discharged_under_retained_foundations=True,
 same_frequency_compact_modulation_and_exact_five_moment_repair_gate_discharged=True,
 complete_FINAL_joint_leading_and_global_physical_scale_joins_proved=False,
 uniform_actual_pulses_complete_covariance_axial_pressure_means_residual_or_blowup_verified=False,
 independent_foundation_review_or_formal_certificate=False)
REMAINING=['join_complete_compact_pressure_flux_and_repaired_state_bounds_to_inner_and_outer_regions at the retained physical scales',
 'verify global pressure moment identities coefficient and reserved-mean supports after the now-fixed cutoff choice',
 'prove uniform actual-pulse size derivative support bounds complete-wave covariance axial-pressure means and the full viscosity-cutoff-interaction residual']
def protocol():
 return dict(parent_commit=PARENT,manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
  correction='C.9 says let zetaL be a smooth nonnegative function with range and plateaus; it does not select a quantitative shape. Earlier actual-retained wording is corrected without editing frozen evidence.',
  choice='phi=exp(-1/t-1/(1-t));K=1-integral(phi)/total;zetaL(v)=K(8(v-2)/deltaL-1);C-infinity flat, decreasing, same exact window',
  derivative='phi^(n)=P_n(u,v)phi,u=1/t,v=1/(1-t);u^a exp(-u)<=a!;integral(phi)>exp(-9/2)/3 and e<3 give inverse<729;B0..4=1,4,80,3600,292320',
  gates='Z4<3000000,Z5<1e9;fourth log gate <9A20<A32;fifth gate divided by G is <1100001<G using G>A32;combine parent R5<A32',
  consequence='under the retained native/regional/root premises, actual target J and B5 now hold for this fixed admissible source choice; same N and exact five-row repair/frequency gate follow',
  scope='analytical source-completion checkpoint;historical compressed-cutoff obstruction remains valid for the class;not FINAL joint/global pulse covariance mean residual or blowup proof')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def module(name,path):
 sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sys.modules[name]=m;sp.loader.exec_module(m);return m
@lru_cache(None)
def history():
 p=HERE.parent/'native_sixth_angular';prior=sys.path[:];sys.path.insert(0,str(p))
 try:parent=module('explicit_cutoff_frozen_parent',p/'audit.py')
 finally:sys.path[:]=prior
 names=sorted([Path(x).parent.name for x in parent.history()['frozen_records']]+['native_sixth_angular'])
 digest=hashlib.sha256(json.dumps(names,separators=(',',':')).encode()).hexdigest()
 if len(names)!=137 or len(set(names))!=137 or digest!=MANIFEST_SHA:raise ValueError('Fixed 137-record manifest changed')
 p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
 try:m=module('explicit_cutoff_historical_bindings',p/'audit.py')
 finally:sys.path[:]=prior
 bindings=m.verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])['bindings']
 if not bindings or not all(bindings.values()):raise ValueError('Historical source binding failed')
 return dict(all_bindings_verified=True,binding_count=len(bindings),fixed_manifest_sha256=digest,
  binding_result_sha256=hashlib.sha256(json.dumps(bindings,sort_keys=True).encode()).hexdigest(),
  frozen_records={str((HERE.parent/n/'evidence.json').relative_to(PROJECT)):sha(HERE.parent/n/'evidence.json') for n in names})
def sources():
 files=[HERE/n for n in ('README.md','cutoff.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
 files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
 files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','native_sixth_angular/audit.py','fourth_loop_shape_gate/loop_fourth.py','fifth_loop_input_gate/moving5.py')]
 return {str(p.relative_to(PROJECT)):sha(p) for p in files}
def generate():
 if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
 r=deepcopy(review());h=history();g=dict(r['selection']['checks']);g['historical_source_bytes']=h['all_bindings_verified']
 g['all_selection_consequences_typed']=all(type(v) is bool for v in r['consequences'].values())
 return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
  inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
 expected=generate()
 if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Explicit cutoff record changed')
 if not all(expected['gates'].values()):raise ValueError('Explicit cutoff check failed')
def main(argv=None):
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
 if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
 r=generate();validate(r)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
 print(f'{STATUS}: {len(r["gates"])} checks; explicit admissible cutoff fixed and Z4/Z5 gates discharged; FINAL/global closure remains open.')
if __name__=='__main__':main()
