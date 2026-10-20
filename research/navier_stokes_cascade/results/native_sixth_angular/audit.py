#!/usr/bin/env python3
"""Actual native sixth angular cover and fifth scalar R5; cutoff still open."""
import argparse,hashlib,importlib.util,json,sys
from copy import deepcopy
from functools import lru_cache
from pathlib import Path
from native6 import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
PARENT='1ebd014108a145266238882a7c3b447946088270'
STATUS='actual-native-sixth-angular-moments-and-fifth-scalar-R5-bound'
MANIFEST_SHA='d1e03294d1fe73722773345d812a3e0ab1052ebdf5fcf6450ee8bf6036106f91'
DEPENDENCIES=('fifth_loop_input_gate','reference_fifth_transfer','actual_fifth_native',
 'native_fifth_mixed_cover','shape_join_pressure','annular_exterior_foundation',
 'prescribed_continuation_foundation','fifth_moment_volterra','final_fifth_moment_interface')
FLAGS=dict(seventh_native_packet_reference_pressure_and_sixth_sources_derived=True,
 actual_PREMOD_sixth_angular_fields_moments_and_same_annular_C6_root_derived=True,
 native_first_radial_fifth_shear_pressure_and_R5_A32_under_retained_foundation_derived=True,
 native_R5_estimate_independent_of_variance_cutoff_and_joint_target_J=True,
 full_axis_Cp_energy_both_squares_and_nonzero_axial_pressure_memory_retained=True,
 all_parameters_widths_cutoffs_supports_roots_gauges_and_history_unchanged=True,
 actual_Z4_Z5_or_joint_target_J_or_reduced_fifth_cutoff_gate_verified=False,
 complete_native_joint_J6_or_FINAL_modulated_sixth_field_bound_derived=False,
 unconditional_global_leading_C4_accepted=False,
 uniform_actual_pulses_complete_covariance_axial_pressure_means_residual_or_blowup_verified=False,
 independent_foundation_review_or_formal_certificate=False)
REMAINING=['enclose_unchanged_actual_Z4_Z5_or_direct_joint_target_J_and_fifth_compositions;native_R5_is_no_longer_missing_under_retained_foundation',
 'join_complete_compact_pressure_flux_maps_to_global_regions_and_retained_physical_scale_inequalities',
 'prove_coefficient_reserved_mean_supports_uniform_actual_pulses_complete_covariance_axial_pressure_means_and_full_viscous_cutoff_interaction_residual']
def protocol():
 return dict(parent_commit=PARENT,manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
  premises='retained exact core coefficient certificate and native_regional_existence_positive_floors_prior_roots;native estimates do not assume variance cutoff J;no manuscript final theorem premise',
  native='D7=11(1+S)(1+rho^-1)^7;Y derivative C6 costs18;N7=1e30(1+S)^10(1+rho^-1)^64<Kold;axisC7<=1650P2;starlogC7<=10LambdaKold9',
  reference='same eta_independent monotone positive taper average through7;T7=1e10 exp(-Lambda)Lambda7Kold64;fullPiC7<2000P2<Kold;complete sourceC6 Sq<=1e4LambdaKold12 Sn<=1e3Kold2;reference p1C6<=1e8X LambdaKold26 nsC6<=1e5Kold3;full incoming data',
  actual='actual slopes not flattened reference physical shear;logCE_C6<=J6=1e11LambdaKold27<Bold/100;g_C6<3e-18;shapeE_C6<=same Eb exp204Bold/C;regular core Cp integral before using Xa floor',
  root='same EARLIER eta_dependent Q,C6 normalizers5_25,target<2.86e-17,quadratic<=100,rootC6<=5.72e-14 inside same C1 uniqueness ball;retain Q_eta through4 and both energy squares;no unit mass change',
  cover='SUM all native pure angular regions E_C6<A2,U_C6<A;first_radial_C5 E_y U_y<A5 using actual slopes and charged fixed width beta_y;no sigma7 and no complete J6 claim;full original moments_and_axis C6<A8 at every partial radius',
  pressure='inverseE_C5<A17,a_bs_C5<A23;complete W<A10,Ns<A12,Hx_inverse<A18,Qs<A29,p1_p2_C5<A31;derivative moment6 cost6,full axis_Cp and nonzero -WU retained;R5<A32,p1 not loop argument',
  gate='under separate unverified J,Lambda5<exp(1100000G+10A20)*Z5*A32;remaining Z5gate logZ5+32logA+10A20+1100000G<=G2;actual Z4_Z5 J and gate unverified;not K_N reselection or global acceptance',
  scope='conditional native analytical enclosure only;PREMOD not FINAL,specified sixth angular plus one radial fifth not complete J6;exact/outward checks diagnostics not continuum proof;no global pulse covariance mean residual blowup')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def module(name,path):
 sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sys.modules[name]=m;sp.loader.exec_module(m);return m
@lru_cache(None)
def history():
 p=HERE.parent/'fifth_loop_input_gate';prior=sys.path[:];sys.path.insert(0,str(p))
 try:parent=module('native6_frozen_parent',p/'audit.py')
 finally:sys.path[:]=prior
 names=sorted([Path(p).parent.name for p in parent.history()['frozen_records']]+['fifth_loop_input_gate'])
 digest=hashlib.sha256(json.dumps(names,separators=(',',':')).encode()).hexdigest()
 if len(names)!=136 or len(set(names))!=136 or digest!=MANIFEST_SHA:raise ValueError('Fixed 136-record manifest changed')
 p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
 try:m=module('native6_historical_bindings',p/'audit.py')
 finally:sys.path[:]=prior
 bindings=m.verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])['bindings']
 if not bindings or not all(bindings.values()):raise ValueError('Historical source binding failed')
 return dict(all_bindings_verified=True,binding_count=len(bindings),fixed_manifest_sha256=digest,
  binding_result_sha256=hashlib.sha256(json.dumps(bindings,sort_keys=True).encode()).hexdigest(),
  frozen_records={str((HERE.parent/n/'evidence.json').relative_to(PROJECT)):sha(HERE.parent/n/'evidence.json') for n in names})
def sources():
 files=[HERE/n for n in ('README.md','native6.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
 files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
 files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','fifth_loop_input_gate/audit.py','fifth_loop_input_gate/moving5.py',
  'reference_fifth_transfer/reference_bounds.py','actual_fifth_native/native_bounds.py',
  'prescribed_continuation_foundation/continuation_checks.py','core_all_order_foundation/coordinate_checks.py')]
 return {str(p.relative_to(PROJECT)):sha(p) for p in files}
def generate():
 if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
 r=deepcopy(review());h=history();g=dict(r['checks']);g['historical_source_bytes']=h['all_bindings_verified']
 return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
 expected=generate()
 if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Native sixth angular record changed')
 if not all(expected['gates'].values()):raise ValueError('Native sixth angular check failed')
def main(argv=None):
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
 if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
 r=generate();validate(r)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
 print(f'{STATUS}: {len(r["gates"])} checks; native R5 enclosed under foundation, actual cutoff/J and FINAL acceptance remain open.')
if __name__=='__main__':main()
