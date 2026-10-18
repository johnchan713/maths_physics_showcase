#!/usr/bin/env python3
"""Explicit compact leading-stress floor, relative to the retained same-profile lemmas."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from floor_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='explicit-three-region-compact-leading-stress-floor'
DEPENDENCIES=('core_growth_activation','compact_pressure_bridge','regional_cone_join',
 'leading_frame_foundation','heat_exterior_audit','global_pressure_moment_foundation','first_seminorm_bound')
FLAGS=dict(positive_compact_leading_stress_floor_derived_from_same_current_regional_maps=True,
 actual_tc_over_sixteen_offset_and_exact_inner_flat_exponent_retained=True,
 compact_modulation_and_exact_repair_D_margin_converted_to_amplitude=True,
 compensated_exterior_finite_radius_and_edited_shear_errors_retained=True,
 exterior_stress_identification_requires_same_exact_functional_five_moments_and_gauge=True,
 free_tau_input_in_first_seminorm_bound_replaced_by_fstar_over_64G10=True,
 pressure_memory_Pi1_Pi2_and_all_previous_first_source_costs_unchanged=True,
 earlier_core_cone_heat_sign_and_seed_results_consumed_not_recounted=True,
 log_floor_kept_without_exponentiating_actual_giant_parameters=True,
 independent_unconditional_exact_core_and_regional_lemmas_proved_here=False,
 positive_unweighted_floor_at_Xa_Xb_or_on_full_open_annulus=False,
 actual_collar_operator_and_jet_envelopes_quantitatively_evaluated=False,
 same_existing_tilt_cone_radius_or_complete_value_gate_discharged=False,
 amplitude_floor_implies_all_derivative_or_full_wave_covariance_bounds=False,
 historical_first_band_old_field_or_old_force_recertified=False,
 unconditional_pulse_full_residual_closure_or_blowup_proved=False,formal_certificate=False)
REMAINING=['quantify_actual_activated_collar_leading_jet_and_global_pressure_source_envelopes_and_bump_cutoff_jet_costs',
 'supply_quantitative_same_existing_tilt_cone_radius_and_check_2h_times_ell_minus_two_at_least_log2_Cstar_over_radius_with_all_other_value_gates',
 'retain_all_scale_pulse_covariance_mean_viscosity_cutoff_interaction_and_full_residual_obligations_without_historical_first_band_promotion']

def protocol():
    return dict(parent_commit='5136b026a5626b21a283fba4f74a9ee424d66060',
       domain='Xinner=Xa exp(tc/16)<=X<=Xc, terminal heat coordinate y=1/2; eta in [-1,1]; viscosity normalized to one',
       cover=['[Xinner,Xminus] protected inner','[Xminus,Xrep] compact modulation plus exact repair','[Xrep,Xc] same-pressure compensated exterior'],
       leading_swirl='f_star=exp(-2Lambda M)/(4C)*exp[-(Ab/2)log(Xb/Xa)]; Ab=4Aenv; Aenv=C^4096',
       inner='s>=1/(160sqrt(log C)); |T0|>=(9/16) f_star C^-25602; G>C^4096 makes inner/common>36 C^15358',
       compact='T0=F(p1-a,p2+b); t=-b/a; D=p1+t p2-v; J=p2-t p1; norm^2=F^2(D^2+J^2)/(1+t^2); D>=kappa/2, R=4G^3, kappa=G^-4',
       exterior='Xrep=Xstar; p1=XQ/L>=hXstar/4000>10^12; a<5; angular stress alone exceeds F; exact global moments identify this same stress',
       floor='tau_star=f_star/(64G^10); log tau_star=log f_star-log64-10C^131072',
       first_relative='A_W<=64G^10/f_star*(Xc Htheta/2+sqrt(Xc) Hz); all previous pressure and jet costs unchanged',
       gate='2h(ell-2)>=log2(Cstar/r_Kc) remains sufficient for the whole 4Q chart; same-tilt radius and other value gates remain inputs',
       scope='conditional leading amplitude floor, not a cone/derivative/full-background/wave/residual bound or actual old-first-band acceptance')

@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('regional_floor_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','floor_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/'harmonic_closure_foundation/audit.py',HERE.parent/'leading_frame_foundation/frame_floors.py',
        HERE.parent/'compact_pressure_bridge/bridge.py',HERE.parent/'heat_exterior_audit/bounds.py',
        HERE.parent/'first_seminorm_bound/seminorms.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
@lru_cache(None)
def reviews():return review()
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(reviews());h=history();g=dict(r['checks'])
    g['historical_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Regional stress-floor record changed')
    if not all(expected['gates'].values()):raise ValueError('Regional stress-floor check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; same-tilt radius, collar envelopes and unconditional theorem remain unverified.')
if __name__=='__main__':main()
