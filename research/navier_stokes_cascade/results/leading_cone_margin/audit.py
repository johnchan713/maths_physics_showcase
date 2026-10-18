#!/usr/bin/env python3
"""Quantitative leading direction margins and a gate for the SAME existing tilt."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from cone_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='global-leading-direction-margins-and-existing-tilt-gate'
DEPENDENCIES=('core_growth_activation','compact_pressure_bridge','regional_cone_join','heat_exterior_audit',
 'outer_edge_stability','background_gap_certificate','regional_stress_floor','leading_frame_foundation')
FLAGS=dict(global_leading_orientation_and_raw_cone_deficit_derived_from_same_regional_maps=True,
 normalized_fourth_gap_denominator_and_finite_p_exterior_terms_retained=True,
 limiting_unit_directions_used_at_both_flat_stress_edges=True,
 same_gauge_functional_moment_and_complete_pressure_inputs_preserved=True,
 exact_existing_tilt_compatibility_and_padded_test_derived=True,
 failed_conservative_tilt_test_not_identified_as_actual_cone_failure=True,
 original_gap_map_receives_full_band_factors_and_squared_discriminant_factor=True,
 existing_tilt_frame_carriers_roots_frequency_width_and_supports_unchanged=True,
 earlier_core_heat_sign_amplitude_and_seed_results_consumed_not_recounted=True,
 actual_existing_tilt_satisfies_new_sufficient_inequality_verified=False,
 actual_activated_collar_jet_source_envelopes_or_combined_scale_gate_evaluated=False,
 independent_unconditional_core_and_all_regional_lemmas_proved_here=False,
 leading_direction_bounds_are_complete_wave_covariance_or_residual_estimates=False,
 historical_first_band_old_field_or_old_force_recertified=False,
 unconditional_pulse_full_residual_closure_or_blowup_proved=False,formal_certificate=False)
REMAINING=['quantify_actual_activated_collar_leading_jet_global_pressure_source_and_cutoff_bump_envelopes',
 'verify_existing_u_squared_greater_than_one_minus_epsilon_star_over_epsilon_star_or_derive_sharper_same_u_raw_deficit_without_rechoosing_tilt',
 'combine_actual_Cstar_and_original_gap_radius_with_all_value_freezing_geometry_and_inverse_gates_before_actual_scale_or_full_residual_claims']

def protocol():
    return dict(parent_commit='1f031a3c84cfafc8461bdb531d05adde904581b9',
       domain='closed [Xa,Xb] x [-1,1], using the smooth limiting unit leading stress at Xa and Xb',
       hierarchy='C>=2; Aenv=C^4096; G=exp(Aenv^32)>C^4096; actual giant parameters not evaluated',
       compact='H0>=kappa/(10R^3)=1/(640G^13); 1-rho^2>=kappa/(64R^6)=1/(2^18 G^22); R=4G^3,kappa=G^-4',
       inner='H0>=9/(32C^4),rho^2<=1/2; global orientation comparison 180G^13/C^4>180C^53244',
       exterior='D/p>1/4,D/p<=14000,|J|/p<2000,fourth/p^2>=4/5 give H0>1/64000 and deficit>=1/490000000',
       terminal='N=(-1,0), |Tz/Ttheta|<=4h^4,c0^2<=h give H0>1/2 and rho^2<=16h^9<1/2',
       global_constants='omega_star=1/(640G^13); epsilon_star=1/(2^18 G^22)',
       same_tilt='delta_u=epsilon_star-(1-epsilon_star)/u^2; delta_u>0 suffices with gamma_u=delta_u/2; padded existing u^2>=2/epsilon_star gives gamma=epsilon_star/4',
       band_inputs='f=fstar/5,Fbar=5fupper,t=alpha*fstar/5,Gbar=5fupper*(Ab+Bb),d=alpha*nu*fstar^2/25; direction margins invariant',
       obstruction='existential admissibility of old u does not imply the conservative explicit sufficient gate; failed gate is uncertified, not actual cone failure',
       scope='conditional leading direction and fixed-tilt input estimates only; no changed tilt or actual old-first-band or unconditional pulse/residual/blowup acceptance')

@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('leading_cone_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','cone_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/'harmonic_closure_foundation/audit.py',HERE.parent/'background_gap_certificate/gap_bounds.py',
        HERE.parent/'compact_pressure_bridge/bridge.py',HERE.parent/'actual_band_pulse_bounds/README.md']
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Leading cone-margin record changed')
    if not all(expected['gates'].values()):raise ValueError('Leading cone-margin check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual unchanged-tilt and collar/scale acceptance remain unverified.')
if __name__=='__main__':main()
