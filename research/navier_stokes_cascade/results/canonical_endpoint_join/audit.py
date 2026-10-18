#!/usr/bin/env python3
"""Canonical punctured endpoint join and sharp exterior Gevrey-two bounds."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from endpoint_bounds import load,universal_bounds,universal_checks,canonical_input_review,closed_jet_review,gevrey_review
from endpoint_checks import obstruction_review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='conditional-canonical-punctured-endpoint-join-and-sharp-exterior-Gevrey-two'
DEPENDENCIES=('canonical_local_realization','canonical_primary_foundation','endpoint_foundation',
    'physical_pulse_foundation','compact_mean_foundation','auxiliary_geometry_audit',
    'core_all_order_foundation','annular_exterior_foundation')
FLAGS=dict(same_canonical_local_field_and_complete_fixed_gauge_representatives=True,
    closed_enlarged_atlas_and_next_time_jet_bound_joined_to_canonical_maps=True,
    compatible_punctured_one_sided_potential_swirl_velocity_pressure_jets_conditional=True,
    own_canonical_primary_entrance_roots_and_current_generated_sources_used=True,
    all_pressure_memory_axial_velocity_and_cutoff_curl_terms_retained=True,
    same_two_diagonal_schedules_labels_carriers_supports_width_and_parameters=True,
    exterior_full_pressure_convolution_and_mixed_jet_majorants_retained=True,
    exterior_heat_and_pressure_terminal_Gevrey_two_order_sharp_for_nonzero_C=True,
    exact_zero_residual_bounded_velocity_counterexample_has_no_terminal_trace=True,
    all_order_estimates_distinct_from_finite_exact_and_numerical_controls=True,
    historical_source_bytes_and_evidence_preserved=True,
    old_endpoint_or_global_field_identified_with_canonical_field=False,
    full_canonical_interior_Gevrey_order_or_analyticity_claimed=False,
    flat_residual_or_bounded_velocity_alone_used_to_infer_terminal_trace=False,
    future_heat_integral_analytic_continuation_or_zero_pressure_gluing_licensed=False,
    actual_infinite_coefficients_or_terminal_canonical_field_evaluated=False,
    regular_velocity_or_pressure_at_singular_origin_inferred=False,
    canonical_localization_global_force_energy_and_comparison_joined=False,
    retained_core_background_and_analytical_premises_discharged=False,
    unconditional_complete_residual_or_blowup_proved=False,formal_certificate=False)
REMAINING=['join_axisymmetric_spatial_time_localization_and_zero_initialization_to_this_same_canonical_complete_field',
    'construct_smooth_compact_force_by_compatible_jets_retaining_all_curls_pressure_viscosity_and_nonlinear_terms_not_analytic_heat_or_zero_gluing',
    'join_same_field_whole_space_energy_dissipation_and_exact_global_comparison',
    'independently_discharge_retained_exact_core_background_and_finite_analytical_defining_existence_inverse_gap_premises']

def protocol():
    return dict(parent_commit='4b9bdd40944f851fc0aebc10fb47bf22f30b4e37',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        same_field='canonical_local_realization; own canonical Br entrance,Hc/ac and all current generated source/mean maps; not old realized fields',
        unchanged='Md64,unit viscosity,original h,N,positive lambda,radius,width C_setup^-20,carriers,labels,supports,qstar,both diagonal schedules,fixed operator gauges',
        representatives='Ac=Aseed+sum chi*Aj; Bc=Bseed+sum chi*bj*e_theta; uc=curl Ac+Bc; pc=p0+sum chi*(pc_j-pc_jminus1)',
        closed_atlas='1/4<=q/Q<=4,T>=0; frozen closed meshes; L=1-2h*eta^2>=1-2h; no 1/(1-eta^2)',
        observation='q_e=min(qstar/16,2^(-ell0-2)); tau<=q_e/16; |z|<=(q_e/16)^D; enlarged 4Q<=16q<=6q_e',
        blocking_inequality='sup_(Eplus x [1-delta,1)) |D_x^alpha partial_t^j(Ac,Bc,pc)|<infinity, including NEXT time derivative',
        trace='uniform Cauchy limits; finite Taylor remainder M_(alpha,k+1)*tau^(k+1)/(k+1)!; potential m+2 for velocity Cauchy',
        exterior='K=C*r^(-1-2h)*H(4*tau/r^2), p=-integral_r^infinity K^2/rho; full Leibniz S_j',
        sharp_bounds='h*(j!)^2/j<=A_j<=h*(j+1)*(j!)^2; 2A_j<=S_j<=3*2^j*(j!)^2; pressure divisor 2+4h+2j retained',
        Gevrey_scope='sharp order 2 at fixed h>0,r>0,C!=0 for EXTERIOR heat and complete pressure; no full-interior Gevrey assertion',
        obstruction='u=sin(1/tau)*e1,p=-x1*cos(1/tau)/tau^2; full residual zero and bounded velocity but unequal terminal subsequences',
        diagnostics='exact factorial/atlas/trace/oscillatory controls and reused independent positive Gamma mixed heat/pressure checks, not actual canonical fields',
        scope='conditional same-canonical-field punctured one-sided C-infinity jets; no origin regularity,future velocity,global force,energy,comparison or unconditional theorem')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'canonical_local_realization';old=sys.path[:];sys.path.insert(0,str(p))
    try:return load('canonical_endpoint_canonical_parent',p/'audit.py')
    finally:sys.path[:]=old

@lru_cache(None)
def parent_record():
    r=json.loads((HERE.parent/'canonical_local_realization/evidence.json').read_text())
    parent_module().validate(r);return r

def history():
    m=parent_module()
    for _ in range(5):m=m.parent_module()
    return m.parent_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])

def sources():
    paths=[HERE/n for n in ('README.md','audit.py','endpoint_bounds.py','endpoint_checks.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    paths += [HERE.parent/'endpoint_foundation'/n for n in ('terminal_bounds.py','terminal_checks.py')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def reviews():
    b=universal_bounds();c=universal_checks()
    return deepcopy(dict(canonical_current_realization_input=canonical_input_review(parent_record()),
        canonical_closed_parameter_jet_join=closed_jet_review(),closed_enlarged_terminal_atlas=b.atlas_review(),
        complete_axis_seed=b.seed_review(),complete_mixed_exterior_majorants=b.mixed_review(),
        sharp_exterior_Gevrey_two=gevrey_review(),exact_residual_endpoint_obstruction=obstruction_review(),
        independent_positive_Gamma_checks=c.integral_review(),finite_axis_Cauchy_and_Taylor_checks=c.exact_review()))

def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=reviews();h=history();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g['historical_recursive_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=h,source_hashes=sources(),remaining_obligations=REMAINING[:])

def validate(r):
    if r.get('status')!=STATUS:raise ValueError('Canonical punctured endpoint scope changed')
    compare_records(protocol(),r.get('protocol'),'protocol');compare_records(FLAGS,r.get('flags'),'flags')
    compare_records(REMAINING,r.get('remaining_obligations'),'limits')
    expected=reviews();compare_records(expected,r.get('reviews'),'reviews')
    g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in expected.items()};g['historical_recursive_source_bytes']=True
    compare_records(g,r.get('gates'),'gates')
    if not all(g.values()):raise ValueError('Canonical endpoint gate failed')
    compare_records(sources(),r.get('source_hashes'),'sources');compare_records(history(),r.get('inherited_provenance'),'history')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:
        f=json.loads(a.verify_record.read_text());validate(f);compare_records(r,f)
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('Conditional punctured jets; sharp Gevrey order for exterior only; global and analytical obligations remain.')
if __name__=='__main__':main()
