#!/usr/bin/env python3
"""Canonical local realization relative to retained analytical premises."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from realization import canonical_input_review,representative_review,comparison_review,scale_review,universal_controls
from gauge_checks import review as gauge_review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='conditional-canonical-single-local-realization-and-all-jet-flat-residual'
DEPENDENCIES=('canonical_finite_cycle_join','realization_foundation','compact_mean_foundation',
    'physical_pulse_foundation','background_summation','finite_induction_foundation')
FLAGS=dict(own_canonical_finite_states_and_operator_representatives_used=True,
    canonical_representative_orders_match_same_domain_all_fixed_jets=True,
    one_canonical_real_smooth_divergence_free_local_field_relative_to_retained_lemmas=True,
    all_fixed_Cartesian_jets_of_complete_residual_flat_relative_to_retained_lemmas=True,
    full_curl_direct_swirl_pressure_viscosity_both_crosses_and_square_retained=True,
    same_pressure_base_tails_axial_memory_and_all_operation_differences_retained=True,
    fixed_operator_potential_and_pressure_gauges_before_stage_cutoffs=True,
    exact_rational_whole_interval_log_majorant_available=True,
    one_cutoff_sequence_for_arbitrary_finite_stage_constants_and_all_fixed_jets=True,
    Gaussian_ledger_used_at_one_fixed_comparison_state_never_forward_inverted=True,
    historical_frozen_fields_and_evidence_preserved=True,
    finite_exact_controls_not_actual_fields_or_continuum_proofs=True,
    old_realized_field_identified_with_canonical_field=False,
    arbitrary_stage_gauges_licensed_by_equal_uncut_states=False,
    actual_infinite_coefficients_or_scales_numerically_evaluated=False,
    raw_series_convergence_uniform_constants_or_arbitrary_j_of_q_substitution_proved=False,
    canonical_closed_time_endpoint_atlas_joined=False,
    canonical_global_force_energy_and_comparison_joined=False,
    retained_core_background_and_analytical_premises_discharged=False,
    unconditional_complete_residual_improvement_or_blowup_proved=False,formal_certificate=False)
REMAINING=['join_closed_time_auxiliary_atlas_and_all_complete_one_sided_potential_velocity_pressure_jets_to_this_canonical_field',
    'then_join_same_field_localization_initialization_admissible_compact_force_energy_and_global_comparison',
    'independently_discharge_retained_exact_core_background_and_analytical_defining_existence_inverse_gap_premises']

def protocol():
    return dict(parent_commit='20bb5a62a4770cc93fb7372705b1991d03905a4e',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        canonical_input='own B0=3/2-2kappa,delta=1/2-kappa,rho=1-kappa; same Hc/ac, current sources and means',
        unchanged='Md64,unit viscosity,original h,N,fixed positive lambda,radius,width C_setup^-20,carriers,supports,qstar',
        representatives='Delta u_j=curl(A_j)+b_j*e_theta,partial_theta b_j=0; p_j=p_c^[j]-p_c^[j-1], fixed operator potential/pressure gauges',
        orders='B_(j-1)=delta*j+rho; wave velocity B-kappa; potential C-kappa; full mean/swirl/pressure C-2kappa; radial C+1-2kappa',
        loss='g_j=h*delta*j; ell_m=2A+2m; L_m=ell_(m+1); H_m=L_m+1; K_m=2A+3/2+2m',
        cut='V_j=curl(chi(a_j*q)*A_j)+chi(a_j*q)*b_j*e_theta; P_j=chi(a_j*q)*p_j',
        gauge='grad chi cross grad phi and c_j(t) grad chi retained; a common c(t) cancels from pressure increments',
        scale='log(a_j)>max(previous+1,log(qstar^-1)+1,1,2*(logC_upper+j)/gamma+P^2/gamma^2); gamma=g_j/2',
        tail='2^-J*q^(h*delta*(J+1)/2-L_m),q<=1/(2*a_J),fixed J>=max(1,m)',
        compare='fix m,N before J; J>=m+2; h*delta*(J+1)/2>=N+1+L_(m+2)+H_(m+2); hB_J-K_m>=N+1',
        residual='e_t-Delta e+grad pi+u_J.grad e+e.grad u_J+e.grad e; Gaussian ledger at ONE fixed finite state',
        scope='one canonical field for t<1,q>0 with complete all-fixed-jet flat residual CONDITIONAL on retained analytic premises; no endpoint/global/unconditional theorem')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'canonical_finite_cycle_join';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('canonical_realization_current_parent',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old

@lru_cache(None)
def parent_record():
    r=json.loads((HERE.parent/'canonical_finite_cycle_join/evidence.json').read_text())
    parent_module().validate(r);return r

def history():
    m=parent_module()
    for _ in range(4):m=m.parent_module()
    return m.parent_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])

def sources():
    paths=[HERE/n for n in ('README.md','audit.py','realization.py','gauge_checks.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    paths += [HERE.parent/'realization_foundation/realized_checks.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def reviews():
    return deepcopy(dict(canonical_finite_input=canonical_input_review(parent_record()),
        canonical_complete_representative_ledger=representative_review(),
        rational_whole_interval_scale=scale_review(),fixed_stage_complete_residual_comparison=comparison_review(),
        potential_and_pressure_gauge_controls=gauge_review(),**universal_controls().reviews()))

def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=reviews();h=history();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g['historical_recursive_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=h,source_hashes=sources(),remaining_obligations=REMAINING[:])

def validate(r):
    if r.get('status')!=STATUS:raise ValueError('Canonical realization scope changed')
    compare_records(protocol(),r.get('protocol'),'protocol');compare_records(FLAGS,r.get('flags'),'flags')
    compare_records(REMAINING,r.get('remaining_obligations'),'limits')
    expected=reviews();compare_records(expected,r.get('reviews'),'reviews')
    g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in expected.items()};g['historical_recursive_source_bytes']=True
    compare_records(g,r.get('gates'),'gates')
    if not all(g.values()):raise ValueError('Canonical realization gate failed')
    compare_records(sources(),r.get('source_hashes'),'sources');compare_records(history(),r.get('inherited_provenance'),'history')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:
        f=json.loads(a.verify_record.read_text());validate(f);compare_records(r,f)
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('Conditional canonical field only: endpoint, force, energy and retained analytical premises remain.')
if __name__=='__main__':main()
