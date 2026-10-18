#!/usr/bin/env python3
"""Canonical pressure recovery, mixed moments and conditional competitor comparison."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from join_bounds import load,universal_bounds,universal_checks,canonical_interface,moment_review,interval_review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='conditional-canonical-pressure-recovery-mixed-moment-and-global-comparison'
DEPENDENCIES=('canonical_force_join','comparison_foundation','global_solution_comparison',
    'canonical_endpoint_join','canonical_local_realization','core_all_order_foundation',
    'core_contraction_audit','global_pressure_moment_foundation')
FLAGS=dict(same_canonical_complete_field_force_zero_datum_and_pressure_used=True,
    ordinary_comparison_hypotheses_checked_on_canonical_force_parent=True,
    full_individual_forced_pressure_and_exact_difference_cancellation_retained=True,
    time_averaged_momentum_recovers_gradient_without_pressure_growth_hypothesis=True,
    compact_mixed_pressure_moment_expansion_with_uniform_radius_four_error=True,
    compact_divergence_free_tensor_witness_has_nonzero_trace_free_moment=True,
    complete_expanding_ball_energy_and_integrated_local_dissipation_decay_radius_two=True,
    canonical_global_competitor_comparison_and_exclusion_relative_to_retained_premises=True,
    canonical_protected_path_uses_same_localization_radius=True,
    all_historical_fields_sources_and_frozen_evidence_preserved=True,
    no_global_competitor_gradient_dissipation_or_energy_inequality_hypothesis_added=True,
    compactness_and_divergence_freedom_assumed_to_force_moment_isotropy=False,
    individual_pressure_replaced_by_unforced_Riesz_formula=False,
    Poisson_alone_used_to_remove_harmonic_pressure_gradient=False,
    old_special_field_identified_with_canonical_field=False,
    tensor_witness_claimed_same_force_NS_flow_pair=False,
    actual_infinite_fields_or_competitor_numerically_evaluated=False,
    comparison_constants_claimed_uniform_through_terminal_time=False,
    future_velocity_or_pressure_constructed=False,
    retained_core_background_and_finite_analytical_premises_discharged=False,
    unconditional_admissibility_residual_or_blowup_proved=False,formal_certificate=False)
REMAINING=['trace_exact_core_pressure_family_to_completed_global_backward_pressure_at_same_parameters_and_prove_first_outstanding_defining_map_identity_or_precise_obstruction',
    'independently_discharge_retained_background_and_finite_defining_existence_inverse_gap_premises_before_unconditional_admissibility_residual_or_blowup_claims']

def protocol():
    return dict(parent_commit='634356dff6d7b6c0b350c714cd11403d52c216ae',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        fixed='Md64,unit viscosity,original h,N,positive lambda,pulse radius,width C_setup^-20,labels,carriers,supports,qstar,both diagonal schedules,fixed operator gauges',
        field='canonical_force_join: THIS u=curl(c*Ac)+c*Bc,p=c*pc,full twelve-term force and zero datum; not historical force or field',
        competitor='same datum and force; smooth on compact sets; finite sup_[0,T] L2 for every T<1; no global pressure gradient dissipation or solution energy inequality premise',
        pressure='g=vv-uu=ww+h,h=wu+uw; -Delta pi=DDg only AFTER same force cancels; individual -Delta p=DD(uu)-div f',
        recovery='pi_star in H^-2; time-averaged momentum puts harmonic gradient difference in H^-3; Fourier weighted L2 function supported at zero is zero',
        multipole='pi_h=k(x):M+error; M=int h; |error|<=384/(4*pi)*L*H_T*R^-4 for R>=max(1,2L); directional derivative <=24/(4*pi)*r^-4',
        witness='V=(-x2,x1,0)*rho(|x|^2) smooth compact radial; u=w=V tensor only; M=2a diag(1,1,0),axis contraction -4a/(4*pi*R^3),a>0',
        cutoff='phi_R^8; local Sobolev only; quadratic main L^(3/2); commutator L^(4/3) with both near cancellation and far cap',
        radius='Young decays 4,14/5; offsets and mixed flux 5/2; complete Laplacian remainder 2',
        energy='0.5 E_R_prime+0.5 D_R^2<=L_T E_R+C_T R^-2; E_R<=2 C_T T exp(2 L_T T) R^-2',
        dissipation='integral_0^t D_R^2<=2 C_T t*(1+2 L_T T exp(2 L_T T))*R^-2',
        order='fix T<1,then R->infinity,then approach protected core; no uniform terminal comparison constant',
        core='tau<=min(tau0/2,rL^2/(8*Xin)); canonical localization radius rL, retained e0,Aexp positive',
        external='Riesz Lp boundedness; Knapp Advanced Real Analysis III.5 Theorem 3.26 and III.6 Example 1, printed pp83,93; rechecked 2026-10-18',
        scope='conditional canonical downstream comparison/exclusion implication only; exact-core pressure identification and retained analytical premises remain')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'canonical_force_join';old=sys.path[:];sys.path.insert(0,str(p))
    try:return load('canonical_comparison_force_parent',p/'audit.py')
    finally:sys.path[:]=old

@lru_cache(None)
def parent_record():
    r=json.loads((HERE.parent/'canonical_force_join/evidence.json').read_text());parent_module().validate(r);return r

@lru_cache(None)
def provenance_tools():
    m=parent_module()
    for _ in range(12):
        if hasattr(m,'parent_tools'):return m.parent_tools()
        m=m.parent_module()
    raise ValueError('Historical provenance helper not found in explicit parent chain')

def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])

def sources():
    paths=[HERE/n for n in ('README.md','audit.py','join_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    paths += [HERE.parent/'comparison_foundation'/n for n in ('compact_bounds.py','compact_checks.py')]
    paths += [HERE.parent/'global_solution_comparison/comparison_bounds.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def reviews():
    b=universal_bounds()
    return deepcopy(dict(canonical_force_comparison_interface=canonical_interface(parent_record()),
        pressure_gradient_and_retained_mixed_memory=b.pressure_review(),
        complete_radius_two_budget=b.radius_review(),
        new_compact_mixed_moment_controls=moment_review(),
        fixed_interval_energy_dissipation_and_canonical_core=interval_review(),
        **universal_checks().reviews()))

def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=reviews();h=history();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g['historical_recursive_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=h,source_hashes=sources(),remaining_obligations=REMAINING[:])

def validate(r):
    if r.get('status')!=STATUS:raise ValueError('Canonical comparison scope changed')
    for expected,actual,name in ((protocol(),r.get('protocol'),'protocol'),(FLAGS,r.get('flags'),'flags'),
        (REMAINING,r.get('remaining_obligations'),'remaining'),(reviews(),r.get('reviews'),'reviews'),
        (sources(),r.get('source_hashes'),'sources'),(history(),r.get('inherited_provenance'),'history')):
        compare_records(expected,actual,name)
    expected={k:bool(v['checks']) and all(v['checks'].values()) for k,v in reviews().items()}
    expected['historical_recursive_source_bytes']=bool(history()['bindings']) and all(history()['bindings'].values())
    compare_records(expected,r.get('gates'),'gates')
    if not all(expected.values()):raise ValueError('Canonical comparison gate failed')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:
        f=json.loads(a.verify_record.read_text());validate(f);compare_records(r,f)
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('Same canonical comparison is conditional; retained analytical premises and unconditional blow-up remain open.')
if __name__=='__main__':main()
