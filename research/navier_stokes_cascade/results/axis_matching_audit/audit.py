#!/usr/bin/env python3
"""Reproduce the conditional annulus, its bounds, and the actual-axis gap."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import mpmath as mp
from bounds import HERE,RESULTS,certificate
from moments import Jet,MomentMap,fixture,zeta,relative
from annulus import Annulus,ideal,radius_requirements,accept_entry_bounds

STATUS='conditional-five-moment-annulus-bounded-axis-entry-unverified'
GATES={
    'inherited_heat_checkpoint','retained_Md_failure','outward_constants',
    'interval_refinement','uniform_inverse','exact_root_smallness',
    'finite_relaxed_cone_bound','explicit_axis_data','actual_axis_gap_retained',
    'fixtures_labelled','all_five_moment_rows','all_five_eta_rows',
    'positive_quadratic_grade','graded_remainder','root_refinement',
    'edit_scaled_profile_refinement','sampled_relaxed_cone',
    'full_cone_not_promoted','endpoint_fields_restored','pressure_preserved',
    'original_source_identities','independent_radial_source_derivatives',
    'independent_physical_moment_ODE','independent_quadratic_integrals',
    'tiny_U_square_retained','independent_eta_derivative',
    'linear_only_failure','four_moments_do_not_preserve_pressure',
    'zero_value_nonzero_derivative','frozen_derivative_failure',
    'old_radius_failure','old_sigma_failure','finite_h_reduction_failure','entry_budget_enforced',
    'amplitude_order_and_core_gap'
}


def encode(x):
    if isinstance(x,mp.mpf):return mp.nstr(x,110)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    if hasattr(x,'item'):return x.item()
    return x


def provenance():
    paths=[HERE/n for n in ('README.md','protocol.json','requirements.txt','bounds.py',
                            'moments.py','annulus.py','audit.py','test_audit.py')]
    paths += [RESULTS/'heat_exterior_audit'/n for n in ('README.md','protocol.json','bounds.py','evidence.json')]
    paths += [RESULTS/'outer_pressure_pilot'/'schedule.py',RESULTS/'axial_stress_audit'/'bounds.py']
    return {str(path.relative_to(RESULTS.parent)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def inherited():
    h=json.loads((RESULTS/'heat_exterior_audit'/'evidence.json').read_text())
    return dict(heat_status=h['status'],heat_gates=h['gates'],
                old_Md_failure=h['inherited']['old_Md_failure'])


def root_record(matrix,data,eta,P,degree):
    solution=matrix.solve(data,eta,P)
    c,ce=mp.matrix(solution['root']),mp.matrix(solution['root_eta'])
    scale=max(*(abs(v) for v in solution['target']),*(abs(v) for v in solution['target_eta']))
    residual=matrix.B*c+matrix.quadratic(c,c,eta,P)-mp.matrix(solution['target'])
    derivative=matrix.jacobian(c,eta,P)*ce+matrix.quadratic_eta(c,eta,P)-mp.matrix(solution['target_eta'])
    graded=matrix.graded(solution['target'],eta,P,degree)
    surrogate=[mp.fsum(row[j] for row in graded['components']) for j in range(5)]
    return dict(eta=eta,fixture=data['name'],is_axis_solution=False,**solution,
        row_errors=[abs(v)/scale for v in residual],eta_row_errors=[abs(v)/scale for v in derivative],
        graded=graded,graded_root_gap=max(abs(a-b) for a,b in zip(surrogate,c)),
        nonlinear_grade_norm=max(abs(v) for v in graded['components'][1]),
        U_square_contribution=zeta(eta,P).v*mp.fsum(matrix.QU[j]*c[j]**2 for j in range(2)))


def scalar_gap(rows_a,rows_b,names):
    return max(relative(x,y) for a,b in zip(rows_a,rows_b)
               for key in names for x,y in zip(a[key],b[key]))


def profile_refinement(a,b):
    gaps={}
    for key,background in (('Q_defect',0),('N_over_P2_defect',0),('a',mp.mpf('.8')),('P_times_bs',0)):
        scale=max(*(abs(v[key]-background) for v in a),*(abs(v[key]-background) for v in b))
        gaps[key]=max(abs(x[key]-y[key]) for x,y in zip(a,b))/scale if scale else mp.mpf(0)
    return gaps


def summarize_states(states):
    """Retain extrema and their scope; the complete grid is fixed in protocol."""
    result=[]
    for label in ('moderate','Md4-diagnostic'):
        for name in ('odd','edge'):
            rows=[v for v in states if v['scale_name']==label and v['fixture']==name]
            worst=min(rows,key=lambda v:v['Pc'])
            result.append(dict(scale_name=label,fixture=name,samples=len(rows),
                minimum_Q=min(v['Q'] for v in rows),minimum_E_over_P=min(v['E_over_P'] for v in rows),
                minimum_G=min(v['G'] for v in rows),minimum_Pc=worst['Pc'],
                minimum_Pc_at=[worst['y'],worst['eta']],maximum_vs=max(v['vs'] for v in rows),
                maximum_N_over_P2=max(abs(v['N_over_P2']) for v in rows),
                maximum_P_times_bs=max(abs(v['P_times_bs']) for v in rows),
                strict_relaxed_cone=all(v['strict_relaxed_cone'] for v in rows),
                full_admissible_cone_claimed=any(v['full_admissible_cone_claimed'] for v in rows)))
    return result


def omitted_pressure_control(matrix,record,eta,P):
    """Solve four equations with four bumps, then measure the missing fifth."""
    target=mp.matrix(record['target'])
    restricted=matrix.B[:4,:4]
    c=mp.matrix(5,1);c[:4,0]=restricted**-1*target[:4,0]
    for _ in range(10):
        residual=matrix.B*c+matrix.quadratic(c,c,eta,P)-target
        change=mp.lu_solve(matrix.jacobian(c,eta,P)[:4,:4],residual[:4,0])
        for j in range(4):c[j]-=change[j]
    residual=matrix.B*c+matrix.quadratic(c,c,eta,P)-target
    scale=max(abs(v) for v in target)
    return dict(four_row_error=max(abs(residual[j]) for j in range(4))/scale,
                pressure_row_gap=abs(residual[4])/scale)


def controls(matrix,annulus,record,P):
    eta=mp.mpf('.5')
    c=matrix.inverse*mp.matrix(record['target'])
    q=matrix.quadratic(c,c,eta,P)
    epsilon=mp.mpf('1e-16')
    d=fixture(0);s=matrix.solve(d,0,P)
    correct=annulus.state('-5',0,d,s)
    frozen=annulus.state('-5',0,d,s,True)
    zero=dict(entry=[Jet(0)]*5,g=Jet(0),epsilon=epsilon)
    none=dict(root=[mp.mpf(0)]*5,root_eta=[mp.mpf(0)]*5)
    small=Annulus(matrix,annulus.h,P,100).state('-8',0,zero,none)
    reference=annulus.state('-7',eta,zero,none)
    # The rejected hand reduction used -h*eta^2/2 instead of +2*h*eta^2.
    correct_ideal=ideal(eta,annulus.h,P,mp.exp(-7))['Q']
    wrong_ideal=correct_ideal-mp.mpf('2.5')*annulus.h*eta**2
    j=mp.mpf('1e-18');old_sigma=mp.mpf('.002')
    return dict(linear_pressure_gap=q[4],linear_pressure_gap_over_coefficient_squared=q[4]/max(abs(v) for v in c)**2,
        omitted_pressure=omitted_pressure_control(matrix,record,eta,P),
        zero_value_roots=s['root'],nonzero_eta_roots=s['root_eta'],
        correct_Q_gap_over_entry=abs(correct['Q_defect'])/epsilon,
        frozen_Q_gap_over_entry=abs(frozen['Q_defect'])/epsilon,
        old_radius_Pc=small['Pc'],old_radius_vs=small['vs'],
        wrong_ideal_Q_gap_over_h=abs(reference['Q']-wrong_ideal)/annulus.h,
        old_sigma_chi_at_zero=j*j/(j*j+old_sigma**2),
        sigma_control_in_small_Z_region=bool(mp.mpf('4.505')*j<j*100**2/100),
        old_j_budget_ratio=mp.mpf('.05')/(epsilon/12),
        rejected_large_entry=not accept_entry_bounds(['1e-14']*5,0),
        early_amplitude=radius_requirements(50,mp.log(16),1000),
        later_amplitude=radius_requirements(120,mp.log(16),1000))


def build():
    p=json.loads((HERE/'protocol.json').read_text())
    mp.mp.dps=400
    cs=[certificate(d) for d in p['interval_digits']]
    interval_gap=max(relative(cs[0][key],cs[1][key]) for key in
        ('Pc_lower','pressure_mass_upper','axis_chi_lower','axis_complex_radius_lower'))
    cases=[]
    for index,order in enumerate(p['orders']):
        mp.mp.dps=p['digits'][index]
        print('Five-moment annulus: order='+str(order),file=sys.stderr,flush=True)
        matrix=MomentMap(order,p['panels'])
        T=mp.exp(p['large_scale_diagnostic_Md'])+10
        scales=[('moderate',mp.mpf(p['moderate_P']),mp.mpf(p['moderate_h'])),
                ('Md4-diagnostic',mp.exp(T+1),mp.exp(-8*T))]
        roots=[];states=[];endpoints=[];algebra_gaps=[]
        for label,P,h in scales:
            annulus=Annulus(matrix,h,P,p['matching_XR_min'])
            for name in p['fixtures']:
                for eta in map(mp.mpf,p['angles']):
                    data=fixture(eta,name,p['entry_C1_bound'])
                    record=root_record(matrix,data,eta,P,p['graded_degree'])
                    record['scale_name']=label
                    roots.append(record)
                    solution={key:record[key] for key in ('root','root_eta')}
                    for y in p['stress_y']:
                        row=annulus.state(y,eta,data,solution)
                        row.update(scale_name=label,fixture=name)
                        states.append(row)
                    endpoint=annulus.fields('-5',eta,data,solution)
                    endpoints.append(dict(scale_name=label,fixture=name,eta=eta,
                        U_gap=abs(endpoint['U'].v-4*eta),relative_E_gap=abs(endpoint['epsilon'].v),
                        slope_gap=abs(endpoint['epsilon_t']),axial_slope_gap=abs(endpoint['Ut']),
                        pressure_defect_over_P2=states[-1]['pressure_defect_over_P2']))
            zero=dict(entry=[Jet(0)]*5,g=Jet(0),epsilon=mp.mpf(p['entry_C1_bound']))
            none=dict(root=[mp.mpf(0)]*5,root_eta=[mp.mpf(0)]*5)
            for eta in map(mp.mpf,p['angles']):
                r=annulus.state('-7',eta,zero,none)
                reference=ideal(eta,h,P,mp.exp(-7))
                algebra_gaps.extend([abs(r['Q']-reference['Q']),abs(r['N_over_P2']-reference['N_over_P2'])])
        cases.append(dict(roots=roots,states=states,endpoints=endpoints,algebra_gap=max(algebra_gaps)))
    mp.mp.dps=max(p['digits'])
    root_gap=scalar_gap(cases[0]['roots'],cases[1]['roots'],('root','root_eta','target','target_eta'))
    stress_gap=profile_refinement(cases[0]['states'],cases[1]['states'])
    P=mp.mpf(p['moderate_P']);eta=mp.mpf('.5')
    data=fixture(eta);annulus=Annulus(matrix,p['moderate_h'],P)
    record=next(x for x in cases[-1]['roots'] if x['scale_name']=='moderate' and x['fixture']=='odd' and x['eta']==eta)
    solution={key:record[key] for key in ('root','root_eta')}
    print('Independent original moment densities and source derivatives',file=sys.stderr,flush=True)
    with mp.workdps(p['independent_digits']):
        ode=[matrix.independent_integrals(record['root'],eta,P,step) for step in p['ode_steps']]
        print('Independent quadratic terms: tanh-sinh',file=sys.stderr,flush=True)
        quadratic=matrix.tanh_sinh_quadratic(record['root'],eta,P)
        source=[]
        for y in ('-7.5','-5.72','-5.48'):
            print('Original source derivative: y='+y,file=sys.stderr,flush=True)
            source.append(annulus.source_check(y,eta,data,solution))
    print('Independent angular derivative and failure controls',file=sys.stderr,flush=True)
    dx=mp.mpf('1e-4')
    nearby=[matrix.solve(fixture(eta+k*dx),eta+k*dx,P)['root'] for k in (-2,-1,1,2)]
    derivative=[(nearby[0][i]-8*nearby[1][i]+8*nearby[2][i]-nearby[3][i])/(12*dx) for i in range(5)]
    eta_gap=max(relative(a,b) for a,b in zip(derivative,record['root_eta']))
    result=encode(dict(status=STATUS,protocol=p,provenance=provenance(),inherited=inherited(),
        certificate=cs[-1],interval_gap=interval_gap,roots=cases[-1]['roots'],
        stress_summary=summarize_states(cases[-1]['states']),endpoints=cases[-1]['endpoints'],
        algebra_gap=cases[-1]['algebra_gap'],root_refinement_gap=root_gap,
        profile_refinement=stress_gap,independent_ODE=ode,independent_quadratic=quadratic,
        source_checks=source,eta_derivative_gap=eta_gap,controls=controls(matrix,annulus,record,P),
        interpretation=dict(actual_axis_entry_verified=False,activation_collar_constructed=False,
            Lambda_chosen=False,full_admissible_cone_realized=False,full_PDE_corrected=False,
            admissible_force_verified=False,blowup_verified=False,
            manufactured_fixtures_are_axis_solutions=False)))
    result['gates']={k:bool(v) for k,v in predicates(result).items()}
    return result


def predicates(r):
    n=mp.mpf;p=r['protocol'];limits={k:n(v) for k,v in p['checks'].items()}
    b=r['certificate'];roots=r['roots'];controls=r['controls']
    nonzero=[v for v in roots if n(v['graded']['majorant_a'])>0]
    final=r['endpoints']
    return dict(
        inherited_heat_checkpoint=r['inherited']['heat_status']=='heat-compensated-outer-reference-bounded-axis-and-global-blowup-unverified'
            and len(r['inherited']['heat_gates'])==33 and all(v is True for v in r['inherited']['heat_gates'].values()),
        retained_Md_failure=n(r['inherited']['old_Md_failure']['Pc_over_ps1'])<0,
        outward_constants=b['status']=='conditional-five-moment-annulus-bounded' and all(v is True for v in b['checks'].values()),
        interval_refinement=n(r['interval_gap'])<limits['interval_relative_gap'],
        uniform_inverse=b['inverse_bound']==1000 and max(n(v) for v in b['E_inverse_row_upper'])<1000,
        exact_root_smallness=0<n(b['coefficient_C1_bound'])<=n('6.0000001e-13') and n(b['contraction_upper'])<1,
        finite_relaxed_cone_bound=n(b['G_lower'])>=n('.99') and n(b['vs_upper'])<1 and b['XR_min']==10000 and n(b['Pc_lower'])>2,
        explicit_axis_data=n(b['axis_j'])==n('1e-18') and n(b['axis_sigma'])==n('1e-20') and n(b['axis_chi_lower'])>n('.99') and n(b['axis_complex_radius_lower'])>0,
        actual_axis_gap_retained=all(v is False for v in r['interpretation'].values()) and all(b[k] is False for k in ('actual_axis_entry_verified','full_admissible_cone_realized','blowup_verified')),
        fixtures_labelled=len(roots)==2*len(p['fixtures'])*len(p['angles']) and all(v['is_axis_solution'] is False for v in roots),
        all_five_moment_rows=all(len(v['row_errors'])==5 and max(n(x) for x in v['row_errors'])<limits['normalized_moment_residual'] for v in roots),
        all_five_eta_rows=all(len(v['eta_row_errors'])==5 and max(n(x) for x in v['eta_row_errors'])<limits['normalized_derivative_residual'] for v in roots),
        positive_quadratic_grade=all(n(v['nonlinear_grade_norm'])>0 for v in nonzero),
        graded_remainder=all(0<n(v['graded']['remainder']) and n(v['graded_root_gap'])<n(v['graded']['remainder']) for v in nonzero),
        root_refinement=n(r['root_refinement_gap'])<limits['quadrature_relative_gap'],
        edit_scaled_profile_refinement=max(n(v) for v in r['profile_refinement'].values())<limits['quadrature_relative_gap'],
        sampled_relaxed_cone=len(r['stress_summary'])==2*len(p['fixtures']) and sum(v['samples'] for v in r['stress_summary'])==len(roots)*len(p['stress_y'])
            and all(v['strict_relaxed_cone'] is True and n(v['maximum_vs'])<1 and n(v['minimum_Pc'])>2 for v in r['stress_summary']),
        full_cone_not_promoted=all(v['full_admissible_cone_claimed'] is False for v in r['stress_summary']),
        endpoint_fields_restored=len(r['endpoints'])==len(roots) and all(all(n(v[k])==0 for k in ('U_gap','relative_E_gap','slope_gap','axial_slope_gap')) for v in r['endpoints']),
        pressure_preserved=all(abs(n(v['pressure_defect_over_P2']))<limits['normalized_moment_residual']*n(p['entry_C1_bound']) for v in final),
        original_source_identities=n(r['algebra_gap'])<limits['source_identity_gap'],
        independent_radial_source_derivatives=all(max(n(x) for x in v['normalized_source_gaps'])<limits['differentiated_moment_gap'] for v in r['source_checks']),
        independent_physical_moment_ODE=len(r['independent_ODE'])==len(p['ode_steps']) and all(max(n(x) for x in v['gaps'])<limits['independent_ODE_gap'] for v in r['independent_ODE']),
        independent_quadratic_integrals=max(n(x) for x in r['independent_quadratic']['gaps'])<limits['independent_tanh_sinh_gap'],
        tiny_U_square_retained=all(n(v['U_square_contribution'])>0 for v in nonzero) and n(r['independent_quadratic']['U_square'])>0,
        independent_eta_derivative=n(r['eta_derivative_gap'])<limits['eta_difference_gap'],
        linear_only_failure=n(controls['linear_pressure_gap'])>0 and n(controls['linear_pressure_gap_over_coefficient_squared'])>n('1e-4'),
        four_moments_do_not_preserve_pressure=n(controls['omitted_pressure']['four_row_error'])<limits['normalized_moment_residual'] and n(controls['omitted_pressure']['pressure_row_gap'])>n('.001'),
        zero_value_nonzero_derivative=all(n(v)==0 for v in controls['zero_value_roots']) and any(n(v)!=0 for v in controls['nonzero_eta_roots']),
        frozen_derivative_failure=n(controls['correct_Q_gap_over_entry'])<limits['normalized_moment_residual'] and n(controls['frozen_Q_gap_over_entry'])>n('.001'),
        old_radius_failure=n(controls['old_radius_Pc'])<2 and n(controls['old_radius_vs'])<1 and n(b['old_XR100_Pc_upper'])<2,
        old_sigma_failure=n(controls['old_sigma_chi_at_zero'])<n('.99') and controls['sigma_control_in_small_Z_region'] is True,
        finite_h_reduction_failure=abs(n(controls['wrong_ideal_Q_gap_over_h'])-n('.625'))<limits['normalized_moment_residual'],
        entry_budget_enforced=controls['rejected_large_entry'] is True and n(controls['old_j_budget_ratio'])>1,
        amplitude_order_and_core_gap=controls['early_amplitude']['transition_before_restoration'] is False and controls['later_amplitude']['transition_before_restoration'] is True and controls['later_amplitude']['core_moments_verified'] is False)


def validate(r):
    if r.get('status')!=STATUS:raise ValueError('Scientific scope changed')
    if r.get('protocol')!=json.loads((HERE/'protocol.json').read_text()):raise ValueError('Protocol changed')
    if r.get('provenance')!=provenance():raise ValueError('Provenance changed')
    if r.get('inherited')!=inherited():raise ValueError('Inherited checkpoint changed')
    if set(r.get('gates',{}))!=GATES or not all(v is True for v in r['gates'].values()):
        raise ValueError('Missing or failed gate')
    with mp.workdps(400):actual=predicates(r)
    if set(actual)!=GATES or not all(v is True for v in actual.values()):
        raise ValueError('Data fail gates: '+str([k for k,v in actual.items() if not v]))


def compare(a,b,path=''):
    if isinstance(a,dict):
        if not isinstance(b,dict) or set(a)!=set(b):raise ValueError('Keys differ: '+path)
        for k in a:compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):raise ValueError('Lengths differ: '+path)
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
    elif isinstance(a,bool):
        if a is not b:raise ValueError('Boolean differs: '+path)
    elif a!=b:
        try:x,y=mp.mpf(a),mp.mpf(b)
        except (ValueError,TypeError):raise ValueError('Value differs: '+path) from None
        if not mp.isfinite(x) or not mp.isfinite(y):raise ValueError('Nonfinite value: '+path)
        if path.endswith('/evaluations'):return
        ode='independent_ODE' in path
        noise=any(k in path for k in ('gap','error','residual','/ledger','pressure_defect'))
        delta=abs(x-y)/(1+max(abs(x),abs(y))) if ode else (abs(x-y) if noise else relative(x,y))
        if delta>mp.mpf('1e-9' if ode else '1e-18'):raise ValueError('Numeric difference: '+path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify-record',type=Path)
    args=parser.parse_args()
    record=build()
    print(json.dumps(dict(status=record['status'],gates=len(record['gates']),
        failed=[k for k,v in record['gates'].items() if not v],
        root_refinement=record['root_refinement_gap'],profile_refinement=record['profile_refinement'])),flush=True)
    validate(record)
    if args.verify_record:
        old=json.loads(args.verify_record.read_text());validate(old);compare(old,record)
    args.output.write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':main()
