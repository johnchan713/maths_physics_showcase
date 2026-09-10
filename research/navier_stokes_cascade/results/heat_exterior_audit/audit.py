#!/usr/bin/env python3
"""Reproduce compensated heat-reference bounds and explicit failure controls."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import mpmath as mp
from bounds import heat_bound,HERE,RESULTS
from construction import (TailMoments,Compensation,Rule,relative,selected_fixture,
                          heat_polynomial,heat_integral,pole_control)
from exterior import Exterior

STATUS = 'heat-compensated-outer-reference-bounded-axis-and-global-blowup-unverified'
PARENTS = ('intermediate_decay_audit','pulse_moment_audit','pulse_stress_audit','post_pulse_stress_audit')
GATES = {
    'inherited_bounds','outward_constants','interval_refinement','exact_moment_scope',
    'positive_correction_budget','finite_radius_cone','finite_shear_excess',
    'positive_stress_error_budget','moderate_moment_rows','moderate_nonlinear_root',
    'graded_root_rows','positive_graded_remainder','positive_taylor_target_errors',
    'selected_target_refinement','moderate_target_refinement','independent_moment_ODE',
    'finite_quadratic_grade_retained','quadratic_rounding_failure',
    'infinite_angular_tail_not_truncated','finite_h_tail_change_retained',
    'pole_value_and_derivative_distinguished','pole_stress_restored','pole_omission_failure',
    'independent_heat_kernel','finite_heat_polynomial_residual',
    'positive_outer_stress','outer_direction_bound','exterior_refinement',
    'independent_flat_integrals','positive_edge_coefficients','smooth_direction_limit',
    'flat_underflow_failure','global_proof_gap_retained'
}


def encode(x):
    if isinstance(x,mp.mpf): return mp.nstr(x,110)
    if isinstance(x,dict): return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [encode(v) for v in x]
    if hasattr(x,'item'): return x.item()
    return x


def provenance():
    paths = [HERE/n for n in ('README.md','protocol.json','requirements.txt','bounds.py','construction.py','exterior.py','audit.py','test_audit.py')]
    paths += [RESULTS/p/n for p in PARENTS for n in ('README.md','bounds.py','protocol.json','evidence.json')]
    paths += [RESULTS/'outer_pressure_pilot'/'schedule.py',RESULTS/'axial_stress_audit'/'bounds.py',RESULTS.parent/'GOAL.md']
    return {str(p.relative_to(RESULTS.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def inherited():
    result = {}
    for name in PARENTS:
        p = json.loads((RESULTS/name/'evidence.json').read_text())
        result[name] = dict(status=p['status'],gates=p['gates'],certificate=p['certificate'])
    post = json.loads((RESULTS/'post_pulse_stress_audit'/'evidence.json').read_text())
    result['old_Md_failure'] = post['inherited']['old_Md_failure']
    result['old_lambda_failure'] = post['inherited']['old_lambda_failure']
    return result


def root_record(tail,matrix,eta,log_xstar,log_rh,degree):
    record = tail.target(eta,log_xstar,log_rh)
    target = mp.matrix(record['target'])
    linear = matrix.inverse*target
    graded = matrix.graded(target,degree)
    root = matrix.newton(target)
    derivative = matrix.derivative(root,record['target_eta'])
    scale = graded['scale']
    residual = matrix.B*root+matrix.quadratic(root,root)-target
    surrogate = [mp.fsum(v[j] for v in graded['components']) for j in range(3)] if scale else list(root)
    nonlinear = max(abs(v) for v in graded['components'][1]) if scale else mp.mpf(0)
    record.update(root=list(root),root_eta=list(derivative),graded=graded,
        moment_row_errors=[abs(v)/(scale or 1) for v in residual],
        surrogate_newton_gap=max(abs(v-w) for v,w in zip(surrogate,root)),
        nonlinear_grade_norm=nonlinear,
        quadratic_grade_lost_on_addition=bool(scale and all(linear[j]+graded['components'][1][j]==linear[j] for j in range(3))),
        truncation_budget=graded['remainder']+100*max(record['taylor_target_C1_error']))
    return record


def row_refinement(coarse,fine):
    gaps = []
    for a,b in zip(coarse,fine):
        for name in ('target','target_eta','root','root_eta'):
            gaps += [relative(x,y) for x,y in zip(a[name],b[name])]
        if a['graded']['scale']:
            gaps += [relative(x,y) for ar,br in zip(a['graded']['components'],b['graded']['components']) for x,y in zip(ar,br)]
    return max(gaps)


def collar_limit_gap(rows,angles,delta_count):
    """Use the declared grid indices, not equality of rounded mpf decimals."""
    if delta_count<1 or len(rows)!=len(angles)*delta_count:
        raise ValueError('Incomplete collar grid')
    selected = rows[delta_count-1::delta_count]
    return max(relative(x['ratio_over_delta6'],x['expected_ratio_over_delta6'])
               for eta,x in zip(angles,selected) if mp.mpf(eta)!=0)


def predicates(r):
    n = mp.mpf
    p,b = r['protocol'],r['certificate']
    limits = {k:n(v) for k,v in p['checks'].items()}
    selected,moderate = r['selected']['rows'],r['moderate']['rows']
    nonzero = [x for x in selected if n(x['graded']['scale'])>0]
    all_nonzero = [x for x in selected+moderate if n(x['graded']['scale'])>0]
    pole = selected[-1]
    pc = r['pole']
    return {
        'inherited_bounds': all(len(r['inherited'][name]['gates'])==count and all(v is True for v in r['inherited'][name]['gates'].values())
            for name,count in zip(PARENTS,(21,34,28,29))),
        'outward_constants': b['status']=='heat-compensated-outer-reference-bounded' and all(v is True for v in b['checks'].values()),
        'interval_refinement': n(r['interval_gap'])<limits['interval_relative_gap'],
        'exact_moment_scope': b['exact_moments']==['M','J','S(infinity)','Cp(infinity)','renormalized I'],
        'positive_correction_budget': 0<n(b['coefficient_universal_upper'])<n('1e-4'),
        'finite_radius_cone': b['XR']==100 and n(b['finite_ps1_lower'])==n('1e12') and n(b['finite_normalized_cone_gap_lower'])>n('.8'),
        'finite_shear_excess': b['a_minus_two_floor']=='h' and 0<n(b['shear_error_over_h_upper'])<n('1e-6000'),
        'positive_stress_error_budget': 0<n(b['absolute_w_error_upper'])<n('1e-6000') and 0<n(b['relative_Q_error_upper'])<n('1e-6000'),
        'moderate_moment_rows': [n(x['eta']) for x in moderate]==[n(v) for v in p['moderate_angles']] and all(max(n(v) for v in x['moment_row_errors'])<limits['graded_row_error'] for x in moderate),
        'moderate_nonlinear_root': all(n(x['surrogate_newton_gap'])<n(x['graded']['remainder']) for x in moderate if n(x['graded']['scale'])>0),
        'graded_root_rows': [n(x['eta']) for x in selected]==[n(v) for v in p['selected_angles']] and len(nonzero)==2 and all(len(x['graded']['row_errors'])==p['graded_degree'] and max(n(v) for v in x['graded']['row_errors'])<limits['graded_row_error'] for x in all_nonzero),
        'positive_graded_remainder': all(0<n(x['graded']['remainder'])<n(x['graded']['scale']) for x in all_nonzero),
        'positive_taylor_target_errors': all(all(n(v)>0 for v in x['taylor_target_C1_error']) for x in selected+moderate),
        'selected_target_refinement': n(r['selected']['refinement_gap'])<limits['quadrature_relative_gap'],
        'moderate_target_refinement': n(r['moderate']['refinement_gap'])<limits['quadrature_relative_gap'],
        'independent_moment_ODE': [x['max_step'] for x in r['independent_ODE']]==p['ode_steps'] and all(0<n(row['quadratic_omission_upper']) and max(n(v) for v in row['gaps'])<limits['ode_normalized_gap'] for row in r['independent_ODE']),
        'finite_quadratic_grade_retained': all(n(x['nonlinear_grade_norm'])>0 for x in nonzero),
        'quadratic_rounding_failure': all(x['quadratic_grade_lost_on_addition'] is True for x in nonzero),
        'infinite_angular_tail_not_truncated': abs(n(nonzero[1]['normalized']['first_order_angular_tail_omitted_at_y1000']))>1,
        'finite_h_tail_change_retained': n(nonzero[1]['normalized']['angular_first_infinite_change_from_h0'])>0,
        'pole_value_and_derivative_distinguished': all(n(v)==0 for v in pole['root']) and any(n(v)!=0 for v in pole['root_eta']),
        'pole_stress_restored': pc['heat_value_edit']==0 and n(pc['ps1_minus_two'])>0 and abs(n(pc['stress_over_F']))<n('1e-300') and n(pc['shear_excess'])>0,
        'pole_omission_failure': pc['omitted_derivative_ps1']==0 and n(pc['omitted_derivative_stress_over_F'])<n('-1.9'),
        'independent_heat_kernel': all(n(x['relative_gap'])<limits['kernel_relative_gap'] and n(x['absolute_gap'])<n(x['normalized_remainder'])*n('1.1') for x in r['kernel']),
        'finite_heat_polynomial_residual': all(abs(n(x['ode_residual']))>0 and n(x['ode_relative_gap'])<limits['kernel_ode_relative_gap'] for x in r['kernel']),
        'positive_outer_stress': len(r['exterior']['rows'])==len(p['exterior_angles'])*len(p['exterior_y']) and all(n(x['theta_over_rho'])>0 and n(x['a_minus_two'])>n(p['exterior_h']) for x in r['exterior']['rows']),
        'outer_direction_bound': all(abs(n(x['ratio']))<n(x['ratio_upper']) and n(x['directional_gap'])>n('1.99') for x in r['exterior']['rows']),
        'exterior_refinement': n(r['exterior']['refinement_gap'])<limits['exterior_refinement_relative_gap'],
        'independent_flat_integrals': max(n(v) for v in r['exterior']['path_gaps'])<limits['exterior_path_relative_gap'],
        'positive_edge_coefficients': all(n(x['btheta'])>0 and n(x['expected_btheta'])>0 for x in r['exterior']['collar']),
        'smooth_direction_limit': n(r['exterior']['limit_gap'])<limits['collar_limit_relative_gap'] and all(x['theta']==x['axial']==0 and x['direction']==[1,0] and n(x['positive_angular_coefficient'])>0 for x in r['exterior']['endpoints']),
        'flat_underflow_failure': any(n(x['positive_flat_factor'])>0 and x['binary64_flat_factor']==0 for x in r['exterior']['collar']),
        'global_proof_gap_retained': r['status']==STATUS and r['selected']['fixture']['old_axial_failure_retained'] is True and n(r['inherited']['old_Md_failure']['Pc_over_ps1'])<0,
    }


def build():
    p = json.loads((HERE/'protocol.json').read_text())
    mp.mp.dps = 160
    bounds = [heat_bound(p['selected_Md'],d,p['XR']) for d in p['interval_digits']]
    interval_gap = max(relative(bounds[0][k],bounds[1][k]) for k in ('absolute_w_error_upper','coefficient_universal_upper','finite_normalized_cone_gap_lower'))
    selected_rows,moderate_rows,ext_rows,collars = [],[],[],[]
    for index,order in enumerate(p['orders']):
        mp.mp.dps = p['selected_digits'][index]
        print('Selected-family heat moments: order='+str(order),file=sys.stderr,flush=True)
        fixture = selected_fixture(order,p['schedule_panels'])
        tail = TailMoments(fixture['h'],fixture['log_xtail'],order,p['panels'],p['taylor_degrees'][index])
        matrix = Compensation(fixture['lam'],Rule(order,p['panels']))
        selected_rows.append([root_record(tail,matrix,eta,fixture['log_xstar'],fixture['log_hpow_ratio'],p['graded_degree']) for eta in p['selected_angles']])
        mp.mp.dps = p['moderate_digits']
        print('Moderate heat compensation: order='+str(order),file=sys.stderr,flush=True)
        h = mp.mpf(p['moderate_h'])
        log_xstar,log_xtail = mp.log(p['moderate_Xstar']),mp.log(p['moderate_Xtail'])
        tail = TailMoments(h,log_xtail,order,p['panels'],p['taylor_degrees'][index])
        matrix = Compensation(p['moderate_lambda'],Rule(order,p['panels']))
        moderate_rows.append([root_record(tail,matrix,eta,log_xstar,4*mp.log(h)-14,p['graded_degree']) for eta in p['moderate_angles']])
        print('Outer stress and flat edge: order='+str(order),file=sys.stderr,flush=True)
        exterior = Exterior(p['exterior_h'],p['exterior_Xtail'],p['exterior_amplitude'],order,p['panels'])
        ext_rows.append([exterior.direct(y,eta) for eta in p['exterior_angles'] for y in p['exterior_y']])
        collars.append([exterior.flat(delta,eta) for eta in p['exterior_angles'] for delta in p['collar_delta']])
    mp.mp.dps = max(p['selected_digits'])
    selected_gap = row_refinement(*selected_rows)
    moderate_gap = row_refinement(*moderate_rows)
    exterior_gap = max(relative(a[k],b[k]) for a,b in zip(*ext_rows) for k in ('theta_over_rho','axial_over_rho','ratio'))
    exterior_gap = max(exterior_gap,max(relative(a[k],b[k]) for a,b in zip(*collars) for k in ('btheta','bz','ratio')))
    ode = [dict(max_step=s,**matrix.independent_moments(moderate_rows[-1][2]['root'],s)) for s in p['ode_steps']]
    kernel = []
    for z in ('.001','.0001'):
        print('Independent heat kernel: Z='+z,file=sys.stderr,flush=True)
        q = heat_polynomial(p['moderate_h'],z,24)
        independent = heat_integral(p['moderate_h'],z,48)
        kernel.append(dict(Z=z,normalized_delta=q['normalized_delta'],independent=independent,
            absolute_gap=abs(q['normalized_delta']-independent),relative_gap=relative(q['normalized_delta'],independent),
            normalized_remainder=q['remainder']/(mp.mpf(p['moderate_h'])*mp.mpf(z)),
            ode_residual=q['ode_residual'],ode_relative_gap=relative(q['ode_direct'],q['ode_residual'])))
    path_gaps = []
    for delta in ('.5','.25'):
        direct = exterior.direct(3-mp.mpf(delta),'.5')
        flat = exterior.flat(delta,'.5')
        path_gaps.append(relative(direct['ratio'],flat['ratio']))
    limit_gap = collar_limit_gap(collars[-1],p['exterior_angles'],len(p['collar_delta']))
    record = encode(dict(status=STATUS,protocol=p,provenance=provenance(),inherited=inherited(),certificate=bounds[-1],
        interval_gap=interval_gap,selected=dict(fixture=fixture,rows=selected_rows[-1],refinement_gap=selected_gap),
        moderate=dict(rows=moderate_rows[-1],refinement_gap=moderate_gap),independent_ODE=ode,
        pole=pole_control(fixture['h']),kernel=kernel,
        exterior=dict(rows=ext_rows[-1],collar=collars[-1],refinement_gap=exterior_gap,path_gaps=path_gaps,
            limit_gap=limit_gap,endpoints=[exterior.endpoint(eta) for eta in p['exterior_angles']])) )
    record['gates'] = {k:bool(v) for k,v in predicates(record).items()}
    return record


def validate(r):
    if r.get('status')!=STATUS: raise ValueError('Scientific scope changed')
    if r.get('protocol')!=json.loads((HERE/'protocol.json').read_text()): raise ValueError('Protocol changed')
    if r.get('provenance')!=provenance(): raise ValueError('Provenance changed')
    if r.get('inherited')!=inherited(): raise ValueError('Inherited evidence changed')
    if set(r.get('gates',{}))!=GATES or not all(v is True for v in r['gates'].values()):
        raise ValueError('Failed or missing gates: '+str([k for k,v in r.get('gates',{}).items() if v is not True]))
    with mp.workdps(500): actual = predicates(r)
    if set(actual)!=GATES or not all(v is True for v in actual.values()):
        raise ValueError('Data fail gates: '+str([k for k,v in actual.items() if v is not True]))


def compare(a,b,path=''):
    if isinstance(a,dict):
        if not isinstance(b,dict) or set(a)!=set(b): raise ValueError('Keys differ: '+path)
        for k in a: compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b): raise ValueError('Rows differ: '+path)
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+'/'+str(i))
    elif isinstance(a,bool):
        if a is not b: raise ValueError('Boolean differs: '+path)
    elif a!=b:
        try: x,y = mp.mpf(a),mp.mpf(b)
        except (ValueError,TypeError): raise ValueError('Value differs: '+path) from None
        if not mp.isfinite(x) or not mp.isfinite(y): raise ValueError('Nonfinite value: '+path)
        ode = 'independent_ODE' in path
        noise = any(k in path for k in ('gap','error','residual','stress_over_F')) and '/certificate/' not in path
        difference = abs(x-y)/(1+max(abs(x),abs(y))) if ode else (abs(x-y) if noise else relative(x,y))
        if difference>mp.mpf('1e-9' if ode else '1e-18'): raise ValueError('Numeric difference: '+path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify-record',type=Path)
    args = parser.parse_args()
    record = build()
    summary = dict(status=record['status'],gates=len(record['gates']),
        failed=[k for k,v in record['gates'].items() if not v],selected_gap=record['selected']['refinement_gap'],
        moderate_gap=record['moderate']['refinement_gap'],exterior_gap=record['exterior']['refinement_gap'])
    print(json.dumps(summary),flush=True)
    validate(record)
    if args.verify_record:
        old = json.loads(args.verify_record.read_text());validate(old);compare(old,record)
    args.output.write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':
    main()
