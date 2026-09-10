#!/usr/bin/env python3
"""Reproduce finite post-pulse bounds, patch pressure and terminal failures."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import mpmath as mp
from bounds import post_bound, HERE, PULSE, MOMENT, OUTER
from diagnostics import fixture, Patch, terminal_control, independent_terminal, exterior_control, Rule, relative

STATUS = 'post-pulse-reference-cone-bounded-heat-collar-and-global-profile-unverified'
GATE_NAMES = {
    'inherited_exact_moments','inherited_pulse_floor','outward_constants',
    'declared_finite_endpoint','interval_precision_refinement','positive_Q_floor',
    'strict_shear_excess','finite_w_bound','finite_radius_cone','actual_angular_edits',
    'actual_patch_moment_rows','differentiated_patch_rows','partial_pressure_not_zero',
    'complete_pressure_restored','positive_quadratic_omission_bound',
    'independent_patch_pressure_ODE','patch_quadrature_refinement',
    'terminal_positive_until_endpoint','terminal_half_floor','independent_terminal_ODE',
    'terminal_refinement','finite_h_difference_retained','binary64_difference_failure',
    'independent_exterior_ramps','positive_wait_and_hold_factors',
    'pure_power_cancellation','omitted_h_energy_failure','infinite_tail_cone_failure',
    'historical_failures_retained'
}


def encode(v):
    if isinstance(v,mp.mpf): return mp.nstr(v,105)
    if isinstance(v,dict): return {k:encode(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)): return [encode(x) for x in v]
    if hasattr(v,'item'): return v.item()
    return v


def provenance():
    paths = [HERE/n for n in ('README.md','protocol.json','requirements.txt','bounds.py','diagnostics.py','audit.py','test_audit.py')]
    paths += [p/n for p in (PULSE,MOMENT) for n in ('README.md','bounds.py','protocol.json','evidence.json')]
    paths += [OUTER/'schedule.py',HERE.parent/'axial_stress_audit'/'bounds.py']
    return {str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def inherited():
    pulse = json.loads((PULSE/'evidence.json').read_text())
    moment = json.loads((MOMENT/'evidence.json').read_text())
    return dict(moment_status=moment['status'],moment_gates=moment['gates'],
        moment_certificate=moment['certificate'],pulse_status=pulse['status'],pulse_gates=pulse['gates'],
        pulse_Q_lower=pulse['certificate']['Q_lower'],old_Md_failure=pulse['inherited']['old_Md_failure'],
        old_lambda_failure=pulse['inherited']['old_lambda_failure'])


def predicates(r):
    n = mp.mpf
    p,b,old = r['protocol'],r['certificate'],r['inherited']
    limits = {k:n(v) for k,v in p['checks'].items()}
    patches = r['patches']
    terminal = r['terminal_controls']
    ext = r['exterior']
    return {
        'inherited_exact_moments': old['moment_status']=='reference-moments-closed-pulse-cone-and-global-profile-unverified'
            and len(old['moment_gates'])==34 and all(v is True for v in old['moment_gates'].values())
            and old['moment_certificate']['Md']==64 and old['moment_certificate']['Tf']==1000,
        'inherited_pulse_floor': old['pulse_status']=='corrected-reference-pulse-cone-bounded-later-stages-and-global-profile-unverified'
            and len(old['pulse_gates'])==28 and all(v is True for v in old['pulse_gates'].values())
            and old['pulse_Q_lower']=='(lambda+eta^2)/4',
        'outward_constants': b['status']=='post-pulse-reference-cone-bounded-through-terminal-half'
            and all(v is True for v in b['checks'].values()),
        'declared_finite_endpoint': p['terminal_cone_endpoint']==b['terminal_endpoint']=='0.5',
        'interval_precision_refinement': n(r['interval_relative_gap'])<limits['interval_relative_gap'],
        'positive_Q_floor': b['Q_floor']=='h/2000' and n(b['terminal_Q_over_h_lower'])>n(1)/2000,
        'strict_shear_excess': n(b['log_a_minus_two_lower_interval'][0])>n(b['log_h_interval'][0]),
        'finite_w_bound': 0<n(b['w_universal_upper'])<n(p['claimed_w_upper'])==n('1e-440')
            and 0<n(b['second_expression_universal_upper'])<n('.01'),
        'finite_radius_cone': b['XR_sufficient_lower']==100 and b['finite_ps1_lower']==100
            and n(b['universal_log_ps1_lower'])>mp.log(100) and n(b['finite_normalized_cone_gap_lower'])>n('1.84'),
        'actual_angular_edits': 0<n(b['angular_relative_edit_upper'])<n('.01')
            and 0<n(b['angular_eta_log_edit_upper'])<1 and n(b['angular_slope_edit_upper'])>0,
        'actual_patch_moment_rows': len(patches)==len(p['patch_angles'])
            and all(max(n(v) for v in x['moment_relative_errors'])<limits['patch_moment_relative_error'] for x in patches),
        'differentiated_patch_rows': all(n(x['derivative_row_relative_error'])<limits['patch_moment_relative_error'] for x in patches),
        'partial_pressure_not_zero': all(abs(n(x['mid_patch_pressure_over_cscale']))>n('.001') for x in patches)
            and all(abs(n(x['gap_pressure_over_cscale']))>n('.001') for x in patches),
        'complete_pressure_restored': all(abs(n(x['complete_pressure_relative_residual']))<limits['patch_moment_relative_error']
            and n(x['post_patch_pressure'])==0 for x in patches),
        'positive_quadratic_omission_bound': all(0<n(x['quadratic_normalized_positive_upper'])<limits['independent_ode_normalized_gap'] for x in patches),
        'independent_patch_pressure_ODE': all(n(row['normalized_gap'])<limits['independent_ode_normalized_gap']
            for x in patches for run in x['independent_runs'] for row in run['rows']),
        'patch_quadrature_refinement': n(r['patch_refinement_relative_gap'])<limits['patch_quadrature_relative_gap'],
        'terminal_positive_until_endpoint': all(all(n(x['Q_over_h'])>0 for x in t['terminal_rows'][:-1])
            and n(t['terminal_rows'][-1]['Q_over_h'])==0 for t in terminal),
        'terminal_half_floor': all(n(t['terminal_rows'][1]['Q_over_h'])>n(1)/2000 for t in terminal),
        'independent_terminal_ODE': all(n(x['normalized_gap'])<limits['independent_ode_normalized_gap'] for run in r['terminal_ODE'] for x in run['rows']),
        'terminal_refinement': n(r['terminal_refinement_relative_gap'])<limits['terminal_refinement_relative_gap'],
        'finite_h_difference_retained': all(n(t['stable_finite_h_difference_over_h'])<0
            and n(t['direct_relative_gap'])<limits['finite_h_difference_relative_gap'] for t in terminal),
        'binary64_difference_failure': all(t['binary64_Qp_over_h_difference']==0 for t in terminal),
        'independent_exterior_ramps': n(ext['maximum_ODE_relative_gap'])<limits['independent_ode_normalized_gap'],
        'positive_wait_and_hold_factors': n(ext['wait'])>0 and n(ext['q_before_wait'])>1 and 0<n(ext['Qp'])<1
            and abs(n(ext['wait_log_residual']))<n('1e-100')
            and all(0<n(ext[a]) and relative(n(ext[a]),n(ext[b]))<n('1e-100') for a,b in (
                ('E_hold_factor','expected_E_hold_factor'),('XE2_hold_factor','expected_XE2_hold_factor'))),
        'pure_power_cancellation': all(abs(n(t['pure_power_N_residual']))<n('1e-100')
            and n(t['pure_power_energy_term'])>n('.9') and n(t['pure_power_pressure_term'])<n('-.9') for t in terminal),
        'omitted_h_energy_failure': all(abs(n(t['omitted_h_energy_N_over_E2']))>n('.9') for t in terminal),
        'infinite_tail_cone_failure': all(t['tail_Q']==0 and t['tail_Pc']==0 and t['tail_cone_pass'] is False for t in terminal),
        'historical_failures_retained': n(old['old_Md_failure']['Pc_over_ps1'])<0 and n(old['old_lambda_failure']['lambda_w_squared'])>1,
    }


def build():
    mp.mp.dps = 280
    p = json.loads((HERE/'protocol.json').read_text())
    mp.mp.dps = p['patch_fixture_digits']
    bounds = [post_bound(p['selected_Md'],d,p['XR_sufficient_lower'],p['terminal_cone_endpoint']) for d in p['interval_digits']]
    interval_gap = max(relative(bounds[0][k],bounds[1][k]) for k in ('w_universal_upper','angular_relative_edit_upper','finite_normalized_cone_gap_lower'))
    old_rows,rows = [],[]
    for order in p['patch_orders']:
        s = fixture(order,p['patch_fixture_lambda'],p['patch_integration_panels'],p['patch_moment_panels'])
        for angle in p['patch_angles']:
            print('Actual angular patch: order='+str(order)+', eta='+angle,file=sys.stderr,flush=True)
            row = Patch(s,angle).record(p['ode_steps'])
            (old_rows if order==p['patch_orders'][0] else rows).append(row)
    gaps = [relative(a[k],b[k]) for a,b in zip(old_rows,rows)
        for k in ('coefficient_scale','mid_patch_pressure_over_cscale','gap_pressure_over_cscale','mid_patch_pressure_eta')]
    for a,b in zip(old_rows,rows):
        gaps.extend(relative(x,y) for x,y in zip(a['coefficient_eta'],b['coefficient_eta']))
    terminals = [terminal_control(d,o,p['terminal_panels']) for d,o in zip(p['terminal_digits'],p['terminal_orders'])]
    terminal_gap = max([relative(terminals[0][k],terminals[1][k]) for k in ('Qp_over_h','stable_finite_h_difference_over_h')]
        +[relative(x['Q_over_h'],y['Q_over_h']) for x,y in zip(terminals[0]['terminal_rows'],terminals[1]['terminal_rows'])])
    with mp.workdps(320):
        r = Rule(p['terminal_orders'][-1],p['terminal_panels'])
        terminal_ODE = [dict(max_step=c,rows=independent_terminal(terminals[-1]['h'],r,c)) for c in p['ode_steps']]
    record = encode(dict(status=STATUS,protocol=p,provenance=provenance(),inherited=inherited(),certificate=bounds[-1],
        interval_relative_gap=interval_gap,patches=rows,patch_refinement_relative_gap=max(gaps),
        terminal_controls=terminals,terminal_refinement_relative_gap=terminal_gap,
        terminal_ODE=terminal_ODE,exterior=exterior_control(s)))
    record['gates'] = {k:bool(v) for k,v in predicates(record).items()}
    return record


def validate(record):
    if record.get('status')!=STATUS: raise ValueError('Scientific scope changed')
    if record.get('protocol')!=json.loads((HERE/'protocol.json').read_text()): raise ValueError('Protocol changed')
    if record.get('provenance')!=provenance(): raise ValueError('Provenance changed')
    if record.get('inherited')!=inherited(): raise ValueError('Inherited evidence changed')
    if set(record.get('gates',{}))!=GATE_NAMES or not all(v is True for v in record['gates'].values()):
        raise ValueError('Missing or failed gates: '+str([k for k,v in record.get('gates',{}).items() if v is not True]))
    with mp.workdps(280): actual = predicates(record)
    if set(actual)!=GATE_NAMES or not all(v is True for v in actual.values()):
        raise ValueError('Recorded data fail predicates: '+str([k for k,v in actual.items() if v is not True]))


def compare_record(a,b,path=''):
    if isinstance(a,dict):
        if not isinstance(b,dict) or set(a)!=set(b): raise ValueError('Keys differ: '+path)
        for k in a: compare_record(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b): raise ValueError('Lengths differ: '+path)
        for i,(x,y) in enumerate(zip(a,b)): compare_record(x,y,path+'/'+str(i))
    elif isinstance(a,bool) or a is None:
        if a is not b: raise ValueError('Value differs: '+path)
    elif a!=b:
        try: x,y = mp.mpf(a),mp.mpf(b)
        except (ValueError,TypeError): raise ValueError('Value differs: '+path) from None
        if not mp.isfinite(x) or not mp.isfinite(y): raise ValueError('Nonfinite value: '+path)
        ode = 'ODE' in path or 'independent_runs' in path
        noise = '/certificate/' not in path and any(k in path for k in ('error','gap','residual'))
        difference = abs(x-y)/(1+max(abs(x),abs(y))) if ode else (abs(x-y) if noise else relative(x,y))
        if difference>mp.mpf('1e-9' if ode else '1e-18'): raise ValueError('Numeric record differs: '+path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify-record',type=Path)
    args = parser.parse_args(argv)
    record = build()
    print(json.dumps(dict(patch_refinement=record['patch_refinement_relative_gap'],
        terminal_refinement=record['terminal_refinement_relative_gap'],
        failed_gates=[k for k,v in record['gates'].items() if v is not True])),file=sys.stderr,flush=True)
    validate(record)
    if args.verify_record:
        old = json.loads(args.verify_record.read_text());validate(old);compare_record(old,record)
    args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(status=record['status'],gates=len(record['gates']),w_upper=record['certificate']['w_universal_upper'],
        patch_refinement=record['patch_refinement_relative_gap'],terminal_refinement=record['terminal_refinement_relative_gap'])))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
