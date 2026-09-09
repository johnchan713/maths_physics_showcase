#!/usr/bin/env python3
"""Reproduce the intermediate continuum bounds and independent controls."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import mpmath as mp
from bounds import decay_bound, finite_cone_bound, axial, HERE, AXIAL
from reference import Reference, OUTER

STATUS = 'axial-and-intermediate-ratios-bounded-global-profile-unverified'
GATE_NAMES = {
    'all_elementary_inequalities','finite_parameter_order','axial_bound_inherited',
    'both_intermediate_ratios_bounded','interval_precision_refinement',
    'positive_Q_in_samples','Q_lower_bounds_in_samples','N_upper_bounds_in_samples',
    'sampled_stress_below_certificate','source_equation_crosscheck','independent_ramp_ODE',
    'integrated_energy_scale','quadrature_refinement','precision_refinement',
    'finite_h_pole_increment','omitted_h_control_detected','binary64_zero_control_detected',
    'future_pressure_error_bounded','old_lambda_failure_retained','omitted_energy_factor_detected',
    'finite_radius_cone_on_covered_stages'
}


def encode(v):
    if isinstance(v,mp.mpf): return mp.nstr(v,90)
    if isinstance(v,dict): return {k:encode(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)): return [encode(x) for x in v]
    return v


def relative(a,b):
    a,b=mp.mpf(a),mp.mpf(b)
    return abs(a-b)/max(abs(a),abs(b)) if a or b else mp.mpf(0)


def provenance():
    files=[HERE/name for name in ('bounds.py','reference.py','audit.py','test_audit.py','protocol.json','README.md','requirements.txt')]
    files += [AXIAL/name for name in ('bounds.py','diagnostics.py','protocol.json','evidence.json')]
    files += [OUTER/name for name in ('schedule.py','protocol.json')]
    return {str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def angle(r,name):
    return {'sqrt_delta':mp.sqrt(r.delta),'-sqrt_delta':-mp.sqrt(r.delta),
            'sqrt_lambda':mp.sqrt(r.lam)}.get(name,mp.mpf(name) if 'sqrt' not in name else None)


def coordinate(r,name):
    return {'T':r.T,'3T':3*r.T,'4T':4*r.T,'Tw/2':r.Tw/2,'Tw':r.Tw}.get(
        name,mp.mpf(name) if 'T' not in name else None)


def pole_control(r):
    start=r.state('ramp',1,1)
    end=r.state('power',r.Tw,1)
    increment=end['v']-start['v']
    source=r.source_power(r.Tw,1)['v_increment']
    dropped=r.state('power',r.Tw,1,h_override=0)['v']-r.state('ramp',1,1,h_override=0)['v']
    return dict(increment=increment,source_increment=source,
                increment_over_h_y=increment/(r.h*r.Tw),
                relative_gap=relative(increment,source),omitted_h_increment=dropped,
                binary64_increment=float(end['v'])-float(start['v']))


def pilot(p):
    with mp.workdps(p['diagnostic_digits'][-1]):
        r=Reference(p['diagnostic_orders'][0])
        rows,ode_rows,source_errors,energy_errors=[],[],[],[]
        q_ratios,n_ratios=[],[]
        # Reassemble log(E^2) by integrating the source slopes from E=P*f.
        first_log_change=2*r.rule.integrate(lambda t:mp.mpf('.6')*(1-axial_step(t))-mp.mpf('.5'))
        for stage,names in (('ramp',p['ramp_samples']),('power',p['power_samples'])):
            for name in names:
                y=coordinate(r,name)
                angles=[(a,angle(r,a)) for a in p['angles']]
                if stage=='power':
                    angles.append(('moving_small_angle',mp.sqrt(r.delta*mp.exp(-r.beta*y)+r.lam)))
                for label,e in angles:
                    row=r.state(stage,y,e)
                    f=r.ramp_factors(y) if stage=='ramp' else r.power_factors(y)
                    expected=mp.exp(2*(r.T+1)+first_log_change-r.T-f['x']-2*f['S'])/(1+e*e)**2
                    energy_errors.append(relative(row['E_squared'],expected))
                    lower=e*e/10+r.delta*mp.exp(-y)/20 if stage=='ramp' else e*e/10+r.lam/20+r.delta*mp.exp(-r.beta*y)/60
                    q_ratios.append(row['Q']/lower)
                    n_ratios.append(abs(row['v'])/(3*(r.T+f['x']+2)*mp.exp(2*f['S'])))
                    if stage=='power':
                        source=r.source_power(y,e)
                        source_errors.extend(relative(row[k],source[k]) for k in ('Q','v'))
                    rows.append(dict(angle_label=label,**row))
        ode_errors=[]
        for label in p['angles']:
            e=angle(r,label)
            for cadence in p['source_ode_max_steps']:
                for row in r.independent_ramp(e,[float(x) for x in p['ramp_samples']],cadence):
                    actual=r.state('ramp',row['x'],e)
                    gaps=[relative(actual['Q']/row['scale'],row['Q_over_scale']),
                          relative(actual['v'],row['v']),relative(r.ramp_factors(row['x'])['Z'],row['Z'])]
                    ode_errors.extend(gaps)
                    ode_rows.append(dict(angle_label=label,max_step=cadence,relative_gap=max(gaps),**row))
        pole=pole_control(r)
        reference_points=[r.state('power',3*r.T,mp.sqrt(r.lam)),r.state('power',0,mp.sqrt(r.delta))]
        fine=Reference(p['diagnostic_orders'][-1])
        fine_points=[fine.state('power',3*fine.T,mp.sqrt(fine.lam)),fine.state('power',0,mp.sqrt(fine.delta))]
        quadrature_gap=max(relative(a[k],b[k]) for a,b in zip(reference_points,fine_points) for k in ('Q','v','lambda_w_squared'))
        energy_control=reference_points[1]
        bad_energy=dict(omitted_factor='exp(-T)',energy_relative_ratio=1/r.delta,
            wrong_lambda_w_squared=energy_control['lambda_w_squared']/r.delta,
            corrected_lambda_w_squared=energy_control['lambda_w_squared'])
        outputs=dict(parameters=dict(Md=4,T=r.T,lambda_value=r.lam,h=r.h,Tw=r.Tw),
                     moment_samples=rows,ode_samples=ode_rows,pole=pole,
                     old_implementation_energy_control=bad_energy,
                     minimum_Q_ratio=min(q_ratios),maximum_N_ratio=max(n_ratios),
                     maximum_sampled_stress=max(x['lambda_w_squared'] for x in rows),
                     maximum_source_error=max(source_errors),maximum_ode_error=max(ode_errors),
                     maximum_energy_error=max(energy_errors),quadrature_gap=quadrature_gap,
                     maximum_tail_v_error=max(x['v_tail_error_upper'] for x in rows))
    with mp.workdps(p['diagnostic_digits'][0]):
        low=Reference(p['diagnostic_orders'][0])
        low_pole=pole_control(low)
        outputs['precision_gap']=relative(low_pole['increment_over_h_y'],pole['increment_over_h_y'])
    return outputs


def axial_step(t):
    # Imported after the local audit namespace is established, avoiding module-name collisions.
    from reference import step
    return step(t)


def old_failure():
    spec=importlib.util.spec_from_file_location('pinned_axial_diagnostics',AXIAL/'diagnostics.py')
    legacy=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(legacy)
    with mp.workdps(80):
        return legacy.Moments().intermediate_start('.5')


def build():
    mp.mp.dps=110
    p=json.loads((HERE/'protocol.json').read_text())
    limits={k:mp.mpf(v) for k,v in p['checks'].items()}
    certificates=[decay_bound(p['selected_Md'],d) for d in p['interval_digits']]
    cert=certificates[-1]
    gap=relative(certificates[0]['log_stress_upper_interval'][1],cert['log_stress_upper_interval'][1])
    preceding=axial.axial_bound(p['selected_Md'],1024,50)
    finite_cone=finite_cone_bound(preceding,p['interval_digits'][-1])
    data=pilot(p)
    old=old_failure()
    logs=cert['log_lambda_interval'],cert['log_h_interval']
    wide_lambda=mp.exp(-4*(mp.exp(p['selected_Md'])+10))
    float_control=dict(lambda_value=math.exp(float(logs[0][1])),h=math.exp(float(logs[1][1])),
                       a_minus_two=2+2*math.exp(float(logs[0][1]))-2,
                       wide_lambda_positive=wide_lambda>0,wide_h_positive=wide_lambda**2>0,
                       mpmath_110_digits_a_minus_two=2+2*wide_lambda-2)
    gates={
        'all_elementary_inequalities':all(cert['elementary_checks'].values()),
        'finite_parameter_order':logs[1][1]<logs[0][0]<mp.log(mp.mpf('.1')) and logs[1][1]<-2*cert['Td_interval'][1] and cert['Tw_interval'][0]>25,
        'axial_bound_inherited':preceding['status']=='axial-ratio-bound-passed' and all(axial.elementary_checks().values()),
        'finite_radius_cone_on_covered_stages':finite_cone['all_pass'],
        'both_intermediate_ratios_bounded':cert['status']=='intermediate-ratios-bounded' and cert['universal_lambda_w_squared_upper']<mp.mpf(p['claimed_lambda_w_squared_upper']),
        'interval_precision_refinement':gap<limits['interval_precision_relative_gap'],
        'positive_Q_in_samples':all(x['Q']>0 for x in data['moment_samples']),
        'Q_lower_bounds_in_samples':data['minimum_Q_ratio']>=1,
        'N_upper_bounds_in_samples':data['maximum_N_ratio']<=1,
        'sampled_stress_below_certificate':data['maximum_sampled_stress']<mp.mpf(p['claimed_lambda_w_squared_upper']),
        'source_equation_crosscheck':data['maximum_source_error']<limits['moment_source_relative_gap'],
        'independent_ramp_ODE':data['maximum_ode_error']<limits['independent_ode_relative_gap'],
        'integrated_energy_scale':data['maximum_energy_error']<limits['energy_scale_relative_gap'],
        'quadrature_refinement':data['quadrature_gap']<limits['quadrature_relative_gap'],
        'precision_refinement':data['precision_gap']<limits['precision_relative_gap'],
        'finite_h_pole_increment':data['pole']['increment']<0 and abs(data['pole']['increment_over_h_y']+2)<mp.mpf('1e-30') and data['pole']['relative_gap']<limits['pole_increment_relative_gap'],
        'omitted_h_control_detected':data['pole']['omitted_h_increment']==0 and data['pole']['binary64_increment']==0,
        'binary64_zero_control_detected':all(float_control[k]==0 for k in ('lambda_value','h','a_minus_two','mpmath_110_digits_a_minus_two')) and float_control['wide_lambda_positive'] and float_control['wide_h_positive'],
        'future_pressure_error_bounded':0<data['maximum_tail_v_error']<limits['maximum_pilot_tail_error'],
        'old_lambda_failure_retained':old['Q']>0 and old['lambda_w_squared']>1,
        'omitted_energy_factor_detected':data['old_implementation_energy_control']['wrong_lambda_w_squared']>mp.mpf(p['claimed_lambda_w_squared_upper']) and data['old_implementation_energy_control']['corrected_lambda_w_squared']<mp.mpf(p['claimed_lambda_w_squared_upper']),
    }
    return encode(dict(status=STATUS,protocol=p,provenance=provenance(),certificate=cert,
        interval_precision_gap=gap,preceding_axial_bound=preceding,finite_cone=finite_cone,pilot=data,
        binary64_control=float_control,old_lambda_failure=old,gates=gates))


def validate(record):
    if record.get('status')!=STATUS: raise ValueError('Scientific scope changed')
    if record.get('protocol')!=json.loads((HERE/'protocol.json').read_text()): raise ValueError('Protocol changed')
    if record.get('provenance')!=provenance(): raise ValueError('Provenance changed')
    gates=record.get('gates',{})
    if set(gates)!=GATE_NAMES or not all(v is True for v in gates.values()):
        raise ValueError('Missing or failed audit gates: '+str(sorted(GATE_NAMES-set(gates)))+' '+str([k for k,v in gates.items() if v is not True]))
    if mp.mpf(record['old_lambda_failure']['lambda_w_squared'])<=1: raise ValueError('Old failure erased')
    if mp.mpf(record['pilot']['pole']['increment'])>=0: raise ValueError('Finite-h increment erased')
    if record['certificate']['claimed_lambda_w_squared_upper']!='1e-42': raise ValueError('Claimed bound changed')


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
        try: x,y=mp.mpf(a),mp.mpf(b)
        except (ValueError,TypeError): raise ValueError('Value differs: '+path) from None
        # Tiny increments use a relative comparison, without a unit floor.
        tolerance=mp.mpf('1e-8') if 'ode' in path else mp.mpf('1e-18')
        if not mp.isfinite(x) or not mp.isfinite(y): raise ValueError('Nonfinite record: '+path)
        if 'error' in path or 'gap' in path:
            changed=abs(x-y)>tolerance
        else: changed=relative(x,y)>tolerance
        if changed: raise ValueError('Numeric record differs: '+path)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify-record',type=Path)
    args=parser.parse_args(argv)
    result=build()
    validate(result)
    if args.verify_record:
        old=json.loads(args.verify_record.read_text())
        validate(old);compare_record(old,result)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],bound=result['certificate']['claimed_lambda_w_squared_upper'],
        ode_error=result['pilot']['maximum_ode_error'],pole_increment=result['pilot']['pole']['increment_over_h_y'],gates=len(result['gates']))))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
