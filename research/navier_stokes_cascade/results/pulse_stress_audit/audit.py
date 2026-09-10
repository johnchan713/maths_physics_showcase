#!/usr/bin/env python3
"""Reproduce the pulse continuum reduction and independent finite controls."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import mpmath as mp
from bounds import pulse_bound, HERE, PULSE, AXIAL
from diagnostics import (SyntheticPulse, algebra_controls, affine_control,
                         finite_cone_controls, angular_equality_control)

STATUS = 'corrected-reference-pulse-cone-bounded-later-stages-and-global-profile-unverified'
GATE_NAMES = {
    'inherited_moment_argument','inherited_parameter_schedule','all_outward_constants',
    'interval_precision_refinement','uniform_positive_Q','all_finite_errors_retained',
    'uniform_w_remainder','first_A24_margin','second_A24_margin','finite_radius_cone',
    'strict_positive_v_excess','source_moment_algebra','finite_cone_algebra',
    'independent_source_ODE','ODE_step_refinement','finite_h_control',
    'angular_derivative_control','energy_control','end_bump_shear_control',
    'insufficient_precision_failure','high_precision_cancellation_check',
    'first_convolution_derivative_required','binary64_excess_failure',
    'actual_family_underflow_control','ratio_only_radius_failure',
    'old_Md_failure_retained','old_lambda_failure_retained',
    'exact_angular_equality_control',
}


def encode(value):
    if isinstance(value,mp.mpf): return mp.nstr(value,105)
    if isinstance(value,dict): return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return [encode(v) for v in value]
    if hasattr(value,'item'): return value.item()
    return value


def relative(a,b):
    a,b = mp.mpf(a),mp.mpf(b)
    return abs(a-b)/max(abs(a),abs(b)) if a or b else mp.mpf(0)


def provenance():
    paths = [HERE/name for name in ('README.md','protocol.json','requirements.txt',
        'bounds.py','diagnostics.py','audit.py','test_audit.py')]
    paths += [PULSE/name for name in ('README.md','protocol.json','bounds.py','pulse.py','evidence.json')]
    paths += [AXIAL/name for name in ('README.md','bounds.py','evidence.json')]
    paths += [HERE.parent/'intermediate_decay_audit'/name for name in ('README.md','bounds.py','evidence.json')]
    paths += [HERE.parent/'outer_pressure_pilot'/'schedule.py']
    return {str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def inherited():
    moment = json.loads((PULSE/'evidence.json').read_text())
    axial = json.loads((AXIAL/'evidence.json').read_text())
    return dict(moment_status=moment['status'],moment_certificate=moment['certificate'],
        moment_gates=moment['gates'],old_Md_failure=axial['old_Md_failure'],
        old_lambda_failure=moment['old_lambda_failure'])


def predicates(r):
    """Recompute gates from recorded measurements, not just saved booleans."""
    n = mp.mpf
    p,b = r['protocol'],r['certificate']
    limits = {k:n(v) for k,v in p['checks'].items()}
    old = r['inherited']
    m = old['moment_certificate']
    algebra,cone = r['algebra'],r['cone_controls']
    rows = [row for run in r['ode_runs'] for row in run['rows']]
    low = [x for x in r['affine_controls'] if x['digits']==p['cancellation_digits'][0]]
    high = [x for x in r['affine_controls'] if x['digits']>p['cancellation_digits'][0]]
    small = cone['small_radius']
    return {
        'inherited_moment_argument': old['moment_status']=='reference-moments-closed-pulse-cone-and-global-profile-unverified'
            and len(old['moment_gates'])==34 and all(v is True for v in old['moment_gates'].values())
            and all(v is True for v in m['checks'].values())
            and m['amplitude_lower']=='1.0100502' and m['amplitude_upper']=='1.0100504'
            and m['amplitude_eta_derivative_bound']=='3000 T lambda'
            and m['axial_c_C1_bound']=='lambda^20 at fixed Amp; 2 lambda^20 after substituting the root',
        'inherited_parameter_schedule': m['Md']==b['Md']==p['selected_Md']==64 and m['Tf']==1000
            and m['co']=='0.001' and n(b['log_h_interval'][1])<n(b['log_lambda_interval'][0])<-256,
        'all_outward_constants': b['status']=='corrected-reference-pulse-cone-bounded'
            and all(v is True for v in b['checks'].values()),
        'interval_precision_refinement': n(r['interval_relative_gap'])<limits['interval_relative_gap'],
        'uniform_positive_Q': b['Q_lower']=='(lambda+eta^2)/4' and 0<n(b['Q_relative_error_universal_upper'])<n('.25'),
        'all_finite_errors_retained': all(n(v)>0 for v in b['w_remainder_parts'].values())
            and n(b['bs_remainder_universal_upper'])>0,
        'uniform_w_remainder': 0<n(b['w_remainder_universal_upper'])<n(p['claimed_w_remainder_upper']),
        'first_A24_margin': n(b['maximum_bsw_upper'])<n(p['claimed_bsw_upper'])
            and n(b['first_A24_margin_lower'])>n('1.48'),
        'second_A24_margin': n(b['maximum_second_expression_upper'])<n(p['claimed_second_expression_upper'])
            and n(b['second_A24_margin_lower'])>n('.82'),
        'finite_radius_cone': b['XR_sufficient_lower']==p['XR_sufficient_lower']==100
            and b['finite_ps1_sufficient_lower']==20000000
            and n(b['universal_log_ps1_lower'])>mp.log(b['finite_ps1_sufficient_lower'])
            and n(b['finite_normalized_quadratic_gap_lower'])>n('.55'),
        'strict_positive_v_excess': n(b['log_v_minus_two_lower_interval'][0])>n(b['log_lambda_interval'][0])
            and all(mp.isfinite(n(x)) for x in b['log_v_minus_two_lower_interval']),
        'source_moment_algebra': n(algebra['maximum_source_identity_error'])<limits['algebra_absolute_gap'],
        'finite_cone_algebra': n(cone['identity_maximum_error'])<limits['algebra_absolute_gap'],
        'independent_source_ODE': len(r['ode_runs'])==2*len(p['synthetic_angles']) and bool(rows)
            and all(n(x['normalized_gap'])<limits['source_moment_normalized_gap'] for x in rows),
        'ODE_step_refinement': n(r['ode_refinement_normalized_gap'])<limits['ode_refinement_normalized_gap'],
        'finite_h_control': n(algebra['omitted_h_maximum_gap'])>limits['omitted_term_detectable'],
        'angular_derivative_control': n(algebra['omitted_eta_maximum_gap'])>limits['omitted_term_detectable'],
        'energy_control': n(algebra['omitted_energy_maximum_gap'])>limits['omitted_term_detectable'],
        'end_bump_shear_control': max(n(x['omitted_bump_bs_gap']) for x in rows)>limits['omitted_term_detectable'],
        'insufficient_precision_failure': len(low)==4 and all(n(x['direct_relative_error'])>n('.9') for x in low),
        'high_precision_cancellation_check': len(high)==8 and all(n(x['direct_relative_error'])<limits['high_precision_relative_gap'] for x in high),
        'first_convolution_derivative_required': n(low[0]['omitted_derivative_absolute_error'])>2,
        'binary64_excess_failure': all(n(x['a_minus_two'])>0 and x['binary64_rounded_a_minus_two']==0 for x in low+high)
            and all(n(x['rounded_a_minus_two'])==0 for x in low),
        'actual_family_underflow_control': r['binary64_actual_lambda']==0 and n(b['log_lambda_interval'][1])<0,
        'ratio_only_radius_failure': small['passes_ratio_test'] is True and small['passes_finite_cone'] is False
            and n(small['first_margin'])>0 and n(small['second_margin'])>0 and n(small['Pc'])<n(small['v']),
        'old_Md_failure_retained': n(old['old_Md_failure']['Pc_over_ps1'])<0,
        'old_lambda_failure_retained': n(old['old_lambda_failure']['lambda_w_squared'])>1,
        'exact_angular_equality_control': r['angular_equality']['exact_rational_checks'] is True
            and r['angular_equality']['floating_comparison_passed'] is False
            and 0<n(r['angular_equality']['relative_rounding_excess'])<n('1e-110'),
    }


def build():
    mp.mp.dps = 120
    p = json.loads((HERE/'protocol.json').read_text())
    certificates = [pulse_bound(p['selected_Md'],digits,p['XR_sufficient_lower']) for digits in p['interval_digits']]
    b = certificates[-1]
    fields = ('w_remainder_universal_upper','maximum_bsw_upper','maximum_second_expression_upper',
              'Q_relative_error_universal_upper','bs_remainder_universal_upper')
    interval_gap = max(relative(certificates[0][k],b[k]) for k in fields)
    ode_runs,refinement = [],[]
    for eta in p['synthetic_angles']:
        print('Independent pulse source check: eta='+eta,file=sys.stderr,flush=True)
        pair = []
        for bulk,patch in zip(p['synthetic_bulk_steps'],p['synthetic_patch_steps']):
            rows = SyntheticPulse(eta,float(p['synthetic_lambda'])).integrate(bulk,float(patch))
            ode_runs.append(dict(eta=eta,bulk_step=bulk,patch_step=patch,rows=rows))
            pair.append(rows)
        for x,y in zip(*pair):
            refinement.extend(abs(x[k]-y[k])/(1+abs(y[k])) for k in ('Q_source','N_over_E_source'))
    record = encode(dict(status=STATUS,protocol=p,provenance=provenance(),
        certificate=b,interval_relative_gap=interval_gap,inherited=inherited(),
        algebra=algebra_controls(),cone_controls=finite_cone_controls(),
        angular_equality=angular_equality_control(),ode_runs=ode_runs,
        ode_refinement_normalized_gap=max(refinement),
        affine_controls=[affine_control(digits,angle) for digits in p['cancellation_digits'] for angle in ('0','.5','2','10')],
        binary64_actual_lambda=math.exp(float(b['log_lambda_interval'][1]))))
    record['gates'] = {k:bool(v) for k,v in predicates(record).items()}
    return record


def validate(record):
    if record.get('status')!=STATUS: raise ValueError('Scientific scope changed')
    if record.get('protocol')!=json.loads((HERE/'protocol.json').read_text()): raise ValueError('Protocol changed')
    if record.get('provenance')!=provenance(): raise ValueError('Provenance changed')
    if record.get('inherited')!=inherited(): raise ValueError('Inherited evidence changed')
    gates = record.get('gates',{})
    if set(gates)!=GATE_NAMES or not all(v is True for v in gates.values()):
        raise ValueError('Missing or failed audit gates')
    with mp.workdps(120):
        recomputed = predicates(record)
    if set(recomputed)!=GATE_NAMES or not all(v is True for v in recomputed.values()):
        raise ValueError('Recorded measurements fail scientific predicates: '+str([k for k,v in recomputed.items() if v is not True]))


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
        noise = '/certificate/' not in path and any(k in path for k in ('error','gap','residual'))
        tolerance = mp.mpf('1e-9') if '/ode_runs/' in path or 'ode_refinement' in path else mp.mpf('1e-18')
        # Manufactured ODE observations use their declared normalized metric.
        # All physical/certificate tiny quantities retain relative comparison.
        difference = abs(x-y)/(1+max(abs(x),abs(y))) if '/ode_runs/' in path else (abs(x-y) if noise else relative(x,y))
        if difference>tolerance:
            raise ValueError('Numeric record differs: '+path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify-record',type=Path)
    args = parser.parse_args(argv)
    record = build()
    validate(record)
    if args.verify_record:
        old = json.loads(args.verify_record.read_text())
        validate(old)
        compare_record(old,record)
    args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(status=record['status'],gates=len(record['gates']),
        first_margin=record['certificate']['first_A24_margin_lower'],
        second_margin=record['certificate']['second_A24_margin_lower'],
        maximum_ode_error=max(x['normalized_gap'] for run in record['ode_runs'] for x in run['rows']))))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
