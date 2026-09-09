#!/usr/bin/env python3
"""Reproduce finite pulse moment bounds, numerical fixtures and failures."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import mpmath as mp
from bounds import closure_bound, HERE, AXIAL
from pulse import Prefix, Continuation, MomentSolver, family_increment, pulse_constant, relative, OUTER, incoming_integrals, Rule

STATUS = 'reference-moments-closed-pulse-cone-and-global-profile-unverified'
GATE_NAMES = {
    'all_interval_constants','nested_shape_enclosures','finite_parameter_order',
    'amplitude_bracket','unique_root_derivative','uniform_finite_remainder',
    'angular_correction_contraction','axial_inverse_loss_bounded',
    'actual_fixture_M_J_rows','fixture_S_quadratic','independent_end_patch_ODE',
    'main_moment_quadrature_refinement','coefficient_quadrature_refinement',
    'fixture_amplitude_refinement','actual_angular_rows','positive_full_tail_mass',
    'nonzero_bump_energy_retained','nonzero_affine_prefix_retained',
    'main_prefix_omission_bounded','only_M_correction_detected',
    'missing_end_corrections_detected','leading_root_not_exact',
    'large_lambda_bracket_failure','binary64_singular_matrix_detected',
    'family_positive_amplitude_increment','family_increment_refinement',
    'family_direct_increment_crosscheck','family_omitted_bump_bound_positive',
    'binary64_amplitude_increment_lost','family_post_swirl_retained',
    'old_parameter_failure_retained','K_diagnostics_inside_enclosure',
    'incoming_moment_refinement','coarse_incoming_resolution_failure'
}


def encode(value):
    if isinstance(value,mp.mpf): return mp.nstr(value,90)
    if isinstance(value,dict): return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [encode(v) for v in value]
    return value


def provenance():
    paths=[HERE/name for name in ('bounds.py','pulse.py','audit.py','test_audit.py','README.md','protocol.json','requirements.txt')]
    paths += [OUTER/'schedule.py', AXIAL/'bounds.py',
              HERE.parent/'intermediate_decay_audit'/'protocol.json',
              HERE.parent/'intermediate_decay_audit'/'evidence.json']
    return {str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def fixture(lam, order, digits, angles, ode_steps=(), incoming_panels=16):
    with mp.workdps(digits):
        prefix=Prefix(lam,order,incoming_panels=incoming_panels)
        solver=MomentSolver(prefix,order)
        continuation=Continuation(prefix,order)
        rows=[solver.root(eta,continuation) for eta in angles]
        ode=[dict(eta=mp.mpf(row['eta']),**solver.source_check(row['eta'],row['amplitude'],step))
             for row in rows for step in ode_steps]
        return dict(lambda_value=prefix.lam,Md=4,h=continuation.schedule.h,
            digits=digits,order=order,T=prefix.T,K=solver.K,
            incoming_moments=dict(Ka=prefix.Ka,K2a=prefix.K2a,Ja=prefix.Ja,panels=incoming_panels),
            main_moments=solver.main,coefficient_log_scale=solver.log_scale,
            normalized_matrix=solver.matrix,rows=rows,ode_checks=ode)


def compare_fixtures(coarse, fine):
    log_gap=max(abs(a['log_value']-b['log_value']) for a,b in zip(coarse['main_moments'],fine['main_moments']))
    rescale=mp.exp(coarse['coefficient_log_scale']-fine['coefficient_log_scale'])
    coefficient_gap=max(relative(rescale*a,b)
        for ra,rb in zip(coarse['rows'],fine['rows'])
        for key in ('scaled_coefficients','affine_a_scaled','affine_b_scaled')
        for a,b in zip(ra[key],rb[key]))
    amplitude_gap=max(relative(a['amplitude'],b['amplitude']) for a,b in zip(coarse['rows'],fine['rows']))
    incoming_gap=max(relative(coarse['incoming_moments'][k],fine['incoming_moments'][k]) for k in ('Ka','K2a','Ja'))
    return dict(main_log_gap=log_gap,coefficient_relative_gap=coefficient_gap,
                amplitude_relative_gap=amplitude_gap,incoming_relative_gap=incoming_gap)


def numeric_controls(p):
    lam=mp.mpf('1e-20')
    f=float(lam)
    naive=math.exp(2*(.5-2*f))-math.exp(2*(.5-f))
    stable=mp.exp(1-2*lam)*mp.expm1(-2*lam)
    # The larger lambda is a deliberate failure; no passing tolerance is changed.
    prefix=Prefix(p['large_lambda_control'],p['fixture_orders'][0])
    K=pulse_constant(prefix.rule)
    S=prefix.moments(0)['s0']
    C=-mp.expm1(-26)/4
    # Post-swirl energy has a nonpositive sign in F. End-bump energy has a
    # positive 40 lambda^41 majorant here (the absorption inequality holds).
    F_upper=mp.mpf('1.2')**2*K-C+prefix.lam*S+40*prefix.lam**41
    resolutions=[]
    for panels in p['failed_incoming_panels']:
        a,b=[incoming_integrals(4,Rule(order,panels)) for order in p['fixture_orders']]
        resolutions.append(dict(panels=panels,Ka_relative_gap=relative(a[0],b[0]),
                                K2a_relative_gap=relative(a[1],b[1])))
    return dict(incoming_resolution_failures=resolutions,determinant=dict(lambda_value=lam,binary64_difference=naive,
        stable_difference=stable,stable_difference_over_lambda=stable/lam),
        large_lambda=dict(lambda_value=prefix.lam,prefix_s0=S,
            prefix_only_amplitude=mp.sqrt((C-prefix.lam*S)/K),
            F_at_1p2_upper_diagnostic=F_upper,
            scope='Numerical bracket-failure control, not an interval proof or globally admissible fixture'))


def build():
    mp.mp.dps=110
    p=json.loads((HERE/'protocol.json').read_text())
    limits={k:mp.mpf(v) for k,v in p['checks'].items()}
    print('Checking continuum interval bounds',file=sys.stderr,flush=True)
    certificates=[closure_bound(n,d,p['selected_Md']) for n,d in zip(p['shape_boxes'],p['interval_digits'])]
    cert=certificates[-1]
    fixtures=[]
    for lam in p['fixture_lambdas']:
        print('Checking pulse fixture lambda='+lam+' at both quadrature orders',file=sys.stderr,flush=True)
        coarse=fixture(lam,p['fixture_orders'][0],p['fixture_digits'][0],p['fixture_angles'],incoming_panels=p['incoming_moment_panels'])
        fine=fixture(lam,p['fixture_orders'][-1],p['fixture_digits'][-1],p['fixture_angles'],p['ode_max_steps'],p['incoming_moment_panels'])
        fine['refinement']=compare_fixtures(coarse,fine)
        fixtures.append(fine)
    family=[]
    for digits,order in zip(p['family_diagnostic_digits'],p['family_diagnostic_orders']):
        print('Checking finite family increment at '+str(digits)+' digits',file=sys.stderr,flush=True)
        with mp.workdps(digits):
            family.append(family_increment(order,p['incoming_moment_panels']))
    family_gap=max(relative(a['amplitude_increment'],b['amplitude_increment'])
                   for a,b in zip(family[0]['rows'],family[-1]['rows']))
    finite=family[-1]
    controls=numeric_controls(p)
    preceding=json.loads((HERE.parent/'intermediate_decay_audit'/'evidence.json').read_text())
    old_failure=preceding['old_lambda_failure']
    rows=[r for f in fixtures for r in f['rows']]
    small_rows=finite['rows']
    gates={
        'all_interval_constants':all(cert['checks'].values()),
        'nested_shape_enclosures':certificates[0]['shape']['K_lower']<cert['shape']['K_lower']<cert['shape']['K_upper']<certificates[0]['shape']['K_upper'],
        'finite_parameter_order':cert['log_h_interval'][1]<cert['log_lambda_interval'][0]<-256 and cert['Tf']==1000,
        'amplitude_bracket':cert['left_residual_upper']<0<cert['right_residual_lower'] and [cert['amplitude_lower'],cert['amplitude_upper']]==p['amplitude_bracket'],
        'unique_root_derivative':cert['S_amplitude_derivative_lower']>mp.mpf('.35'),
        'uniform_finite_remainder':0<cert['remainder_universal_upper']<mp.mpf(p['claimed_remainder_bound']),
        'angular_correction_contraction':0<cert['angular_contraction_upper']<1,
        'axial_inverse_loss_bounded':cert['axial_inverse_times_lambda_upper']<20,
        'actual_fixture_M_J_rows':max(v for r in rows for v in r['linear_relative_errors'])<limits['linear_moment_relative_error'],
        'fixture_S_quadratic':max(r['quadratic_relative_residual'] for r in rows)<limits['quadratic_relative_residual'],
        'independent_end_patch_ODE':max(v['maximum_error'] for f in fixtures for v in f['ode_checks'])<limits['independent_source_error'],
        'main_moment_quadrature_refinement':max(f['refinement']['main_log_gap'] for f in fixtures)<limits['quadrature_relative_gap'],
        'coefficient_quadrature_refinement':max(f['refinement']['coefficient_relative_gap'] for f in fixtures)<limits['quadrature_relative_gap'],
        'incoming_moment_refinement':max(f['refinement']['incoming_relative_gap'] for f in fixtures)<limits['quadrature_relative_gap'],
        'coarse_incoming_resolution_failure':all(v['Ka_relative_gap']>limits['quadrature_relative_gap'] and v['K2a_relative_gap']>limits['quadrature_relative_gap'] for v in controls['incoming_resolution_failures']),
        'fixture_amplitude_refinement':max(f['refinement']['amplitude_relative_gap'] for f in fixtures)<limits['quadrature_relative_gap'],
        'actual_angular_rows':max(v for r in rows for v in r['post_mass']['angular_relative_errors'])<limits['linear_moment_relative_error'],
        'positive_full_tail_mass':all(x['mass']>0 for r in rows+small_rows for x in r['post_mass']['parts']),
        'nonzero_bump_energy_retained':all(r['S_parts']['bump_quadratic_increment']>0 for r in rows),
        'nonzero_affine_prefix_retained':all((r['eta']==0 and r['S_parts']['bump_linear_coefficient']==0) or (r['eta']!=0 and r['S_parts']['bump_linear_coefficient']!=0) for r in rows),
        'main_prefix_omission_bounded':all(0<v['prefix_relative_error_upper']<mp.mpf('1e-100') for f in fixtures for v in f['main_moments']),
        'only_M_correction_detected':all(r['only_M_correction_J_relative_error']>mp.mpf('.8') for r in rows),
        'missing_end_corrections_detected':all(r['omitted_end_corrections_relative_error']==1 for r in rows),
        'leading_root_not_exact':all(abs(r['leading_root_S_error'])>mp.mpf('1e-4') for r in rows),
        'large_lambda_bracket_failure':controls['large_lambda']['F_at_1p2_upper_diagnostic']<mp.mpf('-.1'),
        'binary64_singular_matrix_detected':controls['determinant']['binary64_difference']==0 and controls['determinant']['stable_difference']<0 and abs(controls['determinant']['stable_difference_over_lambda']+2*mp.e)<mp.mpf('1e-18'),
        'family_positive_amplitude_increment':all(r['amplitude_increment']>0 for r in small_rows),
        'family_increment_refinement':family_gap<limits['family_increment_relative_gap'],
        'family_direct_increment_crosscheck':max(r['direct_relative_error'] for r in small_rows)<limits['family_increment_relative_gap'],
        'family_omitted_bump_bound_positive':all(0<r['omitted_bump_amplitude_error_upper']<r['amplitude_increment']*mp.mpf('1e-100') for r in small_rows),
        'binary64_amplitude_increment_lost':all(r['binary64_increment']==0 for r in small_rows),
        'family_post_swirl_retained':all(r['post_mass']['total']>0 and any(v!=0 for v in r['post_mass']['angular_energy_increments']) for r in small_rows),
        'old_parameter_failure_retained':mp.mpf(old_failure['lambda_w_squared'])>1,
        'K_diagnostics_inside_enclosure':all(cert['shape']['K_lower']<v<cert['shape']['K_upper'] for v in [finite['K']]+[f['K'] for f in fixtures]),
    }
    print('Encoding the completed evidence and controls',file=sys.stderr,flush=True)
    return encode(dict(status=STATUS,protocol=p,provenance=provenance(),certificate=cert,
        coarser_shape=certificates[0]['shape'],fixtures=fixtures,
        family_diagnostic=finite,family_refinement_relative_gap=family_gap,
        controls=controls,old_lambda_failure=old_failure,gates=gates))


def validate(record):
    if record.get('status')!=STATUS: raise ValueError('Scientific scope changed')
    if record.get('protocol')!=json.loads((HERE/'protocol.json').read_text()): raise ValueError('Protocol changed')
    if record.get('provenance')!=provenance(): raise ValueError('Source provenance changed')
    gates=record.get('gates',{})
    if set(gates)!=GATE_NAMES or not all(v is True for v in gates.values()):
        raise ValueError('Missing or failed gates: '+str([k for k in GATE_NAMES if gates.get(k) is not True]))
    c=record['certificate']
    if c['claimed_remainder_upper']!='1e-100' or c['status']!='reference-moments-closed-pulse-cone-unverified':
        raise ValueError('Certificate scope or bound changed')
    if not mp.mpf(c['left_residual_upper'])<0<mp.mpf(c['right_residual_lower']): raise ValueError('Root bracket erased')
    if mp.mpf(record['old_lambda_failure']['lambda_w_squared'])<=1: raise ValueError('Old failure erased')
    if mp.mpf(record['controls']['large_lambda']['F_at_1p2_upper_diagnostic'])>=0: raise ValueError('Large-lambda failure erased')
    if any(mp.mpf(r['Ka_relative_gap'])<=mp.mpf('1e-20') for r in record['controls']['incoming_resolution_failures']):
        raise ValueError('Under-resolved incoming moment failure erased')
    for r in record['family_diagnostic']['rows']:
        if not 0<mp.mpf(r['omitted_bump_amplitude_error_upper'])<mp.mpf(r['amplitude_increment']):
            raise ValueError('Finite increment or positive remainder erased')
    def finite_tree(v):
        if isinstance(v,dict):
            for x in v.values(): finite_tree(x)
        elif isinstance(v,list):
            for x in v: finite_tree(x)
        elif isinstance(v,(str,float,int)) and not isinstance(v,bool):
            try: number=mp.mpf(v)
            except (ValueError,TypeError): return
            if not mp.isfinite(number): raise ValueError('Nonfinite evidence')
    finite_tree(record)


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
        tolerance=mp.mpf('1e-9') if 'ode_checks' in path else mp.mpf('1e-18')
        if not mp.isfinite(x) or not mp.isfinite(y): raise ValueError('Nonfinite record: '+path)
        # Recomputed error estimates vary near rounding noise. Physical tiny
        # coefficients, increments and positive bounds always compare relatively.
        diagnostic=('relative_error' in path or 'relative_residual' in path or 'gap' in path or 'ode_checks' in path)
        logarithmic=('log_' in path)
        gap=abs(x-y) if diagnostic or logarithmic else relative(x,y)
        if gap>tolerance: raise ValueError('Numeric record differs: '+path)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify-record',type=Path)
    args=parser.parse_args(argv)
    result=build()
    validate(result)
    if args.verify_record:
        old=json.loads(args.verify_record.read_text())
        validate(old)
        compare_record(old,result)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],gates=len(result['gates']),
        amplitude_interval=[result['certificate']['amplitude_lower'],result['certificate']['amplitude_upper']],
        remainder_upper=result['certificate']['remainder_universal_upper'])))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
