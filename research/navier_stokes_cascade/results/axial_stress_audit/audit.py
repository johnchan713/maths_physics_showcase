#!/usr/bin/env python3
"""Reproduce bounded axial progress and retain both finite-parameter failures."""
import argparse
import hashlib
import json
from pathlib import Path
import mpmath as mp
from bounds import axial_bound, elementary_checks
from diagnostics import Moments, independent_ode, HERE, PREVIOUS

STATUS = 'axial-ratios-bounded-global-profile-unverified'


def encode(value):
    if isinstance(value, mp.mpf):
        return mp.nstr(value, 60)
    if isinstance(value, dict):
        return {k:encode(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


def relative(a, b):
    a, b = mp.mpf(a), mp.mpf(b)
    return abs(a-b)/max(abs(a),abs(b)) if a or b else mp.mpf(0)


def provenance():
    files = [HERE/name for name in ('bounds.py','diagnostics.py','audit.py','test_audit.py','protocol.json','README.md','requirements.txt')]
    files += [PREVIOUS/name for name in ('schedule.py','protocol.json','evidence.json')]
    return {str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def build():
    mp.mp.dps = 80
    p = json.loads((HERE/'protocol.json').read_text())
    limits = {k:mp.mpf(v) for k,v in p['checks'].items()}
    grid = [axial_bound(m, p['interval_boxes'][0], 50) for m in p['Md_candidates']]
    refinement = [axial_bound(p['selected_Md'], n, 50) for n in p['interval_boxes']]
    selected = refinement[-1]
    low = axial_bound(p['selected_Md'], p['interval_boxes'][-1], 30)
    precision_gap = abs(low['maximum_bsw_absolute_upper']-selected['maximum_bsw_absolute_upper'])
    m = Moments(p['moment_orders'][0])
    midpoint = mp.expm1(m.s.md/2)
    ys = [midpoint if y=='midpoint' else mp.mpf(y) for y in p['ode_y']]
    rows, ode_rows, errors, identity_errors = [], [], [], []
    for eta in p['ode_eta']:
        exact = [m.axial(y,eta) for y in ys]
        rows.extend(exact)
        for row in exact:
            identity_errors.append(relative(row['Q'], row['Q_from_positive_identity']))
        for step_size in p['ode_max_steps']:
            ode = independent_ode(m, eta, [float(y) for y in ys], step_size)
            for actual, expected in zip(ode, exact):
                gap = max(relative(actual[k],expected[k]) for k in ('Q','n'))
                errors.append(gap)
                ode_rows.append(dict(max_step=step_size, relative_error=gap, **actual))
    failure = m.axial(midpoint, '.5')
    later_failure = m.intermediate_start('.5')
    fine = Moments(p['moment_orders'][1])
    fine_failure = fine.axial(midpoint, '.5')
    fine_later = fine.intermediate_start('.5')
    quadrature_gap = max([relative(failure[k],fine_failure[k]) for k in ('Q','n','bsw')]
                         +[relative(later_failure[k],fine_later[k]) for k in ('Q','n','lambda_w_squared')])
    first_passing = next(r['Md'] for r in grid if r['status']=='axial-ratio-bound-passed')
    gates = elementary_checks()
    gates.update({
        'selected_candidate_is_first_passing_bound': first_passing==p['selected_Md'],
        'all_axial_bsw_upper_below_limit': selected['maximum_bsw_absolute_upper']<limits['maximum_bsw_upper'],
        'all_axial_first_margin_positive': selected['first_A24_margin_lower']>0,
        'all_axial_second_margin_above_limit': selected['second_A24_margin_lower']>limits['minimum_second_margin'],
        'Q_remains_positive_with_finite_h': 0<selected['relative_Q_loss_upper']<1,
        'finite_shear_remainder_retained': selected['bs_squared_upper']>0,
        'nested_boxes_tighten_bound': all(a['maximum_bsw_absolute_upper']>=b['maximum_bsw_absolute_upper'] for a,b in zip(refinement,refinement[1:])),
        'interval_precision_refinement': precision_gap<limits['interval_precision_gap'],
        'independent_source_ODE': max(errors)<limits['independent_ode_relative_error'],
        'positive_Q_identity': max(identity_errors)<limits['positive_Q_identity_relative_error'],
        'moment_quadrature_refinement': quadrature_gap<limits['moment_quadrature_relative_error'],
        'old_Md_rejection_retained': failure['Q']>0 and failure['Pc_over_ps1']<0,
        'old_lambda_rejection_retained': later_failure['Q']>0 and later_failure['lambda_w_squared']>1,
    })
    return encode(dict(status=STATUS, protocol=p, provenance=provenance(),
                       candidate_bounds=grid, selected=selected,
                       box_refinement=refinement, interval_precision_gap=precision_gap,
                       initial_moment_constants=dict(rI1=m.rI1,energy1=m.energy1),
                       old_Md_failure=failure, old_lambda_failure=later_failure,
                       moment_samples=rows, ode_samples=ode_rows,
                       maximum_ode_error=max(errors), maximum_Q_identity_error=max(identity_errors),
                       maximum_moment_quadrature_error=quadrature_gap, gates=gates))


def validate(record):
    if record.get('status') != STATUS:
        raise ValueError('Incorrect scientific scope/status')
    if record.get('protocol') != json.loads((HERE/'protocol.json').read_text()):
        raise ValueError('Protocol changed')
    if record.get('provenance') != provenance():
        raise ValueError('Source provenance changed')
    gates = record.get('gates', {})
    if len(gates)!=17 or not all(v is True for v in gates.values()):
        raise ValueError('Missing or failed audit gates')
    if mp.mpf(record['old_Md_failure']['Pc_over_ps1']) >= 0:
        raise ValueError('Lost the Md=4 counterexample')
    if mp.mpf(record['old_lambda_failure']['lambda_w_squared']) <= 1:
        raise ValueError('Lost the lambda counterexample')
    if record['selected']['status'] != 'axial-ratio-bound-passed':
        raise ValueError('Selected axial bound did not pass')


def compare_record(expected, actual, path=''):
    if isinstance(expected, dict):
        if not isinstance(actual,dict) or set(expected)!=set(actual):
            raise ValueError(f'Record keys differ at {path}')
        for k,v in expected.items(): compare_record(v,actual[k],path+'/'+k)
    elif isinstance(expected, list):
        if not isinstance(actual,list) or len(expected)!=len(actual):
            raise ValueError(f'Record lengths differ at {path}')
        for i,(a,b) in enumerate(zip(expected,actual)): compare_record(a,b,path+'/'+str(i))
    elif isinstance(expected, bool) or expected is None:
        if expected is not actual: raise ValueError(f'Record differs at {path}')
    elif expected!=actual:
        try:
            a,b = mp.mpf(expected),mp.mpf(actual)
        except (ValueError,TypeError):
            raise ValueError(f'Record differs at {path}') from None
        # ODE roundoff varies across SciPy/platforms. Every gate is recomputed.
        tolerance = mp.mpf('1e-8') if 'ode' in path else mp.mpf('1e-12')
        if not mp.isfinite(a) or not mp.isfinite(b) or abs(a-b)>tolerance*max(1,abs(a),abs(b)):
            raise ValueError(f'Numeric record differs at {path}')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify-record',type=Path)
    args = parser.parse_args(argv)
    result = build()
    validate(result)
    if args.verify_record:
        old = json.loads(args.verify_record.read_text())
        validate(old)
        compare_record(old,result)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','selected','maximum_ode_error','maximum_moment_quadrature_error')}))
    print(f"{len(result['gates'])} gates passed")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
