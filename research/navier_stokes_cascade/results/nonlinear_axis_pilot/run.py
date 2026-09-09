#!/usr/bin/env python3
"""Reproduce a source-bound nonlinear pilot without promoting it to a proof."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time
import mpmath as mp
from inner import Parameters, construct, diagnostics, direct_heat_join
from physical import physical_point, residual, finite_difference, cartesian_finite_difference

HERE = Path(__file__).resolve().parent
SOURCES = ('protocol.json','inner.py','physical.py','run.py','test_pilot.py')
STATUS = 'nonlinear-axis-pilot-passed-direct-join-failed'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def eta_value(label,p):
    if label.startswith('Hzero'):
        suffix = label[len('Hzero'):]
        return p.hzero()+(mp.mpf(suffix) if suffix else 0)
    return mp.mpf(label)


def encode(value):
    if isinstance(value, dict):
        return {k:encode(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)):
        return [encode(v) for v in value]
    if isinstance(value, mp.mpf):
        if not mp.isfinite(value):
            raise ValueError('Nonfinite numerical evidence')
        return mp.nstr(value,50)
    return value


def relative(a,b):
    return abs(a-b)/max(mp.mpf(1),abs(a),abs(b))


def comparison(a,b,Ys):
    errors = []
    for y in Ys:
        for name in ('phi','u'):
            for dy,de in ((0,0),(1,0),(2,0),(0,1),(0,2),(0,3)):
                errors.append(relative(a.value(getattr(a,name),y,radial_order=dy,eta_order=de),
                                       b.value(getattr(b,name),y,radial_order=dy,eta_order=de)))
        # g is nonzero and needs a relative check without a unit floor.
        errors.append(abs(mp.expm1(a.log_g0-b.log_g0)))
    return max(errors)


def moment_quadrature(profile,Y):
    """Factor out tiny g before quadrature, so absolute tolerances cannot hide it."""
    from inner import inner_moments
    p=profile.parameters
    g=mp.exp(profile.log_g0)
    phi=lambda y:profile.value(profile.phi,y)
    U=lambda y:4*profile.eta0+p.j+profile.value(profile.u,y)/p.lam
    refs=[mp.quad(U,[0,Y]),mp.quad(lambda y:y*phi(y),[0,Y]),
          mp.quad(lambda y:y*U(y)*phi(y),[0,Y]),
          mp.quad(lambda y:U(y)**2-g*g*y*phi(y)**2/p.lam,[0,Y]),
          mp.quad(lambda y:phi(y)**2,[0,Y])]
    actual=inner_moments(profile,Y)
    scales=[p.lam,p.lam**2/(2*g),p.lam**2/(2*g),p.lam,p.lam/(g*g)]
    normalized=[a*s for a,s in zip(actual,scales)]
    return dict(normalized_polynomial=normalized,normalized_quadrature=refs,
                maximum_error=max(relative(a,b) for a,b in zip(normalized,refs)))


def audit(protocol,progress=False):
    mp.mp.dps=max(protocol['decimal_precisions'])
    p=Parameters.from_dict(protocol['parameters'])
    Ys=[mp.mpf(v) for v in protocol['radial_samples']]
    rows=[]
    fine_profiles={}
    worst_leading=worst_degree=worst_precision=mp.mpf(0)
    min_phi=mp.inf
    axis_error=mp.mpf(0)
    negative_angular=negative_pressure=mp.mpf(0)
    for label in protocol['eta_samples']:
        start=time.monotonic()
        levels=[]
        for degree in protocol['radial_degrees']:
            with mp.workdps(min(protocol['decimal_precisions'])):
                pp=Parameters.from_dict(protocol['parameters'])
                levels.append(construct(pp,eta_value(label,pp),degree,protocol['eta_derivatives_retained']))
        fine=construct(p,eta_value(label,p),max(protocol['radial_degrees']),protocol['eta_derivatives_retained'])
        fine_profiles[label]=fine
        degree_error=comparison(levels[-2],levels[-1],Ys)
        precision_error=comparison(levels[-1],fine,Ys)
        worst_degree=max(worst_degree,degree_error)
        worst_precision=max(worst_precision,precision_error)
        histories=[]
        for level in levels+[fine]:
            diagnostics_rows=[diagnostics(level,y) for y in Ys]
            histories.append(dict(degree=level.degree,precision=(min(protocol['decimal_precisions']) if level is not fine else max(protocol['decimal_precisions'])),samples=diagnostics_rows))
        samples=histories[-1]['samples']
        for r in samples:
            worst_leading=max(worst_leading,*(r[k] for k in ('angular_residual','axial_residual','pressure_residual')))
            min_phi=min(min_phi,r['Phi'])
            negative_angular=max(negative_angular,r['omitted_eta_transport_gap'])
            negative_pressure=max(negative_pressure,r['omitted_axial_pressure_gap'])
        axis_error=max(axis_error,samples[0]['axis_U_error'],abs(samples[0]['Phi']-1))
        rows.append(dict(eta_label=label,eta=fine.eta0,degree_difference=degree_error,precision_difference=precision_error,
                         refinements=histories,direct_join=direct_heat_join(fine),
                         local_jet_coefficients_first_four={name:[r[:4] for r in getattr(fine,name)] for name in ('phi','u','pressure_increment')}))
        if progress:
            print(json.dumps(dict(eta=label,seconds=round(time.monotonic()-start,3),degree_gap=mp.nstr(degree_error,5))),flush=True)
    physical=[]
    pc=protocol['physical_checks']
    for label in pc['eta']:
        profile=fine_profiles[label]
        point=physical_point(profile,mp.mpf(pc['Y']),mp.mpf(pc['q']))
        r=residual(profile,point)
        physical.append(dict(eta_label=label,point=point,**r))
    fd_profile=fine_profiles[pc['finite_difference_eta']]
    point=physical_point(fd_profile,mp.mpf(pc['Y']),mp.mpf(pc['q']))
    exact=residual(fd_profile,point)
    fd_errors=[]
    for step in pc['relative_steps']:
        approximation=finite_difference(fd_profile,point,mp.mpf(step))
        fd_errors.append(max(abs(a-b)/s for a,b,s in zip(approximation,exact['total'],exact['scales'])))
    rejected=cartesian_finite_difference(fd_profile,point,mp.mpf(pc['relative_steps'][-1]))
    rejected_error=abs(rejected[1]-exact['total'][1])/exact['scales'][1]
    derivatives=[]
    dc=protocol['eta_derivative_check']
    epsilon=mp.mpf(dc['step']);y=mp.mpf(dc['Y'])
    for label in dc['eta']:
        center=fine_profiles[label]
        neighbors=[construct(p,center.eta0+i*epsilon,center.degree,protocol['eta_derivatives_retained']) for i in (-2,-1,1,2)]
        errors=[]
        for name in ('phi','u'):
            a,b,c,d=[v.value(getattr(v,name),y) for v in neighbors]
            middle=center.value(getattr(center,name),y)
            first=(a-8*b+8*c-d)/(12*epsilon)
            second=(-a+16*b-30*middle+16*c-d)/(12*epsilon*epsilon)
            errors.extend([relative(first,center.value(getattr(center,name),y,eta_order=1)),
                           relative(second,center.value(getattr(center,name),y,eta_order=2))])
        derivatives.append(dict(eta_label=label,maximum_error=max(errors)))
    moments=moment_quadrature(fine_profiles['0.5'],mp.mpf(4))
    thresholds=protocol['checks']
    checks={}
    def at_most(name,value,limit):
        checks[name]=dict(value=value,kind='at-most',limit=mp.mpf(limit),passed=bool(0<=value<=mp.mpf(limit)))
    def at_least(name,value,limit):
        checks[name]=dict(value=value,kind='at-least',limit=mp.mpf(limit),passed=bool(value>=mp.mpf(limit)))
    at_most('original-leading-equations',worst_leading,thresholds['maximum_leading_residual'])
    at_most('radial-degree-refinement',worst_degree,thresholds['maximum_degree_difference'])
    at_most('arithmetic-precision-refinement',worst_precision,thresholds['maximum_precision_difference'])
    at_most('axis-data',axis_error,thresholds['maximum_axis_error'])
    at_least('positive-normalized-swirl',min_phi,thresholds['minimum_Phi'])
    at_least('detect-missing-eta-transport',negative_angular,thresholds['minimum_negative_control_gap'])
    at_least('detect-missing-axial-pressure',negative_pressure,thresholds['minimum_negative_control_gap'])
    at_most('cartesian-divergence',max(r['divergence'] for r in physical),pc['maximum_divergence'])
    at_most('cartesian-leading-equations',max(v for r in physical for v in r['normalized_leading']),pc['maximum_leading_error'])
    at_most('physical-finite-difference',fd_errors[-1],pc['maximum_finest_fd_error'])
    at_least('physical-finite-difference-refinement',fd_errors[0]/fd_errors[-1],pc['minimum_fd_reduction'])
    at_most('independently-resolved-eta-derivatives',max(r['maximum_error'] for r in derivatives),dc['maximum_error'])
    at_most('five-moment-quadrature',moments['maximum_error'],'1e-40')
    matching=protocol['matching_checks']
    joins=[dict(eta_label=r['eta_label'],
                passed=bool(r['direct_join']['U_gap']<=mp.mpf(matching['maximum_normalized_U']) and
                            r['direct_join']['M_gap']<=mp.mpf(matching['maximum_normalized_M']) and
                            r['direct_join']['heat_slope_necessary_gap']<=mp.mpf(matching['maximum_heat_angular_slope_gap']))) for r in rows]
    local_pass=all(c['passed'] for c in checks.values())
    return dict(status=STATUS if local_pass and not all(j['passed'] for j in joins) else 'pilot-needs-review',
                full_proof_verified=False,new_candidate=False,matched_paper_profile=False,
                prior_unforced_status='unchanged-resolution-failure',unverified=protocol['unverified'],
                checks=checks,samples=rows,physical=physical,
                finite_difference=dict(steps=pc['relative_steps'],errors=fd_errors),
                rejected_mixed_cartesian_stencil=dict(angular_relative_error=rejected_error,
                    explanation='Absorption of tiny swirl in larger radial samples; retained failure, not a physical divergence.'),
                eta_derivatives=derivatives,moment_quadrature=moments,direct_join_gates=joins,
                pressure_seed_not_paper_schedule=True)


def validate(report,protocol):
    frozen=protocol['checks'];pc=protocol['physical_checks']
    limits={
        'original-leading-equations':('at-most',frozen['maximum_leading_residual']),
        'radial-degree-refinement':('at-most',frozen['maximum_degree_difference']),
        'arithmetic-precision-refinement':('at-most',frozen['maximum_precision_difference']),
        'axis-data':('at-most',frozen['maximum_axis_error']),
        'positive-normalized-swirl':('at-least',frozen['minimum_Phi']),
        'detect-missing-eta-transport':('at-least',frozen['minimum_negative_control_gap']),
        'detect-missing-axial-pressure':('at-least',frozen['minimum_negative_control_gap']),
        'cartesian-divergence':('at-most',pc['maximum_divergence']),
        'cartesian-leading-equations':('at-most',pc['maximum_leading_error']),
        'physical-finite-difference':('at-most',pc['maximum_finest_fd_error']),
        'physical-finite-difference-refinement':('at-least',pc['minimum_fd_reduction']),
        'independently-resolved-eta-derivatives':('at-most',protocol['eta_derivative_check']['maximum_error']),
        'five-moment-quadrature':('at-most','1e-40')}
    expected={'original-leading-equations','radial-degree-refinement','arithmetic-precision-refinement','axis-data',
              'positive-normalized-swirl','detect-missing-eta-transport','detect-missing-axial-pressure','cartesian-divergence',
              'cartesian-leading-equations','physical-finite-difference','physical-finite-difference-refinement',
              'independently-resolved-eta-derivatives','five-moment-quadrature'}
    if set(report['checks']) != expected:
        raise ValueError('Missing or unexpected checks')
    if [r['eta_label'] for r in report['samples']] != protocol['eta_samples']:
        raise ValueError('Changed eta samples')
    def finite(value):
        if isinstance(value,dict):
            for child in value.values():
                finite(child)
        elif isinstance(value,list):
            for child in value:
                finite(child)
        elif isinstance(value,str):
            try:
                numeric=mp.mpf(value)
            except ValueError:
                return
            if not mp.isfinite(numeric):
                raise ValueError('Nonfinite raw evidence')
        elif isinstance(value,float) and not mp.isfinite(value):
            raise ValueError('Nonfinite raw evidence')
    finite(report)
    for r in report['samples']:
        expected_levels=[(n,min(protocol['decimal_precisions'])) for n in protocol['radial_degrees']]+[(max(protocol['radial_degrees']),max(protocol['decimal_precisions']))]
        if [(v['degree'],v['precision']) for v in r['refinements']] != expected_levels:
            raise ValueError('Changed refinement levels')
        for v in r['refinements']:
            if [mp.mpf(s['Y']) for s in v['samples']] != [mp.mpf(y) for y in protocol['radial_samples']]:
                raise ValueError('Changed radial samples')
    for name,check in report['checks'].items():
        value,limit=mp.mpf(check['value']),mp.mpf(check['limit'])
        if (check['kind'],limit)!=(limits[name][0],mp.mpf(limits[name][1])):
            raise ValueError('Changed numerical threshold')
        good=mp.isfinite(value) and mp.isfinite(limit) and value>=0 and (value<=limit if check['kind']=='at-most' else value>=limit)
        if not good or check['passed'] is not True:
            raise ValueError('Numerical pilot gate failed')
    if report['status'] != STATUS or any(report[k] is not False for k in ('full_proof_verified','new_candidate','matched_paper_profile')):
        raise ValueError('Unsupported scientific promotion')
    if report['prior_unforced_status'] != 'unchanged-resolution-failure' or report['unverified'] != protocol['unverified'] or report['pressure_seed_not_paper_schedule'] is not True:
        raise ValueError('Lost scientific limitations')
    if len(report['direct_join_gates']) != len(protocol['eta_samples']) or any(v['passed'] is not False for v in report['direct_join_gates']):
        raise ValueError('Direct join failure not preserved')
    if [v['eta_label'] for v in report['physical']]!=pc['eta'] or [v['eta_label'] for v in report['eta_derivatives']]!=protocol['eta_derivative_check']['eta']:
        raise ValueError('Changed independent samples')
    if report['finite_difference']['steps']!=pc['relative_steps'] or len(report['finite_difference']['errors'])!=len(pc['relative_steps']):
        raise ValueError('Changed finite-difference refinement')
    if mp.mpf(report['rejected_mixed_cartesian_stencil']['angular_relative_error'])<=1:
        raise ValueError('Discovered Cartesian absorption failure lost')


def compare_record(old,new):
    # Check measurements, excluding timing/environment metadata and tiny cancellation errors.
    for group in ('samples','physical','finite_difference','eta_derivatives','moment_quadrature','direct_join_gates'):
        def compare(a,b):
            if isinstance(a,dict):
                if set(a)!=set(b):
                    raise ValueError('Reproduction keys differ')
                for k in a:
                    compare(a[k],b[k])
            elif isinstance(a,list):
                if len(a)!=len(b):
                    raise ValueError('Reproduction lengths differ')
                for x,y in zip(a,b):
                    compare(x,y)
            elif isinstance(a,str):
                try:
                    aa,bb=mp.mpf(a),mp.mpf(b)
                except ValueError:
                    if a!=b:
                        raise ValueError('Reproduction labels differ')
                else:
                    if not mp.isfinite(aa) or not mp.isfinite(bb) or relative(aa,bb)>mp.mpf('1e-35'):
                        raise ValueError('Reproduction measurements differ')
            elif a!=b:
                raise ValueError('Reproduction values differ')
        compare(old[group],new[group])


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify-record',type=Path)
    parser.add_argument('--paper',type=Path)
    args=parser.parse_args(argv)
    if args.output.exists():
        parser.error('Existing evidence is never overwritten')
    protocol=json.loads((HERE/'protocol.json').read_text())
    source_hashes={name:sha(HERE/name) for name in SOURCES}
    if args.paper and sha(args.paper)!=protocol['paper_sha256']:
        raise ValueError('Paper hash mismatch')
    old=None
    if args.verify_record:
        old=json.loads(args.verify_record.read_text())
        if old['source_sha256']!=source_hashes:
            raise ValueError('Source binding changed')
        validate(old,protocol)
    start=time.monotonic()
    report=encode(audit(protocol,progress=True))
    report.update(source_sha256=source_hashes,paper_bytes_verified=bool(args.paper),
                  environment=dict(python=platform.python_version(),mpmath=mp.__version__),elapsed_seconds=time.monotonic()-start)
    with args.output.open('x') as output:
        json.dump(report,output,indent=2,sort_keys=True,allow_nan=False)
        output.write('\n')
    validate(report,protocol)
    if old:
        compare_record(old,report)
    print(json.dumps(dict(status=report['status'],checks=len(report['checks']),seconds=report['elapsed_seconds'])))


if __name__=='__main__':
    main()
