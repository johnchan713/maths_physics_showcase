#!/usr/bin/env python3
"""Reproduce scheduled pressure, actual angular correction, and inner coupling."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time
import mpmath as mp
from schedule import Schedule,step
from coupling import (CoupledParameters,SeedDatum,Parameters,construct,compare,diagnostics,
                      full_residual,physical_finite_difference,PREVIOUS)

HERE=Path(__file__).resolve().parent
STATUS='scheduled-pressure-coupled-global-matching-open'
SOURCES=('protocol.json','schedule.py','coupling.py','audit.py','test_outer.py',
         '../nonlinear_axis_pilot/inner.py','../nonlinear_axis_pilot/physical.py',
         '../nonlinear_axis_pilot/protocol.json')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def encode(x):
    if isinstance(x,dict):
        return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [encode(v) for v in x]
    if isinstance(x,mp.mpf):
        if not mp.isfinite(x):
            raise ValueError('Nonfinite evidence')
        return mp.nstr(x,60)
    return x


def rel(a,b):
    return abs(a-b)/max(mp.mpf(1),abs(a),abs(b))


def eta_value(label,p):
    return p.hzero() if label=='Hzero' else mp.mpf(label)


def independent_schedule_checks(s):
    """Double-precision ODE integration does not reuse the quadrature primitive."""
    import math
    from scipy.integrate import solve_ivp
    def sigma(x):
        if x<=0:
            return 0.
        if x>=1:
            return 1.
        a,b=-1/x**2,-1/(1-x)**2
        m=max(a,b)
        return math.exp(a-m)/(math.exp(a-m)+math.exp(b-m))
    rows=[]
    for stage in s.stages:
        if stage.kind!='ramp':
            continue
        left,right,h=float(stage.left_slope),float(stage.right_slope),float(s.h)
        q0=float(s.q_initial)
        def rhs(t,z):
            slope=left+(right-left)*sigma(t)
            return [slope-.5,math.exp(2*z[0]),-(1+slope)*z[2]-slope-h]
        solution=solve_ivp(rhs,[0,1],[0,0,q0],method='DOP853',rtol=2e-13,atol=2e-15,max_step=.01)
        if not solution.success:
            raise ValueError('Independent ODE integration failed')
        predicted=[s.log_shape(stage,mp.mpf(1),mp.mpf(0)),
                   sum(w for w,t in next(v for v in s.components if v['name']==stage.name)['nodes']),
                   s.propagate_q(mp.mpf(q0),stage.left_slope,stage.right_slope)]
        rows.append(dict(name=stage.name,ode=[mp.mpf(float(v)) for v in solution.y[:,-1]],
                         quadrature=predicted,error=max(rel(a,mp.mpf(float(b))) for a,b in zip(predicted,solution.y[:,-1]))))
    # Tanh-sinh integration cross-checks the non-ramp finite shapes separately.
    for stage in s.stages:
        if stage.kind not in ('interpolation','terminal'):
            continue
        actual=sum(w for w,t in next(v for v in s.components if v['name']==stage.name)['nodes'])
        reference=stage.length*mp.quad(lambda x:mp.exp(2*s.log_shape(stage,stage.length*x,mp.mpf(0))),
                                       [0,mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.75'),1])
        rows.append(dict(name=stage.name,tanh_sinh=reference,quadrature=actual,error=rel(actual,reference)))
    primitive=[]
    for x in map(mp.mpf,['.1','.37','.5','.83']):
        reference=mp.quad(step,[0,x/2,x])
        primitive.append(abs(s.rule.step_integral(x)-reference))
    return dict(rows=rows,primitive_errors=primitive,maximum_error=max(primitive+[v['error'] for v in rows]))


def quadrature_rows(s,degree):
    return [dict(eta=e,parts=s.pressure_parts(mp.mpf(e),degree)) for e in ('0','0.5','1')]


def quadrature_gap(a,b):
    errors=[]
    for x,y in zip(a,b):
        for u,v in zip(x['parts'],y['parts']):
            errors.append(abs(u['log_scale']-v['log_scale']))
            errors.extend(rel(c,d) for c,d in zip(u['jet'],v['jet']))
    return max(errors)


def independent_angular_checks(s,corrections,rejected_panels=4):
    """Integrate the edited moments by tanh-sinh, normalizing before subtraction."""
    from schedule import step_prime,Rule
    beta,rate,width=1-s.lam,1+2*s.lam,mp.mpf('.3')
    B=mp.matrix(2,2);Q=[]
    for i,center in enumerate((mp.mpf(0),mp.mpf(2))):
        y=lambda z:center+width*(z-mp.mpf('.5'))
        B[0,i]=width*mp.quad(lambda z:mp.exp(beta*y(z))*step_prime(z),[0,mp.mpf('.5'),1])
        B[1,i]=2*width*mp.quad(lambda z:mp.exp(-rate*y(z))*step_prime(z),[0,mp.mpf('.5'),1])
        Q.append(width*mp.quad(lambda z:mp.exp(-rate*y(z))*step_prime(z)**2,[0,mp.mpf('.5'),1]))
    rows=[]
    for row in corrections:
        c=row['coefficients'];normal=max(abs(v) for v in c)
        angular=mp.fsum([B[0,i]*(c[i]/normal) for i in range(2)]+[-row['discrepancy']/normal])
        pressure=mp.fsum([B[1,i]*(c[i]/normal) for i in range(2)]+[Q[i]*c[i]*(c[i]/normal) for i in range(2)])
        rows.append(dict(label=row['label'],angular_error=abs(angular),pressure_error=abs(pressure)))
    # A direct integral of H gives an independent actual discrepancy at eta=.5.
    # This center avoids cancellation down to the exponentially smaller endpoint memory.
    eta=mp.mpf('.5');a=mp.log(2/(1+eta*eta))
    integral=mp.quad(lambda y:mp.exp(-beta*(s.tf-y)+a*step(1-y/s.tf)),
                     [0,s.tf/2,s.tf*mp.mpf('.8'),s.tf*mp.mpf('.9'),s.tf])
    # Bound the pre-interpolation memory instead of silently counting it as zero.
    memory_bound=mp.exp(-beta*(s.tw+13/s.lam+s.tf)+a)/beta
    reference=mp.fsum([mp.exp(-beta*s.tf+a)/beta,integral,-1/beta])
    predicted=s.interpolation_moment_discrepancy(eta)
    relative_discrepancy=abs(predicted-reference)/abs(predicted)
    rejected=s.interpolation_moment_discrepancy(eta,Rule(s.rule.order,rejected_panels))
    return dict(rows=rows,discrepancy_relative_error=relative_discrepancy,
                rejected_panel_count=rejected_panels,rejected_relative_error=abs(rejected-reference)/abs(reference),
                omitted_memory_relative_bound=memory_bound/abs(predicted),
                maximum_error=max([relative_discrepancy]+[v[k] for v in rows for k in ('angular_error','pressure_error')]))


def check_limits(protocol):
    c=protocol['checks'];pc=protocol['physical_check']
    return {
        'stage-quadrature-refinement':('at-most',c['quadrature_refinement']),
        'independent-schedule-integrals':('at-most',c['independent_quadrature']),
        'pressure-jet-derivatives':('at-most',c['pressure_derivatives']),
        'pressure-evenness':('at-most',c['pressure_derivatives']),
        'pressure-lower-bound-slack':('at-least','0'),
        'pressure-monotonicity':('at-least',c['negative_control']),
        'angular-moment-and-pressure-restoration':('at-most',c['angular_moment_correction']),
        'independent-edited-moment-integrals':('at-most',c['independent_quadrature']),
        'angular-edit-size':('at-most','1e-8'),
        'angular-correction-smallness-estimate':('at-most','1'),
        'independent-inner-assembly':('at-most','1e-100'),
        'original-leading-equations':('at-most',c['leading_equations']),
        'radial-refinement':('at-most',c['radial_refinement']),
        'precision-refinement-including-amplitude':('at-most',c['precision_refinement']),
        'positive-normalized-swirl':('at-least',c['minimum_Phi']),
        'physical-leading-equations':('at-most',c['leading_equations']),
        'physical-divergence':('at-most','1e-100'),
        'physical-finite-difference':('at-most',pc['maximum_error']),
        'physical-finite-difference-refinement':('at-least',pc['minimum_reduction']),
        'detect-omitted-pressure-constraint':('at-least',c['negative_control']),
        'detect-underresolved-angular-integral':('at-least',c['independent_quadrature']),
        'detect-insufficient-amplitude-precision':('at-least','1')}


def measured_checks(report,protocol):
    p=report['pressure_samples'];corrections=report['angular_corrections'];profiles=report['inner_samples']
    numeric=lambda x:mp.mpf(x)
    values={
        'stage-quadrature-refinement':report['quadrature']['last_gap'],
        'independent-schedule-integrals':report['independent_schedule']['maximum_error'],
        'pressure-jet-derivatives':max(numeric(r['derivative_error']) for r in p),
        'pressure-evenness':max(numeric(r['evenness_error']) for r in p),
        'pressure-lower-bound-slack':min(numeric(r['lower_bound_slack']) for r in p),
        'pressure-monotonicity':min(numeric(r['eta_times_derivative']) for r in p if numeric(r['eta'])!=0),
        'angular-moment-and-pressure-restoration':max(numeric(v) for r in corrections for v in r['relative_errors']),
        'independent-edited-moment-integrals':report['independent_angular']['maximum_error'],
        'angular-edit-size':max(numeric(r['maximum_relative_edit_bound']) for r in corrections),
        'angular-correction-smallness-estimate':max(numeric(r['smallness_bound']) for r in corrections),
        'independent-inner-assembly':report['seed_solver_agreement'],
        'original-leading-equations':max(numeric(r[k]) for a in profiles for r in a['refinements'][-1]['samples']
                                         for k in ('angular_residual','axial_residual','pressure_residual')),
        'radial-refinement':max(numeric(r['degree_gap']) for r in profiles),
        'precision-refinement-including-amplitude':max(numeric(r['precision_gap']) for r in profiles),
        'positive-normalized-swirl':min(numeric(r['Phi']) for a in profiles for r in a['refinements'][-1]['samples']),
        'physical-leading-equations':max(map(numeric,report['physical']['normalized_leading'])),
        'physical-divergence':report['physical']['divergence'],
        'physical-finite-difference':report['physical']['fd_errors'][-1],
        'physical-finite-difference-refinement':numeric(report['physical']['fd_errors'][0])/numeric(report['physical']['fd_errors'][-1]),
        'detect-omitted-pressure-constraint':min(numeric(r['negative_control_pressure_gap']) for r in corrections),
        'detect-underresolved-angular-integral':report['independent_angular']['rejected_relative_error'],
        'detect-insufficient-amplitude-precision':report['rejected_precision']['absolute_log_amplitude_gap']}
    out={}
    for name,(kind,limit) in check_limits(protocol).items():
        v,l=numeric(values[name]),numeric(limit)
        passed=mp.isfinite(v) and v>=0 and (v<=l if kind=='at-most' else v>=l)
        out[name]=dict(value=v,kind=kind,limit=l,passed=bool(passed))
    return out


def audit(protocol,progress=False):
    started=time.monotonic()
    def note(**data):
        if progress:
            print(json.dumps(data),flush=True)
    mp.mp.dps=max(protocol['decimal_precisions'])
    histories=[]
    for order in protocol['quadrature_orders']:
        s=Schedule(protocol['parameters'],order,protocol['quadrature_panels'],protocol['angular_moment_panels'])
        histories.append(dict(order=order,rows=quadrature_rows(s,protocol['pressure_jet_order'])))
        note(phase='schedule',order=order,seconds=round(time.monotonic()-started,2))
    qgap=quadrature_gap(histories[-2]['rows'],histories[-1]['rows'])
    independent=independent_schedule_checks(s)
    p=CoupledParameters.from_schedule(s,protocol['parameters'])
    pressure_rows=[];corrections=[]
    for label in protocol['eta_samples']:
        eta=eta_value(label,p)
        jet=s.pressure_jet(eta,8,normalized=True)
        derivatives=[rel(mp.factorial(k)*jet[k],mp.diff(lambda e:s.pressure(e,normalized=True),eta,k)) for k in (1,2,3,5,8)]
        pressure_rows.append(dict(label=label,eta=eta,normalized_pressure=jet[0],jet=jet,
            derivative_error=max(derivatives),evenness_error=rel(jet[0],s.pressure(-eta,True)),
            lower_bound_slack=-jet[0]-mp.mpf('2.5')/(1+eta*eta)**2,eta_times_derivative=eta*jet[1]))
        corrections.append(dict(label=label,eta=eta,**s.angular_correction(eta)))
    independent_angular=independent_angular_checks(s,corrections,protocol['rejected_angular_moment_panels'])
    # Separate the tiny non-f^2 terms; their magnitude is not compared with 1.
    parts=s.pressure_parts(mp.mpf(0),0)
    first_late=next(i for i,v in enumerate(parts) if v['name']=='parameter-interpolation')
    core=-mp.fsum(mp.exp(v['log_scale'])*v['jet'][0] for v in parts[:first_late])
    tail=-mp.fsum(mp.exp(v['log_scale'])*v['jet'][0] for v in parts[first_late:])
    note(phase='pressure-and-angular-correction',quadrature_gap=mp.nstr(qgap,8))
    from inner import construct as seed_construct
    seed_data=json.loads((PREVIOUS/'protocol.json').read_text())['parameters']
    seed=Parameters.from_dict(seed_data)
    coupled_seed=CoupledParameters(seed.h,seed.j,seed.sigma,seed.lam,mp.mpf(0),SeedDatum(seed.pressure))
    seed_agreement=max(compare(seed_construct(seed,e,18),construct(coupled_seed,e,18),[mp.mpf(4)])
                       for e in (mp.mpf('.5'),seed.hzero()))
    with mp.workdps(min(protocol['decimal_precisions'])):
        low=Schedule(protocol['parameters'],max(protocol['quadrature_orders']),protocol['quadrature_panels'],protocol['angular_moment_panels'])
        low_p=CoupledParameters.from_schedule(low,protocol['parameters'])
    ys=list(map(mp.mpf,protocol['radial_samples']))
    inner_rows=[];fine_profiles={}
    for label in protocol['eta_samples']:
        levels=[]
        with mp.workdps(min(protocol['decimal_precisions'])):
            for degree in protocol['radial_degrees']:
                levels.append(construct(low_p,eta_value(label,low_p),degree,protocol['eta_derivatives_retained']))
        fine=construct(p,eta_value(label,p),max(protocol['radial_degrees']),protocol['eta_derivatives_retained'])
        fine_profiles[label]=fine
        degree_gap=compare(levels[-2],levels[-1],ys)
        amplitude_gap=abs(mp.expm1(levels[-1].log_g0-fine.log_g0))
        precision_gap=max(compare(levels[-1],fine,ys),amplitude_gap)
        samples=[dict(degree=v.degree,precision=min(protocol['decimal_precisions']) if v is not fine else max(protocol['decimal_precisions']),
                      samples=[diagnostics(v,y) for y in ys]) for v in levels+[fine]]
        inner_rows.append(dict(label=label,eta=fine.eta0,degree_gap=degree_gap,precision_gap=precision_gap,
                               relative_amplitude_gap=amplitude_gap,refinements=samples))
        note(phase='inner',eta=label,degree_gap=mp.nstr(degree_gap,8),precision_gap=mp.nstr(precision_gap,8))
    with mp.workdps(protocol['rejected_precision']):
        rejected=Schedule(protocol['parameters'],min(protocol['quadrature_orders']),protocol['quadrature_panels'],protocol['angular_moment_panels'])
        rejected_p=CoupledParameters.from_schedule(rejected,protocol['parameters'])
        rejected_log=rejected_p.log_g(mp.mpf('.5'))
    rejected_gap=abs(rejected_log-fine_profiles['0.5'].log_g0)
    pc=protocol['physical_check']
    from physical import physical_point
    with mp.workdps(pc['digits']):
        profile=fine_profiles[pc['eta']]
        point=physical_point(profile,mp.mpf(pc['Y']),mp.mpf(pc['q']))
        physical=full_residual(profile,point)
        fd_errors=[]
        for step in pc['steps']:
            approximate=physical_finite_difference(profile,point,mp.mpf(step))
            fd_errors.append(max(abs(x-y)/scale for x,y,scale in zip(approximate,physical['total'],physical['scales'])))
        physical.update(point=point,fd_steps=pc['steps'],fd_errors=fd_errors,digits=pc['digits'])
    report=dict(status=STATUS,full_proof_verified=False,new_candidate=False,matched_paper_profile=False,
        global_thresholds_verified=False,prior_unforced_status='unchanged-resolution-failure',unverified=protocol['unverified'],
        derived=dict(Td=s.td,log_P=s.log_p,h=s.h,Lambda=p.lam,lambda_outer=s.lam,Tw=s.tw,Tf=s.tf,
                     terminal_Qp=s.qp,terminal_wait=s.wait,log_C=p.log_c()),
        stages=[dict(name=v.name,kind=v.kind,start=v.start,length='infinite' if mp.isinf(v.length) else v.length,
                     log_amplitude_relative_to_P=v.log_amplitude,theta_at_start=v.theta) for v in s.stages],
        quadrature=dict(histories=histories,last_gap=qgap),independent_schedule=independent,
        pressure_samples=pressure_rows,angular_corrections=corrections,independent_angular=independent_angular,
        pressure_decomposition=dict(coefficient_of_negative_P_squared_f_squared=core,
            positive_tail_mass_normalized=tail,uniform_relative_tail_bound=4*tail/core,
            log10_uniform_relative_tail_bound=mp.log10(4*tail/core)),
        seed_solver_agreement=seed_agreement,inner_samples=inner_rows,physical=physical,
        rejected_precision=dict(digits=protocol['rejected_precision'],absolute_log_amplitude_gap=rejected_gap,
            explanation='Normalized fields alone conceal the loss of absolute log-amplitude accuracy; this precision is rejected.'))
    report['checks']=measured_checks(report,protocol)
    if not all(v['passed'] for v in report['checks'].values()):
        report['status']='pilot-needs-review'
    return report


def validate(report,protocol):
    def finite(x):
        if isinstance(x,dict):
            for v in x.values():
                finite(v)
        elif isinstance(x,list):
            for v in x:
                finite(v)
        elif isinstance(x,str):
            try:
                v=mp.mpf(x)
            except ValueError:
                return
            if not mp.isfinite(v):
                raise ValueError('Nonfinite evidence')
        elif isinstance(x,float) and not mp.isfinite(x):
            raise ValueError('Nonfinite evidence')
    finite(report)
    # Scientific limits are part of the contract, even when numerical gates pass.
    if report['status']!=STATUS or any(report[k] is not False for k in
        ('full_proof_verified','new_candidate','matched_paper_profile','global_thresholds_verified')):
        raise ValueError('Unsupported scientific status')
    if report['unverified']!=protocol['unverified'] or report['prior_unforced_status']!='unchanged-resolution-failure':
        raise ValueError('Lost scientific limitations')
    if [r['label'] for r in report['inner_samples']]!=protocol['eta_samples']:
        raise ValueError('Changed inner sampling')
    for name in ('pressure_samples','angular_corrections'):
        if [r['label'] for r in report[name]]!=protocol['eta_samples']:
            raise ValueError('Changed parameter sampling')
    if [r['order'] for r in report['quadrature']['histories']]!=protocol['quadrature_orders']:
        raise ValueError('Changed quadrature refinement')
    for row in report['inner_samples']:
        levels=[(n,min(protocol['decimal_precisions'])) for n in protocol['radial_degrees']]+[(max(protocol['radial_degrees']),max(protocol['decimal_precisions']))]
        if [(r['degree'],r['precision']) for r in row['refinements']]!=levels:
            raise ValueError('Changed inner refinements')
        for r in row['refinements']:
            if [str(v['Y']) for v in r['samples']]!=[str(v) for v in map(lambda x:mp.nstr(mp.mpf(x),60),protocol['radial_samples'])]:
                raise ValueError('Changed radial samples')
    expected=encode(measured_checks(report,protocol))
    if set(report['checks'])!=set(expected):
        raise ValueError('Missing or extra numerical gate')
    for name,value in expected.items():
        actual=report['checks'][name]
        if actual['kind']!=value['kind'] or mp.mpf(actual['limit'])!=mp.mpf(value['limit']):
            raise ValueError('Changed acceptance limit')
        if not value['passed'] or actual['passed'] is not True or rel(mp.mpf(actual['value']),mp.mpf(value['value']))>mp.mpf('1e-55'):
            raise ValueError('Numerical gate inconsistent with measurements: '+name)
    if report['physical']['fd_steps']!=protocol['physical_check']['steps'] or report['physical']['digits']!=protocol['physical_check']['digits']:
        raise ValueError('Changed physical check')
    if mp.mpf(report['pressure_decomposition']['positive_tail_mass_normalized'])<=0:
        raise ValueError('Discarded pressure tail')


def compare_record(a,b):
    """Compare signed stage data, including tiny tail masses, without a unit floor."""
    def visit(x,y,path=''):
        if isinstance(x,dict):
            if set(x)!=set(y):
                raise ValueError('Reproduction keys differ: '+path)
            for k in x:
                visit(x[k],y[k],path+'/'+k)
        elif isinstance(x,list):
            if len(x)!=len(y):
                raise ValueError('Reproduction lengths differ')
            for i,(u,v) in enumerate(zip(x,y)):
                visit(u,v,path+'/'+str(i))
        elif isinstance(x,str):
            try:
                u,v=mp.mpf(x),mp.mpf(y)
            except ValueError:
                if x!=y:
                    raise ValueError('Reproduction labels differ')
                return
            if not mp.isfinite(u) or not mp.isfinite(v):
                raise ValueError('Nonfinite reproduction')
            # Near-zero residuals may vary with library arithmetic. Measured
            # gate limits are always rechecked; profile/tail values use relative error.
            residual=any(k in path for k in ('error','gap','residual','divergence','agreement'))
            denominator=max(1,abs(u),abs(v)) if residual else max(abs(u),abs(v))
            if denominator and abs(u-v)/denominator>mp.mpf('1e-11'):
                raise ValueError('Reproduction value differs: '+path)
        elif x!=y:
            raise ValueError('Reproduction metadata differs: '+path)
    for group in ('derived','stages','quadrature','independent_schedule','pressure_samples','angular_corrections',
                  'pressure_decomposition','independent_angular','seed_solver_agreement','inner_samples','physical','rejected_precision'):
        visit(a[group],b[group],group)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--paper',type=Path)
    parser.add_argument('--verify-record',type=Path)
    args=parser.parse_args(argv)
    if args.output.exists():
        parser.error('Existing evidence is not overwritten')
    protocol=json.loads((HERE/'protocol.json').read_text())
    mp.mp.dps=max(protocol['decimal_precisions'])
    hashes={name:sha(HERE/name) for name in SOURCES}
    if args.paper and sha(args.paper)!=protocol['paper_sha256']:
        raise ValueError('Paper hash mismatch')
    old=None
    if args.verify_record:
        old=json.loads(args.verify_record.read_text())
        if old['source_sha256']!=hashes:
            raise ValueError('Source binding changed')
        validate(old,protocol)
    started=time.monotonic()
    result=encode(audit(protocol,progress=True))
    import scipy
    result.update(source_sha256=hashes,paper_bytes_verified=bool(args.paper),
                  environment=dict(python=platform.python_version(),mpmath=mp.__version__,scipy=scipy.__version__),
                  elapsed_seconds=time.monotonic()-started)
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False)
        stream.write('\n')
    validate(result,protocol)
    if old:
        compare_record(old,result)
    print(json.dumps(dict(status=result['status'],checks=len(result['checks']),seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
