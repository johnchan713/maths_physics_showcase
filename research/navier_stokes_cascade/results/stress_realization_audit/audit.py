#!/usr/bin/env python3
"""Reproduce the conditional stress mechanism and its unresolved input ledger."""
import argparse
import hashlib
import json
from pathlib import Path
import mpmath as mp
from loop import (point,lower,upper,variance,variance_iv,root_bracket,
                  shear_ratio,zero_pressure_primitives)
from repair import Repair,divided_weight,precondition_target
from bounds import (repair_certificate,fixture_cap_certificate,active_gap_bounds,uniform_mu_log_bound,
                    cone_tolerance,frequency_requirement)

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent
STATUS='conditional-stress-realization-bounded-global-frequency-uninstantiated'
FLAGS=dict(conditional_shear_loop_argument_checked=True,
           conditional_five_moment_restoration_bounded=True,
           actual_compact_input_bounds_instantiated=False,
           actual_joined_profile_frequency_selected=False,
           full_admissible_stress_realized=False,
           full_PDE_corrections_verified=False,
           smooth_force_verified=False,
           formal_proof_assistant_certificate=False,
           independent_peer_review_completed=False,
           blowup_verified=False)


def validate_protocol(protocol):
    expected=dict(status=STATUS,
        source_pdf_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        source_pages=[29,31,38,127,128,129,158,159,160,161,162,163,164,165],
        interval_digits=[80,110],fixture_digits=80,fixture_a='0.8',fixture_ts='0.2',
        fixture_p1='10',fixture_p2=['-1','-0.25','-1e-30','0','1e-30','0.25','1'],
        fixture_d0='1',fixture_mu_cap='64',fixture_delta='1e-8',repair_lambda='0.01',
        quadrature_rules=[[32,8],[48,12]],actual_joined_profile_frequency_selected=False)
    if protocol!=expected or protocol.get('actual_joined_profile_frequency_selected') is not False:
        raise ValueError('Protocol differs from the implemented certificate and scope')


def source_hashes():
    paths=[HERE/name for name in ('README.md','protocol.json','requirements.txt',
                                 'loop.py','repair.py','bounds.py','audit.py','test_audit.py')]
    paths += [RESULTS/name/'evidence.json' for name in
              ('axis_core_attachment','axis_matching_audit','heat_exterior_audit',
               'intermediate_decay_audit')]
    paths += [RESULTS/'outer_pressure_pilot'/'schedule.py']
    return {str(p.relative_to(RESULTS.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths}


def inherited_provenance():
    checks={}
    for folder,key in (('axis_core_attachment','source_hashes'),
                       ('heat_exterior_audit','provenance')):
        data=json.loads((RESULTS/folder/'evidence.json').read_text())
        checks[folder]=all(hashlib.sha256((RESULTS.parent/p).read_bytes()).hexdigest()==sha
                           for p,sha in data[key].items())
    return checks


def loop_fixture(protocol,p2):
    a,ts,p1,d0,cap,delta=[mp.mpf(protocol['fixture_'+key]) for key in
                          ('a','ts','p1','d0','mu_cap','delta')]
    p2=mp.mpf(p2)
    target=(2+delta/2-a*(1+ts*ts))/a
    left,right=root_bracket(target,p2,d0,cap)
    mu=(left+right)/2
    v=a*(1+ts*ts+variance(mu,p2,d0))
    def t(theta):return shear_ratio(theta,mu,p2,ts,d0)
    def avg(f):return mp.quad(f,[0,mp.pi/2,mp.pi,3*mp.pi/2,2*mp.pi])/(2*mp.pi)
    phase_mass=avg(lambda th:a*(1+t(th)**2)/v)
    # Changed-variable cancellation is exact pointwise in the a row.
    mean_a=avg(lambda th:v/(1+t(th)**2)*a*(1+t(th)**2)/v)
    mean_minus_b=avg(lambda th:v*t(th)/(1+t(th)**2)*a*(1+t(th)**2)/v)
    unweighted_a=avg(lambda th:v/(1+t(th)**2))
    p=point
    target_iv=(2+p(protocol['fixture_delta'])/2
               -p(protocol['fixture_a'])*(1+p(protocol['fixture_ts'])**2))/p(protocol['fixture_a'])
    return dict(p2=p2,mu_bracket=[left,right],target_variance=target_iv,
                root_left_variance=variance_iv(left,p2,d0),
                root_right_variance=variance_iv(right,p2,d0),
                root_signs_certified=(upper(variance_iv(left,p2,d0))<lower(target_iv)
                                     and lower(variance_iv(right,p2,d0))>upper(target_iv)),
                cap_reaches_target=lower(variance_iv(cap,p2,d0))>3/a,
                phase_mass_error=abs(phase_mass-1),
                weighted_a_error=abs(mean_a-a),weighted_minus_b_error=abs(mean_minus_b-a*ts),
                unweighted_a_error=abs(unweighted_a-a),
                manufactured_input=True,integrated_p_from_profile_verified=False)


def modulation_controls():
    a,s,E=mp.mpf('.8'),mp.mpf('.1'),mp.mpf('1.3')
    rows=[]
    for theta in (mp.mpf('.3'),mp.mpf('1.1'),mp.mpf('2.2')):
        for N in (10,100,1000):
            f=zero_pressure_primitives(theta,a,E=E)
            At=mp.diff(lambda th:zero_pressure_primitives(th,a,E=E)['A'],theta)
            Bt=mp.diff(lambda th:zero_pressure_primitives(th,a,E=E)['B'],theta)
            ap=1-2*(s+At/f['phase_derivative'])
            en=E*mp.exp(f['A']/N)
            bp=2*(s*f['B']/N+Bt/f['phase_derivative'])/en
            wanted=mp.exp(-f['A']/N)*(f['b']+2*s*f['B']/(N*E))
            omitted=f['b']+2*s*f['B']/(N*E)
            rows.append(dict(N=N,theta=theta,a_identity_error=abs(ap-f['a']),
                             b_identity_error=abs(bp-wanted),
                             omitted_exponential_error=abs(bp-omitted)))
    derivatives=[]
    for N in (10,100,1000):
        first=mp.diff(lambda y:mp.sin(2*mp.pi*N*y)/N,0)
        second=mp.diff(lambda y:mp.sin(2*mp.pi*N*y)/N,mp.mpf(1)/(4*N),2)
        eta=mp.diff(lambda e:mp.sin(2*mp.pi*N*e)/N,0)
        derivatives.append(dict(N=N,radial_first=first,radial_second_over_N=second/N,
                                eta_dependent_phase_derivative=eta))
    return dict(cases=rows, derivative_controls=derivatives,
                maximum_identity_error=max(max(r['a_identity_error'],r['b_identity_error']) for r in rows),
                omitted_exponential_detected=min(r['omitted_exponential_error'] for r in rows)>mp.mpf('1e-8'),
                radial_derivative_is_not_small=all(r['radial_first']>6 and abs(r['radial_second_over_N'])>39 for r in derivatives),
                eta_dependent_phase_can_destroy_C1_smallness=all(r['eta_dependent_phase_derivative']>6 for r in derivatives),
                manufactured_input=True)


def repair_fixture(lam,order,panels):
    repair=Repair(lam,order,panels)
    target=mp.matrix([mp.mpf('1e-16')*x for x in (1,-2,3,-1,2)])
    coefficients=repair.solve(target)
    residual=max(abs(x) for x in repair.value(coefficients)-target)
    direct=max(abs(x) for x in repair.direct_changes(coefficients)-repair.ordinary_changes(coefficients))
    missing_cp=repair.solve([0,0,0,0,mp.mpf('1e-16')])
    remaining=repair.value(missing_cp)
    return dict(order=order,panels=panels,coefficients=list(coefficients),
                numerical_residual=residual,direct_delta_integrand_error=direct,
                inverse_max_row_sum=max(mp.fsum(abs(repair.inverse[i,j]) for j in range(5)) for i in range(5)),
                bump_mass_error=max(abs(repair.integral(j,lambda y:1)-1) for j in range(5)),
                omitted_pressure_control=dict(first_four_moment_error=max(abs(remaining[j]) for j in range(4)),
                    remaining_pressure_increment=remaining[4],rejected_as_exact_restoration=True),
                quadrature_is_an_exact_profile=False)


def cancellation_controls():
    log_lam=-4*(mp.exp(64)+10)
    lam=mp.exp(log_lam)
    y=mp.mpf('1.5')
    stable=divided_weight(lam,y)
    naive=(1-mp.exp(-lam*y))/lam
    rounded_incoming=precondition_target([1,0,1-lam,0,0],lam)[1]
    moderate=precondition_target([mp.mpf('1e-16'),0,0,0,0],mp.mpf('1e-5'))
    return dict(selected_log_lambda=log_lam,selected_lambda_remains_positive=lam>0,
                stable_divided_weight=stable,naive_divided_weight=naive,
                cancellation_detected=naive==0 and stable>1,
                discrepancy_preconditioning_loss=moderate[1]/moderate[0],
                original_loss_log=-log_lam,zero_lambda_not_a_valid_repair=True,
                rounded_incoming_difference_lost=rounded_incoming==0,
                incoming_difference_needs_separate_stable_ledger=True,
                stable_limit_does_not_remove_discrepancy_loss=True)


def generate():
    protocol=json.loads((HERE/'protocol.json').read_text())
    validate_protocol(protocol)
    with mp.workdps(125):
        certificates=[repair_certificate(d) for d in protocol['interval_digits']]
        cap_certificates=[fixture_cap_certificate(d) for d in protocol['interval_digits']]
        active=active_gap_bounds('.8','.2',10,1,'7.8',1,64,'1e-8')
        uniform=uniform_mu_log_bound('.8',1,1)
    with mp.workdps(protocol['fixture_digits']):
        loops=[loop_fixture(protocol,p) for p in protocol['fixture_p2']]
        modulation=modulation_controls()
        repairs=[repair_fixture(protocol['repair_lambda'],*rule) for rule in protocol['quadrature_rules']]
        refinement=max(abs(a-b) for a,b in zip(repairs[0]['coefficients'],repairs[1]['coefficients']))
        cancellation=cancellation_controls()
        epsilon=cone_tolerance('.1',20,'1e-4')
        frequency=frequency_requirement(epsilon=epsilon,state_constant=100,
            correction_constant=100,discrepancy_constant=1000,coefficient_tolerance='1e-8')
    inherited=json.loads((RESULTS/'axis_core_attachment'/'evidence.json').read_text())
    provenance=inherited_provenance()
    gates=dict(
        inherited_source_hashes_unchanged=all(provenance.values()),
        inherited_core_checks_pass=all(inherited['gates'].values()),
        Md4_failure_retained=mp.mpf(inherited['inherited_Md4_outer_failure']['Pc_over_ps1'])<0,
        analytic_inverse_and_elementary_bounds=all(all(c['checks'].values()) for c in certificates),
        active_all_phase_cone=all(active['checks'].values()) and all(lower(g)>0 for g in active['gaps']),
        exact_variance_root_brackets=all(r['root_signs_certified'] for r in loops),
        individual_fixture_caps_reach=all(r['cap_reaches_target'] for r in loops),
        full_fixture_pressure_interval_covered=all(r['full_pressure_range_covered'] for r in cap_certificates),
        phase_reparametrization_has_unit_mass=max(r['phase_mass_error'] for r in loops)<mp.mpf('1e-55'),
        both_prescribed_shear_means=max(max(r['weighted_a_error'],r['weighted_minus_b_error']) for r in loops)<mp.mpf('1e-55'),
        unweighted_average_control=min(r['unweighted_a_error'] for r in loops)>mp.mpf('.01'),
        exact_modulation_identities=modulation['maximum_identity_error']<mp.mpf('1e-60'),
        exponential_denominator_control=modulation['omitted_exponential_detected'],
        angular_and_radial_derivative_controls=(modulation['radial_derivative_is_not_small'] and modulation['eta_dependent_phase_can_destroy_C1_smallness']),
        five_numerical_moments_restored=max(r['numerical_residual'] for r in repairs)<mp.mpf('1e-70'),
        independent_delta_integrands=max(r['direct_delta_integrand_error'] for r in repairs)<mp.mpf('1e-70'),
        independent_quadrature_refinement=refinement<mp.mpf('1e-40'),
        positive_bump_mass_check=max(r['bump_mass_error'] for r in repairs)<mp.mpf('1e-25'),
        missing_pressure_moment_detected=all(r['omitted_pressure_control']['remaining_pressure_increment']>mp.mpf('9e-17') for r in repairs),
        extreme_lambda_cancellation_retained=cancellation['cancellation_detected'],
        lambda_loss_not_erased=(cancellation['discrepancy_preconditioning_loss']==100000
            and cancellation['rounded_incoming_difference_lost']
            and cancellation['incoming_difference_needs_separate_stable_ledger']),
        actual_global_inputs_not_promoted=not any(FLAGS[k] for k in
            ('actual_compact_input_bounds_instantiated','actual_joined_profile_frequency_selected','full_admissible_stress_realized','blowup_verified')),
    )
    return dict(status=STATUS,flags=FLAGS,protocol=protocol,certificates=certificates,
                fixture_cap_certificates=cap_certificates,
                active_fixture_all_phase_bounds=active,universal_cap_example=uniform,
                loop_fixtures=loops,modulation_controls=modulation,
                repair_fixtures=repairs,repair_refinement=refinement,
                cancellation_controls=cancellation,
                manufactured_frequency_budget=dict(epsilon=epsilon,**frequency),
                unresolved_actual_inputs=[
                    'Certified compact interval endpoints and retained analytic/edge collars',
                    'Bounds on a, bs, ps, E and their required mixed/eta derivatives',
                    'Variance cap and positive loop plus unmodified-patch cone margins',
                    'C1 normalized discrepancy constant including lambda^-1',
                    'State/pressure comparison and correction constants',
                    'One finite N selected from those actual inputs before physical q',
                ],
                inherited_provenance_checks=provenance,
                inherited_core_status=inherited['status'],
                inherited_Md4_outer_failure=inherited['inherited_Md4_outer_failure'],
                source_hashes=source_hashes(),gates=gates)


def encode(value):
    if hasattr(value,'_mpi_'):
        with mp.workdps(125):
            return dict(lower=mp.nstr(lower(value),120),upper=mp.nstr(upper(value),120))
    if isinstance(value,mp.mpf):return mp.nstr(value,100)
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encode(v) for v in value]
    return value


def compare(expected,actual,path='root'):
    if type(expected) is not type(actual):
        raise ValueError('Record type changed at '+path)
    if isinstance(actual,dict):
        if not isinstance(expected,dict) or set(expected)!=set(actual):
            raise ValueError('Record keys changed at '+path)
        for key in actual:compare(expected[key],actual[key],path+'.'+key)
    elif isinstance(actual,list):
        if not isinstance(expected,list) or len(expected)!=len(actual):
            raise ValueError('Record length changed at '+path)
        for i,v in enumerate(actual):compare(expected[i],v,path+'.'+str(i))
    elif isinstance(actual,str) and actual!=expected:
        if any(k in path for k in ('hash','status','protocol','unresolved','flags')):
            raise ValueError('Frozen metadata changed at '+path)
        try:
            with mp.workdps(125):
                a,b=mp.mpf(actual),mp.mpf(expected)
                if not mp.isfinite(a) or not mp.isfinite(b) or abs(a-b)>mp.mpf('1e-30')*max(abs(a),abs(b),mp.mpf('1e-65')):
                    raise ValueError('Numeric record changed at '+path)
        except (TypeError,ValueError):
            raise ValueError('Record changed at '+path) from None
    elif actual!=expected:raise ValueError('Record changed at '+path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-record',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    record=encode(generate())
    failures=[k for k,v in record['gates'].items() if not v]
    if failures:raise SystemExit('Failed gates: '+str(failures))
    if args.verify_record:compare(json.loads(args.verify_record.read_text()),record)
    args.output.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    print(STATUS)
    print(str(len(record['gates']))+' audit gates passed')
    print('Actual-profile frequency, full PDE corrections and blow-up remain unverified.')


if __name__=='__main__':main()
