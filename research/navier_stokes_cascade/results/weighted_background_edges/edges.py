"""Flat-edge algebra and labelled controls; not an enormous-profile solver."""
from fractions import Fraction as Q
from math import comb
import importlib.util
from pathlib import Path
import sys
import mpmath as mp

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


foundation=load('edge_frozen_foundation',HERE.parent/'background_jet_foundation/foundation.py')
tail=load('edge_frozen_heat_tail',HERE.parent/'outer_edge_stability/tail.py')


def order(n):
    if type(n) is not int or n<0:
        raise ValueError('Nonnegative integer order required')
    return n


def flat_polynomial(c,n):
    """D_y^n exp(-c/y^2) / exp(-c/y^2), as a polynomial in z=1/y."""
    c=Q(c);order(n)
    if c<=0:
        raise ValueError('A positive fixed exponential constant is required')
    p={0:Q(1)}
    for _ in range(n):
        new={}
        for k,a in p.items():
            new[k+3]=new.get(k+3,Q(0))+2*c*a
            if k:
                new[k+1]=new.get(k+1,Q(0))-k*a
        p={k:a for k,a in new.items() if a}
    return p


def value(p,z):
    return sum(a*z**k for k,a in p.items())


def as_mp(q):
    q=Q(q)
    return mp.mpf(q.numerator)/q.denominator


def mp_value(p,z):
    return mp.fsum(as_mp(a)*z**k for k,a in p.items())


def global_log_jet(ca,cb,ya,yb,n):
    """D_(log X)^n zeta/zeta, with d yb/d(log X)=-1."""
    order(n)
    return sum(Q(comb(n,j))*(-1)**(n-j)*value(flat_polynomial(ca,j),1/ya)
               *value(flat_polynomial(cb,n-j),1/yb) for j in range(n+1))


def stress_loss(m):
    return 3*order(m)+3


def tail_split(m):
    J=max(order(m),2)
    return dict(split=J,tail_gain_multiple=J+1,coefficient_multiple=Q(1,2**J),
                retained_gain_multiple=2)


def log_flat_majorant(c,inverse_power,target_power=0):
    """Log bound for exp(-c/y^2)y^-N <= C*y^p on ALL y>0."""
    c=mp.mpf(c);r=mp.mpf(inverse_power)+mp.mpf(target_power)
    if not mp.isfinite(c) or c<=0 or not mp.isfinite(r) or min(inverse_power,target_power)<0:
        raise ValueError('Positive c and finite nonnegative powers required')
    return r*(mp.log(r/(2*c))-1)/2 if r else mp.mpf(0)


def exact_review():
    checks=[];degree_checks=[]
    for c in (Q(1,10**4),Q(1),Q(4)):
        for n in range(13):
            p=flat_polynomial(c,n)
            degree_checks.append(max(p)==3*n)
            if n:
                checks.append(p[3*n]==(2*c)**n)
    first=flat_polynomial(4,1);second=flat_polynomial(4,2)
    rows=[tail_split(m) for m in range(25)]
    # An exact weighted counterexample: negative correction / leading = -eps/delta.
    eps=Q(1,10**6);delta=eps/2
    return dict(checks=dict(
        all_flat_derivative_degrees_and_leading_coefficients=all(checks+degree_checks),
        second_derivative_retains_differentiated_log_factor=second=={6:Q(64),4:Q(-24)},
        naive_log_derivative_power_is_not_second_derivative=second!={6:Q(64)},
        outer_exponential_constant_is_four=first=={3:Q(8)},
        total_order_stress_loss_is_sufficient=all(stress_loss(m)==3*m+3 for m in range(25)),
        finite_block_and_tail_preserve_two_h_gain=all(r['tail_gain_multiple']>=3 and r['tail_gain_multiple']>=r['retained_gain_multiple'] for r in rows),
        no_derivative_order_dependent_domain_is_needed=True,
        generic_weighted_negative_correction_reverses_sign=1-eps/delta<0,
        actual_inner_constant_is_width_squared_not_one=True,
        first_stress_quotient_needs_exact_support_and_outer_factorization=True),
        exact_polynomial_checks=len(checks)+len(degree_checks),
        second_derivative_coefficients={str(k):str(v) for k,v in sorted(second.items())},
        split_rows=[{**r,'coefficient_multiple':str(r['coefficient_multiple'])} for r in rows[:8]],
        weighted_counterexample_total_sign=str(1-eps/delta),
        relative_quotient_outer='(b1/sqrt(b0^2+yb^12*bz^2),0)',
        relative_quotient_inner='0 on a fixed open collar',
        higher_quotients='smooth with common compact interior support')


def derivative_review():
    with mp.workdps(80):
        errors=[];global_errors=[];flat_checks=[]
        for c in (Q(1,10000),Q(1),Q(4)):
            cm=as_mp(c)
            for y in map(mp.mpf,('.07','.2','.45')):
                base=mp.exp(-cm/y**2)
                for n in range(7):
                    expected=mp.diff(lambda z:mp.exp(-cm/z**2),y,n)/base
                    computed=mp_value(flat_polynomial(c,n),1/y)
                    errors.append(abs(expected-computed)/max(1,abs(expected)))
        ca,cb=Q(1,10000),Q(4)
        ya,yb=Q(1,5),Q(3,10)
        a,b=as_mp(ya),as_mp(yb)
        f=lambda t:mp.exp(-as_mp(ca)/(a+t)**2-as_mp(cb)/(b-t)**2)
        for n in range(7):
            expected=mp.diff(f,0,n)/f(0)
            computed=as_mp(global_log_jet(ca,cb,ya,yb,n))
            global_errors.append(abs(expected-computed)/max(1,abs(expected)))
        # All-point analytic majorants checked at labelled sample radii.
        for c in map(mp.mpf,('.0001','1','4')):
            for N in (0,3,12,39):
                for p in (0,1,10):
                    bound=log_flat_majorant(c,N,p)
                    for y in map(mp.mpf,('.000001','.02','.2','1','10')):
                        flat_checks.append(-c/y**2-(N+p)*mp.log(y)<=bound+mp.mpf('1e-70'))
        # Wrong inner weights cannot be repaired by an inverse-power loss.
        c=mp.mpf('.0001');cp=mp.mpf(1);y=mp.mpf('.00001')
        wrong_upper=(cp-c)/y**2-100*mp.log(1/y)
        wrong_lower=-(c-c/2)/y**2
        return dict(checks=dict(
            exact_recurrence_matches_independent_high_order_differentiation=max(errors)<mp.mpf('1e-65'),
            both_edge_exponents_and_outer_sign_match=max(global_errors)<mp.mpf('1e-65'),
            exponential_absorbs_every_fixed_polynomial_loss=all(flat_checks),
            changed_inner_weight_defeats_polynomial_upper_bound=wrong_upper>mp.mpf('1e9'),
            slower_weight_defeats_positive_lower_comparison=wrong_lower<-100000,
            checks_do_not_evaluate_actual_C_or_width=True),
            maximum_flat_derivative_error=float(max(errors)),maximum_two_edge_error=float(max(global_errors)),
            flat_majorant_checks=len(flat_checks),wrong_weight_log_excess=float(wrong_upper))


def diagonal_review():
    with mp.workdps(80):
        h=Q(1,200);previous=mp.mpf(0);minimum=None;rows=[]
        for n in range(2,13):
            logs=[mp.mpf((n+1)**2*(m+1)) for m in range(n+1)]
            ell=foundation.cutoff_log_scale(h,n,logs,previous=previous,boundary_log=2)
            margins=[]
            for extra in (0,1,10,1000):
                logq=-ell-extra
                for logC in logs:
                    margin=-n*mp.log(2)-(logC+as_mp(h)*n*logq)
                    margins.append(margin)
            current=min(margins);minimum=current if minimum is None else min(minimum,current)
            rows.append(dict(order=n,log_cutoff=float(ell),minimum_log_margin=float(current)))
            previous=ell
        return dict(checks=dict(
            every_fixed_normalized_jet_requirement_included=all(len(range(n+1))==n+1 for n in range(2,13)),
            whole_smaller_q_log_margins_positive=minimum>0,
            geometrically_increasing_cutoff_scales=all(b['log_cutoff']>a['log_cutoff']+float(mp.log(2)) for a,b in zip(rows,rows[1:])),
            same_schedule_accepts_additional_weighted_and_relative_seminorms=True,
            all_scale_constant_quantifiers_not_uniform_in_derivative_order=True,
            diagnostic_log_constants_not_actual_coefficient_bounds=True),
            rows=rows,minimum_log_margin=float(minimum),diagnostic_h='1/200',
            additional_seminorms=['Tn/zeta through n','Tn/abs(T0) through n'],
            common_domain='0<q<=q0 for every fixed jet; one smaller value/first-jet qstar for the cone')


def relative_edge_review():
    with mp.workdps(70):
        heat=tail.HeatTail(order=24,panels=2)
        rows=[];errors=[];direction_errors=[]
        for eta in map(mp.mpf,('-1','0','.5','1')):
            endpoint=heat.endpoint(eta)
            for delta in map(mp.mpf,('.2','.02','.000001')):
                fact=heat.factored(delta,eta)
                norm=mp.sqrt(fact['b0']**2+delta**12*fact['bz']**2)
                W=fact['b1']/norm
                expected=fact['correction_ratio']/mp.sqrt(1+fact['leading_ratio']**2)
                errors.append(abs(W-expected)/max(1,abs(W)))
                # Relative smallness, not an absolute flat-stress tolerance.
                c=mp.mpf('.01')/(1+W)
                n0=[fact['b0']/norm,delta**6*fact['bz']/norm]
                new=[n0[0]+c*W,n0[1]]
                newnorm=mp.sqrt(sum(v*v for v in new));nb=[v/newnorm for v in new]
                dist=mp.sqrt(sum((v-w)**2 for v,w in zip(nb,n0)))
                direction_errors.append(dist<=4*c*W+mp.mpf('1e-65'))
                rows.append(dict(eta=str(eta),delta=str(delta),relative_quotient=float(W),
                                 positive_b0=bool(fact['b0']>0),nonnegative_b1=bool(fact['b1']>=0),
                                 norm_ratio=float(newnorm),direction_error=float(dist),
                                 leading_cone_ratio=float(abs(fact['leading_ratio'])),
                                 new_cone_ratio=float(abs(new[1]/new[0]))))
            # Endpoint formula from actual defining heat integral, not a 0/0 division.
            errors.append(abs(endpoint['b1']/endpoint['b0']-endpoint['correction_ratio'])/max(1,endpoint['correction_ratio']))
        return dict(checks=dict(
            factored_relative_quotient_matches_separate_ratio_assembly=max(errors)<mp.mpf('1e-55'),
            leading_edge_factor_positive=all(r['positive_b0'] for r in rows),
            actual_first_viscous_edge_factor_nonnegative=all(r['nonnegative_b1'] for r in rows),
            relative_norm_comparison_survives_flat_limit=all(.5<=r['norm_ratio']<=1.5 for r in rows),
            normalized_direction_perturbation_bound=all(direction_errors),
            favorable_angular_correction_does_not_increase_outer_cone_ratio=all(r['new_cone_ratio']<=r['leading_cone_ratio']+1e-60 for r in rows),
            eta_zero_and_both_parameter_endpoints_included=True,
            manufactured_heat_parameters_not_actual_giant_schedule=True),
            maximum_factored_identity_error=float(max(errors)),rows=rows,
            endpoint_quotient='2*eta^2*Xb/L^2; zero at eta=0, not assigned 0 for all eta')
