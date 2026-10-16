"""Manufactured closed-time q/chart jets and independent heat moment controls."""
from functools import lru_cache
from math import factorial,pi,sqrt,log
from pathlib import Path
import sys
import mpmath as mp
import numpy as np
from endpoint import heat_A,normalized_heat_time_jet,normalized_pressure_time_jet

HERE=Path(__file__).resolve().parent
sys.path.append(str(HERE.parent/'physical_wave_interface'))
from jets import Jet,INDICES,bump,variables

H=.005


def q_value(z,t,h=H):
    tau=1-t
    if not 0<h<.01 or tau<0 or (tau==0 and z==0):
        raise ValueError('Preterminal positive branch or punctured terminal point required')
    if tau==0:return abs(z)**(1/(.5-h))
    lo=max(tau,abs(z)**(1/(.5-h)));hi=2*tau+(2*z*z)**(1/(1-2*h))
    for _ in range(100):
        mid=(lo+hi)/2
        if mid-tau-z*z*mid**(2*h)>0:hi=mid
        else:lo=mid
    return (lo+hi)/2


def q_jet(z,t,h=H):
    q=Jet(q_value(z.value,t.value,h))
    for _ in range(4):q=q-(q-(1-t)-z*z*q**(2*h))/(1-2*h*z*z*q**(2*h-1))
    return q


def curl(A):
    return [A[2].diff(1)-A[1].diff(2),A[0].diff(2)-A[2].diff(0),A[1].diff(0)-A[0].diff(1)]


def field(coords,Q=None,wrong_bump=False):
    x,y,z,t=coords;r=(x*x+y*y)**.5;q=q_jet(z,t);A=.5+H;D=.5-H
    beta=sqrt(2)-1;rho=log(4-sqrt(2))/log(4+sqrt(2));dr=2*((1+H)*rho-H*1e-5)
    Y0=r**dr+beta*t;Y1=-beta*r**dr+t;theta=(y/x).atan()
    if Q is None:
        X=r*r/(2*q);eta=z/q**D;pa=q**(.5-A);pb=q**(-A);pp=q**(-2*A)
    else:
        R=r/Q**.5;Z=z/Q**D;s=q/Q
        X=R*R/(2*s);eta=Z/s**D
        pa=Q**(.5-A)*s**(.5-A);pb=Q**(-A)*s**(-A);pp=Q**(-2*A)*s**(-2*A)
        if wrong_bump:X=R*R/2 # deliberately chart-fixed instead of physical q-based
    envelope=bump((X-1.2)/.8)
    profile=envelope*(1+.2*eta+.1*(2*pi*Y0).sin()+.07*(2*pi*Y1).cos())
    psi=pa*profile
    potential=[-y/r*psi,x/r*psi,.13*psi*(7*theta+2*pi*Y0).sin()]
    swirl=pb*profile
    direct=[-y/r*swirl,x/r*swirl,Jet(0)]
    pressure=pp*envelope*(1+.2*eta+.11*(2*pi*Y0).sin()+.17*(7*theta).cos())
    velocity=[a+b for a,b in zip(curl(potential),direct)]
    return dict(q=q,potential=potential,direct=direct,pressure=pressure,velocity=velocity)


def relative(a,b):return float(np.linalg.norm(a-b)/(1+np.linalg.norm(a)+np.linalg.norm(b)))


@lru_cache(None)
def chart_stats():
    implicit=chart=div=0.;wrong=jump=0.;rows=0;time_sign=slow=wrong_slow=0.
    for sign in (-1,1):
        z=sign*.08**(.5-H)
        for tau in (0.,2**-12,2**-10):
            coords=variables((.4,.1,z,1-tau));physical=field(coords);q=physical['q']
            defect=q-(1-coords[3])-coords[2]*coords[2]*q**(2*H)
            implicit=max(implicit,float(np.linalg.norm(defect.c)/(1+np.linalg.norm(q.c))))
            L=1-2*H*z*z*q.value**(2*H-1)
            time_sign=max(time_sign,abs(q.diff(3).value+1/L))
            for Q in (1/8,1/16):
                candidate=field(coords,Q)
                s=q/Q;Z=coords[2]/Q**(.5-H);T=(1-coords[3])/Q
                defect=s-Z*Z*s**(2*H)-T
                slow=max(slow,float(np.linalg.norm(defect.c)/(1+np.linalg.norm(s.c))))
                bad=s-Z*Z*s**(2*H)-(1-coords[3])/Q**(1+H)
                wrong_slow=max(wrong_slow,abs(bad.diff(3).value))
                for name in ('potential','direct','pressure','velocity'):
                    one=physical[name] if isinstance(physical[name],list) else [physical[name]]
                    two=candidate[name] if isinstance(candidate[name],list) else [candidate[name]]
                    mask=[sum(a)<=2 if name=='velocity' else True for a in INDICES]
                    chart=max(chart,*(relative(a.c[mask],b.c[mask]) for a,b in zip(one,two)))
                rows+=1
            v=physical['velocity'];divergence=sum(v[i].diff(i) for i in range(3))
            mask=[sum(a)<=1 for a in INDICES]
            div=max(div,float(np.linalg.norm(divergence.c[mask])/(1+sum(np.linalg.norm(vv.c) for vv in v))))
            wrong=max(wrong,relative(field(coords,1/8,True)['pressure'].c,field(coords,1/16,True)['pressure'].c))
            if tau==0:jump=max(jump,abs(physical['pressure'].value))
    return dict(implicit_equation_error=implicit,chart_jet_error=chart,divergence_error=div,
        wrong_fixed_chart_bump_gap=wrong,zero_time_extension_pressure_jump=jump,
        q_time_sign_error=time_sign,slow_time_chart_error=slow,wrong_slow_time_power_gap=wrong_slow,chart_comparisons=rows)


@lru_cache(None)
def heat_stats():
    error=pressure_error=0.;comparisons=0
    with mp.workdps(65):
        for numerator,denominator in ((1,200),(9,1000)):
            h=mp.mpf(numerator)/denominator;hf=__import__('fractions').Fraction(numerator,denominator)
            moments=[]
            for j in range(6):
                integral=mp.quad(lambda v:mp.exp(-v)*v**(h+j),[0,1,mp.inf])
                actual=mp.rf(h,j)*integral/mp.gamma(1+h);moments.append(actual)
                exact=heat_A(hf,j);expected=mp.mpf(exact.numerator)/exact.denominator
                error=max(error,float(abs(actual-expected)/(1+abs(expected))));comparisons+=1
            for j in range(6):
                from math import comb
                actual=-4**j*sum(comb(j,k)*moments[k]*moments[j-k] for k in range(j+1))/(2+4*h+2*j)
                exact=normalized_pressure_time_jet(hf,j);expected=mp.mpf(exact.numerator)/exact.denominator
                pressure_error=max(pressure_error,float(abs(actual-expected)/(1+abs(expected))));comparisons+=1
    return dict(heat_moment_error=error,complete_pressure_convolution_error=pressure_error,comparisons=comparisons)


@lru_cache(None)
def Cauchy_stats():
    worst=0.;wrong=0.;gaps=[]
    with mp.workdps(65):
        G=lambda t:mp.exp(t)+mp.sin(3*t)/5
        for k in range(3):
            derivatives=[mp.diff(G,mp.mpf(1),j) for j in range(k+1)]
            for tau in (mp.mpf('0.01'),mp.mpf('0.0001')):
                approximation=sum((-tau)**j*derivatives[j]/mp.factorial(j) for j in range(k+1))
                bound=(mp.e+mp.mpf(3)**(k+1)/5)*tau**(k+1)/mp.factorial(k+1)
                worst=max(worst,float(abs(G(1-tau)-approximation)/bound))
                if k==1:
                    bad=derivatives[0]+tau*derivatives[1]
                    wrong=max(wrong,float(abs(G(1-tau)-bad)/bound))
        for n in (1,10,100):
            a=1/(mp.pi/2+2*mp.pi*n);b=1/(3*mp.pi/2+2*mp.pi*n)
            gaps.append(float(abs(mp.sin(1/a)-mp.sin(1/b))))
    return dict(maximum_finite_Taylor_remainder_ratio=worst,reversed_time_Taylor_gap=wrong,
        bounded_value_without_next_derivative_limit_gaps=gaps)


def diagnostic_review():
    a,b,c=chart_stats(),heat_stats(),Cauchy_stats()
    return dict(manufactured_not_actual_terminal_fields=True,implicit_q_derivative_order=3,
        chart_comparisons=a['chart_comparisons'],heat_and_pressure_moment_comparisons=b['comparisons'],
        tolerances=dict(chart_relative=2e-10,implicit_relative=2e-12,heat_moment_relative=1e-45),
        raw_roundoff_scores_not_refrozen_under_an_ancestor_comparator=True,
        checks=dict(punctured_terminal_q_jet_satisfies_implicit_equation=a['implicit_equation_error']<2e-12,
            correct_q_time_sign_at_both_terminal_poles=a['q_time_sign_error']<2e-12,
            slow_chart_time_is_tau_over_Q_not_tau_over_Q_to_one_plus_h=a['slow_time_chart_error']<2e-12 and a['wrong_slow_time_power_gap']>1e-3,
            complete_physical_representative_and_velocity_jets_agree_across_charts=a['chart_jet_error']<2e-10,
            terminal_curl_plus_direct_swirl_is_divergence_free=a['divergence_error']<2e-12,
            independent_fixed_R_pressure_bumps_do_not_glue=a['wrong_fixed_chart_bump_gap']>1e-3,
            extending_nonzero_terminal_pressure_by_zero_jumps=a['zero_time_extension_pressure_jump']>1e-3,
            independent_positive_Gamma_integrals_match_exact_heat_jets=b['heat_moment_error']<1e-45,
            complete_pressure_Leibniz_sum_matches_independent_moments=b['complete_pressure_convolution_error']<1e-45,
            finite_order_terminal_Taylor_remainder_respects_analytic_bound=c['maximum_finite_Taylor_remainder_ratio']<1,
            wrong_tau_time_sign_fails_Taylor_bound=c['reversed_time_Taylor_gap']>100,
            bounded_values_and_origin_flatness_alone_do_not_imply_terminal_limits=min(c['bounded_value_without_next_derivative_limit_gaps'])>1.99))
