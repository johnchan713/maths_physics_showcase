"""Independent rational Taylor flows and a separate raw second-jet stress check.

All fixtures are manufactured, not the enormous pressure-family core.
"""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from math import factorial
from pathlib import Path
import sys
import mpmath as mp

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m


@lru_cache(None)
def modules():
    return (load('continuation_independent_rational_jets',HERE.parent/'core_all_order_foundation/coordinate_checks.py'),
            load('continuation_old_raw_second_jets',HERE.parent/'compact_jet_envelope/jets.py'))


def jet(v=0,order=5):return modules()[0].Jet(v,n=2,order=order)


def exponential(z):
    if z.value()!=0:raise ValueError('Exact exponential requires zero constant term')
    out=term=jet(1,z.order)
    for i in range(1,z.order+1):term=term*z/i;out=out+term
    return out


def primitive(z):
    return jet({(i+1,j):v/Q(i+1) for (i,j),v in z.terms.items()},z.order)


def coordinates(X,root,E,U,mom,axis,h=Q(1,200)):
    M,I,J,S,Cp=mom;eta=modules()[0].Jet.variable(1,coordinates.center,n=2,order=X.order)
    d=1-eta*eta;L=1-2*h*eta*eta;A=Q(1,2)+h;D=Q(1,2)-h
    W=1-(2*D*eta*M+d*M.derivative(1))/X
    Qs=-W+((1-h)*I-D*eta*I.derivative(1)-d*J.derivative(1)+2*(h-D)*eta*J)/(X*root*E)
    Pi=axis+Cp
    Ns=-W*U+(D*(M-eta*M.derivative(1))+4*h*eta*S-d*S.derivative(1))/X+4*A*eta*Pi-d*Pi.derivative(1)
    return dict(W=W,Qs=Qs,Ns=Ns,Pi=Pi,p1=X*Qs/L,ns=Ns/L,p2=X*Ns/(L*E))


def densities(X,root,E,U):return [X*U,X*root*E,X*root*U*E,X*(U*U-E*E/2),E*E/2]


def reference(center=Q(1,3),order=5,drop_axis=False):
    J=modules()[0].Jet;y=J.variable(0,0,n=2,order=order);eta=J.variable(1,center,n=2,order=order)
    X=2*exponential(y);root=2*exponential(y/2)
    e=2+eta/10+eta*eta/20;G=4*eta+Q(1,10)
    E=e*exponential(y/2);U=G
    mom=[2*G*exponential(y),2*e*exponential(2*y),2*G*e*exponential(2*y),
         2*G*G*exponential(y)-e*e*exponential(2*y)/2,e*e*exponential(y)/2]
    axis=jet(0,order) if drop_axis else -2+eta/7+eta*eta/11
    coordinates.center=center
    state=coordinates(X,root,E,U,mom,axis)
    return dict(X=X,root=root,E=E,U=U,mom=mom,axis=axis,state=state,eta=eta)


def prescribed(center=Q(1,3),kappa=Q(1,10),order=5,drop_axis=False):
    r=reference(center,order,drop_axis);st=r['state'];y=modules()[0].Jet.variable(0,0,n=2,order=order)
    # Integral jets have zero incoming value but retain all incoming eta jets.
    e=r['E']*exponential(-kappa*primitive(st['p1'])/2)
    u=r['U']-kappa*primitive(r['X']*st['ns'])/2
    incoming=[jet({(0,j):v for (i,j),v in a.terms.items() if i==0},order) for a in r['mom']]
    mom=[a+primitive(v) for a,v in zip(incoming,densities(r['X'],r['root'],e,u))]
    coordinates.center=center
    return dict(reference=r,E=e,U=u,mom=mom,axis=r['axis'],
        state=coordinates(r['X'],r['root'],e,u,mom,r['axis']))


def zero_through(z,order):return not any(v for a,v in z.terms.items() if sum(a)<=order)


def mpq(q):return mp.mpf(q.numerator)/q.denominator


def raw_second(z):
    return modules()[1].SecondJet(*(mpq(z.partial(a)) for a in [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]))


@lru_cache(None)
def stats():
    counts=dict(shear_coefficients=0,moment_coefficients=0,pressure_coefficients=0,old_coordinate_components=0)
    all_flows=True;all_memory=True;all_wrong=True;axis_term=True;maximum=mp.mpf(0)
    with mp.workdps(70):
        for center in (Q(-1),Q(0),Q(1,3),Q(1)):
            r=reference(center);n=reference(center,drop_axis=True);eta=r['eta'];h=Q(1,200)
            expected=4*(Q(1,2)+h)*eta*r['axis']-(1-eta*eta)*r['axis'].derivative(1)
            axis_term &= zero_through(r['state']['Ns']-n['state']['Ns']-expected,3)
            for kappa in (Q(0),Q(1,10),Q(1,2)):
                a=prescribed(center,kappa);r=a['reference']
                slopes=[a['E'].derivative(0)-(1-kappa*r['state']['p1'])*a['E']/2,
                        a['U'].derivative(0)+kappa*r['X']*r['state']['ns']/2]
                all_flows &= all(zero_through(v,3) for v in slopes)
                counts['shear_coefficients']+=2*10 # ten multiindices of total order <=3
                for v,rhs in zip(a['mom'],densities(r['X'],r['root'],a['E'],a['U'])):
                    all_flows &= zero_through(v.derivative(0)-rhs,3);counts['moment_coefficients']+=10
                memory=a['state']['Pi']-a['mom'][4]-a['axis']
                all_memory &= zero_through(memory,3);counts['pressure_coefficients']+=10
                all_wrong &= bool(r['state']['p1'].value()) and zero_through(1-2*r['E'].derivative(0)/r['E'],3)
                if kappa:
                    all_wrong &= bool((a['E'].derivative(0)-r['E'].derivative(0)).value())
                else:all_wrong &= zero_through(a['E']-r['E'],3) and zero_through(a['U']-r['U'],3)
                old=modules()[1].pressure_from_moments(mp.log(2),mpq(center),mpq(Q(1,200)),
                    raw_second(a['E']),raw_second(a['U']),[raw_second(v) for v in a['mom']],raw_second(a['axis']))
                for name in ('W','Pi','Qs','Ns','p1','p2'):
                    for attr,alpha in [('v',(0,0)),('y',(1,0)),('e',(0,1))]:
                        err=abs(getattr(old[name],attr)-mpq(a['state'][name].partial(alpha)))
                        maximum=max(maximum,err);counts['old_coordinate_components']+=1
    return dict(counts=counts,all_flows=bool(all_flows),all_memory=bool(all_memory),
        shortcut_detected=bool(all_wrong),axis_term=bool(axis_term),maximum_old_coordinate_error=maximum)


def diagnostic_review():
    s=stats()
    return dict(fixture='smooth-axis-consistent frozen F_r,U_r manufactured reference at X=2, not actual core',
        centers=['-1','0','1/3','1'],kappa_values=['0','1/10','1/2'],
        rational_input_order=5,checked_flow_total_order=3,counts=dict(s['counts']),
        independent_implementation='compact_jet_envelope/jets.py raw second jets; first coordinate outputs only',
        tolerance='1e-55 at 70 decimal digits',raw_roundoff_scores_tested_not_refrozen=True,
        checks=dict(exact_prescribed_shear_and_five_moment_flows=s['all_flows'],
            full_Pi_minus_Cp_equals_original_eta_dependent_axis_datum=s['all_memory'],
            extending_native_shear_shortcut_past_reference_taper_is_detected=s['shortcut_detected'],
            pressure_datum_changes_reference_Ns_by_full_exact_term=s['axis_term'],
            independent_raw_second_jet_coordinates_agree=s['maximum_old_coordinate_error']<mp.mpf('1e-55')))
