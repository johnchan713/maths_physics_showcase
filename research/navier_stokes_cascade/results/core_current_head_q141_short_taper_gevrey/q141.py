"""Value-based Gevrey-3 recurrence on Q72's actual short taper."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,HERE.parent/path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

previous=load('q141_previous','core_current_head_q140_four_to_six_mixed_binding/q140.py')
exact=previous.exact;order=previous.order;lg=previous.lg
IDENTITY='same-Q72-short-taper-reference-and-actual-maps-with-independent-full-memory'

def env(P,A=1,B=1,stage='mixed-Gevrey3'):
    P,A,B=map(exact,(P,A,B))
    if P<0 or min(A,B)<1:raise ValueError('Nonnegative prefactor and bases at least one required')
    return dict(stage=stage,P=P,A=A,B=B)
def need(r):
    if not isinstance(r,dict) or r.get('stage') not in ('mixed-Gevrey3','constant-in-y-Gevrey3'):raise ValueError('Mixed Gevrey-3 envelope required')
    return r
def add(*rs):
    for r in rs:need(r)
    return env(sum(r['P'] for r in rs),max(r['A'] for r in rs),max(r['B'] for r in rs))
def scale(r,c):
    need(r);c=abs(exact(c));return env(c*r['P'],r['A'],r['B'])
def mul(r,s):
    need(r);need(s);return env(r['P']*s['P'],2*max(r['A'],s['A']),2*max(r['B'],s['B']))
def dy(r):
    need(r);return env(r['P']*r['A'],8*r['A'],r['B'])
def de(r):
    need(r);return env(r['P']*r['B'],r['A'],8*r['B'])
def reciprocal(r,floor):
    need(r);floor=exact(floor)
    if floor<=0:raise ValueError('Positive real reciprocal floor required')
    L=8*(1+2*r['P']/floor)
    return env(2/floor,L*r['A'],L*r['B'])
def integrate(q,trace,length):
    need(q);need(trace);length=exact(length)
    if length<0 or trace.get('stage')!='constant-in-y-Gevrey3':raise ValueError('Fixed nonnegative interval and constant-in-y trace required')
    return env(trace['P']+length*q['P'],q['A'],max(q['B'],trace['B']))
def exp_comp(L,value_cap):
    need(L);value_cap=exact(value_cap)
    if value_cap<1:raise ValueError('Outward actual-value exponential cap required')
    f=8*(1+L['P'])
    return env(value_cap,f*L['A'],f*L['B'])

def seed(L0,U0,h0,b,F0,invF0,identity):
    L0,U0,h0,b,F0,invF0=map(exact,(L0,U0,h0,b,F0,invF0))
    if identity!=IDENTITY or min(L0,U0,h0,b,F0,invF0)<=0:raise ValueError('Same Q72 core polydisks and trace caps required')
    return dict(stage='Q72-analytic-seed',
        log_slope=env(L0/h0,max(Q(1),2/h0),max(Q(1),1/b)),
        U_slope=env(U0/h0,max(Q(1),2/h0),max(Q(1),1/b)),
        F_trace=env(F0,1,max(Q(1),1/b),'constant-in-y-Gevrey3'),
        invF_trace=env(invF0,1,max(Q(1),1/b),'constant-in-y-Gevrey3'),
        U_trace=env(U0,1,max(Q(1),1/b),'constant-in-y-Gevrey3'),b=b)
def switch(t,identity):
    t=exact(t)
    if identity!=IDENTITY or t<=0:raise ValueError('Same positive C_setup^-20 short width required')
    return env(162,max(Q(1),3260/t),1)
def fields(seed_recipe,log_source,U_source,length,exp_cap,identity):
    if identity!=IDENTITY or seed_recipe.get('stage')!='Q72-analytic-seed':raise ValueError('Same Q72 seed required')
    need(log_source);need(U_source)
    zero=env(0,1,max(log_source['B'],U_source['B']),'constant-in-y-Gevrey3')
    L=integrate(log_source,zero,length);eL=exp_comp(L,exp_cap)
    F=mul(seed_recipe['F_trace'],eL);invF=mul(seed_recipe['invF_trace'],eL)
    U=integrate(U_source,seed_recipe['U_trace'],length)
    return dict(stage='Q72-value-based-field-recurrence',F=F,invF=invF,U=U,log_increment=L,
        exp_value_cap=exact(exp_cap),order_dependent_norm_exponentiated=False)

def geometry(Xcap,invXcap,Rcap,invRcap):
    Xcap,invXcap,Rcap,invRcap=map(exact,(Xcap,invXcap,Rcap,invRcap))
    if min(Xcap,invXcap,Rcap,invRcap)<=0:raise ValueError('Positive Q72 geometry caps required')
    return dict(stage='same-Q72-log-geometry',X=env(Xcap),invX=env(invXcap),R=env(Rcap),invR=env(invRcap))
def trace_rows(T,B):
    T,B=map(exact,(T,B))
    if T<=0 or B<1:raise ValueError('Positive common five-trace cap and angular base required')
    r=env(T,1,B,'constant-in-y-Gevrey3')
    return dict(stage='same-five-Q72-traces',**{k:r for k in ('M','I','J','S','Pi')})

def moments(F,U,traces,g,length):
    need(F);need(U)
    if traces.get('stage')!='same-five-Q72-traces' or g.get('stage')!='same-Q72-log-geometry':raise ValueError('Same five traces and geometry required')
    X=g['X'];X2=mul(X,X)
    dens=dict(M=mul(X,U),I=scale(mul(X2,F),2),J=scale(mul(mul(X2,F),U),2),
              S=add(mul(X,mul(U,U)),mul(X2,mul(F,F))),Pi=mul(X,mul(F,F)))
    return {k:integrate(dens[k],traces[k],length) for k in dens}

def angular_factors(h):
    h=exact(h)
    if not 0<h<Q(1,2):raise ValueError('Original h range required')
    one=env(1);eta=env(2,1,2);d=env(2,1,2);L=env(1+2*h,1,2)
    return one,eta,d,reciprocal(L,1-2*h),Q(1,2)+h,Q(1,2)-h

def sources(field,traces,h,g,length,identity):
    if identity!=IDENTITY or field.get('stage')!='Q72-value-based-field-recurrence':raise ValueError('Same Q72 reference field required')
    F,invF,U=field['F'],field['invF'],field['U'];P=moments(F,U,traces,g,length)
    one,eta,d,invL,Aphys,D=angular_factors(h);X,invX=g['X'],g['invX']
    W=add(one,scale(mul(invX,add(scale(mul(eta,P['M']),2*D),mul(d,de(P['M'])))),1))
    NQ=add(scale(P['I'],1-exact(h)),scale(mul(eta,de(P['I'])),D),mul(d,de(P['J'])),scale(mul(eta,P['J']),2*abs(exact(h)-D)))
    p1=mul(invL,add(mul(X,W),scale(mul(mul(invX,invF),NQ),Q(1,2))))
    inner=add(scale(P['M'],D),scale(mul(eta,de(P['M'])),D),scale(mul(eta,P['S']),4*exact(h)),mul(d,de(P['S'])))
    ns=mul(invL,add(mul(W,U),mul(invX,inner),scale(mul(eta,P['Pi']),4*Aphys),mul(d,de(P['Pi']))))
    return dict(stage='same-Q72-reference-full-memory-sources',p1=p1,Xns=mul(X,ns),moments=P,
        all_five_incoming_traces_retained=True,Pi_eta_retained=True)

def state(field,traces,h,g,length,identity):
    if identity!=IDENTITY or field.get('stage')!='Q72-value-based-field-recurrence':raise ValueError('Same Q72 actual field required')
    F,invF,U=field['F'],field['invF'],field['U'];P=moments(F,U,traces,g,length)
    E=mul(g['R'],F);invE=mul(g['invR'],invF)
    one,eta,d,invL,Aphys,D=angular_factors(h);X=g['X']
    a=add(one,scale(mul(dy(E),invE),2));b=scale(mul(dy(U),invE),2)
    w=add(X,scale(mul(eta,P['M']),2*D),mul(d,de(P['M'])))
    num=add(mul(w,U),scale(P['M'],D),scale(mul(eta,de(P['M'])),D),scale(mul(eta,P['S']),4*exact(h)),mul(d,de(P['S'])),
            scale(mul(X,mul(eta,P['Pi'])),4*Aphys),mul(X,mul(d,de(P['Pi']))))
    p=mul(mul(num,invL),invE)
    return dict(stage='same-Q72-actual-four-PREMOD-mixed-fields',a=a,b=b,p=p,E=E,invE=invE,U=U,moments=P,
        H=a['P']+b['P']+p['P']+E['P'],Ay=max(a['A'],b['A'],p['A'],E['A']),Be=max(a['B'],b['B'],p['B'],E['B']),
        full_p2_including_Pi_and_Pi_eta=True,actual_value_floors_and_outward_exp_caps_certified=False)

def short_taper(seed_recipe,step,traces,g,h,t,reference_exp_cap,actual_exp_cap,identity):
    if identity!=IDENTITY or step.get('stage')!='mixed-Gevrey3':raise ValueError('Same Q72 short-taper data required')
    t=exact(t);length=2*t
    ref=fields(seed_recipe,mul(step,seed_recipe['log_slope']),mul(step,seed_recipe['U_slope']),length,reference_exp_cap,identity)
    src=sources(ref,traces,h,g,length,identity)
    actual=fields(seed_recipe,scale(mul(step,src['p1']),Q(1,2)),scale(mul(step,src['Xns']),Q(1,2)),length,actual_exp_cap,identity)
    out=state(actual,traces,h,g,length,identity)
    return dict(stage='Q72-closed-value-based-short-taper-Gevrey3',reference=ref,sources=src,actual=actual,state=out,
        inherited_exp_c_Tm_cap_used=False,actual_order_zero_caps_and_floors_certified=False)

def bind(recipe,afloor,delta,G,Eexp,d0,disk,V0,N,identity):
    if identity!=IDENTITY or recipe.get('stage')!='Q72-closed-value-based-short-taper-Gevrey3':raise ValueError('Same closed Q72 short taper required')
    s=recipe['state'];four=previous.four_to_six(s['H'],s['Ay'],s['Be'],s['a']['P'],s['b']['P'],afloor,delta,G,Eexp,d0,previous.IDENTITY)
    joined=previous.bind(four,disk,V0,N,previous.IDENTITY)
    return dict(stage='conditional-Q72-short-taper-to-Ny-pullback',four_to_six=four,joined=joined,
        actual_order_zero_caps_floors_frequency_and_phase_gates_certified=False)

def bound(r,i,k,binary=False):
    need(r);order(i);order(k)
    if binary:
        v=lg(max(Q(1),r['P']))+i*lg(r['A'])+k*lg(r['B'])
        v+=3*i*lg(Q(i)) if i else 0;v+=3*k*lg(Q(k)) if k else 0
    else:
        if max(i,k)>12:raise ValueError('Use binary budget at high order')
        v=r['P']*r['A']**i*r['B']**k*(factorial(i)*factorial(k))**3
    return dict(y_order=i,eta_order=k,upper=v,binary=binary,continuum_input_gates_certified=False)

def record():
    return dict(status='Q72-value-based-short-taper-mixed-Gevrey3-recurrence',
        correction='integrate fixed Gevrey sources and compose exp at an outward order-zero value cap; never exponentiate T_m',
        reference='same analytic core slopes and same psi build F_r,U_r; all five Xa traces build p1r and X*nsr',
        actual='same kappa,p1r,X*nsr build F_a,U_a; moments are reintegrated independently from the original Xa traces',
        pressure='p2 retains complete M,S,Pi memory and both 4*Aphys*eta*Pi and -(1-eta^2)*Pi_eta terms; p1 retains I,J',
        output='conditional common mixed Gevrey-3 H,Ay,Be for actual short-taper a,b,p2,E, passed to Q140/Q139/Q138',
        premises='Q72 order-zero geometry, trace, exponential value and real floor/value/frequency gates remain to be numerically bound on the selected profile',
        next='verify those order-zero caps, assemble the remaining regional four-field constants, then join Q137/Q121 on the whole hull',
        actual_complete_profile_pulse_covariance_residual_or_blowup_certified=False,historical_frozen_records_modified=False)
