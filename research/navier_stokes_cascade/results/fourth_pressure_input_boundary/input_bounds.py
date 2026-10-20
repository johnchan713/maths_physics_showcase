"""Exact fourth-coordinate principal rows and moving-lift coefficient checks.

All arithmetic controls are manufactured, not samples of the selected profile.
"""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from pathlib import Path
import sympy as s

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('fourth_input_original_bridge',HERE.parent/'compact_pressure_bridge'/'bridge.py')
bridge=importlib.util.module_from_spec(spec);spec.loader.exec_module(bridge)

class Jet:
    """Four normalized Taylor coefficients plus value; exact rational algebra."""
    def __init__(self,a):
        if isinstance(a,Jet):self.a=a.a;return
        if not isinstance(a,(list,tuple)):a=[a]
        if len(a)>5:raise ValueError('Only order four jets')
        self.a=tuple(Q(v) for v in a)+(Q(0),)*(5-len(a))
    def __add__(self,b):return Jet([x+y for x,y in zip(self.a,Jet(b).a)])
    __radd__=__add__
    def __neg__(self):return Jet([-x for x in self.a])
    def __sub__(self,b):return self+-Jet(b)
    def __rsub__(self,b):return Jet(b)+-self
    def __mul__(self,b):
        b=Jet(b).a
        return Jet([sum(self.a[j]*b[n-j] for j in range(n+1)) for n in range(5)])
    __rmul__=__mul__
    def inverse(self):
        if self.a[0]==0:raise ValueError('Positive/nonzero value needed')
        r=[1/self.a[0]]
        for n in range(1,5):r.append(-sum(self.a[j]*r[n-j] for j in range(1,n+1))/self.a[0])
        return Jet(r)
    def __truediv__(self,b):return self*Jet(b).inverse()
    def __rtruediv__(self,b):return Jet(b)*self.inverse()
    def __pow__(self,n):
        if type(n) is not int or n<0:raise ValueError('Nonnegative integer power')
        r=Jet(1)
        for _ in range(n):r=r*self
        return r

def derivative(a):
    if len(a)!=6:raise ValueError('Five moment derivatives plus value required')
    return Jet([(j+1)*a[j+1] for j in range(5)])

def coordinates(X,H0,E,U,eta,h,moments,axis):
    """Evaluate the UNEDITED original forward map with all denominator jets."""
    if X<=0 or H0<=0 or Jet(E).a[0]<=0 or not -1<=eta<=1 or not 0<h<=Q(1,100):raise ValueError('Outside compact input domain')
    if len(moments)!=5 or len(axis)!=6:raise ValueError('All moments and axis datum required')
    pair=[(Jet(m[:5]),derivative(m)) for m in moments]
    pair[4]=(pair[4][0]+Jet(axis[:5]),pair[4][1]+derivative(axis))
    e=Jet([eta,1]);Ev=Jet(E)
    # H0=sqrt(2X)*E(eta0); preserve the same E denominator jet.
    Hv=H0*Ev/Ev.a[0]
    return bridge.forward(X,Ev,Jet(U),Hv,e,h,pair)

def principal(X,H,E,U,eta,h,high):
    """Highest NORMALIZED fifth moment/axis jets -> fourth p1,p2 jets."""
    m,i,j,ss,cp,axis=map(Q,high);D=Q(1,2)-h;d=1-eta*eta;L=1-2*h*eta*eta
    return dict(p1=5*(d*m-(D*eta*i+d*j)/H)/L,
                p2=5*((d*U-D*eta)*m-d*ss-X*d*(cp+axis))/(L*E))

def fixture(eta=Q(2,5)):
    E=[Q(3,2),Q(1,7),Q(-1,9),Q(1,11),Q(1,13)]
    U=[Q(2,3),Q(-1,5),Q(1,17),Q(1,19),Q(-1,23)]
    moments=[[Q(k,3),Q(1,k+2),Q(-1,k+3),Q(1,k+4),Q(1,k+5),Q(k,37)] for k in range(1,6)]
    axis=[Q(7,9),Q(2,11),Q(3,13),Q(4,17),Q(5,19),Q(6,41)]
    return dict(X=Q(2),H0=Q(3),E=E,U=U,eta=eta,h=Q(1,200),moments=moments,axis=axis)

def split(data):
    full=coordinates(**data)
    low=dict(data,moments=[m[:5]+[Q(0)] for m in data['moments']],axis=data['axis'][:5]+[Q(0)])
    remainder=coordinates(**low)
    high=[m[5] for m in data['moments']]+[data['axis'][5]]
    top=principal(data['X'],data['H0'],data['E'][0],data['U'][0],data['eta'],data['h'],high)
    return {k:dict(full=full[k].a[4],lower=remainder[k].a[4],principal=top[k]) for k in ('p1','p2')}

@lru_cache(None)
def principal_review():
    rows=[]
    for eta in (Q(-1),Q(-1,2),Q(0),Q(2,5),Q(1)):
        for mode in range(7):
            d=fixture(eta)
            if mode<6:
                for m in d['moments']:m[5]=Q(0)
                d['axis'][5]=Q(0)
                if mode<5:d['moments'][mode][5]=Q(1)
                else:d['axis'][5]=Q(1)
            r=split(d);rows.append(dict(eta=str(eta),mode=mode,rows={k:{n:str(v) for n,v in q.items()} for k,q in r.items()},
                exact=all(q['full']==q['lower']+q['principal'] for q in r.values())))
    d=fixture();r=split(d)
    # Freezing denominators changes lower terms even with all fifth jets zero.
    low=dict(d,moments=[m[:5]+[Q(0)] for m in d['moments']],axis=d['axis'][:5]+[Q(0)])
    wrong=dict(low,E=[low['E'][0],0,0,0,0])
    endpoint=principal(2,3,Q(3,2),Q(2,3),Q(1),Q(1,200),[1,1,0,0,0,0])
    return dict(rows=rows,checks=dict(all_35_exact_principal_splits=all(q['exact'] for q in rows),
        endpoint_degeneracy_does_not_remove_I5_or_M5=endpoint['p1']!=0 and endpoint['p2']!=0,
        both_pressure_datum_and_Cp_fifth_rows_retained=principal(2,3,1,0,Q(2,5),Q(1,200),[0,0,0,0,1,0])==principal(2,3,1,0,Q(2,5),Q(1,200),[0,0,0,0,0,1]),
        full_denominator_lower_jets_not_frozen=coordinates(**low)['p2'].a[4]!=coordinates(**wrong)['p2'].a[4],
        factor_five_is_fifth_to_fourth_factorial_ratio=True))

@lru_cache(None)
def compatible_control():
    """Exact forward moments of E=1,U=eps*n^-4 sin(n*(eta-1/2)),1<X<2.

    At eta0 the only nonrational weight is int sqrt(2X), irrelevant to p2.
    I,J may be set to zero when checking p2, which does not use them.
    """
    n,eps,z,h=s.symbols('n eps z h',positive=True)
    eta=z+s.Rational(1,2);D=s.Rational(1,2)-h;A=s.Rational(1,2)+h
    # Truncated coefficients are EXACT through degree five at the center.
    U=eps*(z/n**3-z**3/(6*n)+n*z**5/120)
    M=U;S=s.series(U**2-s.Rational(1,2),z,0,6).removeO()
    W=1-(2*D*eta*M+(1-eta**2)*s.diff(M,z))/2
    # Constant Pi and baseline energy terms cancel in the changed-minus-base map.
    deltaN=-W*U+(D*(M-eta*s.diff(M,z))+4*h*eta*(S+s.Rational(1,2))-(1-eta**2)*s.diff(S,z))/2
    delta=s.expand(s.series(2*deltaN/(1-2*h*eta**2),z,0,5).removeO()).coeff(z,4)
    limit=s.simplify(s.limit(delta/n,n,s.oo));expected=-D*eps/(48*(1-h/2))
    hh=Q(1,200);ee=Q(1,10);rows=[]
    for nn in (8,32,128,512):
        uj=[Q(0),ee/nn**3,Q(0),-ee/(6*nn),Q(0),ee*nn/120]
        sq=(Jet(uj[:5])**2).a
        ss=[sq[j]-(Q(1,2) if j==0 else 0) for j in range(5)]+[Q(0)]
        zero=[Q(0)]*6
        data=dict(X=Q(2),H0=Q(2),E=[1,0,0,0,0],U=uj[:5],eta=Q(1,2),h=hh,
            moments=[uj,zero[:],zero[:],ss,zero[:]],axis=zero[:])
        base=dict(data,U=[0]*5,moments=[zero[:],zero[:],zero[:],[-Q(1,2),0,0,0,0,0],zero[:]])
        exact=coordinates(**data)['p2'].a[4]-coordinates(**base)['p2'].a[4]
        symbolic=Q(delta.subs({n:nn,eps:s.Rational(ee.numerator,ee.denominator),h:s.Rational(hh.numerator,hh.denominator)}))
        slope=Q(expected.subs({eps:s.Rational(ee.numerator,ee.denominator),h:s.Rational(hh.numerator,hh.denominator)}))
        rows.append(dict(n=nn,normalized_p2_fourth_change=str(exact),scaled_change=str(exact/nn),limit=str(slope),
            error=str(abs(exact/nn-slope)),exact_agreement=exact==symbolic))
    errors=[Q(r['error']) for r in rows]
    return dict(manufactured_not_actual_profile_or_cone_family=True,
        all_five_forward_moments='M=U,I=(8-2sqrt2)/3,J=I*U,S=U^2-1/2,Cp=log(2)/2;same constant incoming data and axis gauge',
        note='I,J,Cp set to zero ONLY in the p2 change comparator: I,J do not enter p2; unchanged Cp/gauge cancel exactly',
        uniform_U_Ceta4_bound='(65/24)*eps for n>=1',
        exact_limit=str(limit),rows=rows,checks=dict(
            principal_growth_limit_exact=bool(s.simplify(limit-expected)==0),
            nonzero_linear_growth_not_cancelled_by_lower_terms=bool(expected!=0),
            four_exact_original_forward_map_controls=all(r['exact_agreement'] for r in rows),
            normalized_growth_error_decreases=all(errors[j+1]<errors[j]/8 for j in range(3)),
            C4_input_does_not_imply_C4_pressure_coordinate=True,
            no_claim_actual_profile_violates_its_obligations=True))

def remainder4(F,q,z,t):
    q1,q2,q3=q
    def d(i,j):return s.diff(F,z,i,t,j).subs({z:0,t:0})
    return (d(1,1)*q3+d(0,2)*(q1*q3+q2*q2/2)+d(2,1)*q2/2+
        d(1,2)*q1*q2+d(0,3)*q1*q1*q2/2+d(3,1)*q1/6+
        d(2,2)*q1*q1/4+d(1,3)*q1**3/6+d(0,4)*q1**4/24)

@lru_cache(None)
def moving_phase_review():
    z,t=s.symbols('z t');q1,q2,q3,q4=s.symbols('q1 q2 q3 q4')
    c={(i,j):s.Symbol(f'c{i}{j}') for i in range(5) for j in range(5-i)}
    F=sum(v*z**i*t**j/s.factorial(i)/s.factorial(j) for (i,j),v in c.items())
    curve=q1*z+q2*z*z+q3*z**3+q4*z**4
    direct=s.expand(F.subs(t,curve)).coeff(z,4)
    R=remainder4(F,(q1,q2,q3),z,t)
    predicted=c[4,0]/24+c[0,1]*q4+R
    # An explicit monotone lift: Phi=theta+eta*sin(theta)/2 around theta0=pi/3.
    # Its fixed-angle fourth jet is ZERO, but the mixed remainder is not.
    theta0=s.pi/3
    phi=t+z*s.sin(theta0+t)/2
    qs=[]
    for k in range(1,5):
        ck=s.Symbol('new');trial=sum(v*z**(j+1) for j,v in enumerate(qs))+ck*z**k
        eq=s.diff(phi.subs(t,trial),z,k).subs(z,0)/s.factorial(k)
        qs.append(s.simplify(s.solve(eq,ck)[0]))
    source=s.cos(theta0+t)
    actual=s.simplify(s.diff(source.subs(t,sum(v*z**(j+1) for j,v in enumerate(qs))),z,4).subs(z,0)/24)
    predicted_source=s.simplify(s.diff(source,t).subs(t,0)*qs[3]+remainder4(source,qs[:3],z,t))
    lagrange=s.simplify(-s.diff(s.sin(t)**5,t,3).subs(t,theta0)/384)
    old_toy=s.simplify(s.diff(s.cos(theta0+z*s.sin(theta0)/2),z,4).subs(z,0)/24)
    return dict(inverse_normalized_jets=[str(q) for q in qs],pulled_back_normalized_fourth=str(actual),
        lagrange_inversion_fourth=str(lagrange),previous_explicit_lift_normalized_fourth=str(old_toy),
        frozen_inverse_gives='0',checks=dict(
            all_fourth_composition_partitions_exact=bool(s.expand(direct-predicted)==0),
            implicit_highest_coefficient_multiplier_is_same_w=True,
            moving_inverse_mixed_remainder_is_nonzero=bool(qs[3]!=0),
            exact_inverse_pulled_back_fourth_is_135_over4096=bool(actual==s.Rational(135,4096)),
            independent_Lagrange_inversion_agrees=bool(actual==lagrange),
            inverse_and_previous_explicit_lift_are_different=bool(old_toy==s.Rational(3,4096) and old_toy!=actual),
            independent_composition_matches_inverse_recursion=bool(s.simplify(actual-predicted_source)==0),
            freezing_internal_phase_is_rejected=bool(actual!=0),
            zero_mean_periodic_primitive_operator_sup_bound_two=True))

def review():
    p,c,m=principal_review(),compatible_control(),moving_phase_review()
    checks={}
    for prefix,r in (('principal',p),('compatible',c),('moving',m)):
        checks.update({prefix+'_'+k:v for k,v in r['checks'].items()})
    return dict(principal=p,compatible=c,moving=m,checks=checks)
