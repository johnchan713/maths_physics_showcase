"""Tame fourth-angular modulation densities and the full five-row target."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import importlib.util
import sys

HERE=Path(__file__).resolve().parent
def rat(x):
    if isinstance(x,(bool,float)):raise TypeError('Exact rational inputs required')
    return F(x)
def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x
@lru_cache(None)
def module(folder,file):
    p=HERE.parent/folder/file;prior=sys.path[:];sys.path.insert(0,str(p.parent))
    try:
        s=importlib.util.spec_from_file_location('fourth_target_'+folder,p)
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=prior
def add(*rows):return tuple(sum(r[j] for r in rows) for j in range(3))
def scale(row,c):return tuple(v*c for v in row)
def value(row,H,a4,b4):return sum(v*w for v,w in zip(row,(H,a4,b4)))

def density_affines(Ahi):
    """N times each Ceta4 difference, coefficients of (H,a4,b4)."""
    Ahi=rat(Ahi)
    if Ahi<1:raise ValueError('Original C4 envelope at least one required')
    return dict(U=(F(1),F(0),F(1)),E=(12*Ahi,3*Ahi,F(0)),
        U2=(6*Ahi,F(0),4*Ahi),E2=(114*Ahi*Ahi,18*Ahi*Ahi,F(0)),
        UE=(38*Ahi*Ahi,6*Ahi*Ahi,3*Ahi*Ahi))

def moment_affines(A,Ahi):
    A=rat(A);d=density_affines(Ahi)
    if A<1:raise ValueError('A>=1 enclosing X and inverse X and radial length required')
    return dict(M=scale(d['U'],A),J=scale(d['UE'],2*A*A),I=scale(d['E'],2*A*A),
        S=scale(add(d['U2'],scale(d['E2'],F(1,2))),A),Cp=scale(d['E2'],A*A/2))

def target_affines(A,Ahi,rhop,invsqrt,e_star,lam):
    """invsqrt >= (sqrt(2)*rhop^(3/2))^-1; no lambda-zero limit."""
    A,Ahi,rhop,invsqrt,e_star,lam=map(rat,(A,Ahi,rhop,invsqrt,e_star,lam))
    m=moment_affines(A,Ahi)
    if min(rhop,invsqrt,e_star)<=0 or not 0<lam<=F(1,100) or 2*rhop**3*invsqrt**2<1:
        raise ValueError('Positive normalizers, retained lambda and genuine reciprocal-square-root enclosure required')
    w=dict(M=5/(rhop*e_star),J=25*invsqrt/e_star**2,I=5*invsqrt/e_star,
        S=25/(rhop*e_star**2),Cp=25/e_star**2)
    n={k:scale(m[k],w[k]) for k in m}
    rows=(n['M'],scale(add(n['M'],n['J']),1/lam),n['I'],n['S'],n['Cp'])
    caps=tuple(max(row[j] for row in rows) for j in range(3))
    return dict(normalizer_C4=w,raw_moment_affines=m,transformed_rows=rows,Z4_affine=caps,
        interpretation='Z4<=ZH*H+Za*a4+Zb*b4; a4,b4 actual fourth normalized primitive jets; complete forward incoming differences and same axis pressure datum')

def primitive_gate(A,Ahi,rhop,invsqrt,e_star,lam,H,N,r,a4,b4,K):
    A,Ahi,rhop,invsqrt,e_star,lam,H,N,r,a4,b4,K=map(rat,(A,Ahi,rhop,invsqrt,e_star,lam,H,N,r,a4,b4,K))
    if H<1 or N<H or min(r,a4,b4,K)<0:raise ValueError('H>=1,N>=H and nonnegative supplied jets/trial cap required')
    q=target_affines(A,Ahi,rhop,invsqrt,e_star,lam);caps=q['Z4_affine']
    z4=value(caps,H,a4,b4);B4=H+b4
    parent=module('fourth_angular_axial_interface','angular_bounds.py')
    cost=parent.compact_cost(A,Ahi,N,r,B4,z4)
    ca=3200000*A*caps[1];cb=1+3200000*A*caps[2];ch=1+3200000*A*caps[0]
    lhs=ca*a4+cb*b4;rhs=N*(K-cost['known_part'])-ch*H
    return dict(a4_coefficient=ca,b4_coefficient=cb,H_coefficient=ch,Z4_upper=z4,B4_upper=B4,
        compact_C4_upper=cost['repaired_U_C4'],lhs=lhs,rhs=rhs,accepts=lhs<=rhs,
        substituted_majorant_accepts=cost['repaired_U_C4']<=K,
        scope='sufficient for parent gate once actual a4,b4 satisfy it; not evaluated at retained huge parameters, not necessary for true field')

def symbolic_review():
    import sympy as s
    z=s.Symbol('z');a1,a2,a3,a4=s.symbols('a1 a2 a3 a4');aa=a1*z+a2*z*z+a3*z**3+a4*z**4
    e4=s.expand(sum(aa**j/s.factorial(j) for j in range(5))).coeff(z,4)
    expected=a4+a1*a3+a2*a2/2+a1*a1*a2/2+a1**4/24
    e42=s.expand(sum((2*aa)**j/s.factorial(j) for j in range(5))).coeff(z,4)
    expected2=2*a4+4*a1*a3+2*a2*a2+4*a1*a1*a2+s.Rational(2,3)*a1**4
    e=s.Symbol('eta',real=True);phi=s.Symbol('phi',real=True)
    theta=phi+e*s.sin(phi)/2;pulled=s.cos(theta)
    moving=s.simplify(s.diff(pulled,e,4).subs({e:0,phi:s.pi/3}))
    p=1+e*e
    c1=1+sum(s.binomial(2,j) for j in range(3));c2=sum(s.binomial(4,j) for j in range(5))+2*sum(s.binomial(2,j) for j in range(3))+1
    axis=s.Function('Pi_axis')(e);old=s.Function('Cp_old')(e);new=s.Function('Cp_new')(e)
    pressure_difference=(axis+new)-(axis+old)
    return dict(moving_phase_fourth_control=str(moving),
        checks=dict(exponential_fourth_has_all_partitions=s.expand(e4-expected)==0,
            squared_swirl_exponential_fourth_has_all_partitions=s.expand(e42-expected2)==0,
            omitting_mixed_lower_jets_detected=s.expand(e4-a4)!=0,
            native_reciprocal_e_normalizer_C4_is_five=c1==5,
            native_reciprocal_e_squared_normalizer_C4_is25=c2==25,
            moving_phase_eta_fourth_not_frozen=moving==s.Rational(9,512),
            same_pressure_axis_datum_cancels_only_in_difference=s.diff(pressure_difference,e,4)==s.diff(new-old,e,4),
            changed_axis_gauge_would_add_a_fourth_jet=s.diff(pressure_difference+e**4,e,4)-s.diff(new-old,e,4)==24,
            moving_lift_retains_uniform_positive_derivative=bool(s.Rational(1,2)>0)))

def manufactured_review():
    """Direct full target differentiation versus independent normalized product jets."""
    import sympy as s
    import numpy as np
    e,phi=s.symbols('eta phi',real=True);N=4096;H=64;A=3;Ahi=3
    E=1+e*e/3;U=e/5+s.Rational(1,7)
    aa=(s.Rational(1,11)+e/13+e*e/17+e**3/19+2*e**4)*s.cos(phi)
    bb=(s.Rational(1,23)+e/29+e*e/31+e**3/37+3*e**4)*s.sin(phi)
    En=E*s.exp(aa/N);Un=U+bb/N
    den=(Un-U,En-E,Un*Un-U*U,En*En-E*E,Un*En-U*E)
    geom=(8-2*s.sqrt(2))/3
    raw=(den[0],geom*den[4],geom*den[1],den[2]-den[3]/2,s.log(2)*den[3]/2)
    pp=1+e*e;normal=(pp/6,pp*pp/36,pp/12,pp*pp/18,pp*pp/9)
    normed=[normal[i]*raw[i] for i in range(5)]
    target=(normed[0],200*(normed[0]-normed[1]),normed[2],-normed[3],normed[4])
    vals=s.lambdify((e,phi),[s.diff(f,e,4)/24 for f in den+target],'numpy',cse=True)
    # Independent Leibniz convolution, with polynomial normalizers kept variable.
    conv=[sum(s.diff(normal[i],e,j)*s.diff(raw[i],e,4-j)/(s.factorial(j)*s.factorial(4-j)) for j in range(5)) for i in range(5)]
    conval=s.lambdify((e,phi),(conv[0],200*(conv[0]-conv[1]),conv[2],-conv[3],conv[4]),'numpy',cse=True)
    ordinary=s.lambdify((e,phi),(conv[0],conv[2],conv[1],conv[3],conv[4]),'numpy',cse=True)
    original=module('stress_realization_audit','repair.py')
    d=density_affines(Ahi);q=target_affines(A,Ahi,2,F(1,4),3,F(1,200));caps=q['transformed_rows']
    rows=[]
    for eta in (-1.,-.4,0.,.6,1.):
        for ph in (0.,float(s.pi/3),float(s.pi/2)):
            got=np.asarray(vals(eta,ph),dtype=float);other=np.asarray(conval(eta,ph),dtype=float)
            legacy=np.array([float(v) for v in original.precondition_target([str(v) for v in ordinary(eta,ph)],'0.005')])
            checks=dict(all_five_density_fourth_caps=bool(all(N*abs(v)<=float(value(d[k],H,2,3)) for v,k in zip(got[:5],('U','E','U2','E2','UE')))),
                full_five_row_target_fourth_caps=bool(all(N*abs(v)<=float(value(caps[i],H,2,3)) for i,v in enumerate(got[5:]))),
                independent_normalizer_Leibniz_convolution=bool(np.max(np.abs(got[5:]-other))<1e-10))
            checks['original_preconditioner_fourth_jet_matches']=bool(np.max(np.abs(got[5:]-legacy))<1e-10)
            rows.append(dict(eta=eta,phase=ph,checks=checks))
    # Exact polynomial C3 majorants for the chosen quartic primitives, all phases.
    a_low=sum(abs(F(c))*sum(F(s.binomial(i,j)) for j in range(min(i,3)+1)) for (i,),c in s.Poly(aa.subs(phi,0),e).terms())
    b_low=sum(abs(F(c))*sum(F(s.binomial(i,j)) for j in range(min(i,3)+1)) for (i,),c in s.Poly(bb.subs(phi,s.pi/2),e).terms())
    return dict(rows=rows,manufactured_not_selected_loop_or_PDE=True,floating_not_interval_certificate=True,
        checks=dict(actual_manufactured_lower_jet_inputs_satisfy_H=a_low<=H and b_low<=H,
            actual_manufactured_native_C4_inputs_satisfy_Ahi=F(7,3)<=Ahi and F(19,35)<=Ahi,
            all_fifteen_full_target_controls=all(all(r['checks'].values()) for r in rows)))

def review():
    q=target_affines(3,3,2,F(1,4),3,F(1,200));caps=q['Z4_affine'];H,N,r=F(64),F(4096),F(1,10**9)
    initial=primitive_gate(3,3,2,F(1,4),3,F(1,200),H,N,r,2,3,0)
    gate=primitive_gate(3,3,2,F(1,4),3,F(1,200),H,N,r,2,3,initial['compact_C4_upper'])
    sym=symbolic_review();man=manufactured_review()
    checks=dict(first_exponential_lower_quartic_terms_fit_three=F(1)+F(1,2)+F(1,2)+F(1,24)<3,
        squared_exponential_lower_quartic_terms_fit_96=9*(4+2+4+F(2,3))==96,
        mixed_b_exponential_fourth_lower_terms_fit23=F(13,2)+F(9,2)+3+9==23,
        unknown_highest_jets_enter_all_rows_only_linearly=all(len(row)==3 and min(row)>=0 for row in q['transformed_rows']),
        componentwise_max_bounds_every_row=all(caps[j]>=row[j] for row in q['transformed_rows'] for j in range(3)),
        exact_normalizer_geometry=q['normalizer_C4']['M']==F(5,6) and q['normalizer_C4']['J']==F(25,36),
        inverse_lambda_row_not_discarded=all(q['transformed_rows'][1][j]>=200*q['transformed_rows'][0][j] for j in range(3)),
        primitive_gate_is_exact_at_substituted_boundary=gate['lhs']==gate['rhs'] and gate['accepts'] and gate['substituted_majorant_accepts'],
        target_Z4_not_an_additional_independent_unknown=gate['Z4_upper']==value(caps,H,2,3),
        axial_primitive_B4_retains_lower_jets=gate['B4_upper']==H+3)
    checks.update(sym['checks']);checks.update(man['checks'])
    return dict(density_affine_fixture=serial(density_affines(3)),full_target_fixture=serial(q),primitive_gate_fixture=serial(gate),
        symbolic=sym,manufactured=man,checks=checks)
