"""Explicit mixed shape premise and conditional same-root cone-state bridge."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util,sys
import mpmath as mp
import sympy as sp
HERE=Path(__file__).resolve().parent
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
high=module('mixed_shape_native_bounds',HERE.parent/'higher_jet_budget/higher_bounds.py')
fourth=module('mixed_shape_fourth_scalar',HERE.parent/'fourth_loop_shape_gate/loop_fourth.py')
repair=module('mixed_shape_original_repair',HERE.parent/'repair_state_bounds/bounds.py')

def shape_majorant(A,Z3,delta_inverse):
    if type(A) is not int or A<2**16:raise ValueError('Exact retained A>=2^16')
    if type(Z3) not in (int,Q) or type(delta_inverse) not in (int,Q) or min(Z3,delta_inverse)<1:raise ValueError('Exact positive shape and inverse width bounds')
    B=A**256;di=Q(delta_inverse);z=Q(Z3)
    zeta=729*z*B**3*di**3;sqrtw=250*B**3*di**3;ainv=4*A**4*B**3
    return dict(raw_product=zeta*sqrtw*ainv,majorant=10**6*z*A**2308*di**6,
        checks=dict(constants_fit=729*250*4<10**6,native_power_exact=4+9*256==2308,
            width_six_kept=3+3==6))

def shape_gate(A,logZ,order):
    if type(A) is not int or A<2**16 or type(logZ) not in (int,Q) or logZ<0 or type(order) is not int or order not in (3,4):raise ValueError('Exact auxiliary A,logZ and order 3 or 4')
    with mp.workdps(90):
        a=mp.mpf(A);z=mp.mpf(logZ.numerator)/logZ.denominator if isinstance(logZ,Q) else mp.mpf(logZ)
        cost=z+(6*mp.log(10)+2308*mp.log(a)+6*a**20 if order==3 else 12*mp.log(10)+20673*mp.log(a)+8*a**20)
        return dict(auxiliary_sufficient_comparison=bool(cost<=a**32),actual_Z_enclosed=False,
            actual_S_or_compact_cone_accepted=False)

def exponent_review():
    native=high.compact_exponents();mu=fourth.implicit_budget(64,10,[80,240,400,560])[:3]
    theta=fourth.implicit_budget(10010,1,[11000,33000,55000,77000])[:3]
    colors=[r['lower_composition_terms']*3**r['order'] for r in theta]
    checks=dict(native_mixed_inputs_fit_Bn=all(v<4096*256 for v in native['exponents'].values()),
        native_mixed_ledger_checks=all(native['checks'].values()),
        mu_joint_three_orders_have_headroom=all(r['accepted'] for r in mu),
        theta_three_external_colors_fit_absorbed_count=max(colors)<10000,
        theta_joint_three_orders_have_headroom=all(r['accepted'] for r in theta),
        fixed_angle_joint3_t_cost=5+3*(600+1+1)+600+1<4000,
        full_three_variable_moving_cost=20000+3*(80000+1)<350000,
        full_phase_primitive_operator_three_fits=350000+1<400000,
        stronger_shape4_implies_mixed3_constant_gap=12-6==6 and 20673-2308==18365 and 8-6==2)
    return dict(native_exponents=native['exponents'],mu=mu,theta=theta,color_counts=colors,
        primitive_joint3_exponent=400000,primitive_joint_operator_bound=3,
        checks={k:bool(v) for k,v in checks.items()})

def universal_review(digits):
    if type(digits) is not int or digits<60:raise ValueError('At least 60 interval digits')
    mp.iv.dps=digits;iv=mp.iv.mpf;log=mp.iv.log
    lo=lambda x:mp.make_mpf(x._mpi_[0]);hi=lambda x:mp.make_mpf(x._mpi_[1])
    a=g=iv(2**16)
    original=repair.scalar_checks(digits)
    checks=dict(native_Bn_below_G=hi(256*log(a)-a**32)<0,
        native_comparison_increases_for_all_A=lo(32-1/log(a))>0,
        colored_count_absorbed_without_materializing_expG=hi(log(10000))<lo(g),
        full_primitive_exp400000G_below_H=lo(g**3)>400000,
        state_prefactor_constant=hi(log(4002))<lo(g),
        state_threshold_193G_below16G4=lo(16*g**3)>193,
        cone_epsilon_tolerance=lo(g**30)>512*4**10,
        joint_shape_gate_strictly_weaker_than_sufficient_shape4=lo(6*log(10)+18365*log(a)+2*a**20)>0,
        original_repair_neighborhood_inverse_and_threshold_checks=all(original['checks'].values()))
    return dict(interval_digits=digits,checks={k:bool(v) for k,v in checks.items()},
        actual_G_H_N_or_cutoff_materialized=False)

def add(*polys):
    out={}
    for p in polys:
        for k,v in p.items():out[k]=out.get(k,0)+v
    return {k:sp.expand(v) for k,v in out.items() if v!=0}
def scale(p,c):return {k:sp.expand(c*v) for k,v in p.items()}
def mul(p,q,degree=3):
    out={}
    for a,v in p.items():
        for b,w in q.items():
            k=tuple(x+y for x,y in zip(a,b))
            if sum(k)<=degree:out[k]=out.get(k,0)+v*w
    return {k:sp.expand(v) for k,v in out.items() if v!=0}
def power(p,n):
    out={(0,0,0):sp.Integer(1)}
    for _ in range(n):out=mul(out,p)
    return out
def collapse(p):return [sum(abs(v) for k,v in p.items() if sum(k)==n) for n in range(4)]

@lru_cache(None)
def colored_review():
    zero=(0,0,0);y={(1,0,0):sp.Integer(1)};e={(0,1,0):sp.Integer(1)};z={(0,0,1):sp.Integer(1)}
    x=sp.Symbol('x');theta0=sp.pi/3
    def trig(name,delta,mult=1):
        fn=sp.sin(mult*x) if name=='sin' else sp.cos(mult*x)
        return add(*[scale(power(delta,j),sp.diff(fn,x,j).subs(x,theta0)/factorial(j)) for j in range(4)])
    def residual(delta):return add(delta,scale(mul(y,trig('sin',delta)),sp.Rational(1,7)),
        scale(mul(e,trig('cos',delta)),sp.Rational(1,9)),
        scale(mul(mul(y,e),trig('sin',delta,2)),sp.Rational(1,11)),scale(z,-1))
    delta={}
    for n in range(1,4):
        res=residual(delta);new={k:-sp.simplify(v) for k,v in res.items() if sum(k)==n}
        delta=add(delta,new)
    res=residual(delta);source=trig('cos',delta);mixed=sp.simplify(source.get((1,1,1),0))
    P={(0,0,0):sp.Rational(1,3),(1,0,0):-2,(0,1,1):3,(2,0,0):-5}
    Qp={(0,0,0):-2,(0,0,1):sp.Rational(1,7),(1,0,0):-6,(1,1,0):-4}
    pc,qc,actual=collapse(P),collapse(Qp),collapse(mul(P,Qp))
    bounds=[sum(pc[j]*qc[n-j] for j in range(n+1)) for n in range(4)]
    checks=dict(all_three_variable_inverse_coefficients_satisfy_identity=all(sp.simplify(v)==0 for v in res.values()),
        external_phase_derivative_one=delta[(0,0,1)]==1,
        mixed_y_eta_phase_composition_nonzero=mixed!=0,
        freezing_phase_would_give_wrong_zero=mixed!=0,
        all_colored_product_degrees_bounded=all(sp.simplify(a-b)<=0 for a,b in zip(actual,bounds)),
        cancellation_makes_majorant_an_inequality=any(a<b for a,b in zip(actual,bounds)))
    return dict(checks={k:bool(v) for k,v in checks.items()},mixed_y_eta_phase_coefficient=str(mixed),
        positive_inverse_jet_monomials=len(delta),manufactured_not_actual_phase=True)

@lru_cache(None)
def state_memory_review():
    E,Ey,Uy,Ay,By,N,Lop,aL,bL=sp.symbols('E Ey Uy Ay By N Lop aL bL',nonzero=True)
    a=1-2*Ey/E;b=2*Uy/E
    Aphi=-(aL-a)/2;Bphi=E*(bL-b)/2
    an=1-2*(Ey/E+Aphi+Ay/N);bn=sp.exp(-Lop/N)*(b+2*Bphi/E+2*By/(N*E))
    eta,X,U,h=sp.symbols('eta X U h');D=sp.Rational(1,2)-h;Ap=sp.Rational(1,2)+h;d=1-eta**2;L=1-2*h*eta**2
    m,i,j,ss,c=[sp.Function(n)(eta) for n in ('m','i','j','S','Cp')]
    DW=-(2*D*eta*m+d*sp.diff(m,eta))/X
    B=(1-h)*i-D*eta*sp.diff(i,eta)-d*sp.diff(j,eta)+2*(h-D)*eta*j
    DQ=-DW+B/(X*sp.Symbol('Hx'))
    DN=-U*DW+(D*(m-eta*sp.diff(m,eta))+4*h*eta*ss-d*sp.diff(ss,eta))/X+4*Ap*eta*c-d*sp.diff(c,eta)
    allzero={v:0 for f in (m,i,j,ss,c) for v in (f,sp.diff(f,eta))}
    eps=sp.Symbol('eps',positive=True)
    cp_only={v:0 for f in (m,i,j,ss) for v in (f,sp.diff(f,eta))};cp_only.update({c:eps*(1+eta),sp.diff(c,eta):eps})
    pressure_only=sp.simplify((X*DN/(L*E)).subs(cp_only).subs(eta,0))
    dropU=sp.simplify((-U*DW).subs({m:1+eta,sp.diff(m,eta):1,eta:0,U:1}))
    checks=dict(fast_radial_a_shear_exact=sp.simplify(an-(aL-2*Ay/N))==0,
        exponential_denominator_b_shear_exact=sp.simplify(bn-sp.exp(-Lop/N)*(bL+2*By/(N*E)))==0,
        dropping_exponential_changes_shear=sp.simplify(bn-(bL+2*By/(N*E)))!=0,
        all_five_exact_moments_restore_both_pressure_coordinates=sp.simplify(DQ.subs(allzero))==0 and sp.simplify(DN.subs(allzero))==0,
        pressure_only_memory_is_not_erased=pressure_only==-eps*X/E,
        nonzero_U_DeltaW_memory_required=dropU!=0)
    return dict(checks={k:bool(v) for k,v in checks.items()},Cp_only_p2_at_equator=str(pressure_only),
        full_actual_axis_datum_unchanged=True)

def auxiliary_acceptance(A=2**16,H_power=600):
    if type(A) is not int or A<2**16 or type(H_power) is not int or H_power<600:raise ValueError('Exact auxiliary domain')
    G=A;H=A**H_power;N=H**32+1;s=Q(2000*H**16,N);rb=repair.repair_budget(A,s)
    eps=Q(1,G**64);kappa=Q(1,G**4);R=4*G**3
    raw=Q(1024*A**9*H,N);state=Q((1+2000*A**128)*H**16,N)
    checks=dict(original_repair_radius_accepts=rb['coefficient_norm']<=rb['radius']/2,
        original_repair_preserves_swirl=rb['relative_swirl_error']<Q(1,2),
        raw_modulation_state_fits_H16=1024*A**9*H<H**16,
        original_raw_repair_state_below_Ccorr=rb['state_error']<=A**128*s,
        full_state_bound_below_epsilon_half=state<eps/2,
        actual_raw_positive_errors_below_state_bound=raw+rb['state_error']<=state,
        complete_five_row_target_with_inverse_lambda_fits=64*A**6*H<H**16,
        exact_C1_quadratic_contraction=Q(8*1000**2*10**5*H**16,N)<1,
        same_frequency_above_4AH=N>4*A*H,
        cone_gap_loss_below_half=256*R**10*state<kappa/2)
    return dict(A_aux=A,G_aux=G,H_aux_power=H_power,N_aux_power=32,
        checks={k:bool(v) for k,v in checks.items()},auxiliary_exact_rational_not_actual_H_N_or_profile=True)

@lru_cache(None)
def review():
    ex=exponent_review();co=colored_review();st=state_memory_review();ac=auxiliary_acceptance();un=[universal_review(d) for d in (80,110)]
    sh=shape_majorant(2**16,7,Q(3,2));A=2**16
    checks=dict(ex['checks']);checks.update(co['checks']);checks.update(st['checks']);checks.update(ac['checks']);checks.update(sh['checks'])
    checks['positive_mixed_shape_majorant']=sh['raw_product']<=sh['majorant']
    checks['outward_universal_and_original_repair_checks']=all(all(r['checks'].values()) for r in un)
    checks['mixed_gate_not_a_reverse_fourth_shape_implication']=shape_gate(A,A**32-7*A**20,3)['auxiliary_sufficient_comparison'] and not shape_gate(A,A**32-7*A**20,4)['auxiliary_sufficient_comparison']
    return dict(exponents=ex,colored=co,state_memory=st,auxiliary=ac,universal=un,checks=checks)
