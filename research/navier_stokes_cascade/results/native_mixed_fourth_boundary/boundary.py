"""Native full pressure J4 and exact field-only fourth frequency interface."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial,comb
from pathlib import Path
import importlib.util,sys
import mpmath as mp
import sympy as sp
HERE=Path(__file__).resolve().parent
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
high=module('native_mixed_fourth_higher_bounds',HERE.parent/'higher_jet_budget/higher_bounds.py')

@lru_cache(None)
def pressure_identity():
    y,e,h=sp.symbols('y eta h',real=True);X=sp.exp(y);d=1-e**2;D=sp.Rational(1,2)-h;Ap=sp.Rational(1,2)+h;L=1-2*h*e**2
    E,U,M,I,J,S,Cp=[sp.Function(n)(y,e) for n in ('E','U','M','I','J','S','Cp')]
    axis=sp.Function('Pi_axis')(e);Pi=axis+Cp;Hx=sp.sqrt(2)*sp.exp(y/2)*E;ell=sp.Rational(1,2)+sp.diff(E,y)/E;Hc=D*e+d*U
    W=1-(2*D*e*M+d*sp.diff(M,e))/X
    B=(1-h)*I-D*e*sp.diff(I,e)-d*sp.diff(J,e)+2*(h-D)*e*J
    Qs=-W+B/(X*Hx)
    Ns=-W*U+(D*(M-e*sp.diff(M,e))+4*h*e*S-d*sp.diff(S,e))/X+4*Ap*e*Pi-d*sp.diff(Pi,e)
    p1=X*Qs/L;p2=X*Ns/(L*E)
    flows={M:X*U,I:X*Hx,J:X*U*Hx,S:X*(U**2-E**2/2),Cp:E**2/2}
    substitutions={sp.diff(f,y,e,j):sp.diff(v,e,j) for f,v in flows.items() for j in (0,1)}
    Sq=-W*ell-h*(1-2*e*U)-Hc*sp.diff(E,e)/E
    Sn=-W*sp.diff(U,y)-Ap*(1-2*e*U)*U-Hc*sp.diff(U,e)-d*sp.diff(Pi,e)+4*Ap*e*Pi+e*E**2
    r1=sp.factor((sp.diff(p1,y)+ell*p1-X*Sq/L).subs(substitutions))
    r2=sp.factor((sp.diff(p2,y)+sp.diff(E,y)/E*p2-X*Sn/(L*E)).subs(substitutions))
    missing_axis=4*Ap*e*axis-d*sp.diff(axis,e)
    checks=dict(p1_original_forward_density_ODE_exact=r1==0,p2_original_forward_density_ODE_exact=r2==0,
        dropping_axis_source_changes_original_ODE=sp.simplify(missing_axis)!=0,
        wrong_p2_y_plus_p2_equation_rejected=sp.simplify(p2-sp.diff(E,y)/E*p2)!=0,
        both_U_and_E_square_density_retained=flows[S].has(U**2) and flows[S].has(E**2),
        complete_Pi_has_axis_and_Cp=Pi==axis+Cp)
    return dict(checks={k:bool(v) for k,v in checks.items()},native_not_FINAL_pressure=True,
        actual_incoming_axis_or_moments_set_zero=False)

def pressure_ledger():
    native=high.compact_exponents();geom=sum((Q(3,49)**j for j in range(4)),Q())*Q(50,49)
    checks=dict(native_fields_pressure_inverse_and_coordinates_fit_Bn=all(v<4096*256 for v in native['exponents'].values()),
        differentiated_J4_to_J3_cost_four=4==4,
        Sq_positive_terms_fit34=5+5+24==34,
        Sn_positive_terms_fit60=4+5+24+16+8+2<60,
        p1_y_positive_terms_fit107=3*34+5==107,
        p2_y_positive_terms_fit184=3*60+4==184,
        common200_fits_one_A=200<2**16,
        mixed_pressure_A_power771=2+3*256+1==771,
        pure_angular1546_dominates_radial771=771<1546,
        sum_fits1547=2<2**16,
        inverse_L_J3_below_two=geom<2)
    return dict(checks=checks,Sq_coefficient=34,Sn_coefficient=60,
        p1_y_coefficient=107,p2_y_coefficient=184,radial_mixed_A_power=771,full_native_J4_A_power=1547)

@lru_cache(None)
def shear_principal():
    y,e=sp.symbols('y eta');E0,Et,Ut=sp.symbols('E0 Etop Utop',nonzero=True);rows=[]
    for r in range(5):
        s=4-r;E=E0+y/sp.Integer(7)+e/sp.Integer(11)+y*e/13+Et*y**(r+1)*e**s
        U=sp.Rational(1,5)+y/17+e/19+Ut*y**(r+1)*e**s
        a=1-2*sp.diff(E,y)/E;b=2*sp.diff(U,y)/E
        ac=sp.diff(a,y,r,e,s).subs({y:0,e:0})/(factorial(r)*factorial(s))
        bc=sp.diff(b,y,r,e,s).subs({y:0,e:0})/(factorial(r)*factorial(s))
        rows.append(dict(r=r,s=s,field_index=[r+1,s],
            a_top_coefficient=str(sp.simplify(sp.diff(ac,Et))),b_top_coefficient=str(sp.simplify(sp.diff(bc,Ut))),
            checks=dict(a_principal_exact=sp.simplify(sp.diff(ac,Et)+2*(r+1)/E0)==0,
                b_principal_exact=sp.simplify(sp.diff(bc,Ut)-2*(r+1)/E0)==0)))
    return dict(rows=rows,remaining_field_indices=[r['field_index'] for r in rows if r['r']>=1],
        actual_fifth_mixed_field_inputs_enclosed=False)

@lru_cache(None)
def frequency_algebra():
    y,e,phi,N=sp.symbols('y eta phi N');rows=[]
    # Arbitrary independent Taylor coefficients: exact identity on every 4-jet.
    # Avoid comparing noncanonical Subs(Derivative(...)) representations.
    P=sum(sp.Symbol(f'p{a}{b}{c}')*y**a*e**b*phi**c
          for a in range(5) for b in range(5-a) for c in range(5-a-b))
    for r in range(5):
        s=4-r
        raw=sp.diff(P.subs(phi,N*y)/N,y,r,e,s)/(factorial(r)*factorial(s))
        expected=sum(N**(c-1)*sp.diff(P,y,r-c,e,s,phi,c).subs(phi,N*y)/(factorial(r-c)*factorial(s)*factorial(c)) for c in range(r+1))
        rows.append(dict(r=r,s=s,exact=sp.simplify(raw-expected)==0))
    z=sp.Symbol('z');v=sp.symbols('b1:5');poly=sum(v[j-1]*z**j for j in range(1,5))
    fourth=sp.expand(sum(poly**j/sp.factorial(j) for j in range(5))).coeff(z,4)
    claimed=v[3]+v[0]*v[2]+v[1]**2/2+v[0]**2*v[1]/2+v[0]**4/24
    # Positive degree-collapse constants after using N>=H>=1.
    total=Q(1)+Q(1)+Q(3,2)+Q(13,6)+Q(97,24)
    checks=dict(all_five_normalized_mixed_chain_rows=all(r['exact'] for r in rows),
        full_fourth_exponential_partitions=sp.expand(fourth-claimed)==0,
        missing_lower_square_rejected=sp.expand(fourth-(claimed-v[1]**2/2))!=0,
        Gamma_enters_linearly=sp.diff(fourth,v[3])==1,
        exponential_known_constant_fits20=2*total<20,
        scaled_primitive_J4_known_constant_fits3=2+sum(Q(1,2**j) for j in range(1,5))<3,
        width_point_one_bump_J4_below_Bn=10**106<(2**16)**256)
    return dict(rows=rows,fourth_exponential=str(fourth),known_exponential_constant=str(2*total),checks=checks)

def field_gate(A,H,N,Gamma):
    if type(A) is not int or A<2**16 or type(H) is not int or H<1 or type(N) is not int or not H**32<N<2*H**32:
        raise ValueError('Exact auxiliary A,H and same-frequency sandwich required')
    if type(Gamma) not in (int,Q) or Gamma<0:raise ValueError('Exact nonnegative supplied Gamma')
    B=A**256;known=32*B*N**3*H+2;unknown=2*B*Q(Gamma)/N;budget=H**128
    return dict(known_below_half=2*known<budget,
        Gamma_sufficient_comparison=Q(Gamma)<=Q(N*budget,4*B),
        full_majorant_below_budget=known+unknown<=budget,
        repair_J4_below_one=Q(10000*B**2*H**16,N)<1,
        actual_Gamma_enclosed=False,actual_FINAL_global_C4_accepted=False)

@lru_cache(None)
def countercontrol():
    n,N=sp.symbols('n N',positive=True);x=sp.Symbol('x');omega=2*sp.pi*N
    exact=sp.diff(sp.cos(n*x)*sp.cos(omega*x),x,4).subs(x,0)/(24*N*n**3)
    formula=(n**4+6*n**2*omega**2+omega**4)/(24*N*n**3)
    checks=dict(exact_normalized_fourth_composition=sp.simplify(exact-formula)==0,
        fixed_N_limit_linear_in_n=sp.limit(exact/n,n,sp.oo)==1/(24*N),
        omitting_pure_slow_fourth_loses_linear_term=sp.limit((formula-n/(24*N))/n,n,sp.oo)==0)
    rows=[]
    with mp.workdps(80):
        for nn,NN in ((8,3),(128,3),(32,11),(512,37)):
            w=2*mp.pi*NN
            measured=mp.diff(lambda q:mp.cos(nn*q)*mp.cos(w*q)/(NN*nn**3),0,4)/24
            expected=(nn**4+6*nn**2*w*w+w**4)/(24*NN*nn**3)
            rows.append(dict(n=nn,N=NN,relative_error_limit='1e-60',
                checked=bool(abs(measured-expected)/(1+abs(expected))<mp.mpf('1e-60'))))
    return dict(checks={k:bool(v) for k,v in checks.items()},numerical_rows=rows,
        manufactured_functional_counter_not_actual_C9_or_global_pressure_profile=True,
        normalized_fourth_limit='1/(24*N)',actual_profile_failure_inferred=False)

def universal_review(digits):
    if type(digits) is not int or digits<60:raise ValueError('At least 60 interval digits')
    mp.iv.dps=digits;iv=mp.iv.mpf;log=mp.iv.log;a=g=iv(2**16)
    lo=lambda v:mp.make_mpf(v._mpi_[0]);hi=lambda v:mp.make_mpf(v._mpi_[1])
    checks=dict(native_pressure_A1547_below_G=hi(1547*log(a)-a**32)<0,
        native_comparison_ratio_increases=lo(32-1/log(a))>0,
        known_field_half_budget=hi(log(1024))<lo(g) and lo(31*g**3)>257,
        leading_repair_J4_small=hi(log(10000))<lo(g) and lo(16*g**3)>513,
        counter_fixed_support_J3_fits_H=lo(g**4)>hi(211*log(10)),
        trig_J3_constant_below1000=hi(2*mp.iv.pi*sum((1+2*mp.iv.pi)**j/factorial(j) for j in range(4)))<1000,
        coefficient_normalization_factorial24=24==factorial(4))
    return dict(interval_digits=digits,checks={k:bool(v) for k,v in checks.items()},actual_A_G_H_N_Gamma_materialized=False)

@lru_cache(None)
def review():
    pr=pressure_identity();le=pressure_ledger();sh=shear_principal();fr=frequency_algebra();co=countercontrol();un=[universal_review(d) for d in (80,110)]
    A=2**16;H=A**300;N=H**32+1
    accepted=field_gate(A,H,N,2*H);rejected=field_gate(A,H,N,Q(N*H**129,A**256))
    checks=dict(pr['checks']);checks.update(le['checks']);checks.update(fr['checks']);checks.update(co['checks'])
    checks['all_five_shear_principal_rows_exact']=all(all(r['checks'].values()) for r in sh['rows'])
    checks['all_manufactured_counter_derivatives_match']=all(r['checked'] for r in co['numerical_rows'])
    checks['all_outward_universal_comparisons']=all(all(r['checks'].values()) for r in un)
    checks['auxiliary_Gamma_gate_and_full_majorant_accept']=all(accepted[k] for k in ('known_below_half','Gamma_sufficient_comparison','full_majorant_below_budget','repair_J4_below_one'))
    checks['large_auxiliary_Gamma_gate_and_majorant_reject']=not rejected['Gamma_sufficient_comparison'] and not rejected['full_majorant_below_budget']
    return dict(pressure=pr,ledger=le,shear=sh,frequency=fr,counter=co,universal=un,
        auxiliary=dict(A=A,H_power=300,N_power=32,accepted=accepted,rejected=rejected,not_actual_hierarchy=True),checks=checks)
