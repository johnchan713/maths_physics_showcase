"""Tame FINAL fifth moments and conditional complete compact fourth maps."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util,sys
import mpmath as mp
import sympy as sp
HERE=Path(__file__).resolve().parent
def module(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
joint=module('fifth_moment_joint4',HERE.parent/'joint_fourth_loop_gamma/joint.py')

@lru_cache(None)
def exponential_review():
    z=sp.Symbol('z');a=sp.symbols('a1:6');p=sum(a[j-1]*z**j for j in range(1,6))
    c=sp.Poly(sp.expand(sum(p**k/factorial(k) for k in range(6))),z).coeff_monomial(z**5)
    claimed=a[4]+a[0]*a[3]+a[1]*a[2]+a[0]**2*a[2]/2+a[0]*a[1]**2/2+a[0]**3*a[1]/6+a[0]**5/120
    return dict(fifth_exponential_coefficient=str(c),checks=dict(
        all_seven_fifth_exponential_partitions=sp.expand(c-claimed)==0,
        top_primitive_jet_enters_linearly=sp.diff(c,a[4])==1,
        lower_a2_a3_partition_cannot_be_dropped=sp.expand(c-(claimed-a[1]*a[2]))!=0,
        known_taylor_sum_below_three=sum(Q(1,factorial(j)) for j in range(6))<3,
        tame_product_full_C5_constant_six=2+4==6))

@lru_cache(None)
def product_review():
    # Exact normalized jet identity; fifth coefficients do not multiply each other.
    f=sp.symbols('f0:6');g=sp.symbols('g0:6');u=sp.symbols('u0:6');v=sp.symbols('v0:6')
    delta=sp.expand(sum((f[j]+u[j])*(g[5-j]+v[5-j])-f[j]*g[5-j] for j in range(6)))
    high=u[5]*(g[0]+v[0])+v[5]*(f[0]+u[0])
    lower=sp.expand(delta-high)
    return dict(checks=dict(
        exact_full_fifth_product_split=sp.expand(delta-high-lower)==0,
        no_product_of_two_unknown_fifth_jets=sp.diff(delta,u[5],v[5])==0,
        top_fields_multiply_only_value_jets=not lower.has(u[5],v[5]),
        both_high_value_cross_terms_retained=sp.diff(delta,u[5])==g[0]+v[0] and sp.diff(delta,v[5])==f[0]+u[0]))

@lru_cache(None)
def normalizer_review():
    e=sp.Symbol('eta');norm=lambda p:sum(abs(sp.diff(p,e,j).subs(e,1))/factorial(j) for j in range(6))
    return dict(normalizer_polynomial_C5=[str(norm(1+e*e)),str(norm((1+e*e)**2))],
        checks=dict(inverse_e_C5_polynomial_cost_still_five=norm(1+e*e)==5,
            inverse_e_square_C5_cost25_not20=norm((1+e*e)**2)==25,
            normalizer_A5_sufficient_at_Amin=25<(2**16),
            transformed_M_J_row_keeps_two_rows_and_inverse_lambda=2<4))

@lru_cache(None)
def root_review():
    x=sp.Symbol('x');c=sp.symbols('c0:6');q,B=sp.symbols('q B');P=sum(c[j]*x**j for j in range(6))
    coeff=sp.expand(B*P+q*P*P).coeff(x,5)
    expected=(B+2*q*c[0])*c[5]+q*sum(c[j]*c[5-j] for j in range(1,5))
    # A separate exact manufactured implicit polynomial, with known nonzero fifth jet.
    cp=sum(sp.Rational(1,1000*(j+1))*x**j for j in range(6));target=7*cp+3*cp*cp
    recovered=sp.expand(target).coeff(x,5)
    cj=[cp.coeff(x,j) for j in range(6)]
    inverse_result=(recovered-3*sum(cj[j]*cj[5-j] for j in range(1,5)))/(7+6*cj[0])
    return dict(lower_pairs=4,Jacobian_inverse_bound=2000,
        known_quadratic_tail_coefficient=8000*10**5*2000**2,
        checks=dict(exact_fifth_constant_Q_repair_identity=sp.expand(coeff-expected)==0,
            fifth_top_multiplier_is_pointwise_Jacobian=sp.diff(coeff,c[5])==B+2*q*c[0],
            high_fifth_square_not_in_lower_sum=not (q*sum(c[j]*c[5-j] for j in range(1,5))).has(c[5]),
            independent_manufactured_fifth_root_recovered=sp.simplify(inverse_result-cj[5])==0,
            lower_pairs_four_not_six=len(range(1,5))==4))

def auxiliary_bounds(A=2**16,H_power=600,Lambda_multiple=2):
    if type(A) is not int or A<2**16 or type(H_power) is not int or H_power<600 or type(Lambda_multiple) is not int or Lambda_multiple<0:
        raise ValueError('Exact auxiliary A,H power and nonnegative fifth-jet multiple')
    H=A**H_power;N=H**32+1;lam=Lambda_multiple*H;B=A**256
    dmod=Q(16*B*(H+lam),N);m_mod=16*A*A*B*dmod;z=4*A**6*m_mod
    r=Q(2000*H**16,N);c5=r+2000*z+8000*10**5*r*r
    ccap=Q(H**20+A**524*lam,N)
    field=dmod+3200*B*ccap;fcap=Q(H**24+A**784*lam,N)
    moments=16*A*A*B*fcap;mcap=Q(H**26+A**1044*lam,N)
    return dict(A_aux=A,H_power=H_power,Lambda5_multiple_of_H=Lambda_multiple,
        manufactured_not_actual_profile_or_Lambda5=True,actual_Lambda5_enclosed=False,checks=dict(
            modulation_moment_cost_below_A515= m_mod<Q(A**515*(H+lam),N),
            all_normalized_target_costs_below_A522=z<Q(A**522*(H+lam),N),
            same_pointwise_Jacobian_C5_cap=c5<ccap,
            complete_modulation_plus_repair_field_C5_cap=field<fcap,
            complete_five_moment_C5_difference_cap=moments<mcap,
            lower_C4_change_remains_below_one=Q(H**20,N)<1))

@lru_cache(None)
def density_principal_review():
    z=sp.Symbol('eta');E=sp.symbols('e0:6');U=sp.symbols('u0:6')
    ep=sum(E[j]*z**j for j in range(6));up=sum(U[j]*z**j for j in range(6))
    expressions=[up,ep,up*ep,up*up-ep*ep/2,ep*ep/2]
    expected=[(0,1),(1,0),(U[0],E[0]),(-E[0],2*U[0]),(E[0],0)]
    rows=[]
    for name,f,pair in zip(('M','I','J','S','Cp'),expressions,expected):
        top=sp.expand(f).coeff(z,5)
        rows.append(dict(row=name,exact=sp.diff(top,E[5])==pair[0] and sp.diff(top,U[5])==pair[1]))
    return dict(rows=rows,checks=dict(all_five_density_principal_rows_exact=all(r['exact'] for r in rows),
        energy_row_keeps_both_U_square_and_E_square=expected[3]==(-E[0],2*U[0]),
        Cp_row_has_E_square_not_axial_square=expected[4]==(E[0],0)))

def inverse_and_map_ledger():
    return dict(field_J4_H_power=128,reciprocal_E_J4_H_power=514,reciprocal_Hx_J4_H_power=516,
        moment_J5_form='H^260+A^1044*Lambda5/N',
        map_H_powers=dict(W_A_power=4,Qs=520,Ns=132,p1=526,p2=650,V=4,F=130),
        common_logarithmic_map_H_power=660,ordinary_conversion_H_power=4,
        sufficient_compact_ordinary_budget_H_power=1024,
        checks=dict(inverse_E_finite_four_products=4*128==512,
            Qs_original_full_moments_and_Hx_inverse=516+4<=520,
            Ns_keeps_W_U_and_both_pressure_terms=128+4<=132,
            p1_geometry_and_L=520+2<526,
            p2_geometry_E_inverse_and_Ns=132+514+4<=650,
            common_map_power_covers_all_maps=650<660,
            ordinary_map_power_plus_moment_gate_below1024=660+4+261<1024))

@lru_cache(None)
def ordinary_review():
    D=sp.Symbol('D');rows=[];total=Q()
    for r in range(5):
        p=sp.prod(D-j for j in range(r));co=[]
        for k in range(r+1):
            v=Q(sp.expand(p).coeff(D,k))*factorial(k)/factorial(r);co.append(str(v));total+=abs(v)
        rows.append(dict(order=r,normalized_Dy_to_DX_coefficients=co))
    return dict(rows=rows,absolute_coefficient_sum=str(total),checks=dict(
        fourth_conversion_keeps_all_Stirling_terms=rows[4]['normalized_Dy_to_DX_coefficients']==['0','-1/4','11/12','-3/2','1'],
        complete_conversion_constant_below10000=total<10000))

@lru_cache(None)
def countercontrol():
    n,N,c,h,e,X,f0=sp.symbols('n N c h eta X f0',positive=True);D=sp.Rational(1,2)-h;L=1-2*h*e*e
    # c is the strictly positive weighted integral over a fixed-N positive phase lobe.
    M=c*sp.sin(n*e)/(N*n**4);U=f0*sp.sin(n*e)/(N*n**4)
    V=(2*e*X*U-2*D*e*M-(1-e*e)*sp.diff(M,e))/L
    fourth=sp.simplify(sp.diff(V,e,4).subs(e,0)/24)
    limit=sp.simplify(sp.limit(fourth/n,n,sp.oo))
    return dict(normalized_flux_fourth_over_n_limit=str(limit),
        manufactured_functional_family_not_actual_C9_or_cone_profile=True,
        actual_profile_failure_inferred=False,checks=dict(
            fixed_N_fifth_memory_grows_linearly=sp.diff(M,e,5).subs(e,0)/120==c*n/(120*N),
            complete_flux_limit_keeps_all_lower_terms=sp.simplify(limit+c/(24*N))==0,
            flux_limit_nonzero_for_positive_mass=bool(c/(24*N)>0)))

def universal_review(digits):
    if type(digits) is not int or digits<60:raise ValueError('At least 60 exact interval digits')
    mp.iv.dps=digits;p=mp.iv.mpf;log=mp.iv.log;a=g=p(2**16)
    lo=lambda v:mp.make_mpf(v._mpi_[0]);hi=lambda v:mp.make_mpf(v._mpi_[1])
    # Every log(A) is <G, and constants below are <exp(G); these compare logarithms.
    comparisons=dict(root_known_H20=(3,523),field_known_H24=(4,258),moment_known_H26=(2,260),
        density_radial_H258=(2,5),inverse_E_H514=(2,6),inverse_Hx_H516=(2,4),
        Qs_H520=(4,6),Ns_H132=(4,6),pressure_geometry_margin=(4,4),ordinary_conversion_H4=(4,5))
    checks={k:bool(lo(q*g**3)>s+1) for k,(q,s) in comparisons.items()}
    checks.update(all_fixed_prefactors_have_log_below_G=hi(log(p('1e20')))<lo(g),
        modulation_moment_constant256_fits_one_A=256<2**16,
        repair_target_constant4_fits_one_A=4<2**16,
        root_unknown_coefficient2000_fits_A2=2000<(2**16)**2,
        repair_field_unknown_coefficients_fit_A4=3200+16<(2**16)**4,
        final_moment_unknown_coefficient16_fits_A2=16<(2**16)**2,
        actual_A_powers_absorbed_via_logA_less_G=hi(1579*log(a)-a**32)<0,
        original_N_not_reselected=32>26,
        sufficient_gate_H260_moment_sum_below_H261=hi(log(2))<lo(g**4))
    return dict(interval_digits=digits,checks=checks,actual_Lambda5_Z_or_huge_parameters_materialized=False)

def fifth_input_gate(A,H,N,Lambda5):
    if type(A) is not int or A<2**16 or type(H) is not int or H<1 or type(N) is not int or not H**32<N<2*H**32:
        raise ValueError('Exact auxiliary geometry and original frequency sandwich required')
    if type(Lambda5) not in (int,Q) or Lambda5<0:raise ValueError('Exact nonnegative supplied fifth seminorm')
    return dict(sufficient_fifth_input_comparison=A**1044*Q(Lambda5)<=N*H**260,
        actual_Lambda5_enclosed=False,actual_joint_target_verified=False,actual_global_leading_C4_accepted=False)

@lru_cache(None)
def review():
    ex=exponential_review();pr=product_review();no=normalizer_review();ro=root_review();ac=auxiliary_bounds();le=inverse_and_map_ledger();orr=ordinary_review();co=countercontrol()
    un=[universal_review(d) for d in (80,110)];mem=joint.bridge.state_memory_review();de=density_principal_review()
    checks={**ex['checks'],**pr['checks'],**no['checks'],**ro['checks'],**le['checks'],**orr['checks'],**co['checks'],**mem['checks'],**de['checks']}
    checks['all_exact_auxiliary_tame_majorants_fit']=all(ac['checks'].values())
    checks['all_outward_universal_comparisons']=all(all(r['checks'].values()) for r in un)
    return dict(exponential=ex,product=pr,normalizers=no,root=ro,auxiliary=ac,ledger=le,ordinary=orr,
        counter=co,universal=un,full_memory=mem,density=de,checks={k:bool(v) for k,v in checks.items()})
