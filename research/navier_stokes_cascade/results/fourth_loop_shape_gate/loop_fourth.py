"""Conditional fourth angular loop/moving-inverse bounds, with explicit shape gate."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial
from pathlib import Path
import importlib.util,sys
import mpmath as mp
import sympy as s

HERE=Path(__file__).resolve().parent
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m

def product(a,b,n):
    out={}
    for (ka,ea),ca in a.items():
        for (kb,eb),cb in b.items():
            if ka+kb<=n:out[ka+kb,ea+eb]=out.get((ka+kb,ea+eb),0)+ca*cb
    return out
def power(a,k,n):
    out={(0,0):1}
    for _ in range(k):out=product(out,a,n)
    return out

def implicit_budget(partial_power,inverse_power,selected):
    """Exponents in exp(c*G); explicit normalized two-input Taylor counts."""
    if type(partial_power) is not int or type(inverse_power) is not int or min(partial_power,inverse_power)<0:
        raise ValueError('Nonnegative integer exponents required')
    if len(selected)!=4 or any(type(v) is not int or v<1 for v in selected):raise ValueError('Four positive selected exponents required')
    rows=[]
    for n,chosen in enumerate(selected,1):
        u={(j,selected[j-1]):1 for j in range(1,n)}
        p={(j,1):1 for j in range(1,n+1)}
        terms={}
        for i in range(n+1):
            for j in range(n+1-i):
                if i+j==0:continue
                for (k,e),c in product(power(u,i,n),power(p,j,n),n).items():
                    if k==n:terms[e]=terms.get(e,0)+c
        count=sum(terms.values())+1 # the target coefficient
        largest=max(terms)
        sufficient=inverse_power+max(1,partial_power+largest)+1
        rows.append(dict(order=n,lower_composition_terms=count,largest_input_power=largest,
            sufficient_power=sufficient,selected_power=chosen,accepted=chosen>=sufficient))
    return rows

@lru_cache(None)
def exponent_review():
    mu=implicit_budget(64,10,[80,240,400,560])
    theta=implicit_budget(10010,1,[11000,33000,55000,77000])
    constants=dict(q_joint5=4,removable_R_joint4=5,R2_C4=10,sqrt_R2_C4=60,
        signed_F_partial_C4=64,inverse_Fmu=10,mu_C4=600,
        fixed_angle_t_joint4=4000,w_joint4=10000,Phi_joint4=10010,
        moving_Theta_C4=80000,fixed_angle_sources_joint4=20000,
        pulled_sources_C4=350000,primitive_C4=400000)
    return dict(mu=mu,theta=theta,constants=constants,checks=dict(
        all_normalized_signed_root_orders_have_headroom=all(v['accepted'] for v in mu),
        all_normalized_inverse_phase_orders_have_headroom=all(v['accepted'] for v in theta),
        all_Taylor_partition_counts_below10000=max(v['lower_composition_terms'] for v in mu+theta)<10000,
        sqrt_finite_Taylor_cost_fits60=5+4*(10+3)+1<60,
        signed_F_partial_composition_fits64=60+1+1<64,
        mu_norm_sum_fits600=max(v['selected_power'] for v in mu)+1<600,
        fixed_angle_removable_t_composition_fits4000=5+4*(600+1+1)+600+1<4000,
        forward_lift_after_w_fits10010=10000+1<10010,
        inverse_phase_norm_sum_fits80000=max(v['selected_power'] for v in theta)+1<80000,
        rational_loop_and_E_source_cost_fits20000=4*4000+1+1+1+1<20000,
        complete_fourth_moving_composition_fits350000=20000+4*(80000+1)<350000,
        periodic_primitive_operator_two_fits400000=350000+1<400000))

def universal_gate(A,G,digits=80):
    if type(A) is not int or A<2**16 or type(G) is not int or G<2**16:raise ValueError('Exact auxiliary A,G>=2^16 required')
    if type(digits) is not int or digits<60:raise ValueError('At least 60 interval digits required')
    mp.iv.dps=digits;iv=mp.iv.mpf;log=mp.iv.log
    upper=lambda v:mp.make_mpf(v._mpi_[1]);lower=lambda v:mp.make_mpf(v._mpi_[0])
    a,g=iv(A),iv(G)
    absorption=12*log(10)+10*log(g)-g
    # These prove all-larger-parameter inequalities by positive derivatives.
    checks=dict(polynomial_absorption_at_minimum=upper(absorption)<0,
        absorption_log_ratio_decreases=lower(1-10/g)>0,
        A10000_below_exp_A32=upper(10000*log(a)-a**32)<0,
        A_power_comparison_ratio_decreases=lower(32-1/log(a))>0,
        primitive_400000G_below_logH_G4=lower(g**3)>400000,
        compact_known_part_constant_fits_A259=lower(a**3)>1+6400+960000000000,
        normalizer_gate_coefficients_fit_A525=lower(a*a)>1+3200000,
        fixed_frequency_root_radius_below_one=lower(g**16)>2000,
        same_C1_small_ball_uniqueness=lower(g**16)>8*1000**2*10**5,
        conditional_repair_target_and_C4_contraction_small=lower(g**29)>8*1000**2*10**5*3,
        logH_exceeds_log_G2=lower(g**4)>upper(2*log(g)),
        conditional_high_cost_G_plus_one_below_H_by_elementary_G2=lower(g*g)>upper(g+1))
    return dict(interval_digits=digits,A_auxiliary=A,G_auxiliary=G,checks=checks,
        actual_variance_shape_or_gate_accepted=False,enormous_actual_G_H_N_materialized=False)

@lru_cache(None)
def q5_review():
    z,h=s.symbols('z h');q=s.series(s.exp(z*h)/(1+z*z/4+z**4/64),z,0,6).removeO()
    coeff=[s.expand(q).coeff(z,j) for j in range(6)]
    R=sum(coeff[j+1]*z**j for j in range(5))
    def average(poly):
        return s.expand(sum(c*(s.Rational(comb(2*(k//2),k//2),4**(k//2)) if k%2==0 else 0)
            for (k,),c in s.Poly(s.expand(poly),h).terms()))
    R2=s.series(average(R*R),z,0,5).removeO()
    without=s.series(average((R-coeff[5]*z**4)**2),z,0,5).removeO()
    mu,p,d=s.symbols('mu p d',real=True)
    signed=d*mu*s.series(s.sqrt(R2),z,0,5).removeO().subs(z,mu*p)
    # Mean-one lift is an exact consequence of the complete variance equation.
    a,ts,r2=s.symbols('a ts r2');v=a*(1+ts*ts)+a*d*d*mu*mu*r2
    mass=s.expand(a*(1+ts*ts+d*d*mu*mu*r2)-v)
    return dict(q5=str(coeff[5]),R2_through_four=str(R2),
        checks=dict(exact_fifth_q_coefficient=coeff[5]==h**5/120-h**3/24+3*h/64,
            exact_removable_R2_zero_and_fourth=R2==s.Rational(1,2)-3*z*z/32+11*z**4/576,
            omitting_q5_changes_fourth_R2=s.expand(R2-without).coeff(z,4)!=0,
            signed_Fmu_at_zero=d/s.sqrt(2)==s.diff(signed,mu).subs(mu,0),
            signed_Fp_at_zero=s.diff(signed,p).subs(mu,0)==0,
            unsigned_variance_mu_derivative_vanishes_at_zero=s.diff(signed*signed,mu).subs(mu,0)==0,
            exact_variance_implies_full_mean_one_lift=mass==0))

@lru_cache(None)
def moving_review():
    h=s.Symbol('h');theta0=s.pi/3;rows=[]
    for eta0 in (-1,0,1):
        jets=[]
        for n in range(1,5):
            delta=sum(jets[j-1]*h**j for j in range(1,n))
            sine=sum(s.diff(s.sin(s.Symbol('t')),s.Symbol('t'),j).subs(s.Symbol('t'),theta0)*delta**j/s.factorial(j) for j in range(5))
            residual=s.expand(delta+(eta0+h)*sine/2-eta0*s.sin(theta0)/2)
            jets.append(s.simplify(-residual.coeff(h,n)/(1+s.Rational(eta0,2)*s.cos(theta0))))
        delta=sum(jets[j-1]*h**j for j in range(1,5))
        sinjet=sum(s.diff(s.sin(s.Symbol('t')),s.Symbol('t'),j).subs(s.Symbol('t'),theta0)*delta**j/s.factorial(j) for j in range(5))
        res=s.expand(delta+(eta0+h)*sinjet/2-eta0*s.sin(theta0)/2)
        cosine=sum(s.diff(s.cos(s.Symbol('t')),s.Symbol('t'),j).subs(s.Symbol('t'),theta0)*delta**j/s.factorial(j) for j in range(5))
        fourth=s.simplify(s.expand(cosine).coeff(h,4))
        rows.append(dict(eta=eta0,implicit_residual_zero=all(s.simplify(res.coeff(h,k))==0 for k in range(1,5)),
            fourth=str(fourth),frozen_phase_gives_wrong_zero=fourth!=0,
            positive_w=bool(1+s.Rational(eta0,2)*s.cos(theta0)>0)))
    return dict(manufactured_not_actual_phase=True,rows=rows,checks=dict(
        all_four_inverse_coefficients_satisfy_identity=all(v['implicit_residual_zero'] for v in rows),
        fixed_origin_and_poles_keep_positive_inverse=all(bool(v['positive_w']) for v in rows),
        moving_fourth_coefficient_135_over4096=rows[1]['fourth']=='135/4096',
        freezing_phase_fails_at_all_three_eta_centers=all(v['frozen_phase_gives_wrong_zero'] for v in rows)))

@lru_cache(None)
def root_diagnostic():
    scalar=module('shape_gate_retained_scalar',HERE.parent/'higher_jet_budget/higher_jets.py')
    def mul(a,b):return [sum(a[j]*b[k-j] for j in range(k+1)) for k in range(5)]
    def powr(a,k):
        v=[mp.mpf(1)]+[mp.mpf(0)]*4
        for _ in range(k):v=mul(v,a)
        return v
    offsets=list(range(-4,5));weights=s.Matrix([[s.Integer(i)**k for k in range(9)] for i in offsets]).inv().row(4)
    rows=[]
    with mp.workdps(90):
        for name,ss,pp in [('zero_pressure',['.02','.003','.004','.0002','.0003'],['0','.2','.03','0','.01']),
            ('zero_variance',['0','.01','.003','.002','.001'],['.4','.02','.01','0','.001'])]:
            target=list(map(mp.mpf,ss));pressure=list(map(mp.mpf,pp));d0=mp.mpf('.1')
            F=lambda m,p:scalar.signed_F(m,p,d0)
            m0=mp.findroot(lambda m:F(m,pressure[0])-target[0],target[0]*mp.sqrt(2)/d0)
            partial={(i,j):mp.diff(F,(m0,pressure[0]),(i,j))/(factorial(i)*factorial(j)) for i in range(5) for j in range(5-i)}
            mu=[m0]
            for n in range(1,5):
                dm=[mp.mpf(0)]+mu[1:]+[mp.mpf(0)]*(5-len(mu))
                dp=[mp.mpf(0)]+pressure[1:]
                rhs=sum(c*mul(powr(dm,i),powr(dp,j))[n] for (i,j),c in partial.items() if i+j)
                mu.append((target[n]-rhs)/partial[1,0])
            step=mp.mpf('1e-4');samples=[]
            for offset in offsets:
                x=step*offset;sv=sum(target[j]*x**j for j in range(5));pv=sum(pressure[j]*x**j for j in range(5))
                samples.append(mp.findroot(lambda m:F(m,pv)-sv,(m0-mp.mpf('.01'),m0+mp.mpf('.01'))))
            interpolated=sum(mp.mpf(w.p)/int(w.q)*v for w,v in zip(weights,samples))/step**4
            error=abs(interpolated-mu[4])/max(1,abs(mu[4]))
            # Freeze the tolerance and outcomes, not version-dependent roundoff strings.
            rows.append(dict(case=name,relative_error_limit='1e-18',solved_values=9,
                independent_solved_root_interpolation=error<mp.mpf('1e-18'),Fmu_positive=partial[1,0]>0))
    return dict(manufactured_floating_not_interval_or_actual_profile=True,rows=rows,checks=dict(
        fourth_implicit_root_matches_independent_nine_solved_values=all(v['independent_solved_root_interpolation'] for v in rows),
        zero_pressure_and_zero_variance_keep_positive_Fmu=all(v['Fmu_positive'] for v in rows)))

def shape_gate(A,Z,logZ=None):
    if type(A) is not int or A<2**16 or type(Z) not in (int,Q) or Z<1:raise ValueError('Exact A>=2^16,Z>=1 required')
    if logZ is not None and (type(logZ) not in (int,Q) or logZ<0):raise ValueError('Nonnegative exact auxiliary logZ required')
    with mp.workdps(80):
        a=mp.mpf(A);z=mp.mpf(Z.numerator)/Z.denominator if isinstance(Z,Q) else mp.mpf(Z)
        lz=mp.log(z) if logZ is None else mp.mpf(logZ.numerator)/logZ.denominator if isinstance(logZ,Q) else mp.mpf(logZ)
        lhs=lz+12*mp.log(10)+20673*mp.log(a)+8*a**20
        return dict(sufficient=lhs<=a**32,actual_shape_budget_supplied=False,actual_compact_gate_accepted=False)

@lru_cache(None)
def target_gate_review():
    target=module('shape_gate_complete_target',HERE.parent/'fourth_target_normalization/target_bounds.py')
    A=2**16;Ahi=A**256
    # Positive-coefficient worst normalizer enclosures, not a new physical patch.
    affine=target.target_affines(A,Ahi,Q(1,A),A*A,Q(1,A),Q(1,A))
    caps=affine['Z4_affine']
    # An auxiliary algebraic H,N only; no actual exponential hierarchy is evaluated.
    H=A**600;N=H**32+1;r=Q(2000)*H**16/N
    gate=target.primitive_gate(A,Ahi,Q(1,A),A*A,Q(1,A),Q(1,A),H,N,r,H,H,H**64)
    return dict(auxiliary_algebraic_parameters_not_actual_H_N=True,
        normalizer_common_A_power=5,moment_coefficient_common_A_power=515,
        transformed_target_coefficient_A_power=522,gate_coefficient_A_power=525,
        compact_known_part_A_power=259,checks=dict(
            all_five_original_C4_normalizers_fit_A5=all(v<=A**5 for v in affine['normalizer_C4'].values()),
            all_three_actual_target_affine_coefficients_fit_A522=all(v<=A**522 for v in caps),
            inverse_lambda_transformed_row_is_retained=all(affine['transformed_rows'][1][j]>=A*affine['transformed_rows'][0][j] for j in range(3)),
            original_substituted_compact_inequality_accepts_auxiliary_bound=gate['accepts'],
            original_positive_compact_majorant_also_accepts=gate['substituted_majorant_accepts'],
            actual_shape_or_unconditional_gate_not_inferred=True))

def review():
    e=exponent_review();q=q5_review();m=moving_review();r=root_diagnostic()
    t=target_gate_review()
    u=[universal_gate(2**16,2**16,d) for d in (80,110)]
    return dict(exponents=e,q5=q,moving=m,root=r,target=t,universal=u,
        checks=dict(**e['checks'],**q['checks'],**m['checks'],**r['checks'],**t['checks'],
            two_precision_universal_monotone_parameter_comparisons=all(all(v['checks'].values()) for v in u),
            auxiliary_Z1e100_satisfies_shape_gate_not_actual_enclosure=shape_gate(2**16,10**100)['sufficient']))
