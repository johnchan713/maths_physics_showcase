"""Global C_eta^3 axial envelope: exact budgets and manufactured root diagnostics."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import importlib.util
from math import comb
from pathlib import Path

def rat(x):
    if isinstance(x,(bool,float)):raise TypeError('Exact rational inputs required')
    return F(x)

def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x

def root_jet_caps(epsilon):
    """Normalized jets of k*a^2+B*a+D=0, |k|<=1, inverse Jacobian<=3."""
    epsilon=rat(epsilon)
    if not 0<=epsilon<=F(1,2**20):raise ValueError('Need 0<=epsilon<=2^-20')
    u=[F(2)]
    for n in range(1,4):
        u.append(3*(sum(u[j]*u[n-j] for j in range(1,n))+
                    epsilon*sum(u[n-j] for j in range(1,n+1))+epsilon))
    return u

def repaired_source_bounds(h,Xv,K_axis,Ahi):
    h,Xv,K_axis,Ahi=map(rat,(h,Xv,K_axis,Ahi))
    if min(K_axis,Ahi)<0:raise ValueError('Nonnegative region bounds required')
    p=Path(__file__).resolve().parent.parent/'signed_pressure_source_moments/moments.py'
    s=importlib.util.spec_from_file_location('global_axial_signed_interface',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    KU=max(K_axis,2*Ahi,F(352));return dict(KU=KU,**m.from_global_axial_norm(h,1,Xv,KU))

def exact_budget():
    u=root_jet_caps(F(1,2**20));caps=(9,279,15939)
    return dict(constants=dict(angular_discrepancy_C3=2**21,angular_target_C3_coefficient=2**29,
        angular_root_C3_coefficient=2**32,post_energy_C3_coefficient=8192,
        strength_normalized_jet_caps=list(caps),strength_C3_upper='2',end_coefficient_C3_multiplier=4,
        f_C3_upper=16,pulse_axial_C3_upper=352),
        binary_upper_exponents=dict(angular_contraction=41-256*30,angular_C3_radius=32-256*30,
            epsilon=13+6-256),
        checks=dict(discrepancy_enclosure=18*2**16<2**21,
            target_decay_conversion=3+1<=4 and 4**4==2**8,
            angular_contraction_strict=F(8*9*32*2**29,2**(256*30))<1,
            radius_enclosure=6*2**29<2**32,
            old_root_in_new_pointwise_uniqueness_ball=4200<2**32,
            pointwise_same_root_lipschitz=F(192*2**32,2**(256*30))<1,
            angular_radius_tiny=32-256*30<-40,
            normalized_interpolation_power_norm=1+4+F(12,2)+F(56,6)<32,
            interpolation_energy_at_T64=32*1000<=512*64,
            complete_unedited_post_energy_bound=6*512+512<4096,
            angular_edit_energy_less_than_one=F(29,2**40)+F(960,2**80)<1,
            complete_energy_norm_bound=4096*64+1<=8192*64,
            s0_polynomial_C3=4+2*15+42==76,
            s0_norm_at_T64=241*64+F(17,6)+19<512*64,
            linear_coefficient_bound=160<8192*64,
            constant_coefficient_jet_budget=512+4096+1<8192,
            quadratic_coefficient_at_most_one=F(1,4)+F(40,2**(256*41))<1,
            epsilon_selection=13+6-256<-20,
            inverse_Jacobian_padding=F(100,49)<3,
            jet_caps_cover_recursion=all(u[n]<=caps[n-1]*F(1,2**20) for n in range(1,4)),
            strength_full_C3_below_two=F(6,5)+sum(caps)*F(1,2**20)<2,
            affine_prefix_polynomial_C3=2+4+3+1==10,
            prefix_C3_from_C1=F(10,6)<2,
            composed_axial_end_coefficients=2+2==4,
            reciprocal_f_factorial_C3=1+2+4+8<16,
            disjoint_main_shape_budget=11*2==22,
            disjoint_end_bump_budget=F(64*4,2**(256*20))<22,
            pulse_factor=16*22==352))

def symbolic_review():
    import sympy as s
    e,theta=s.symbols('eta theta',real=True);p=1+e**2;g=(p/2)**theta
    forms=[2*theta*e*g/p,
        g*(2*theta/p+4*theta*(theta-1)*e**2/p**2),
        g*(12*theta*(theta-1)*e/p**2+8*theta*(theta-1)*(theta-2)*e**3/p**3)]
    identities=[s.simplify(s.diff(g,e,n)-forms[n-1])==0 for n in range(1,4)]
    z=s.symbols('z',positive=True);a=s.Function('a')(z);B=s.Function('B')(z);D=s.Function('D')(z);k=s.symbols('k')
    energy=k*a*a+B*a+D;want=(2*k*a+B)*s.diff(a,z,3)+6*k*s.diff(a,z)*s.diff(a,z,2)+3*s.diff(B,z)*s.diff(a,z,2)+3*s.diff(B,z,2)*s.diff(a,z)+s.diff(B,z,3)*a+s.diff(D,z,3)
    cross=s.simplify(s.diff(energy,z,3)-want)==0
    return dict(checks=dict(interpolation_power_derivatives=all(identities),third_root_jet_all_product_terms=cross,
        omitting_linear_coefficient_derivative_is_detected=s.simplify(3*s.diff(B,z)*s.diff(a,z,2)+3*s.diff(B,z,2)*s.diff(a,z)+s.diff(B,z,3)*a)!=0))

def manufactured_review():
    """Differentiate the explicit positive quadratic root, not selected giant data."""
    import sympy as s
    e=s.symbols('eta');epsilon=s.Rational(1,2**20);k=s.Rational(1,4)
    B=epsilon*(e+e**3)/10;D=-s.Rational(1,4)+epsilon*(e**2+e**4)/20
    alpha=(-B+s.sqrt(B*B-4*k*D))/(2*k);caps=(9,279,15939);rows=[]
    for angle in (s.Rational(-1),s.Rational(-2,5),s.Rational(0),s.Rational(3,5),s.Rational(1)):
        jets=[s.diff(alpha,e,n).subs(e,angle)/s.factorial(n) for n in range(4)]
        bs=[s.diff(B,e,n).subs(e,angle)/s.factorial(n) for n in range(4)]
        ds=[s.diff(D,e,n).subs(e,angle)/s.factorial(n) for n in range(4)]
        checks=[]
        for n in range(1,4):
            lhs=(2*k*jets[0]+bs[0])*jets[n]
            rhs=-k*sum(jets[j]*jets[n-j] for j in range(1,n))-sum(bs[j]*jets[n-j] for j in range(1,n+1))-ds[n]
            checks.append(s.simplify(lhs-rhs)==0)
        got=[abs(complex(s.N(v,40))) for v in jets[1:]]
        rows.append(dict(eta=str(angle),checks=dict(exact_same_root_recursion=all(checks),
            normalized_jets_below_caps=all(v<float(c*epsilon) for v,c in zip(got,caps)))))
    return dict(manufactured_not_selected_parameters=True,rows=rows,
        checks=dict(all_five_root_and_first_three_jet_controls=all(all(r['checks'].values()) for r in rows)))

def regional_review():
    # Domain glue uses C>=1e9, hence log C>1 and (tc/2)/t<1/20.
    # The axis Cauchy norm has factorial cancellation, not extra k! factors.
    d=F(1,12800);K=F(7);axis=K*sum(d**(-j) for j in range(4))
    r=repaired_source_bounds(F(1,200),11,axis,13)
    return dict(fixture=serial(r),axis_fixture=str(axis),
        checks=dict(common_cover_offsets=F(1,4)<F(1,2),
            larger_axis_activation_still_before_modulation=F(1,20)/2<F(1,20),
            all_regions_covered_by_maximum=r['KU']>=axis and r['KU']>=26 and r['KU']>=352,
            regular_radial_average_cost=r['v2']==36*r['KU'],
            both_full_signed_source_bounds_positive=r['I0']>0 and r['I1']>0))

@lru_cache(None)
def review():
    parts=dict(exact=exact_budget(),symbolic=symbolic_review(),manufactured=manufactured_review(),regional=regional_review())
    return serial(dict(reviews=parts,checks={n:all(r['checks'].values()) for n,r in parts.items()}))
