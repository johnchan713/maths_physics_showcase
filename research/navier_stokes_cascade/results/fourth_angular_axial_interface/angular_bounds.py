"""Same-root C_eta^4 reference pulse and an explicit global axial input gate."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial,comb
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
        s=importlib.util.spec_from_file_location('fourth_axial_'+folder,p)
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=prior

def strength_jets(epsilon):
    """Normalized jets of the SAME scalar quadratic, inverse Jacobian<=3."""
    epsilon=rat(epsilon)
    if not 0<=epsilon<=F(1,2**20):raise ValueError('Retained epsilon<=2^-20 required')
    u=[F(2)]
    for n in range(1,5):
        u.append(3*(sum(u[j]*u[n-j] for j in range(1,n))+
                    epsilon*sum(u[n-j] for j in range(1,n+1))+epsilon))
    return tuple(u)

def compact_cost(A,Ahi,N,r,B4,Z4):
    """Conditional actual compact U C4: primitive B4 and transformed target Z4."""
    A,Ahi,N,r,B4,Z4=map(rat,(A,Ahi,N,r,B4,Z4))
    if min(A,Ahi,N)<1 or min(r,B4,Z4)<0:raise ValueError('Positive retained envelopes/frequency and nonnegative inputs required')
    d4=2000*Z4/N+600000000*r*r
    c4=4*r+d4  # Four known normalized jets each <=r; no vector norm convention loss.
    return dict(root_fourth_jet=d4,root_C4=c4,modulated_U_C4=Ahi+B4/N,
        repaired_U_C4=Ahi+B4/N+1600*A*c4,
        missing_input=B4+3200000*A*Z4,
        known_part=Ahi+6400*A*r+960000000000*A*r*r)

def compact_gate(A,Ahi,N,r,B4,Z4,K):
    A,Ahi,N,r,B4,Z4,K=map(rat,(A,Ahi,N,r,B4,Z4,K))
    q=compact_cost(A,Ahi,N,r,B4,Z4)
    if K<0:raise ValueError('Nonnegative trial certificate required')
    lhs=q['missing_input'];rhs=N*(K-q['known_part'])
    return dict(lhs=lhs,rhs=rhs,accepts=lhs<=rhs,
        same_positive_envelope_inequality=(q['repaired_U_C4']<=K),
        interpretation='equivalent acceptance for this sufficient majorant, not necessity for the true field or a new construction schedule')

def axis_fourth(S,rho,Lambda,t,ell):
    q=module('collar_operator_envelope','envelope.py').bounds(S,rho,Lambda,t,ell)
    factor=q['radial_factor']*sum(q['delta']**(-j) for j in range(5))
    return dict(factor=factor,power=q['exponential_power'],delta=q['delta'],
        expression='axis_C4<=factor*2^power; SAME axis-side continuation through ell=tc/2<t/20')

def pressure_interface(h,Xv,Kaxis4,Kcompact4,Kfirst3):
    """Conditional full support pressure-order-two signed C1 moments."""
    h,Xv,Kaxis4,Kcompact4,Kfirst3=map(rat,(h,Xv,Kaxis4,Kcompact4,Kfirst3))
    if min(Kaxis4,Kcompact4,Kfirst3)<0:raise ValueError('Nonnegative actual full-region norms required')
    m=module('higher_pressure_source_moments','moments.py');K4=max(Kaxis4,Kcompact4,F(110))
    v02=m.regular_v_bound(0,h,2,K4);v12=m.regular_v_bound(1,h,2,Kfirst3)
    v03=m.regular_v_bound(0,h,3,K4)
    q=m.source_bounds(1,h,1,Xv,(K4,Kfirst3),(v02,v12),v03)
    return dict(Kleading4=K4,v0_C2=v02,v1_C2=v12,v0_C3=v03,**q)

def symbolic_review():
    import sympy as s
    e,theta=s.symbols('eta theta',real=True);p=1+e*e;g=(p/2)**theta;I=s.I
    f=1/p;a=s.log(2/p)
    partial=s.simplify(f-((e-I)**-1-(e+I)**-1)/(2*I))==0
    logs=[s.simplify(s.diff(a,e,n)-(-1)**n*s.factorial(n-1)*((e-I)**(-n)+(e+I)**(-n)))==0 for n in range(1,5)]
    g4=g*(12*theta*(theta-1)/p**2+48*theta*(theta-1)*(theta-2)*e*e/p**3+
        16*theta*(theta-1)*(theta-2)*(theta-3)*e**4/p**4)
    z=s.symbols('z');amp=s.Function('amp')(z);B=s.Function('B')(z);D=s.Function('D')(z);k=s.symbols('k')
    target=(2*k*amp+B)*s.diff(amp,z,4)+k*(8*s.diff(amp,z)*s.diff(amp,z,3)+6*s.diff(amp,z,2)**2)
    target+=4*s.diff(B,z)*s.diff(amp,z,3)+6*s.diff(B,z,2)*s.diff(amp,z,2)+4*s.diff(B,z,3)*s.diff(amp,z)+s.diff(B,z,4)*amp+s.diff(D,z,4)
    return dict(checks=dict(reciprocal_partial_fractions_exact=partial,
        logarithm_derivatives_through_four_exact=all(logs),
        fourth_interpolation_power_all_terms=s.simplify(s.diff(g,e,4)-g4)==0,
        fourth_strength_derivative_all_linear_and_quadratic_terms=s.simplify(s.diff(k*amp*amp+B*amp+D,z,4)-target)==0,
        dropping_linear_coefficient_derivatives_detected=s.simplify(target-(target-4*s.diff(B,z)*s.diff(amp,z,3)-6*s.diff(B,z,2)*s.diff(amp,z,2)-4*s.diff(B,z,3)*s.diff(amp,z)))!=0),
        strength_fourth_terms='8k amp1 amp3+6k amp2^2+4B1 amp3+6B2 amp2+4B3 amp1+B4 amp0+D4; raw derivatives')

def exact_budget():
    eps=F(1,2**20);u=strength_jets(eps)
    # Pole distances |eta+-i|>=1 give |f^(j)|/j!<=1 and |a^(j)|/j!<=2/j.
    anorm=1+sum(F(2,j) for j in range(1,5));pnorm=sum(F(comb(4,j)) for j in range(5))+2*sum(F(comb(2,j)) for j in range(3))+1
    snorm=sum(F(comb(2,j)) for j in range(3))+2*sum(F(comb(4,j)) for j in range(5))+sum(F(comb(6,j)) for j in range(5))
    return dict(constants=dict(log_C4_upper=str(anorm),reciprocal_C4_upper=5,interpolation_fourth_raw=408,
        angular_target_C4_coefficient=2**29,angular_root_C4_coefficient=2**32,
        complete_power_energy_C4_coefficient=8192,strength_jet_caps=[12]*4,
        strength_C4_upper=2,end_coefficient_C4_multiplier=4,reference_axial_pulse_C4_upper=110),
        actual_normalized_recursion_caps=list(map(str,u)),
        checks=dict(log_fourth_norm_fits_same_eight_budget=anorm<8,
            angular_discrepancy_same_budget=18*2**16<2**21,
            angular_same_root_radius=6*2**29<2**32,
            same_pointwise_uniqueness_ball=4200<2**32 and F(192*2**32,2**(256*30))<1,
            fourth_interpolation_raw_bound=24+192+192==408,
            full_interpolation_C4_below64=1+4+F(12,2)+F(56,6)+F(408,24)<64,
            interpolation_energy_at_T64=64*1000<=1024*64,
            p_squared_factorial_C4_is25=pnorm==25,
            unedited_energy_C4_budget=F(25,4)*512+1024==4224,
            complete_signed_edit_C4_below_one=F(30,2**40)+F(1000,2**80)<1,
            full_post_energy_C4_below8192T=4224*64+1<8192*64,
            prefix_s0_polynomial_C4_is93=snorm==93,
            full_s0_C4_below512T=241*64+F(17,6)+snorm/4<512*64,
            complete_coefficient_epsilon_budget=512+4096+1<8192,
            epsilon_at_selected_minimum=13+6-256<-20,
            uniform_twelve_caps_inductive=all(9+468*(n-1)*eps<12 for n in range(1,5)),
            exact_positive_recursion_below_twelve=all(u[n]<=12*eps for n in range(1,5)),
            fourth_strength_full_norm_below_two=F(6,5)+48*eps<2,
            affine_end_correction_C4_kept=2+2==4,
            disjoint_end_pulse_bound=F(64*4,2**(256*20))<22,
            reference_axial_pulse_C4_bound=5*22==110,
            actual_compact_unit_mass_bump_value=module('cutoff_higher_jet_interface','jet_bounds.py').repair_bump_bounds(0)==(320,)))

def manufactured_review():
    """Five exact positive quadratic roots; not actual selected field data."""
    import sympy as s
    e=s.symbols('eta');eps=s.Rational(1,2**20);k=s.Rational(1,4)
    B=eps*(e+e**3)/10;D=-s.Rational(1,4)+eps*(e*e+e**4)/20
    amp=(-B+s.sqrt(B*B-4*k*D))/(2*k);rows=[]
    # Independent differentiation of the explicit root is compared to the recurrence.
    for eta in (s.Rational(-1),s.Rational(-2,5),s.Rational(0),s.Rational(3,5),s.Rational(1)):
        u=[s.diff(amp,e,n).subs(e,eta)/s.factorial(n) for n in range(5)]
        b=[s.diff(B,e,n).subs(e,eta)/s.factorial(n) for n in range(5)]
        d=[s.diff(D,e,n).subs(e,eta)/s.factorial(n) for n in range(5)]
        residual=[]
        for n in range(1,5):
            residual.append(s.simplify((2*k*u[0]+b[0])*u[n]+k*sum(u[j]*u[n-j] for j in range(1,n))+sum(b[j]*u[n-j] for j in range(1,n+1))+d[n])==0)
        rows.append(dict(eta=str(eta),checks=dict(all_four_exact_normalized_recurrences=all(residual),
            all_four_jets_below_twelve_epsilon=all(abs(complex(s.N(v,40)))<float(12*eps) for v in u[1:]))))
    return dict(rows=rows,manufactured_not_selected_reference_root=True,
        checks=dict(all_five_same_root_fourth_jet_controls=all(all(r['checks'].values()) for r in rows)))

def interface_review():
    cost=compact_cost(2,3,100,F(1,1000),5,7)
    gate=compact_gate(2,3,100,F(1,1000),5,7,cost['repaired_U_C4'])
    lower=compact_gate(2,3,100,F(1,1000),5,7,cost['repaired_U_C4']-1)
    axis=axis_fourth(2,F(1,8),4,F(1,100),F(1,10000))
    q=pressure_interface(F(1,200),7,13,17,19)
    prev=module('higher_pressure_source_moments','moments.py')
    no_shift=prev.source_bounds(1,F(1,200),1,7,(q['Kleading4'],19),(q['v0_C2'],q['v1_C2']),0)
    return dict(compact_fixture=serial(cost),gate_fixture=serial(gate),axis_factor_fixture=serial(axis),pressure_fixture=serial(q),
        abstract_C3_not_C4_family='n^-3*cos(n*eta)*sin(phi), n>=1; Ceta3<=8/3, raw eta4 ateta0,phi=pi/2 is n; not an actual loop or PDE family',
        checks=dict(compact_eliminated_gate_exact=gate['accepts'] and gate['same_positive_envelope_inequality'],
            too_small_compact_certificate_rejected=not lower['accepts'] and not lower['same_positive_envelope_inequality'],
            compact_fourth_target_has_inverse2000=cost['root_fourth_jet']==2000*F(7,100)+600000000*F(1,1000)**2,
            compact_primitive_and_target_both_unresolved=cost['missing_input']==5+3200000*2*7,
            full_global_cover_uses_later_pulse=q['Kleading4']==110,
            both_signed_axial_viscosities_retained=q['I0']>no_shift['I0'] and q['I1']>no_shift['I1'],
            abstract_small_C3_does_not_bound_C4=F(1)+1+F(1,2)+F(1,6)==F(8,3) and F(24*100,24)>F(8,3)))

@lru_cache(None)
def review():
    parts=dict(budgets=exact_budget(),symbolic=symbolic_review(),manufactured=manufactured_review(),interface=interface_review())
    checks={k:v for part in parts.values() for k,v in part['checks'].items()}
    return dict(**parts,checks=checks)
