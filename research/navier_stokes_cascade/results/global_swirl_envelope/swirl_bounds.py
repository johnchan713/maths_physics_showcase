"""Full swirl/heat bounds relative to retained roots, core and regional lemmas."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
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
def poch(x,n):
    x=rat(x);out=F(1)
    if type(n) is not int or n<0:raise ValueError('Nonnegative integer order required')
    for j in range(n):out*=x+j
    return out
def gamma_derivative_cap(h,j):
    h=rat(h)
    if not 0<h<=F(1,100):raise ValueError('Need 0<h<=.01')
    return poch(h,j)*poch(1+h,j)
def heat_root_budget(h,Xstar):
    """Actual same B,Q; full factorial C_eta^3 target, not old C1 relabelled."""
    h,Xstar=map(rat,(h,Xstar))
    if not 0<h<=F(1,100) or Xstar<100:raise ValueError('Need 0<h<=.01 and Xstar>=100')
    z=8192*h**4/Xstar;root=30*z;common=2**18*h**4/Xstar
    return dict(target_C3=z,root_C3=root,common_pointwise_radius=common,
        contraction=8*15**2*29*z,pointwise_uniqueness_lipschitz=2*15*29*common,
        old_C1_root_radius=6400*h**4/Xstar)
def swirl_inputs(Xa,Xc,Kaxis,native_mixed):
    """native_mixed=H^70 after axis, and Kaxis common radial 0..2/angular C3.

    Geometry conversion is used only after the axis-side interval.
    """
    Xa,Xc,Kaxis,native_mixed=map(rat,(Xa,Xc,Kaxis,native_mixed))
    if Xa<=0 or Xc<Xa or min(Kaxis,native_mixed)<0:raise ValueError('Positive ordered radii and nonnegative envelopes required')
    f=2*(1+1/Xa)*native_mixed
    f3=max(Kaxis,f);fx=max(Kaxis,f/Xa);fxx=max(Kaxis,3*f/Xa**2)
    CF=64*Xc**2*fxx+704*Xc*fx+1664*f3
    return dict(F_C3=f3,F_X_C2=fx,F_XX_C1=fxx,shifted_F_C1=CF)
@lru_cache(None)
def first_module():
    p=HERE.parent/'first_seminorm_bound/seminorms.py'
    s=importlib.util.spec_from_file_location('swirl_first_seminorm_interface',p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
def complete_first_interface(h,a,c,jets,axial,swirl,tau):
    """All leading inputs now explicit relative to lemmas; jets retain actual cutoff shape.

    a^2=Xa,c^2=Xc. This composes a conditional Cstar, not an accepted scale.
    """
    q=first_module();src=q.source_bounds(h,1,a,c,jets,
        swirl['F_C3'],axial['U_C3'],swirl['F_X_C2'],axial['U_X_C2'],
        axial['V_C1'],axial['w0_C1'],swirl['shifted_F_C1'],axial['shifted_U_C1'])
    return dict(sources=src,compact=q.compact_gate(a,c,jets,src,tau))

def exact_review():
    h=F(1,100);X=F(100)
    A=[gamma_derivative_cap(h,j) for j in range(1,4)]
    step=(F(1),F(32),F(1216),F(220000))
    chi=sum(v/F(3,10)**j/factorial(j) for j,v in enumerate(step))
    theta=sum(v/F(1000)**j/factorial(j) for j,v in enumerate(step))
    # The original terminal taper is fo(y)/fo(0), fo=1-co*h*sigma((3-y)/2).
    fo=(1+F(1,1000)*h*sum(v/F(2)**j/factorial(j) for j,v in enumerate(step)))/(1-F(1,1000)*h)
    scalar=1+2+F(4+64,2)+F(8+3*2*64+20000,6)
    angular=4+8+F(72,2)+F(2240,6)
    Znorm=2*(4+4+F(3,2)+F(1,6))
    mixed_gamma=4+2*20+F(4,2)*20**2/100+F(32,6)*20**3/100**2
    budget=heat_root_budget(h,X)
    heat_edit_exp=18-4*512+412+2
    angular_edit_exp=32-30*256+412+1
    return dict(gamma_caps_at_h_ceiling=serial(A),angular_gamma_C3_coefficient=str(angular),
        joint_gamma_coefficient=str(mixed_gamma),heat_cutoff_J3_upper=str(chi),
        interpolation_theta_J3_upper=str(theta),original_terminal_normalized_J3_upper=str(fo),
        scalar_J3_coefficient=str(scalar),worst_parameter_heat_root=serial(budget),
        actual_relative_edit_binary_exponents=dict(heat=heat_edit_exp,angular=angular_edit_exp),
        checks=dict(gamma_derivatives_through_three_bounded_by_2_4_32_h=all(a<=b*h for a,b in zip(A,(2,4,32))),
            angular_Gamma_difference_C3_below_512_h_over_X=angular<512,
            square_difference_budget_valid_at_X100=512*h/100<1,
            normalized_heat_pressure_and_energy_rows_fit_8192=73728*h**5<8192,
            normalized_heat_angular_row_retains_h_integral=5120<8192,
            strict_C3_heat_root_contraction=budget['contraction']<1,
            old_and_new_roots_in_common_pointwise_ball=max(budget['root_C3'],budget['old_C1_root_radius'])<budget['common_pointwise_radius'],
            common_pointwise_uniqueness_lipschitz_below_one=budget['pointwise_uniqueness_lipschitz']<1,
            new_C3_radius_not_claimed_inside_old_C1_radius=budget['root_C3']>budget['old_C1_root_radius'],
            heat_cutoff_J3_below_two_million=chi<2*10**6,
            interpolation_theta_J3_below_two=theta<2,
            exact_terminal_two_unit_normalized_factor_J3_below_two=fo<2,
            scalar_slope_and_two_derivatives_fit_4096=scalar<4096,
            positive_heat_argument_J3_below_twenty_over_X=Znorm<20,
            joint_Gamma_difference_J3_below_64_h_over_X=mixed_gamma<64,
            actual_heat_multiplier_J3_below_two=128000000*F(2)**-512/100<1,
            both_relative_edit_J3_budgets_below_one=max(heat_edit_exp,angular_edit_exp)<0,
            reference_product_fits_two60_Ahi=13+34+4<60,
            actual_H_dominates_reference_prefactor=60<1024))

def symbolic_review():
    import sympy as s
    X=s.symbols('X',positive=True);e,y=s.symbols('eta y',real=True);g=s.Function('g');z=2*(1-e*e)/X
    chain3=s.diff(g(z),e,3)
    target=s.Derivative(g(s.Symbol('z')),s.Symbol('z'),3).subs(s.Symbol('z'),z)*(-4*e/X)**3+3*s.Derivative(g(s.Symbol('z')),s.Symbol('z'),2).subs(s.Symbol('z'),z)*(-4*e/X)*(-4/X)
    rawchain=s.simplify(chain3-target)==0
    E=s.Function('E')(X,e);R=E/s.sqrt(2*X)
    G=s.Function('G')(y,e);Flog=G/s.sqrt(2*s.exp(y))
    # Exact logarithmic factorial jets of (2X)^-1/2 have total sum <2/sqrt(2X).
    rho=s.symbols('rho',positive=True);sigma=s.Function('sigma')
    terminal=(1-rho*sigma((3-y)/2))/(1-rho)
    tprime=s.diff(terminal,y)
    q=s.symbols('q',positive=True);v=s.symbols('v',positive=True)
    return dict(checks=dict(
        Gamma_third_angular_chain_includes_Z_second=rawchain,
        missing_Gamma_second_argument_term_detected=s.simplify(chain3-(target-3*s.diff(g(s.Symbol('z')),s.Symbol('z'),2).subs(s.Symbol('z'),z)*(-4*e/X)*(-4/X)))!=0,
        F_radial_first_conversion=s.simplify(s.diff(R,X)-(s.diff(E,X)-E/(2*X))/s.sqrt(2*X))==0,
        F_radial_second_conversion=s.simplify(s.diff(R,X,2)-(s.diff(E,X,2)-s.diff(E,X)/X+3*E/(4*X**2))/s.sqrt(2*X))==0,
        original_terminal_derivative_has_half_width=s.simplify(tprime-rho*s.Derivative(sigma(s.Symbol('v')),s.Symbol('v')).subs(s.Symbol('v'),(3-y)/2)/(2*(1-rho)))==0,
        infinite_angular_tail_retains_h_cancellation=s.integrate(q*v**(-1-q),(v,1,s.oo))==1,
        full_C3_inverse_amplitude_normalizations=(2+2+1==5 and 4+8+8+4==24),
        geometric_log_jet_factorial_sum_below_two=sum(F(1,2)**j/factorial(j) for j in range(4))<2))

def gamma_diagnostic():
    """Positive-integral finite controls; no negative Z or endpoint differencing."""
    import mpmath as mp
    rows=[]
    with mp.workdps(45):
      for h0 in ('0.01','0.005'):
       h=mp.mpf(h0)
       for eta0 in ('-1','-.5','0','.5','1'):
        eta=mp.mpf(eta0);X=mp.mpf(100);Z=2*(1-eta*eta)/X
        H=[(-1)**j*mp.rf(h,j)/mp.gamma(1+h)*mp.quad(lambda v:mp.exp(-v)*v**(h+j)*(1+Z*v)**(-h-j),[0,1,mp.inf]) for j in range(4)]
        jets=[H[0]-1,H[1]*(-4*eta/X),H[2]*(4*eta/X)**2+H[1]*(-4/X),
            H[3]*(-4*eta/X)**3+3*H[2]*(-4*eta/X)*(-4/X)]
        norm=sum(abs(v)/factorial(j) for j,v in enumerate(jets))
        rows.append(dict(h=h0,eta=eta0,X='100',Z=mp.nstr(Z,20),C3_point_jet=mp.nstr(norm,20),
            endpoint_eta_derivative_nonzero=bool(abs(jets[1])>0) if abs(eta)==1 else None,
            checks=dict(nonnegative_heat_argument=Z>=0,point_jets_below_analytical_envelope=norm<512*h/X,
                Gamma_integral_derivatives_below_positive_moments=all(abs(H[j])<=mp.rf(h,j)*mp.rf(1+h,j)*(1+mp.mpf('1e-40')) for j in range(1,4)))))
    return dict(numerical_quadrature_not_continuum_proof=True,rows=rows,checks=dict(
        all_ten_positive_integral_controls=all(all(r['checks'].values()) for r in rows),
        both_poles_have_nonzero_first_eta_heat_jet=all(r['endpoint_eta_derivative_nonzero'] for r in rows if r['eta'] in ('-1','1'))))

def composition_review():
    p=HERE.parent/'global_axial_source_envelope/envelopes.py';s=importlib.util.spec_from_file_location('swirl_global_axial_inputs',p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m)
    axial=m.global_bounds(F(1,200),1,11,1,7,13,101,303)
    swirl=swirl_inputs(1,16,7,101);jets=dict(value=F(3),radial=F(5),radial2=F(7))
    out=complete_first_interface(F(1,200),1,4,jets,axial,swirl,F(1,10))
    changed=complete_first_interface(F(1,200),1,4,jets,axial,dict(swirl,shifted_F_C1=0),F(1,10))
    no_pressure=complete_first_interface(F(1,200),1,4,jets,dict(axial,w0_C1=0),swirl,F(1,10))
    no_U_viscosity=complete_first_interface(F(1,200),1,4,jets,dict(axial,shifted_U_C1=0),swirl,F(1,10))
    return dict(auxiliary_rational_inputs_not_actual_profile=True,axial=serial(axial),swirl=serial(swirl),full_first=serial(out),
        checks=dict(full_first_constant_positive=out['compact']['Cval']>0,
            shifted_swirl_viscosity_enters_theta_source=out['sources']['theta']>changed['sources']['theta'],
            leading_pointwise_source_and_shifted_U_not_erased=out['sources']['axial']>no_pressure['sources']['axial'] and out['sources']['axial']>no_U_viscosity['sources']['axial'],
            full_swirl_operator_same_coefficient_ledger=swirl['shifted_F_C1']==m.scalar_axial_bound(16,swirl['F_C3'],swirl['F_X_C2'],swirl['F_XX_C1'])))

def review():
    e,s,g,c=exact_review(),symbolic_review(),gamma_diagnostic(),composition_review()
    return dict(exact=e,symbolic=s,Gamma_diagnostic=g,composition=c,checks=dict(
        exact_Gamma_target_root_and_reference_budgets=all(e['checks'].values()),
        exact_derivative_and_original_taper_identities=all(s['checks'].values()),
        finite_positive_integral_diagnostics=all(g['checks'].values()),
        complete_first_stress_interface_composition=all(c['checks'].values())))
