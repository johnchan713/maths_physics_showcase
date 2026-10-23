"""Sharp finite-jet obstruction for the ORIGINAL compact drift operators.

These are scalar comparison sources, never substituted for canonical fields.
No large profile seminorm or current-source norm is materialized here.
"""
from fractions import Fraction as Q
from math import ceil
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

parent=load('compact_sharp_parent',HERE.parent/'core_next_sourced_mean_tail_join/budgets.py')
exact,order=parent.exact,parent.order

def pell(j):
    """(P+Q sqrt2)=(1+sqrt2)^j; k=(P-Q,Q), |k|_1=P."""
    order(j)
    if j<1:raise ValueError('Positive Pell index required')
    P,Qv=1,0
    for _ in range(j):P,Qv=P+2*Qv,P+Qv
    return dict(index=j,P=P,Q=Qv,k=(P-Qv,Qv),norm=P*P-2*Qv*Qv)

def select_mode(M):
    """Works for ANY M>=1, including the unchanged rounded band values.

    The loop is logarithmic in M, not in a gigantic derivative request.
    It gives sqrt(M)<=P<3sqrt(M), |M omega|<=P, Q>=P/2.
    """
    M=exact(M)
    if M<1:raise ValueError('The asymptotic comparison uses M>=1')
    P,Qv,j=1,1,1
    while P*P<M:P,Qv,j=P+2*Qv,P+Qv,j+1
    return dict(M=str(M),index=j,P=P,Q=Qv,k=(P-Qv,Qv),norm=P*P-2*Qv*Qv,
        exact_omega='P-sqrt(2)*Q = norm/(P+sqrt(2)*Q)',
        squared_size_ratio=Q(P*P)/M)

def finite_jet_obstruction(m,J,p):
    """Uniform N_m(A_e f)<=C M^-p N_J(f) needs J>=m+2p.

    Necessity concerns this ISOTROPIC FINITE-JET estimate over all compact
    zero-Haar sources, not all proofs and not the canonical current source.
    """
    order(m);order(J);p=exact(p)
    if p<0:raise ValueError('Nonnegative desired compact gain required')
    gap=Q(m)+2*p-J
    return dict(output_jet=m,source_jet=J,M_gain=p,
        ratio_growth_in_P=gap,uniform_estimate_obstructed=gap>0,
        necessary_source_order=ceil(m+2*p),
        sufficient_source_order=(m+2*int(p)+3 if p.denominator==1 and p>=1 else None),
        fixed_profile_output_factor='abs(k_e(U*) chi_U(U*)) * integral(phi)>0',
        source_raw_norm_upper='16^J A_J, for amplitude P^-J and A_J=N_J(f_e^-1 phi B)',
        canonical_source_or_all_jet_counterexample=False)

def epsilon_target(m,J,a,target,kappa=Q(1,100000)):
    """Finite-jet necessity for epsilon^a -> epsilon^target, up to S powers."""
    order(m);order(J);a,target,k=map(exact,(a,target,kappa))
    if k<=0 or target<a:raise ValueError('Positive kappa and nonnegative target gain required')
    threshold=Q(m)+2*(target-a)/k
    exponent=a+k*(J-m)/2
    return dict(output_jet=m,source_jet=J,source_order=a,target_order=target,
        resonant_output_epsilon_order=exponent,
        epsilon_gap=target-exponent,
        necessary_source_order=max(0,ceil(threshold)),
        finite_jet_route_obstructed=exponent<target,
        equality_case_decided=False,
        S_factor='S^(rho_g*(J-m)/2), within a fixed Tg-dependent factor',
        all_jet_or_actual_source_failure=False)

def compare_request(e,m,*,physical_power=None,h=Q(1,200),kappa=Q(1,100000)):
    """Compare the sharp obstruction with the parent's retained raw rule.

    Physical rows use the parent's SUFFICIENT transfer loss; this is NOT
    a necessary lower bound for every physical/anisotropic estimate.
    """
    if e not in (0,1):raise ValueError('Original A0 or A1 required')
    order(m);r=parent.rates(kappa);a=r['B'] if e==0 else r['E']
    if physical_power is None:
        target=r['C0'];p=r['A0_integrations'] if e==0 else r['A1_integrations']
    else:
        d=parent.compact_physical_orders(e,m,physical_power,h,kappa)
        target=d['normalized_target'];p=d['integration_order']
    J=m+2*p+3;v=epsilon_target(m,J,a,target,kappa)
    return dict(operator='A'+str(e),physical_power=physical_power,
        source_order=a,normalized_target=target,retained_integrations=p,
        retained_raw_source_order=J,necessary_finite_jet_source_order=v['necessary_source_order'],
        excess_over_necessary_order=J-v['necessary_source_order'],
        strict_epsilon_margin=a+p*exact(kappa)-target,
        scope='isotropic scalar compact estimate on the retained transfer route; not a canonical or physical lower bound')

def actual_frontier():
    return dict(parent_commit='df3b93b4f07489e4f2806f359200dc34f5537d84',
        actual_required_inequality='N_J(Qcirc gr_final) <= epsilon^B P_J(S+K_all), and N_J(Qcirc Epostwave_z) <= epsilon^E Q_J(S+K_all), on every original band and complete support hull',
        source_identification='gr_final uses old + actual new fast + fresh slow means and all six current eight-mode covariance entries; Epostwave_z uses full forced and reconstructed pressure',
        profiles='same f_e,k_e,chi,positive hull, raw fixed-profile derivatives and output support gap; constants depend on these fixed profiles',
        finite_jet_improvement_obstructed_by_resonance=True,
        possible_escape='derive additional structure or anisotropic/Fourier decay of the SAME canonical sources, rather than dropping requested isotropic jets',
        all_jet_class_warning='The P^-J fixture has only a uniformly bounded finite jet set; higher derivatives grow with M. It is not in the original all-fixed-jet M_a class with the same a.',
        cutoff_support_and_zero_Haar_proved_for_comparison_sources=True,
        actual_current_source_high_jet_norms_evaluated=False,
        actual_input_validity_residual_decay_or_blowup_certified=False)
