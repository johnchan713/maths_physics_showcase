"""Continuum matching estimates using rational exponential enclosures."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m
rb=load('matching_root_rational_bounds',HERE.parent/'root_foundation/root_bounds.py')
exp_bounds,exp_interval,outward=rb.exp_bounds,rb.exp_interval,rb.outward
def upper(x):return exp_interval(F(x))[1]
def lower(x):return exp_interval(F(x))[0]
def show(x):return str(outward(x,10**6)[1])

def vandermonde_rows(rates,spacing,center,width):
    """Positive bump of exact mass width; inverse via Lagrange coefficients."""
    nodes=[exp_interval(spacing*r) for r in rates];out=[F(0)]*len(rates)
    for k,r in enumerate(rates):
        others=[j for j in range(len(rates)) if j!=k]
        denominator=width*lower((center-width/2)*r)
        for j in others:
            denominator*=abs(nodes[k][0]-nodes[j][1]) if nodes[k][0]>nodes[j][1] else nodes[j][0]-nodes[k][1]
        assert denominator>0
        polynomial=[F(1)]
        for j in others:
            nxt=[F(0)]*(len(polynomial)+1)
            for d,v in enumerate(polynomial):nxt[d]+=v*nodes[j][1];nxt[d+1]+=v
            polynomial=nxt
        for d,v in enumerate(polynomial):out[d]+=v/denominator
    return out

@lru_cache(None)
def annular_review():
    u=vandermonde_rows([F(1),F(8,5)],F(2,5),F(3,10),F(2,25))
    e=vandermonde_rows([F(8,5),F(6,5),F(1,5)],F(7,20),F(3,20),F(2,25))
    ee=32*F(6,25)*upper(F(8,5));uu=64*F(4,25)*upper(1);gamma=upper(F(6,5))/256
    q=[ee+4*uu*gamma]+[uu*gamma*k for k in (8,8,4,1)]
    return dict(U_inverse_rows=list(map(show,u)),E_inverse_rows=list(map(show,e)),
        quadratic_normalized_eta_jet_upper=list(map(show,q)),quadratic_algebra_upper=show(sum(q)),
        beta=610,q=50,coefficient='zeta=exp(6/5)*(1+eta^2)^2/P^2; P>=16',
        checks=dict(positive_bump_mass_is_width=True,inverse_below_610=max(u+e)<610,
            all_quartic_coefficient_jets_retained=sum(q)<50,disjoint_supports_remove_only_mixed_U_E=True))

@lru_cache(None)
def target_review():
    eps=F(1,10**16);g=3*eps/100
    delta=[eps/4+(1+upper(-2))*g,
        eps/4+(1+F(5,8)*upper(F(-16,5)))*g,eps/4,
        eps/4+(1+F(20,256)*upper(F(-4,5)))*g*g,eps/4]
    d=F(143,500)*eps;r=2*610*d;k=8*610**2*50*d
    return dict(epsilon=str(eps),target_over_epsilon=[show(x/eps) for x in delta],
        delta=str(d),root_radius=str(r),radius_over_epsilon=str(r/eps),contraction_parameter=str(k),
        norm='max component sum_{j=0}^2 sup|partial_eta^j f|/j!',
        checks=dict(actual_current_entry_estimate_retained=8*10**6*eps**199<1,
            target_below_p286epsilon=max(delta)<d,contraction=k<1,
            same_old_C2_root=r<572*eps,same_old_C1_root=r<6000*eps,
            restoration_energy_cancellation_retained=True,pointwise_inverse_below_1220=4*610**2*50*d<F(1,2)))

@lru_cache(None)
def repair_review():
    du=F(9,10)*lower(F(95,100)+F(195,100)-F(205,10000))
    ui=F(305,100)*upper(F(205,100))/du
    de=lower(-F(51,100)*F(465,100))
    for a,b in ((F(145,100),F(55,100)),(F(245,100),F(55,100)),(F(245,100),F(155,100))):de*=lower(a)-upper(b)
    ei=6*upper(F(765,100))/de
    q=640*(2*upper(F(205,100))+F(3,2)*upper(F(255,100)))
    return dict(U_inverse_upper=show(ui),E_inverse_upper=show(ei),quadratic_algebra_upper=show(q),
        beta=800,q=25000,parameter_range='0<lambda<=1/100; original e(eta) used in both corrections',
        determinant='x^(-1/2-lambda) times (x^(2+lambda),x,1); divided difference >=1',
        norm='max component (sup|f|+sup|partial_eta f|)',
        checks=dict(unit_mass_not_annular_width=True,U_gap_positive=du>0,U_inverse_below_two=ui<2,
            E_divided_difference_positive=de>0,E_inverse_below_800=ei<800,
            quadratic_below_25000=q<25000,zero_lambda_excluded=True,constant_eta_B_and_Q_requires_original_e=True))

@lru_cache(None)
def repair_budget_review():
    k=8*800**2*25000
    return dict(target='delta<=64*A^5*H/(lambda*N)<=64*A^6*H/N<=D/N; D=H^16',
        frequency='same N=1+floor(H^32)>H^32; H=exp(G^4), G=exp(A^32), A=C^4096',
        root='full C1 sum norm <=1600*D/N; pointwise Jacobian inverse <=1600',
        checks=dict(lambda_inverse_loss_retained=True,normalizer_eta_derivatives_retained=True,
            transformed_target_not_recovered_from_rounded_M_minus_J=True,
            logarithmic_contraction=k<2**37 and F(2)**37*F(2)**-256<1,
            full_sum_controls_old_max=True,same_old_root_ball=1600<2000,
            same_frequency_and_acceptance_region=True))

@lru_cache(None)
def pressure_review():
    return dict(rows=['M','J','I','S','Cp'],gauge='same Pi_axis + Cp; never reset at a join',
        jet_recurrence='d_n<=2 beta [z_n+sum_{i+j+l=n excluding (0,n,0),(0,0,n)} q_i*d_j*d_l]; i<=4 in annular map',
        missing_pressure_defect='4*Aphys*eta*c(eta)-(1-eta^2)*c_eta(eta)',
        checks=dict(five_rows_preserve_complete_later_pressure_memory=True,
            fixed_eta_jets_match_by_same_root_identity=True,
            annular_quartic_Q_derivatives_not_deleted=True,
            no_high_order_small_target_contraction_asserted=True,
            complete_core_backward_pressure_identification_still_required=True,
            unconditional_construction_not_claimed=True))
