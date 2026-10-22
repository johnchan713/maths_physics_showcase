"""Low-order diagonal summation and the complete tangential normal floor."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('tangential_head',HERE.parent/'core_background_finite_head_budget/budgets.py')
head=importlib.util.module_from_spec(spec);spec.loader.exec_module(head)
exact=head.exact


def diagonal_tail(t,m=1,start=1):
    """Sum of majorants, NOT a measured fluid jet. t=q^h in (0,1]."""
    t=exact(t)
    if type(m) is not int or m not in (0,1):raise ValueError('Only the head-free orders zero and one are certified')
    if type(start) is not int or start<1 or not 0<t<=1:raise ValueError('Positive tail start and 0<t<=1 required')
    return (t/2)**start/(1-t/2)


def comparison(Lead1,Rmin,Rmax,U,lambda_max,gmin,Llow,Lhigh,Dbox):
    """No first coefficient, K2, leading second jet, or schedule is supplied anew.

    The SAME complete Euler diagonal for EVERY n>=1 remains a premise.
    Leading first tuple norm and physical chart/value hypotheses remain premises.
    """
    lead,lo,hi,U,lam,g,Llo,Lhi,D=map(exact,(Lead1,Rmin,Rmax,U,lambda_max,gmin,Llow,Lhigh,Dbox))
    if lead<0 or min(lam,g,Llo,Lhi)<=0 or min(U,D)<0 or Lhi<Llo:raise ValueError('Valid leading and unchanged whole-path geometry bounds required')
    chart=head.chart_constants(1,lo,hi);fg=chart['F']+chart['G']
    carriers=head.phase.carrier_bounds(lam,U,g,hi,Llo)
    p,pz=carriers['theta'],carriers['axial'];B=max(Q(1),lam)
    groups=dict(tilt=2*B*U/(Llo*g),angular_rounding=1/lo,
                angular_box=(p-1)*D/lo**2,axial_transport=Lhi*(p+pz)*fg*(lead+1))
    return dict(Ctan=sum(groups.values(),Q(0)),groups=groups,carriers=carriers,
                chart=chart,leading_FG=fg*lead,full_FG=fg*(lead+1),error_FG_at_epsilon=2*fg,
                lambda_max=lam,tilt_upper=U,
                independent_inputs_removed=['first coefficient-one tuple norm','Borel cutoff K2','leading second jet'],
                retained_premise='same complete Euler diagonal C_(n,m)*a_n^(-nh)<=2^-n for every n>=1 and m<=1; same actual chart and leading gaps',
                actual_infinite_schedule_verified=False,actual_normal_floor_verified=False)


def tangential_error(result,S,epsilon,epsilon2,k_inverse):
    S,e,e2,ki=map(exact,(S,epsilon,epsilon2,k_inverse))
    if S<1 or min(e,e2,ki)<0:raise ValueError('Whole-path scale and nonnegative rounding data required')
    g=result['groups']
    # Component-sum majorant retains BOTH tilt errors and full H_Z transport.
    return g['tilt']/S+g['angular_rounding']*ki+g['angular_box']/S**3+g['axial_transport']*S*e


def normal_gate(result,lambda_min,S,epsilon,epsilon2,k_inverse):
    lam,S,e,e2,ki=map(exact,(lambda_min,S,epsilon,epsilon2,k_inverse))
    if lam<=0 or lam>result['lambda_max']:raise ValueError('Consistent positive leading growth floor required')
    err=tangential_error(result,S,e,e2,ki);U=result['tilt_upper'];C=result['Ctan']
    beta=lam/(4*(1+U*U)**2)
    slow=S*S*(e+e2+ki)<=1
    return dict(beta=beta,inverse_normal_square=4/beta,error_upper=err,
                gates=dict(same_whole_path_slow_gate=slow,tangential_margin=4*C*C<=beta*S*S),
                error_le_Ctan_over_S_if_slow=not slow or err<=C/S,
                existing_scale=S,scales_or_length_selected_or_changed=False,
                actual_normal_floor_verified=False)


def length_compatibility(result,lambda_min,radius,Lreq):
    lam,r,L=map(exact,(lambda_min,radius,Lreq))
    if min(lam,r,L)<=0 or lam>result['lambda_max']:raise ValueError('Consistent positive original selector data required')
    beta=lam/(4*(1+result['tilt_upper']**2)**2)
    margin=beta*L*L-16*r*r*result['Ctan']**2
    return dict(beta=beta,margin=margin,sufficient_gate=margin>=0,
                actual_existing_selector_verified=False,scales_or_length_selected_or_changed=False)
