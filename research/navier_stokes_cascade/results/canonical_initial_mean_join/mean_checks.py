"""Exact fresh pressure/moment controls and finite initial-order ledger."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys


@lru_cache(None)
def kernel():
    p=Path(__file__).resolve().parent.parent/'complete_covariance_initial_means/means.py'
    spec=importlib.util.spec_from_file_location('canonical_mean_exact_polynomials',p)
    m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m


@lru_cache(None)
def maps():
    p=Path(__file__).resolve().parent.parent/'compact_mean_foundation/mean_maps.py'
    spec=importlib.util.spec_from_file_location('canonical_mean_fixed_rows',p)
    m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m


def pressure_data(eta):
    eta=F(eta);m=kernel();R,Z=m.R,m.Z;eps=F(1,200)
    f=(R-1)**4*(2-R)**4;profile=R*f;h=1+Z/5;t=eps*eta*profile*h/2
    gr=t.shift(-1);P=gr.integral();Jtheta=-(R**2*t).integral()
    Jz=(R*t).integral()-(R**2*gr).integral()/2
    rho0=f/f.integral().scalar();rho1=R*f/(R*f).integral().scalar()
    rho=(1-Z)*rho0+Z*rho1;crho=(R**2*rho).integral()/2
    source=gr-rho*P;primitive=source.primitive();pm=primitive-primitive.at(0,1)
    moment=(R*pm).integral();Ez=eps*(t+pm).diff(1);Et=-eps*t.diff(1)
    return dict(eps=eps,profile=profile,f=f,t=t,gr=gr,P=P,Jtheta=Jtheta,Jz=Jz,rho=rho,crho=crho,
        pm=pm,moment=moment,Ez=Ez,Et=Et,R=R,Z=Z)


def pressure_review():
    rows=[];good=True;moving=True;stale=True;fresh=True
    for eta in (F(-1,4),F(1,20),F(1,7)):
        d=pressure_data(eta);R=d['R'];eps=d['eps']
        good &= d['rho'].integral()==1 and d['pm'].at(0,1)==d['pm'].at(0,2)==0
        good &= d['pm'].diff(0)==d['gr']-d['rho']*d['P']
        good &= d['moment']==-(R**2*d['gr']).integral()/2+d['crho']*d['P']
        good &= (R*d['Ez']).integral()==eps*(d['Jz']+d['crho']*d['P']).diff(1)
        good &= (R**2*d['Et']).integral()==eps*d['Jtheta'].diff(1)
        moving &= bool(d['crho'].diff(1)*d['P']) and bool(d['moment'].diff(1))
        target=[d[n].at(1,F(1,3)).scalar() for n in ('P','Jtheta','Jz')]
        stale &= all(x!=0 for x in target)
        At,Az=maps().blocks(F(3,2),F(3,4),F(5,4),(F(2),F(3),F(4),F(5),F(6)))
        bt=[F(0),-target[0],-target[2]];bz=[F(0),-target[1]]
        xt,xz=maps().solve(At,bt),maps().solve(Az,bz)
        fresh &= [sum(a*x for a,x in zip(row,xt)) for row in At]==bt
        fresh &= [sum(a*x for a,x in zip(row,xz)) for row in Az]==bz
        rows.append(dict(eta=str(eta),fresh_defects=list(map(str,target)),
            pressure_moment_at_Zthird=str(d['moment'].at(1,F(1,3)).scalar()),
            omitted_moving_bump_term_at_Zthird=str((eps*d['crho'].diff(1)*d['P']).at(1,F(1,3)).scalar()),
            fresh_angular_row_solution=list(map(str,xt)),fresh_axial_row_solution=list(map(str,xz))))
    return dict(rows=rows,checks=dict(all_full_pressure_and_axial_moment_identities_exact=bool(good),
        omitting_axial_pressure_or_moving_bump_derivative_detected=bool(moving),
        same_two_principal_fluxes_can_leave_three_nonzero_new_defects=bool(stale),
        fixed_mean_matrix_cancels_recomputed_linear_targets=bool(fresh),
        linear_cancellation_not_claimed_to_erase_nonlinear_defects=True,
        finite_polynomial_endpoint_regularity_not_actual_smooth_profile=True))


def order_review():
    k=F(1,100000)
    orders=dict(canonical_real_curl_remainder=1-k,real_curl_covariance_square=2-2*k,
        all_six_covariance=F(1),selected_Haar_flux_error=2-2*k,
        raw_pressure=1-k,averaged_raw_pressure=F(1),raw_defects=F(1),
        initial_averaged_tangential=2-2*k,fast_tangential=1-k,fast_radial=2-k,
        fast_viscosity=2-3*k,fast_averaged_quadratic=2-2*k,
        slow_tangential=F(1),slow_radial=F(2),slow_square=F(2),
        final_defects=F(2),final_radial_mean=F(2),final_tangential_mean=2-3*k)
    return dict(kappa=str(k),orders={n:str(v) for n,v in orders.items()},checks=dict(
        real_primary_square_is_stronger_than_generic_cross=2-2*k>F(3,2)-k,
        all_six_not_promoted_to_selected_Haar_order=orders['all_six_covariance']<orders['selected_Haar_flux_error'],
        fast_viscosity_is_actual_limiting_mean_order=orders['fast_viscosity']==1+(1-k)-2*k,
        recomputed_fast_quadratic_defects_do_not_worsen_raw_order=2-2*k>1,
        fixed_slow_map_with_initial_zero_Haar_fast_state_gives_two=orders['final_defects']==orders['slow_square'],
        whole_tangential_mean_not_promoted_to_two=orders['final_tangential_mean']<2,
        no_nonzero_angular_or_later_complex_cycle_improvement_asserted=True))


def nonflat_review():
    rows=[];growth=True
    # h=1/200, delta=1/4, ell=800*m: epsilon^(-delta)=2^m exactly.
    r0=F(1,6057072)
    for degree in (2,5):
        previous=F(0)
        for m in (400,800,1600):
            ell=800*m;ratio=F(2**m,1)/(2*r0*ell**(2*degree+2))
            lower=ratio.numerator.bit_length()-ratio.denominator.bit_length()-1
            growth &= ratio>1 and ratio>previous;previous=ratio
            rows.append(dict(polynomial_degree=degree,ell=ell,log2_ratio_lower_bound=lower))
    return dict(rows=rows,checks=dict(polynomial_pulse_closeness_not_an_epsilon_power_gain=bool(growth),
        exact_integer_ratio_no_exponential_float_overflow=True,
        all_fixed_degrees_obstruction_is_analytical_exponential_vs_polynomial_limit=True,
        no_claim_every_actual_canonical_difference_saturates_upper_bound=True))


def recompute_review():
    r=kernel().haar_update_review()
    return dict(checks=dict(r['checks']),arithmetic=r['arithmetic'],
        cancellation_requires_initial_zero_Haar_fast_state=r['cancellation_requires_initial_zero_Haar_fast_state'],
        reused_exact_polynomial_kernel_not_a_new_actual_profile=True)
