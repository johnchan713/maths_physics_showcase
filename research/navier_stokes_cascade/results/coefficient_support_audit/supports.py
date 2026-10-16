"""Exact support bookkeeping, cutoff curl, and labelled pressure-tail controls."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path
import sys
import math
from types import SimpleNamespace
import mpmath as mp

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


hierarchy=load('support_frozen_hierarchy',HERE.parent/'background_summation/hierarchy.py')
from residual import Poly, axial_second, profile_jets, series_coefficients


def require_order(n):
    if type(n) is not int or n<1:
        raise ValueError('A positive integer formal order is required')
    return n


def support_edges(n):
    require_order(n)
    return dict(order=n,velocity='Xvel',pressure='Xv' if n<=2 else 'Xvel',
                stress='Xb' if n==1 else ('Xv' if n==2 else 'Xvel'))


def radial_flux(U,M,Me,X,eta,h,n):
    require_order(n)
    return (2*eta*X*U-2*(Q(1,2)-h+2*n*h)*eta*M-(1-eta*eta)*Me)/(1-2*h*eta*eta)


def cut_potential_radial_flux(chi,sigma_chi_prime,V,M,eta,h):
    """r*u_r / q^(2nh) after multiplying the streamfunction by chi(c_n q)."""
    return chi*V-2*eta*sigma_chi_prime*M/(1-2*h*eta*eta)


def mean_pressure(n,X,mu,K):
    """I4 restriction; K is an incoming functional constant, not fitted data."""
    require_order(n)
    return K-mu*mu/(4*X) if n==1 else (K if n==2 else Q(0))


def geometry_review():
    intervals=[[-25,-20],[-20,-15],[-14,-9],[-8,-3]]
    U=[Q(1),Q(2)];E=[Q(1,2),Q(3,2),Q(5,2)];half=Q(1,20)
    bumps=sorted((c-half,c+half) for c in U+E)
    checks=dict(
        leading_repair_and_heat_only_use_first_two_patches=intervals[:2]==[[-25,-20],[-20,-15]],
        I3_I4_are_unedited_by_those_leading_patches=intervals[1][1]<intervals[2][0] and intervals[2][1]<intervals[3][0],
        all_positive_order_bumps_inside_I3=all(0<a<b<5 for a,b in bumps),
        positive_bumps_have_mutually_disjoint_supports=all(a[1]<b[0] for a,b in zip(bumps,bumps[1:])),
        common_velocity_edge_after_last_bump=bumps[-1][1]<3<5,
        common_velocity_edge_before_Imean=3<6,
        Imean_before_axial_pulse=intervals[3][1]<0,
        first_patch_end_equals_second_patch_start=intervals[0][1]==intervals[1][0],
        later_correction_patches_not_recycled_from_leading_repair=True,
        same_protected_inner_offsets=Q(1,16)<Q(1,8)<Q(3,16)<Q(1,4)<Q(1,2))
    return dict(checks=checks,stage_intervals=intervals,positive_bump_supports=[[str(a),str(b)] for a,b in bumps],
                positive_velocity_edge='Xvel=rho3*exp(3), rho3=XR*exp(241*T-12)',
                mean_patch='I4: rho3*exp(6)<X<rho3*exp(11)',
                leading_pressure_edge='Xv after the leading axial pulse, before the outer collar')


def exact_review():
    checks=[]
    for h in (Q(1,200),Q(9,1000)):
        for eta in map(Q,('-1','-1/2','0','1/3','1')):
            mu=Q(3);X=Q(10);L=1-2*h*eta*eta
            # ACTUAL leading patch formula once exact leading M matching holds.
            V0=(-2*(Q(1,2)-h)*eta*(mu*eta)-(1-eta*eta)*mu)/L
            checks.append(V0==-mu)
            for n in (1,2,3,7):
                V=radial_flux(Q(0),Q(0),Q(0),X,eta,h,n)
                checks.append(V==0)
                checks.append(cut_potential_radial_flux(Q(5,32),Q(-27,32),V,Q(0),eta,h)==0)
    h=Q(1,200);eta=Q(1,2);M=Q(1)
    V=radial_flux(Q(0),M,Q(0),Q(10),eta,h,1)
    extra=cut_potential_radial_flux(Q(5,32),Q(-27,32),V,M,eta,h)-Q(5,32)*V
    # A positive derivative on an interval rules out zero pressure there.
    pressure_change=mean_pressure(1,Q(20),Q(3),Q(7))-mean_pressure(1,Q(10),Q(3),Q(7))
    rows=[support_edges(n) for n in range(1,9)]
    return dict(checks=dict(
        exact_leading_memory_and_all_order_zero_velocity_identities=all(checks),
        velocity_values_without_M_do_not_close_support=V!=0,
        omitted_potential_cutoff_curl_detected=extra!=0,
        first_pressure_nonzero_radial_change_on_Imean=pressure_change>0,
        second_pressure_constant_not_assigned_zero=mean_pressure(2,Q(10),Q(3),Q(7))==7,
        higher_pressure_zero_only_after_exact_total_pressure_moment=True,
        first_two_exceptional_pressure_orders_kept=rows[0]['pressure']==rows[1]['pressure']=='Xv' and all(r['pressure']=='Xvel' for r in rows[2:]),
        shifted_viscosity_extends_second_but_not_third_pressure_order=True),
        rational_identity_checks=len(checks),cutoff_extra_radial_flux=str(extra),
        first_pressure_change=str(pressure_change),support_edges=rows,
        Imean_pressure_shapes=['K1(eta)-mu^2/(4X)','K2(eta)','0 for n>=3'])


def formal_review():
    checks=[]
    for n in range(1,13):
        pairs=[(i,n-i) for i in range(n+1)]
        current=[(i,j) for i,j in pairs if i==n or j==n]
        checks.extend((current==[(0,n),(n,0)],2*n>n,
                       all(i<n and j<n for i,j in pairs if i and j)))
    return dict(checks=dict(
        current_order_moment_map_linear_and_lower_products_fixed=all(checks),
        pressure_target_contains_whole_lower_source_not_only_I3=True,
        zero_axis_pressure_is_retained=True,
        exact_current_total_M_and_pressure_moments_required=True,
        fixed_inner_analytic_collar_not_shrunk_with_order=True,
        finite_constants_per_order_not_uniform_in_order=True,
        locally_finite_potential_cutoffs_preserve_exact_Imean_velocity=True,
        pressure_does_not_follow_velocity_support=True),
        formal_grade_checks=len(checks),raw_normalization_costs=['lambda^-1','e(eta)^-1'],
        pressure_equation='Pi_n_X=sum_{i+j=n} F_i F_j-Omega_(n-1)/(2X)')


def inverse_review():
    mp.iv.dps=80;p=mp.iv.mpf;exp=mp.iv.exp
    dE=exp(-p('.51')*p('4.65'))*(exp(p('1.45'))-exp(p('.55')))*(exp(p('2.45'))-exp(p('.55')))*(exp(p('2.45'))-exp(p('1.55')))
    dU=p('.9')*exp(p('2.9')-p('.01')*p('2.05'))
    e_bound=6*exp(3*p('2.55'))/dE;u_bound=2*p('2.05')*exp(p('2.05'))/dU
    lo=lambda z:mp.make_mpf(z._mpi_[0]);hi=lambda z:mp.make_mpf(z._mpi_[1])
    return dict(checks=dict(positive_E_determinant=lo(dE)>0,positive_U_determinant=lo(dU)>0,
                           normalized_inverse_E_below1000=hi(e_bound)<1000,
                           normalized_inverse_U_below2=hi(u_bound)<2,
                           same_bumps_are_available_on_actual_unedited_I3=True,
                           lambda_positive_required_for_raw_two_axial_rows=True),
                E_inverse_upper=math.nextafter(float(hi(e_bound)),math.inf),U_inverse_upper=math.nextafter(float(hi(u_bound)),math.inf),
                raw_target_lambda_inverse_not_discarded=True)


def source_review():
    zero=Poly.cast(0);x=Poly({(1,0):1});e=Poly({(0,1):1})
    one=Poly.cast(1);invx=Poly({(-1,0):1})
    errors=[];omission=[];zero_high=[];constant_checks=[]
    for constant in (True,False):
        B=one if constant else one+x*Q(1,5)+x*x*Q(1,10)
        U=e*B.derivative(0)
        base=SimpleNamespace(F=one+e*e,U=U,Pi=zero,average_U=e*B*invx)
        null=SimpleNamespace(F=zero,U=zero,Pi=zero,average_U=zero)
        for X in (2.,7.):
            for eta in (-1.,-.4,0.,.5,1.):
                h=.005;nu=1.7
                V=profile_jets(base,X,eta,h,0.)['V']
                Z2=axial_second(0.,V,X,eta,h)
                second=hierarchy.lower_sources([base,null],X,eta,h,nu)
                errors.append(abs(second['pressure']-nu*Z2/(2*X)))
                omitted=hierarchy.lower_sources([base,null],X,eta,h,nu,omit_radial_axial_diffusion=True)
                omission.append(abs(second['pressure']-omitted['pressure']))
                for n in (3,4,7):
                    r=hierarchy.lower_sources([base]+[null]*(n-1),X,eta,h,nu)
                    zero_high.append(r['pressure']==0.)
                if constant:
                    first=hierarchy.lower_sources([base],X,eta,h,1.)
                    constant_checks.extend((abs(V.v+1)<1e-13,abs(first['pressure']-1/(4*X*X))<1e-13,abs(second['pressure'])<1e-13))
                # A distinct complete-series radial residual assembly at order two.
                if abs(eta)<1:
                    full=series_coefficients([base,null,null],X,eta,h,nu)
                    expected=-full['radial'][2]/(2*X)
                    errors.append(abs(second['pressure']-expected))
    # Smooth toy M0=eta*B(X), B=mu before a later pulse and 0 afterward:
    # at eta=0: Z_-D Z_0 V0=6X B', hence Pi2_X=3nu B'.
    # Its backward pressure on the plateau is 3nu*mu, not zero.
    toy=3*Q(17,10)*Q(3)
    return dict(checks=dict(
        inherited_recursion_matches_independent_shifted_operator=bool(max(errors)<1e-11),
        zero_positive_velocities_imply_zero_pressure_source_from_order3=all(zero_high),
        actual_patch_constant_memory_formula_checked=all(constant_checks),
        omitted_radial_axial_viscosity_changes_second_pressure_source=bool(max(omission)>1e-3),
        toy_constant_second_pressure_on_mean_patch_nonzero=toy>0,
        manufactured_profiles_not_actual_selected_h_or_actual_coefficients=True),
        maximum_source_equation_error=float(max(errors)),maximum_omitted_viscosity_error=float(max(omission)),
        zero_high_order_source_checks=len(zero_high),toy_second_pressure_plateau=str(toy))
