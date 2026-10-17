"""Exact annular derivative and endpoint heat majorants at fixed positive parameters."""
from fractions import Fraction as Q
from math import factorial


def order(n):
    if not isinstance(n,int) or isinstance(n,bool) or n<0:raise ValueError('Nonnegative integer order required')


def rising(x,n):
    order(n);out=Q(1)
    for j in range(n):out*=Q(x)+j
    return out


def heat_A(h,n):
    order(n);h=Q(h)
    if not 0<h<1:raise ValueError('Fixed 0<h<1 required')
    return rising(h,n)*rising(1+h,n)


def angular_majorant(h,k):
    order(k)
    if k==0:return 2*heat_A(h,1)
    return sum(Q(factorial(k)*4**(k-2*b)*2**b,factorial(k-2*b)*factorial(b))*heat_A(h,k-b)
        for b in range(k//2+1))


def derivative_bounds(beta,targets,qjets,d0):
    beta,d0=Q(beta),Q(d0)
    if beta<=0 or d0<0 or not targets or any(Q(v)<0 for v in [*targets,*qjets]):raise ValueError('Invalid positive inverse or seminorm data')
    d=[d0]
    for k in range(1,len(targets)):
        lower=Q(0)
        for i in range(min(k+1,len(qjets))):
            for j in range(k-i+1):
                l=k-i-j
                if (i,j,l) in ((0,k,0),(0,0,k)):continue
                lower+=Q(qjets[i])*d[j]*d[l]
        d.append(2*beta*(Q(targets[k])+lower))
    return d


def foundation_review():
    return dict(annular_Q='eta-dependent U-square coefficient exp(6/5)*(1+eta^2)^2/P^2',
        normalized_polynomial_derivative_suprema=[4,8,8,4,1],
        same_root_inverse='J=B+2Q_eta(c,.), ||J^-1||<=2000',
        derivative_recursion='triple i+j+l=k excluding only (0,k,0),(0,0,k); Q derivatives do not vanish until i>4',
        exact_root_registry=['annular beta=1000,q=100','angular beta=5,q=32',
            'axial affine inverse<20/lambda','strength F_A>.35','scalar eta-independent waiting',
            'heat compensation beta=50,q=100'],
        checks=dict(polynomial_all_order_norm_does_not_require_new_Ck_smallness=sum([4,8,8,4,1])==25,
            actual_P_floor_keeps_polynomial_multiplier_below_one=Q(25*4,16**2)<1,
            accepted_annular_C1_root_gives_Neumann_floor=2*1000*100*Q(6,10**13)<Q(1,2),
            annular_Q_eta_derivatives_not_replaced_by_constant_Q_rule=True,
            finite_outer_prefix_uses_eta_independent_flat_radial_endpoints=True,
            smooth_whole_compact_input_supplied_relative_to_exact_core_and_retained_finite_existence_bounds=True,
            same_fixed_N_all_order_map_theorem_applies_without_reselection=True,
            no_scalar_strength_swirl_or_waiting_circularity=True))


def tail_review():
    return dict(heat_derivative='(-1)^j*(h)_j*Gamma_integral/ Gamma(1+h), bound A_j=(h)_j*(1+h)_j',
        angular_chain='k!*sum_b H^(k-b)*(-4eta/X)^(k-2b)*(-2/X)^b/((k-2b)!*b!)',
        field_difference='|Dy^r Deta^k(E-Epow)|<=C_rk,h*c_infinity*X^(-3/2-h)',
        global_identities=['M_inf=0','J_inf=0','S_inf=0','renormalized_R=0','Pi_axis+Cp_inf=0'],
        checks=dict(all_fixed_endpoint_jets_have_Gamma_dominating_integrals=True,
            raw_I_total_not_set_to_zero_or_finite=True,
            fixed_h_positive_kept_before_infinite_integration=True,
            all_fixed_global_moment_eta_derivatives_commute_with_integrals=True,
            global_complex_eta_strip_at_heat_poles_is_obstructed=True,
            analytic_first_collar_not_global_heat_analyticity_is_required=True,
            exact_global_ledger_uses_exact_roots_not_quadrature_coefficients=True))


def support_review():
    intervals=dict(I1=[-25,-20],I2=[-20,-15],I3=[-14,-9],I4=[-8,-3])
    return dict(stage_offsets=intervals,exact_endpoint='Xrep=XR*exp(241T-18)=Xstar',
        reserved_coordinates='I3:0<log(X/rho3)<5; I4:6<log(X/rho3)<11',
        plateau='U0=0, E0=e3(eta)*exp(-(.5+lambda)*y3), M0=mu*eta, V0=-mu, mu>0',
        pressure='K3(eta)-e3(eta)^2*exp(-(1+2lambda)*y3)/(2*(1+2lambda)); incoming K3 not reset',
        checks=dict(first_repair_bumps_end_before_I2=Q(-25)+Q(5,2)+Q(1,20)<-20,
            heat_compensation_bumps_end_before_I3=Q(-20)+Q(9,2)+Q(3,20)<-14,
            heat_compensation_last_bump_precedes_rho3=Q(-18)+Q(9,2)+Q(3,20)<-12,
            I3_and_I4_are_disjoint=intervals['I3'][1]<intervals['I4'][0],
            reserved_intervals_precede_axial_pulse=intervals['I4'][1]<0,
            U_zero_does_not_set_M_or_V_to_zero=True,
            compact_repair_restores_M_and_heat_compensation_changes_only_E=True,
            leading_pressure_keeps_compensation_memory_and_higher_pressure_tails_not_dropped=True,
            exterior_M_V_vanish_only_after_completed_axial_support_and_exact_total_M=True,
            all_order_weighted_T1_over_T0_join_not_inferred_from_smoothness=True))
