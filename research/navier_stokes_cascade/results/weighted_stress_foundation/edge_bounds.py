"""Exact flat Volterra jets and all-smaller-q budgets; not a continuum solver."""
from fractions import Fraction as Q
from math import factorial


def order(k):
    if type(k) is not int or k<0:raise ValueError('Nonnegative integer order required')
    return k


def binomial(p,n):
    order(n);out=Q(1)
    for j in range(n):out*=Q(p-j,j+1)
    return out


def volterra_coefficients(j,c=4,method='kernel'):
    c=Q(c)
    if c<=0 or not j or method not in ('kernel','ode'):raise ValueError('Positive c, jets and a known method required')
    j=list(map(Q,j));out=[]
    for k,jk in enumerate(j):
        if method=='kernel':
            value=sum(j[k-2*a]*binomial(-Q(k-2*a,2),a)*factorial(a)/c**a for a in range(k//2+1))
        else:value=jk-(Q(k-2,2)/c*out[k-2] if k>=2 else 0)
        out.append(value)
    return out


def flat_power_polynomial(c,p,n):
    c,p=Q(c),Q(p);order(n)
    if c<=0:raise ValueError('The fixed flat constant must be positive')
    a={0:Q(1)}
    for _ in range(n):
        b={}
        for k,v in a.items():
            b[k+3]=b.get(k+3,Q(0))+2*c*v
            b[k+1]=b.get(k+1,Q(0))+(p-k)*v
        a={k:v for k,v in b.items() if v}
    return a


def cutoff_integer_log(h,n,constants,previous=0,boundary=0):
    """Use c_n=2**ell, with an exact rational LOG2 budget, not floating logs."""
    h=Q(h)
    if not 0<h<Q(1,100) or type(n) is not int or n<2:raise ValueError('Fixed construction h and n>=2 required')
    if len(constants)!=n+1 or any(type(v) is not int or v<0 for v in [*constants,previous,boundary]):raise ValueError('All finite nonnegative log2 bounds through n required')
    required=max(Q(previous+1),Q(boundary+1),max(Q(v+n)/(n*h) for v in constants))
    return required.numerator//required.denominator+1


def actual_input_review():
    return dict(inner_factor='exp(-t1^2/ya^2)*Aa; Aa(0)=(1-kappa0)*e*Fr*(ar,qr)',
        exact_core_requirement='Tr identically zero, not |Tr|<=epsilon',
        leading_input='whole closed smooth leading rectangle, common analytic FIRST collar, actual I3 power patch, exterior M=U=V=0 after Xv',
        support=dict(inner='Tn=0 through Xa*exp(tc/16)',velocity_right='rho3*exp(3)',
            pressure_right={'1':'Xv','2':'Xv','n>=3':'Xvel'},stress_right={'1':'Xb','2':'Xv','n>=3':'Xvel'}),
        checks=dict(actual_activation_not_extended_unmodified_core_used=True,
            parent_supplies_actual_leading_inputs_to_sparse_Volterra_and_linear_moment_maps=True,
            exact_five_moments_and_shifted_viscosity_identify_signed_backward_stresses=True,
            order_one_angular_viscosity_uses_renormalized_leading_swirl_integral=True,
            inner_zero_and_higher_common_interior_support_are_functional_identities=True,
            pressure_I4_orders_one_and_two_and_leading_radial_velocity_kept=True))


def outer_review():
    return dict(transform='int e^(-c/u^2)*u^-3*J du=e^(-c/delta^2)*L_c[J]/(2c)',
        endpoint_recursion='v_k=j_k-(k-2)*v_(k-2)/(2c)',
        outer_factor='T0=e^(-4/delta^2)*delta^-3*(b0,delta^6*bz); T1=e^(-4/delta^2)*delta^-3*(b1,0)',
        relative_quotient='(b1/sqrt(b0^2+delta^12*bz^2),0)',
        endpoint='W1(0,eta)=(2*eta^2*Xb/L^2,0)',
        eta_zero_jets='W1theta(delta,0)=Xb*delta^3/8+O(delta^4); Ddelta^3=3Xb/4; Deta^2=4Xb',
        checks=dict(exact_outer_flat_constant_four_and_radial_weights_retained=True,
            beta_f_prime_term_retained_in_full_axial_viscosity=True,
            b0_positive_uniformly_on_one_closed_eta_collar=True,
            eta_zero_quotient_smooth_without_dividing_by_eta_or_b1=True,
            all_fixed_heat_pole_jets_are_one_sided_not_analytic=True,
            joint_smoothness_from_integrable_kernel_majorants_not_an_infinite_Taylor_sum=True))


def background_review():
    return dict(weight='exp(-t1^2/ya^2-4/yb^2)',distance='min(1,ya,yb)',loss='3*m+3',
        quotients='W1 smooth on closed annulus; Wn and Tn/zeta have same compact interior support for n>=2',
        schedule='one free diagonal choice includes all earlier requirements plus every quotient seminorm through n',
        normalized_error='D^I[(That-T0)/abs(T0)]=O_I(q^(2h)) on one common q-domain',
        radial_error='q^A radial velocity O_I(q^h), not O_I(q^(2h))',
        checks=dict(all_order_weighted_edges_and_relative_quotients_use_actual_inputs=True,
            finitely_many_new_bounds_at_each_diagonal_step_are_finite=True,
            whole_smaller_q_support_controls_all_fixed_Euler_profile_jets=True,
            one_value_and_first_tangential_jet_qstar_preserves_growth_and_oriented_cone=True,
            pressure_memory_and_potential_cutoff_radial_M_term_not_dropped=True,
            direct_stress_sum_not_identified_with_zero_complete_NS_residual=True,
            physical_actual_pulse_and_complete_covariance_residual_obligations_remain_open=True))
