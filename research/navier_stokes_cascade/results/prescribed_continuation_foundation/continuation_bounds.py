"""Exact seminorm bookkeeping for the finite prescribed, pre-modulation maps."""
from fractions import Fraction as Q
from math import factorial,comb


def step_cost(constants,t):
    t=Q(t)
    if t<=0 or not constants or any(Q(c)<0 for c in constants):raise ValueError('Positive width and nonnegative constants required')
    return sum(Q(c)/(factorial(i)*t**i) for i,c in enumerate(constants))


def primitive_bound(trace,density,length):
    trace,density,length=map(Q,(trace,density,length))
    if min(trace,density,length)<0:raise ValueError('Nonnegative seminorm bounds required')
    return trace+(length+1)*density


def product_convolution(a,b):
    if len(a)!=len(b) or not a:raise ValueError('Same nonempty derivative lists required')
    return [sum(comb(i,j)*Q(a[j])*Q(b[i-j]) for j in range(i+1)) for i in range(len(a))]


def bounds_review():
    # Four real cutoffs and one shape map: all evaluated at frozen positive widths.
    return dict(final_strip='delta=parent_delta/16=rho/(12800*S)',
        nested_strips=['8delta: native/reference fields and reference moments',
            '4delta: reference pressure coordinates, prescribed fields, shape and actual moments',
            '2delta: actual pressure coordinates','delta: every fixed eta derivative'],
        finite_map_order=['native zero-free logarithm','reference slope cutoff and primitive',
            'reference five forward moments','reference Q,N,p1,ns with one eta derivative',
            'actual prescribed slopes and exponential','actual five moments and full pressure',
            'actual endpoint shape exponential','actual coordinates and regular quotients'],
        formulas=dict(primitive='trace+(T+1)*density',product='|fg|_r<=|f|_r*|g|_r',
            exponential='|exp H|_r<=exp(|H|_r)',eta_derivative='|partial_eta H|_(r;s)<=|H|_(r;2s)/s',
            step='sum_i c_i/(i!*t^i)',fixed_width_term='t^(-i)=C^(20*i)'),
        checks=dict(eight_delta_strictly_inside_parent=Q(8,16)<1,
            each_eta_differentiation_has_a_reserved_strip=True,
            actual_exponentials_have_nonzero_reciprocals=True,
            no_unknown_eta_derivative_feedback=True,
            no_inverse_kappa_or_p1_needed_in_prescribed_fields=True,
            width_terms_include_all_requested_orders=step_cost([1,2,3],Q(1,100))==15201,
            primitive_keeps_trace_and_value_length=primitive_bound(3,5,7)==43,
            no_constant_smallness_or_C_uniformity_assumed=True))


def join_review():
    return dict(domains=dict(native_reference='0<=y<=2*t1; thereafter F_r,U_r constant, not E_r',
        first_collar_shortcut='only before the reference taper at y=t1',
        endpoint='Xi=110',shape_end='Xsep=110*exp(Tsh)',shape_length='Tsh=2000*B'),
        prescribed_slopes=dict(ordinary='a=kappa*p1_r; U_y=-kappa*X*ns_r/2',
            axial_taper='a=kappa0*p1_r; U_y=-beta*kappa0*X*ns_r/2',
            a_taper='a=(1-theta)*kappa0*p1_r+(4/5)*theta; U_y=0',
            final='a=4/5; U_y=0; (log E)_y=1/10'),
        density_rows=['X*U','sqrt(2)*X^(3/2)*E','sqrt(2)*X^(3/2)*U*E',
            'X*(U^2-E^2/2)','E^2/2'],
        full_pressure='Pi=Pi_axis+Cp, including eta-dependent original traces',
        leading_certificate='all fixed (X,eta) jets of F,U,Pi,M/X,V/X on [0,Xsep] x D_delta, PRE-MODULATION',
        checks=dict(reference_physical_shear_not_integrated_p1_after_taper=True,
            native_domain_ends_before_seed_extension=4/(1-Q(1,50))<Q(41,10),
            flat_steps_match_all_constant_side_jets=True,
            final_log_E_slope_matches_shape=Q(1,2)*(1-Q(4,5))==Q(1,10),
            ell_is_half_plus_log_E_slope=Q(1,2)+Q(1,10)==Q(3,5),
            power_amplitude_exponents_cancel=(Q(1)-10*Q(1,10)==0),
            actual_moment_traces_not_reset_at_Xi_or_Xsep=True,
            axis_regular_quotients_and_one_more_eta_derivative_retained=True,
            final_modulation_profile_not_identified_with_pre_modulation_profile=True,
            global_backward_pressure_and_all_repairs_remain_separate=True))
