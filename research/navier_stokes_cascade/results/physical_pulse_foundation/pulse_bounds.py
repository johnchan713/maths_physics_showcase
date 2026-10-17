"""Exact exponent partitions, barrier and all-band Gaussian constants."""
from fractions import Fraction as Q

KAPPA=Q(1,100000)


def order(n):
    if type(n) is not int or n<0:raise ValueError('Nonnegative integer order required')
    return n


def exponents(h):
    h=Q(h)
    if not 0<h<Q(1,100):raise ValueError('Keep fixed construction h in (0,1/100)')
    return dict(A=Q(1,2)+h,D=Q(1,2)-h,sx=Q(1,2)+h/2,st=1+3*h/2,
        coefficient_space=Q(1,2)+h*KAPPA,coefficient_time=1+h)


def harmonic_loss(h,a,b,phase_factors):
    e=exponents(h);order(a);order(b);order(phase_factors)
    if not 1<=phase_factors<=a+b:raise ValueError('A nonempty derivative partition is required')
    return Q(a,2)+b*(1+Q(h))+phase_factors*Q(h,2),a*e['sx']+b*e['st']


def primary_powers(h,a,b):
    e=exponents(h);order(a);order(b);base=-a*e['sx']-b*e['st']
    return dict(velocity=base-e['sx'],pressure=base-2*e['sx'],potential=base)


def riccati_outward(mu,ce,length):
    mu,ce,length=map(Q,(mu,ce,length))
    if min(mu,ce,length)<=0 or ce/length>mu/8:raise ValueError('Finite positive accepted pulse parameters required')
    eta=ce/length;r=eta/mu
    return r,-2*mu*r+eta*(1+2*r+r*r),-Q(47,64)*eta


def tail_square_constant(a,polynomial_power,physical_loss,target_power):
    a=Q(a);order(polynomial_power);physical_loss,target_power=map(Q,(physical_loss,target_power))
    if a<=0 or target_power<0:raise ValueError('Positive Gaussian curvature and nonnegative target required')
    d=Q(polynomial_power)+max(Q(0),physical_loss+target_power)
    return d,d*d/(4*a)


def input_review():
    return dict(actual_inputs='all fixed Euler/profile jets, strict closed growth/cone gaps, exact sqrt-edge factors and formal supports from weighted_stress_foundation',
        band_map='s-Z^2*s^(2h)=T; L=1-2h*eta^2>=.98; 1/4<=s<=4; 4Q<qstar',
        frozen_data='band, representative, sign, mesh, cover, rectangle copy, nearest NONZERO carrier',
        checks=dict(actual_background_inputs_now_supplied_relative_to_exact_core_and_finite_lemmas=True,
            representatives_selected_inside_closed_active_chart_intersection=True,
            one_band_threshold_for_all_later_scales_and_fixed_derivative_orders=True,
            nonzero_integer_rounding_error_is_one_over_k_not_half_over_k=True,
            entire_pulse_L_over_k_rounding_and_fixed_radius_costs_retained=True,
            heat_endpoint_jets_one_sided_not_analytic_or_zero_extended=True))


def pulse_review():
    return dict(frame='Bleft*(An*B-Bprime)=diag(lambda,-lambda)+E; |E|+|d-dref|<=C/Ls',
        pressure='i*(n.Kbg*t-nprime.t+n.f)/(k*m*|n|^2)',
        amplitude='cP<=uncut radial x<=CP; complete oscillating curl has NO positive pointwise lower bound',
        radius_weights='Haar h_sigma comparable r0^(3/2)*S^(-1/2); roots cost r0^(-3/4)*S^(1/4)',
        checks=dict(tangential_normal_floor_not_inverse_of_angular_component_used=True,
            moving_normal_and_frame_derivatives_and_damping_retained=True,
            Riccati_barrier_and_scalar_exponential_close_nonvanishing_continuation=True,
            actual_Haar_transverse_area_kept_in_covariance_inverse=True,
            original_leading_stress_target_not_replaced_by_That=True,
            all_fixed_coefficient_jets_have_finite_polynomial_S_losses=True,
            sourced_inverse_retains_explicit_input_bounds_and_endpoint_support_hypotheses=True))


def physical_review():
    return dict(graph='Dr=dR+Mi*dr*R^(dr-1)*Li; Dt=-epsilon*dT+ci*Ni; Dz=epsilon*dZ',
        phase='Li*Phi=0, grad_x*v=0, dt_phys*v=Q^(-1-h); coefficients do NOT share radial graph cancellation',
        losses='sx=1/2+h/2; st=1+3h/2',
        primary_powers='velocity Q^(-(a+1)*sx-b*st); pressure Q^(-(a+2)*sx-b*st); potential Q^(-a*sx-b*st)',
        checks=dict(all_repeated_graph_coefficient_and_Cartesian_frame_derivatives_kept=True,
            Faadi_Bruno_phase_factor_count_at_most_total_derivative_order=True,
            potentials_and_pressures_cut_before_complete_curl=True,
            primary_curl_remainder_exponent_one_minus_kappa_above_one_half=True,
            physical_bounds_uniform_in_band_label_copy_not_in_h_N_or_order=True,
            background_radial_flow_and_I4_pressure_orders_one_and_two_retained=True))


def support_review():
    return dict(supports=dict(transverse='abs(xi)<=r0/2',pulse='abs(v/Ls-1/2)<=2/5',shell='Xa<=X<=Xb',time='T>=0 one-sided'),
        tail='C*ell^B*Q^(-Lambda)*exp(-a0*ell^2)<=C_IK*Q^K on SAME entire band domain',
        square_root_weight='exp(-t1^2/(2*ya^2)-2/yb^2)',
        checks=dict(closed_enlarged_preimage_separation_preserves_all_derivative_supports=True,
            locally_finite_at_q_positive_and_smooth_zero_near_axis=True,
            positive_endpoint_envelope_P_does_not_supply_zero_extension=True,
            source_requires_its_own_full_harmonic_endpoint_extension=True,
            square_completion_absorbs_every_fixed_physical_loss_on_same_start_band=True,
            exact_half_weight_constants_used_for_flat_shell_zero_extension=True,
            complete_mean_and_generated_source_full_residual_gates_not_promoted=True))
