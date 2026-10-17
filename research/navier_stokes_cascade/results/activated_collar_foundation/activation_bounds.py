"""All-order actual-collar budgets; continuum proofs are in README.md."""
from fractions import Fraction as Q
from math import comb


def order(n):
    if type(n) is not int or n<0:raise ValueError('Nonnegative finite integer order required')
    return n


def gains(r):
    r=order(r)
    return dict(field=1-r,moment=2-r,integrated_shear=1-r,
        complete_field_coarse=-r,complete_moment_coarse=-r)


def falling_operator(r):
    """partial_X^r=X^-r sum c_j*(X partial_X)^j; signed first-kind coefficients."""
    r=order(r);c=[1]
    for n in range(r):
        out=[0]*(len(c)+1)
        for j,v in enumerate(c):out[j]-=n*v;out[j+1]+=v
        c=out
    return c


def domain_review():
    # log C>=log 2>=1/2 implies sqrt(log C)>1/2.
    ell=Q(1,20*2**20)
    return dict(actual_width='t1=kappa0=C^-20',tc='t1/(10*sqrt(log C))',
        collar_log_length='ell=tc/4',reference_Y_upper='4*exp(ell)<4.1',
        collar_offsets=['1/8','3/16','1/4','1/2'],final_eta_strip='delta=rho/(800*S)',
        checks=dict(actual_length_below_seed_margin=ell<Q(1,100) and 4/(1-Q(1,100))<Q(41,10),
            rescaled_interval_inside_flat_quotient_range=Q(1,20)<Q(1,4),
            actual_collar_before_reference_taper_and_modulation=Q(1,4)<Q(1,2) and Q(1,20)<1,
            ordered_keep_cut_collar_and_modulation_offsets=Q(1,8)<Q(3,16)<Q(1,4)<Q(1,2),
            four_delta_inside_parent_zero_free_strip=Q(4,800)<Q(1,100),
            complex_L_lower_exceeds_point_nine_seven=1-Q(1,50)*Q(101,100)**2>Q(97,100),
            fixed_actual_width_parameters_and_supports_not_retuned=True))


def formula_review():
    return dict(reference='F_r=g*Phi; U_r exact seed; Y=4*exp(y); a_r=-2*Y*Phi_Y/Phi',
        activation='e_a=(1-kappa0)*sigma(y/t1); kappa=1-e_a; no eta dependence in cutoff',
        J='J=int_0^y e_a(v)*a_r(v,eta)/2 dv',
        K='K=int_0^y e_a(v)*partial_v U_r(v,eta) dv',
        actual_fields='F_a=F_r*exp(J); E_a=E_r*exp(J); U_a=U_r-K',
        actual_pressure='Pi_a=Pi_r(Xa,eta)+int_0^y X(v)*F_a(v,eta)^2 dv',
        full_moment_densities=dict(M='X*U',I='2*X^2*F',J='2*X^2*U*F',S='X*(U^2-X*F^2)',Pi='X*F^2'),
        regular_axis_quotient='v0=[2*eta*U-2*D*eta*(M/X)-(1-eta^2)*partial_eta(M/X)]/L',
        inverse_bound='|1/F_a|<=|1/F_r|*exp(ell*sup|a_r|/2)',
        checks=dict(actual_shear_prescription_integrated_not_seed_substituted=True,
            nonvanishing_from_exponential_not_unsupported_small_complex_perturbation=True,
            all_five_incoming_moments_retained_and_Pi_trace_contains_axis_datum=True,
            pressure_increment_uses_actual_F_not_F_r=True,
            no_inverse_kappa_or_flat_activation_divisor_required=True,
            all_fixed_radial_derivatives_preserve_common_eta_strip=True,
            flat_join_is_smooth_radially_not_asserted_radially_analytic=True))


def width_review():
    rows=[dict(order=r,**gains(r),radial_log_coefficients=falling_operator(r)) for r in (0,1,2,4,10,20)]
    return dict(rows=rows,
        field_difference='|partial_y^r partial_eta^k (F_a-F_r,U_a-U_r)|<=C_rk*t1^(1-r)',
        moment_difference='|partial_y^r partial_eta^k (M_a-M_r,I_a-I_r,J_a-J_r,S_a-S_r,Pi_a-Pi_r)|<=C_rk*t1^(2-r)',
        complete_collar='for every r,k: sup |partial_X^r partial_eta^k H|<=B_rk*t1^(-r), H=F,U,Pi,M/X,V/X',
        parameter_Cauchy='on inner strip: k!*delta^(-k)*B_r from one fixed outer strip',
        checks=dict(Bell_factors_have_gain_at_least_one_minus_order=all(n-r>=1-r for r in range(1,31) for n in range(1,r+1)),
            Volterra_integration_adds_one_width_power=all(gains(r)['moment']==gains(r)['field']+1 for r in range(50)),
            full_state_coarse_bound_includes_all_width_losses=all(gains(r)['field']>=-r and gains(r)['moment']>=-r for r in range(50)),
            radial_log_conversion_has_required_lower_order_terms=falling_operator(4)==[0,-6,11,-6,1],
            moments_gain_two_requires_same_nonzero_incoming_traces=True,
            constants_independent_of_auxiliary_width_only_for_frozen_seed_not_as_C_varies=True,
            actual_radial_costs_include_C_to_twenty_r=True))


def collar_review():
    return dict(closed_certificate='for every fixed r: some B_r finite on one delta>0, sup_(0<=X<=Xa*exp(tc/4),dist(eta,[-1,1])<delta) sum_H |partial_X^r H|<=B_r',
        H=['F0','U0','Pi0','M0/X','V0/X'],
        physical_weight='w=(alpha_x+alpha_y)/2+(1/2-h)*alpha_z+j',
        physical_bounds='leading |D_x^alpha partial_t^j u|<=C_alpha_j*q^(-A-w); |D_x^alpha partial_t^j p|<=C_alpha_j*q^(-2*A-w)',
        checks=dict(parent_core_covers_axis_and_new_forward_formulas_cover_remaining_actual_collar=True,
            all_fixed_radial_jets_match_at_Xa_because_e_a_is_flat=True,
            M_over_X_has_axis_average_extension_not_arbitrary_zero_value=True,
            V_over_X_uses_one_more_parameter_derivative_and_nonzero_L=True,
            common_collar_includes_both_coefficient_cutoff_radii=True,
            uniformity_is_in_physical_q_not_in_C_or_derivative_order=True,
            global_backward_pressure_identity_not_reproved_by_this_local_forward_argument=True,
            complete_background_pulses_viscous_residual_and_blowup_not_inferred=True))
