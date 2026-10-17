"""Exact all-order power/implicit-derivative bookkeeping, not numerical PDE bounds."""
from fractions import Fraction as Q


def frequency_powers(r):
    if not isinstance(r,int) or isinstance(r,bool) or r<0:raise ValueError('Nonnegative integer radial order required')
    return dict(field=r-1,moment=max(r-2,-1),angular_state=-1)


def coefficient_bounds(beta,q,targets):
    beta,q=Q(beta),Q(q)
    if beta<=0 or q<0 or not targets or any(Q(v)<0 for v in targets):raise ValueError('Positive inverse and nonnegative seminorms required')
    d=[2*beta*Q(targets[0])]
    for k in range(1,len(targets)):
        d.append(2*beta*(Q(targets[k])+q*sum(d[j]*d[k-j] for j in range(1,k))))
    return d


def map_review():
    return dict(root='F=d0*mu*sqrt(R2(mu*p)); R2(0)=1/2; signed root',
        coefficient_recursion='mu_alpha=-[h^alpha](F(mu_<n,p)-s)/F_mu(base)',
        phase_recursion='same homogeneous inversion with denominator w=phi_theta>0',
        repair='Jc=B+2Q(c,.), inverse<=2beta; SAME accepted C1 root',
        repair_recurrence='d_k<=2beta*(z_k+q*sum_{j=1}^{k-1}d_j*d_(k-j))',
        checks=dict(zero_variance_and_zero_pressure_have_nonzero_signed_root_derivative=True,
            unsigned_sqrt_or_division_by_V_mu_at_zero_not_used=True,
            moving_inverse_phase_and_all_denominator_derivatives_retained=True,
            inverse_lift_value_bounded_only_on_fundamental_period=True,
            no_simultaneous_contraction_in_all_Ck_norms=True,
            eta_independent_B_and_Q_use_same_original_e_and_fixed_lambda=True,
            pointwise_Jacobian_not_high_derivative_smallness_controls_inverse=True,
            same_low_order_accepted_frequency_and_coefficients_preserved=True,
            common_complex_eta_analyticity_not_asserted_on_modulation=True))


def transfer_review():
    return dict(orders={str(r):frequency_powers(r) for r in (0,1,2,3,8,100)},
        full_pressure='Pi_N=Pi_axis+Cp_N; same original eta-dependent datum, actual fields in density',
        exterior_pressure_defect='c(eta); Ns defect=4*Aphys*eta*c-(1-eta^2)*c_eta',
        input='C-infinity PRE-MODULATION profile on WHOLE compact rectangle; positive floors and accepted finite-order bounds',
        scope='all-order MAP theorem relative to smooth input; parent covers inner/shape only',
        checks=dict(radial_field_differences_not_assumed_small=True,
            moments_gain_one_radial_order_over_fields=True,
            phase_N_log_X_independent_of_eta=True,
            pressure_coordinate_angular_order_k_uses_input_order_k_plus_one=True,
            incoming_moments_normalization_derivatives_and_inverse_lambda_retained=True,
            five_exact_rows_match_as_functions_with_all_fixed_eta_derivatives=True,
            C0_cone_acceptance_is_same_not_new_all_order_smallness=True,
            whole_annular_exterior_smooth_input_not_supplied_by_inner_parent=True))
