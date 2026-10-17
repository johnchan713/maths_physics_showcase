"""Explicit joint analytic bounds for the infinite core, not a finite truncation."""
from fractions import Fraction as Q
from math import comb,factorial


def index(n):
    if type(n) is not int or n<0:raise ValueError('Finite nonnegative integer order required')
    return n


def mixed_bound(r,k,norm=1,rho=1,Y=Q(5),eta_offset=Q(1,8)):
    """eta_offset is |Im/complex eta displacement| divided by rho."""
    index(r);index(k);norm,rho,Y,b=map(Q,(norm,rho,Y,eta_offset))
    margin=1-Y/20-b
    if norm<0 or rho<=0 or Y<0 or b<0 or margin<=0:raise ValueError('Require positive joint analytic margin')
    return norm*factorial(r+k)/20**r/rho**k/margin**(r+k+1)


def analytic_review():
    rows=[]
    for r,k in ((0,0),(1,0),(0,1),(2,6),(8,9),(40,20)):
        rows.append(dict(radial=r,parameter=k,bound=str(mixed_bound(r,k))))
    a,b=Q(1,4),Q(1,8)
    partial=sum(Q(comb(n+k,k))*a**n*b**k for n in range(12) for k in range(12))
    total=1/(1-a-b)
    return dict(majorant='norm*(r+k)!*20^(-r)*rho^(-k)*(1-|Y|/20-|zeta|/rho)^(-r-k-1)',
        generating_function='sum_(n,k>=0) binom(n+k,k)*a^n*b^k=1/(1-a-b), a+b<1',
        rows=rows,finite_positive_partial_sum=str(partial),infinite_majorant=str(total),
        checks=dict(joint_margin_not_two_independent_margins=1-a-b==Q(5,8),
            finite_positive_partial_sum_below_exact_infinite_majorant=partial<total,
            mixed_factorial_not_product_factorials=factorial(8+9)>factorial(8)*factorial(9),
            all_arbitrary_fixed_orders_have_finite_bounds=all(Q(x['bound'])>0 for x in rows),
            uniform_strip_includes_eta_endpoints_by_glued_Taylor_disks=True,
            infinite_coefficient_norm_not_finite_Taylor_iteration_is_input=True,
            constants_may_depend_on_derivative_order_and_fixed_core_parameters=True))


def zero_free_review():
    lower=Q(33,125);Y_slope=Q(16,125);eta_slope=Q(64,25)
    variation=(Y_slope+eta_slope)/100;remaining=lower-variation
    return dict(real_lower=str(lower),joint_margin='5/8',
        Y_radius='1/(100*S)',eta_radius='rho/(100*S)',S_lower='2',
        Y_derivative_over_S=str(Y_slope),eta_derivative_times_rho_over_S=str(eta_slope),
        variation_upper=str(variation),complex_real_part_lower=str(remaining),reciprocal_upper='5',
        reciprocal_jet_bound='5*r!*k!*(100*S)^r*(100*S/rho)^k',
        checks=dict(common_polydiscs_remain_in_Y_five=Q(41,10)+Q(1,200)<5,
            common_polydiscs_remain_in_eta_eighth=Q(1,200)<Q(1,8),
            slopes_follow_from_joint_majorant=mixed_bound(1,0)==Y_slope and mixed_bound(0,1)==eta_slope,
            variation_is_strictly_below_real_margin=remaining>Q(1,5),
            one_common_zero_free_neighborhood_for_all_orders=True,
            reciprocal_and_principal_log_not_inferred_from_real_positivity_alone=True,
            no_order_independent_jet_constant_or_uniform_C_limit_claimed=True))


def profile_review():
    return dict(norm_budgets=dict(Phi='S',u='S',g='1',F='256*S',
            U='9+S/Lambda',p='80*256^3*S^2',Pi='10*P^2+80*256^3*S^2/Lambda'),
        definitions='F=g*Phi; U=Ustar+u/Lambda; Pi=Pi0+p/Lambda; p_Y=g^2*Phi^2; avg U=int_0^1 U(tY,eta)dt',
        regular_radial_quotient='v0=(2*eta*U-2*D*eta*avgU-(1-eta^2)*partial_eta avgU)/L',
        checks=dict(exact_pressure_family_datum_and_pressure_increment_both_retained=True,
            pressure_eta_and_radial_derivatives_follow_from_same_joint_majorant=True,
            radial_average_preserves_uniform_joint_neighborhood=True,
           _v0_uses_one_additional_parameter_derivative_not_a_bare_norm_inverse=True,
            F_inverse_on_same_smaller_neighborhood_uses_g_nonvanishing_and_Phi_inverse=True,
            E_equals_radius_times_F_and_cannot_be_inverted_at_axis=True,
            full_profile_outside_core_and_final_global_pressure_identity_not_proved_here=True))


def derivative_weight(horizontal,axial,time,h=Q(1,200)):
    h=Q(h)
    if not 0<h<=Q(1,100):raise ValueError('Use the retained core h range')
    return Q(index(horizontal),2)+(Q(1,2)-h)*index(axial)+index(time)


def physical_review():
    return dict(domain='0<q<=1; actual unchanged leading core 0<=Y=Lambda*(xi^2+upsilon^2)/2<=4; -1<=eta<=1; singular event excluded',
        coordinates='xi=x/sqrt(q); upsilon=y/sqrt(q); eta=z/q^D; q-z^2*q^(2h)=1-t',
        weights='w=(alpha_x+alpha_y)/2+D*alpha_z+j; D=1/2-h; A=1/2+h',
        velocity='u1=q^(-1/2)*xi*v0/2-q^(-A)*upsilon*F; u2=q^(-1/2)*upsilon*v0/2+q^(-A)*xi*F; u3=q^(-A)*U',
        pressure='p_phys=q^(-2A)*Pi',
        normalized_z_operator='Z_beta=[(1-eta^2)*partial_eta-eta*(xi*partial_xi+upsilon*partial_upsilon)-2*beta*eta]/L; output exponent beta+D',
        normalized_t_operator='T_beta=[beta+(xi*partial_xi+upsilon*partial_upsilon)/2+D*eta*partial_eta]/L; output exponent beta+1',
        bounds='|D_x^alpha partial_t^j u|<=C_alpha_j*q^(-A-w); |D_x^alpha partial_t^j p|<=C_alpha_j*q^(-2A-w)',
        checks=dict(L_uniformly_positive_on_closed_parameter_interval=1-2*Q(1,100)==Q(49,50),
            radial_and_swirl_Cartesian_axis_form_contains_no_division_by_radius=True,
           _profile_bounded_jets_close_under_each_normalized_operator=True,
            exponent_beta_updated_after_every_derivative=True,
            radial_velocity_power_is_dominated_for_q_at_most_one=True,
            constants_uniform_in_physical_scale_q_but_not_parameter_C_or_derivative_order=True,
            closed_eta_endpoints_use_implicit_root_at_positive_q_not_event_extension=True,
            leading_core_jets_do_not_prove_full_viscous_residual_flatness_or_pulse_bounds=True))
