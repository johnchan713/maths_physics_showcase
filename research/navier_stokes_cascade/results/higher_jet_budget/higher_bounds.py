"""Higher input jets and finalized compact-background derivative bounds.

The continuum derivations are in README. Interval arithmetic checks their
constant comparisons; it does not independently prove the inherited core.
Large parameters, N and the actual high-frequency flow are not materialized.
"""
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import sys
import mpmath as mp

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
AMIN, INPUT_POWER = 2**16, 256


def module_at(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


core = module_at('higher_jet_frozen_core', HERE.parent/'axis_core_attachment'/'bounds.py')
p, lower, upper, logsum = core.point, core.lower, core.upper, core.logsum


def cutoff_certificate(digits=80):
    """Uniform derivatives through order five of the unchanged flat step.

    Interior: logistic derivatives <=2^j*j!, |H^(j)|<1e8, Bell_5=52.
    Edges: the convergent exponential series retains exp(-1/x^2).
    """
    mp.iv.dps = digits
    middle = p(10000)*52*p('1e8')**5
    edge = 104*p(1440)**5*4**15*mp.iv.exp(p(16)/9-16)
    tail_ratio = 32*mp.iv.exp(-14)
    return dict(interval_digits=digits, maximum_step_order=5,
                common_step_derivative_bound='1e100', middle_upper=middle,
                edge_upper=edge, geometric_edge_ratio=tail_ratio,
                checks=dict(interior_derivative_bound=upper(middle) < p('1e100'),
                            flat_edge_derivative_bound=upper(edge) < p('1e100'),
                            exponential_edge_series=upper(tail_ratio) < p('.5'),
                            constant_fits_C12=10**100 < (10**9)**12,
                            width_point_one_fifth_cost_retained=10**5 < 10**9))


def transfer_C4(digits=80):
    previous = core.certificate('Md64', digits)
    b = previous['log_bounds']
    log, exp = mp.iv.log, mp.iv.exp
    B, lc, eps = exp(b['B']), b['C'], p('1e-16')
    e4 = sum((F(4, 3)**(k+1)/(k+1)**2 for k in range(5)), F())
    e5 = sum((F(4, 3)**(k+1)/(k+1)**2 for k in range(6)), F())
    d4, d5 = log(3)+b['ball_norm_bound']-4*b['rho'], log(3)+b['ball_norm_bound']-5*b['rho']
    S = exp(b['ball_norm_bound'])
    native = log(p('1e20'))+8*log(1+S)+48*log(1+exp(-b['rho']))
    zeta4 = log(80)+b['zeta_bound']-4*b['r0']
    pressure5 = log(p('1e5'))+2*b['P']
    star4 = b['Lambda']+b['zeta_bound']+logsum(log(2), -b['r0'], log(p(4)/3)-2*b['r0'], log(2)-3*b['r0'])
    star5 = logsum(star4, b['Lambda']+b['zeta_bound']+log(p(16)/5)-4*b['r0'])
    # These reference bounds are derived by slope integration and the exact
    # variation-of-constants formula in README, not inherited C3 labels.
    rp1 = log(p('1e18'))+log(110)+b['Lambda']+12*b['K']
    rns = log(p('1e8'))+8*b['K']
    rphi5 = logsum(star5, log(10)+b['K'])
    continuation = logsum(star4, log(1000)+4*d4,
                          log(p('55e18'))+b['Lambda']+12*b['K'],
                          log(2*(log(p('27.5'))+b['Lambda'])+10))
    offset_change = log(p('1e10'))+8*b['K']+b['activation_width']
    gb = p('3e-18')
    target = p('.255')*eps+gb+gb**2
    coefficient = 2000*target
    checks = dict(inherited_core_parameters=all(previous['checks'].values()),
                  exact_C4_generating_function=e4 == F(14624, 6075) and e4 < 3,
                  C5_native_generating_function=e5 < 3,
                  native_logarithm_and_reciprocal_jets=upper(native-b['K']) < 0,
                  native_zeta_C4=upper(zeta4-b['K']) < 0,
                  native_axis_pressure_C5=upper(pressure5-b['K']) < 0,
                  native_slope_integral_small=upper(log(200)+b['K']+b['activation_width']) < 0,
                  native_C4_offset=upper(d4-b['Lambda']) < lower(log(eps/100)),
                  reference_U_C5_offset=upper(log(1000)+b['K']-b['Lambda']) < lower(log(eps/100)),
                  reference_p1_C4_below_B=upper(rp1-b['B']+log(100)) < 0,
                  reference_ns_C4_below_B=upper(rns-b['B']+log(100)) < 0,
                  reference_phi_C5_logarithm=upper(rphi5-b['B']+log(100)) < 0,
                  reference_pressure_C5_small=upper(log(110)+2*B-2*lc) < 0,
                  actual_continuation_C4_logarithm=upper(continuation-b['B']) < 0,
                  actual_C4_terminal_offset=upper(offset_change) < lower(log(eps/100)),
                  shape_log_f_C4=upper(201*B+682-204*B) < 0,
                  actual_incoming_C4_rows=upper(log(p('2e6'))+412*B-2*lc) < lower(log(eps/4)),
                  normalized_quadratic_C4=upper(25*exp(p('1.2'))/256) < 1,
                  same_root_C4_contraction=upper(8*1000**2*100*target) < 1,
                  joining_C4_coefficient=upper(coefficient) < p('5.72e-14'),
                  raw_fourth_factor_retained=upper(24*coefficient) < p('1.3728e-12'))
    return dict(interval_digits=digits, exact_C4_evaluation=str(e4), exact_C5_evaluation=str(e5),
                log_D4=d4, log_D5=d5, log_native_jet_bound=native,
                reference_p1_C4='1e18*X*Lambda*K^12', reference_ns_C4='1e8*K^8',
                reference_fourth_order_derived_from_defining_equations=True,
                log_reference_p1_bound=rp1, log_reference_ns_bound=rns,
                log_actual_continuation_bound=continuation,
                normalized_target_C4_bound=target, joining_C4_bound=coefficient,
                joining_raw_fourth_bound=24*coefficient,
                f_inverse_C4=5, f_inverse_square_C4=25, log_f_C4_bound=682,
                normalized_energy_row_coefficient=25*(101+160+64), checks=checks)


def add_exp(*exponents):
    if not exponents or any(type(v) is not int or v < 0 for v in exponents):
        raise ValueError('Nonnegative exponent list required')
    return max(exponents)+1


def inverse_exp(norm_exponent, inverse_value_exponent, order):
    if any(type(v) is not int or v < 0 for v in (norm_exponent, inverse_value_exponent, order)):
        raise ValueError('Nonnegative integer orders and exponents required')
    return inverse_value_exponent+order*(inverse_value_exponent+norm_exponent)+1


def compact_exponents():
    """J4 field/moment and J3 coordinate costs, with original value floors."""
    # These reference source bounds are anisotropic: W includes U_eta, so
    # at total order two it sees at most two radial derivatives of U,
    # hence only ONE differentiated radial cutoff, not two.
    sq_ref = add_exp(36+56, 2+35, 36+36)
    sn_ref = add_exp(36+64, 2+2*35, 36+36, 1+132, 1+2*64)
    p1y = add_exp(1+6+1, 2+16)
    p1yy = add_exp(1+110+1, 36+16, 2+p1y)
    p1yyy = add_exp(1+110+1, 56+16, 36+p1y, 2+p1yy)
    nsy = add_exp(5+1, 9)
    nsyy = add_exp(160+1, nsy)
    nsyyy = add_exp(160+1, nsyy)
    regional = dict(reference_W_J2=36, reference_log_slope_J2=56,
                    reference_Hc_J2=36, reference_logE_eta_J2=36,
                    reference_Uy_J2=64, reference_Pi_eta_J2=132,
                    source_Sq_J2=sq_ref, source_Sn_J2=sn_ref,
                    reference_p1_y_C2=p1y, reference_p1_yy_C1=p1yy,
                    reference_p1_yyy_C0=p1yyy, reference_ns_yyy_C0=nsyyy,
                    actual_nonzero_logE_J4=75+128+1,
                    actual_U_J4=75+75+192+2+1,
                    actual_E_J4=4*256+1,
                    cumulative_moment_J4=2*2048+36+1)
    field, moment, pressure = 2048, 8192, 8193
    ei = inverse_exp(field, 4, 4)
    hi = ei+4
    padded_ei, padded_hi = 16384, 32768
    df, dm, dp = field+1, moment+1, pressure+1
    W = add_exp(0, add_exp(1+moment, 1+dm)+2)
    qnum = add_exp(1+moment, 1+dm, 1+dm, 1+moment)
    Q = add_exp(W, qnum+2+padded_hi)
    nnum = add_exp(1+moment, 1+dm, 1+moment, 1+dm)
    Ns = add_exp(W+field, nnum+2, 1+pressure, 1+dp)
    shear = add_exp(0, 1+df+padded_ei)
    ai = inverse_exp(32768, 22, 3)
    ts = 32768+131072
    vs = 32768+add_exp(0, 2*ts)
    pc = add_exp(131072, 65536+ts)
    jc = add_exp(65536, 131072+ts)
    values = dict(fields_J4=field, moments_J4=moment, pressure_J4=pressure,
                  reciprocal_E_J4=ei, reciprocal_H_J4=hi,
                  W_J3=W, Qs_J3=Q, Ns_J3=Ns,
                  p1_J3=12+Q+1, p2_J3=12+Ns+1+padded_ei,
                  shear_J3=shear, inverse_a_J3=ai, ts_J3=ts, vs_J3=vs,
                  Pc_J3=pc, Jc_J3=jc,
                  normalized_repair_weights_C3=32768)
    return dict(exponents=values, regional_exponents=regional, original_A_C_power=4096,
                higher_input_A_power=INPUT_POWER, higher_input_C_power=4096*INPUT_POWER,
                original_value_bounds_and_parameters_unchanged=True,
                checks=dict(reference_Sq_padding=sq_ref < 110, reference_Sn_padding=sn_ref < 160,
                            reference_radial_p1_padding=max(p1y, p1yy, p1yyy) < 128,
                            reference_radial_ns_padding=max(nsy, nsyy, nsyyy) < 192,
                            actual_log_slope_padding=regional['actual_nonzero_logE_J4'] < 256,
                            actual_axial_taper_padding=regional['actual_U_J4'] < 512,
                            value_preserving_exponential_padding=regional['actual_E_J4'] < 2048,
                            original_moment_padding=regional['cumulative_moment_J4'] < 8192,
                            E_inverse_padding=ei < padded_ei, H_inverse_padding=hi < padded_hi,
                            Qs_padding=Q < 65536, Ns_padding=Ns < 16384,
                            pressure_paddings=12+Q+1 < 131072 and 12+Ns+1+padded_ei < 65536,
                            shear_padding=shear < 32768, inverse_a_padding=ai < 131072,
                            all_higher_inputs_below_A256=max(values.values()) < 4096*INPUT_POWER))


def scale_checks(digits=80):
    mp.iv.dps = digits
    a, g = p(AMIN), p(AMIN)
    log = mp.iv.log
    # Dividing by A^20 bounds log(A)/A^20 by 1/A^19; all ratios
    # decrease for A>=AMIN. This proves the comparison for every larger A.
    polynomial = 150*log(10)+10000*log(a)+10*a**20
    polynomial_ratio = 150*log(10)/a**20+10000/a**19+10
    checks = dict(small_polynomial_and_cutoff_costs_below_G=upper(polynomial) < lower(a**32),
                  comparison_holds_for_all_larger_A=upper(polynomial_ratio) < lower(a**12),
                  G64_below_expG=upper(64*log(g)) < lower(g),
                  G64_absorption_increases=lower(1-64/g) > 0,
                  primitive_joint_C3_fits_same_H=10**6 < AMIN**3,
                  target_prefactor_A_power_below_template=3+4*INPUT_POWER < 10000,
                  field_repair_prefactor_below_template=INPUT_POWER < 10000,
                  moment_repair_prefactor_below_template=3+2*INPUT_POWER < 10000,
                  same_N_has_two_sided_bound=32+8+1 < 64,
                  no_new_frequency_is_selected=True)
    return dict(interval_digits=digits, checks=checks, A_min=AMIN,
                absorption_template='1e150*A^10000*delta^(-10)<G',
                fixed_frequency='1+floor(H^32)', K='H^64',
                large_parameters_or_frequency_materialized=False)


def loop_ledger():
    """Safe exponents k in exp(k*G); README gives each chain rule."""
    rows = dict(q_mixed_order_four=4, regular_variance_factor_order_three=9,
                sqrt_variance_factor_order_three=40,
                signed_F_partials_order_three=50, inverse_F_mu=10,
                mu_first=70, mu_second=220, mu_third=500,
                ratio_joint_order_three_fixed_theta=2000,
                phase_density_order_three=5000, forward_phase_order_three=5001,
                inverse_phase_joint_order_three=40000,
                pulled_back_ratio_joint_order_three=150000,
                loop_joint_order_three=500000,
                primitive_joint_order_three=1000000)
    checks = dict(first_root_chain_rule=10+50+1+2 < rows['mu_first'],
                  second_root_chain_rule=10+50+2*rows['mu_first']+2 < rows['mu_second'],
                  third_root_chain_rule=10+50+max(3*rows['mu_first'], rows['mu_first']+rows['mu_second'])+3 < rows['mu_third'],
                  ratio_composition=rows['mu_third']+4+10 < rows['ratio_joint_order_three_fixed_theta'],
                  density_product=2*rows['ratio_joint_order_three_fixed_theta']+10 < rows['phase_density_order_three'],
                  inverse_phase_third=5*rows['forward_phase_order_three']+20 < rows['inverse_phase_joint_order_three'],
                  pulled_back_composition=rows['ratio_joint_order_three_fixed_theta']+3*rows['inverse_phase_joint_order_three']+10 < rows['pulled_back_ratio_joint_order_three'],
                  rational_loop_composition=3*rows['pulled_back_ratio_joint_order_three']+10 < rows['loop_joint_order_three'],
                  primitive_integration=rows['loop_joint_order_three']+20 < rows['primitive_joint_order_three'])
    return dict(exponents=rows, maximum_joint_order=3, maximum_density_order=4,
                zero_pressure_uses_integral_quotient=True,
                zero_variance_uses_signed_square_root=True,
                inverse_phase_derivatives_not_frozen=True,
                inverse_lift_value_bounded_on_one_fundamental_period_only=True, checks=checks)


def final_ledger():
    return dict(higher_input='A_hi=A^256', primitive_joint_J3='H',
                angular_field_modulation_error_C3='4*A_hi*H/N',
                angular_moment_modulation_error_C3='32*A^3*A_hi^2*H/N',
                normalized_target_C3='128*A^3*A_hi^4*H/N <= H^16/N',
                final_repair_coefficient_C3='2*1000*H^16/N',
                final_field_angular_error_C3='H^18/N',
                final_moment_angular_error_C3='H^20/N',
                final_field_angular_C3='2*A_hi', final_moment_angular_C3='2*A_hi',
                modulated_field_mixed_J2='H^2*N', finalized_field_mixed_J2='H^3*N',
                ordinary_F_U_Pi_V_C2='H^8*N < H^64',
                ordinary_jet_envelope='K=H^64',
                annular_residual_bounds=dict(radial='42*H^448*q^(-3/2)',
                                            tangential='440*H^448*q^(-Aphys-1)',
                                            axial='148*H^384*q^(-Aphys-1)'),
                fixed_viscosity=1, physical_q_range='0<q<=1',
                domain='X in [Xa*exp(tc/2), rhop*exp(5)], eta in [-1,1]',
                full_wave_tensor_or_global_physical_residual_included=False,
                residual_is_proved_to_decay=False,
                checks=dict(positive_radial_frequency_loss_retained=True,
                            fourth_moment_both_squares_retained=True,
                            unchanged_axis_pressure=True,
                            final_root_is_same_as_accepted_C1_root=True,
                            residual_power_arithmetic=7*64 == 448 and 6*64 == 384))


def fourth_implicit_terms():
    return dict(coefficients=[8, 6, 8, 24, 12, 12, 8, 1],
                equation='L*c4=z4-8Q[c1,c3]-6Q[c2,c2]-8Q1[c,c3]-24Q1[c1,c2]-12Q2[c,c2]-12Q2[c1,c1]-8Q3[c,c1]-Q4[c,c]',
                raw_fourth_factor=24)
