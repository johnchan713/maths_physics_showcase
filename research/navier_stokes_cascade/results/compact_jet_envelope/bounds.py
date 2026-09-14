"""Compact-input bounds: logarithmic comparisons and a positive exponent ledger.

The continuum arguments justifying the inputs are in README. No sampling of
the exact core, materialization of C, or certified frequency is implied.
"""
import importlib.util
from pathlib import Path
import sys
import mpmath as mp

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


c2 = load_module('compact_previous_c2', RESULTS/'angular_c2_transfer'/'bounds.py')
core, previous = c2.core, c2.previous
axial = load_module('compact_axial', RESULTS/'axial_stress_audit'/'bounds.py')
intermediate = load_module('compact_intermediate', RESULTS/'intermediate_decay_audit'/'bounds.py')
p, lower, upper = core.point, core.lower, core.upper


def sum_exp(*exponents):
    """A nonempty sum of at most C terms, C>=10^9: C^max times count <= C^(max+1)."""
    if not exponents or len(exponents) > 10**9:
        raise ValueError('A nonempty bounded term list is required')
    if any(type(v) is not int or v < 0 for v in exponents):
        raise ValueError('Nonnegative integer upper exponents required')
    return max(exponents)+1


def reciprocal_first(norm_exp, inverse_value_exp):
    """||1/f||J1 <= 1/m + ||f||J1/m^2, with positive m>=C^-inverse_value_exp."""
    return sum_exp(inverse_value_exp, 2*inverse_value_exp+norm_exp)


def exponent_ledger():
    field, moment, pressure = 64, 256, 257
    de, dm, dpi = field+1, moment+1, pressure+1
    xi, x, li = 2, 13, 1
    ei = reciprocal_first(field, 4)
    hi = ei+2
    W = sum_exp(0, sum_exp(1+moment, 1+dm)+xi)
    qnum = sum_exp(1+moment, 1+dm, 1+dm, 1+moment)
    Q = sum_exp(W, qnum+xi+hi)
    nnum = sum_exp(1+moment, 1+dm, 1+moment, 1+dm)
    N = sum_exp(W+field, nnum+xi, 1+pressure, 1+dpi)
    p1, p2 = x+Q+li, x+N+li+ei
    a, bs = sum_exp(0, 1+de+ei), 1+de+ei
    # Pad the common bounds BEFORE forming inverse a and the loop coordinates.
    shear, inv_a = 256, reciprocal_first(256, 22)
    ts = shear+inv_a
    vs = shear+sum_exp(0, 2*ts)
    Pc = sum_exp(512, 512+ts)
    states = dict(E_U_J2=field, moments_J2=moment, Pi_J2=pressure,
                  E_inverse_J1=ei, H_inverse_J1=hi, E_inverse_J2=144,
                  H_inverse_J2=146, X_inverse_J2=xi, L_inverse_J2=li,
                  W_J1=W, Qs_J1=Q, Ns_J1=N, p1_J1=p1, p2_J1=p2,
                  a_J1=a, bs_J1=bs, a_inverse_J1=inv_a, ts_J1=ts,
                  vs_J1=vs, Pc_J1=Pc, Jc_J1=Pc,
                  repair_normalization_operator_C2=290,
                  normalized_actual_moments_C2=moment+290,
                  original_required_gap_inverse_C0=1024)
    return dict(base_C_lower=10**9, exponents=states, envelope_exponent=4096,
                checks=dict(positive_denominator_cost_retained=ei == 73 and hi == 75,
                            angular_differentiation_cost_retained=dm == 257 and dpi == 258,
                            original_pressure_within_C512=max(p1, p2) < 512,
                            shear_within_C256=max(a, bs) < shear,
                            all_input_bounds_within_envelope=max(states.values()) < 4096),
                loop_derivatives_bounded=False, actual_frequency_certified=False)


def scalar_bounds(digits=80):
    inherited = core.certificate('Md64', digits)
    b = inherited['log_bounds']
    log, exp = mp.iv.log, mp.iv.exp
    B, lc, td = exp(b['B']), b['C'], b['P']-1
    L0 = log(p('27.5'))+b['Lambda']
    xa = log(4)-b['Lambda']
    rho = b['XR']+241*td-23
    xrep = rho+5
    e_outer = b['P']-log(2)-p('.5')-td/2-p('.51')-p('.51')*240*td
    bump_J2 = 5*p('5.72e-14')*(64+250000+10**8)
    restoration_J2 = p('3e-18')*(1+64+10000)
    linv_J2 = 1/p('.98')+p('.04')/p('.98')**2 \
        +(p('.04')/p('.98')**2+2*p('.04')**2/p('.98')**3)/2
    checks = dict(
        inherited_core_constants=all(inherited['checks'].values()),
        C_above_absorbed_constants=lower(lc) > upper(log(10**9)),
        P_K_Lambda_B_below_C=max(upper(b[k]) for k in ('P','K','Lambda','B')) < lower(lc),
        native_mixed_U_bound=upper(log(50)+b['K']-b['Lambda']) < lower(log(p('.01'))),
        reference_radial_slope_bound=upper(log(100)+b['K']-lc) < 0,
        reference_angular_field_below_one=upper(B-lc) < 0,
        core_axis_pressure_below_one=upper(412*B-2*lc) < 0,
        radial_length_L0_below_C=upper(log(L0)-lc) < 0,
        shape_angular_raw_derivatives_below_C=upper(log(408)+b['B']-lc) < 0,
        shape_first_y_bound=upper(p('.1')+64*(B+1)/(2000*B)) < 1,
        shape_mixed_y_eta_bound=upper(64*(2*B+2)/(2000*B)) < 1,
        shape_second_y_bound=upper(20000*(B+1)/(2000*B)**2) < 1,
        earlier_join_mixed_correction_below_one=upper(bump_J2) < 1,
        earlier_join_factor_above_half=upper(5*p('5.72e-14')*64) < p('.5'),
        restoration_mixed_correction_below_one=upper(restoration_J2) < 1,
        axial_second_radial_derivative_below_24=upper(4*(p(20000)/64**2+1)) < 24,
        outer_field_J2_below_C2=upper(log(120)+b['P']-2*lc) < 0,
        outer_E_above_exp_minus123Td=lower(e_outer+123*td) > 0,
        outer_E_above_inverse_C=lower(lc-123*td) > 0,
        continuation_E_above_inverse_C_squared=lower(lc-B) > 0,
        continuation_a_above_inverse_C22=lower(lc-log(10)) > 0,
        compact_positive_left=lower(xa+lc) > 0,
        repair_scale_above_one=lower(rho) > 0,
        compact_right_below_C12=upper(xrep-12*lc) < 0,
        full_log_length_below_13logC=upper(xrep-xa-13*lc) < 0,
        full_log_length_below_C=upper(log(13*lc)-lc) < 0,
        axis_pressure_C2_below_C4=upper(log(40)+2*b['P']-4*lc) < 0,
        L_inverse_C2_below_two=upper(linv_J2) < 2,
        inverse_lambda_below_C=upper(4*td-lc) < 0,
        boundary_power_excess_above_inverse_C=lower(log(2)-4*td+lc) > 0,
    )
    return dict(interval_digits=digits, checks=checks,
                log_bounds=dict(C=lc, Xa=xa, repair_scale=rho, Xrep=xrep,
                                outer_E_lower=e_outer, full_radial_length=xrep-xa),
                earlier_join_correction_J2_upper=bump_J2,
                restoration_J2_upper=restoration_J2, L_inverse_J2_upper=linv_J2,
                raw_C_or_activation_width_materialized=False)


def margin_bounds(digits=80):
    base = previous.geometry(digits)
    ax = axial.axial_bound(64, 1024, digits)
    finite = intermediate.finite_cone_bound(ax, digits)
    D0 = p(ax['bs_squared_upper'])/2
    pa = p(finite['axial_ps1_lower'])
    A0 = 1-(2+D0)/pa
    # D0 is an UPPER bound, not the actual bs^2/2. Drop the actual positive
    # D/2 term when lower-bounding this factor; adding D0/2 would be invalid.
    factor_one = A0
    factor_two = p(finite['axial_transformed_gap_lower'])
    third = p(ax['Pc_over_ps1_lower'])*pa-(2+D0)
    fourth = 2*pa**2*factor_one*factor_two
    intermediate_fourth = 2*p(12)**2*p(finite['intermediate_quadratic_gap_lower'])
    checks = dict(inherited_geometry_and_inner_margins=all(base['checks'].values()),
                  inherited_finite_cone=finite['all_pass'],
                  axial_Pc_above_two=lower(p(ax['Pc_over_ps1_lower'])*pa) > p('2.1'),
                  axial_third_above_one=lower(third) > 1,
                  axial_first_factor_above_p94=lower(factor_one) > p('.94'),
                  axial_second_factor_above_p417=lower(factor_two) > p('.417'),
                  axial_fourth_above_one=lower(fourth) > 1,
                  intermediate_p1_above_twelve=min(finite['ramp_ps1_lower'],finite['power_ps1_lower']) > 12,
                  intermediate_third_above_one=12-p('2.2') > 1,
                  intermediate_fourth_above_one=lower(intermediate_fourth) > 1)
    return dict(interval_digits=digits, checks=checks,
                axial_third_lower=third, axial_fourth_lower=fourth,
                axial_first_factor_lower=factor_one, axial_second_factor_lower=factor_two,
                intermediate_fourth_lower=intermediate_fourth,
                inner_required_gap_lower='C^(-1024)', Pc_minus_two_lower='0.1',
                left_collar_vs_minus_two_lower='0.29',
                right_collar_vs_minus_two_exact='2*lambda > C^(-1)',
                positive_original_gaps_required_only_where_vs_at_least_two=True)
