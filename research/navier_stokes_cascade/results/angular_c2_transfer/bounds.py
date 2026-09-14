"""Outward C2 transfer constants for the fixed, implicitly defined profile.

The analytic reasons these quantities bound the actual fields are in README.
No enormous amplitude, tiny width, Taylor surrogate, or sampled moment target
is substituted for the exact profile. Large quantities are kept as logarithms.
"""
from fractions import Fraction
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


core = load_module('c2_core_bounds', RESULTS/'axis_core_attachment'/'bounds.py')
matching = load_module('c2_matching_bounds', RESULTS/'axis_matching_audit'/'bounds.py')
previous = load_module('c2_previous_stress', RESULTS/'joined_stress_construction'/'construction.py')
p, lower, upper, logsum = core.point, core.lower, core.upper, core.logsum


def core_transfer(digits=80):
    old = core.certificate('Md64', digits)
    b = old['log_bounds']
    log, exp = mp.iv.log, mp.iv.exp
    B, lc = exp(b['B']), b['C']
    eps = p('1e-16')

    # The all-index generating-function bound gives ||Phi||2, ||u||2 <= D2.
    ld2 = log(3)+b['ball_norm_bound']-2*b['rho']
    # Cauchy on a disk of radius r0/2 gives ||log(phi_star)||2 <= Lambda M(2+1/r0).
    lstar = b['Lambda']+b['zeta_bound']+logsum(log(2), -b['r0'])
    # In the positive core, ||log(Phi)||2 <= 20 D2^2.
    llogphi = log(20)+2*ld2
    L0 = log(p('27.5'))+b['Lambda']
    # Integrate the reference p1/X bound with dX, not its sup with d(log X).
    lcontinuation = logsum(lstar, llogphi,
                           log(p('55e6'))+b['Lambda']+8*b['K'], log(2*L0+10))
    # log of the bound on ||CE/sqrt(X)||2 in the core (not CE itself).
    laxis_amplitude = log(2)/2+exp(lstar)+ld2
    native = ld2-b['Lambda']
    axial_change = log(p('1e8'))+8*b['K']+b['activation_width']
    log_length = log(L0+2000*B)
    # Cp <= C^-2 exp(408B) [L/2+2] <= C^-2 exp(412B).
    cp_length_factor = logsum(log_length-log(2), log(2))
    core_rows = log(p('2e6'))+412*B-2*lc
    checks = dict(
        inherited_parameter_checks=all(old['checks'].values()),
        all_index_evaluation_constant=Fraction(496, 243) < 3,
        rho_below_one=upper(b['rho']) < 0,
        Phi_lower_allows_reciprocal=Fraction(264, 1000) > Fraction(1, 4),
        Lambda_above_four=lower(b['Lambda']) > upper(log(4)),
        core_angular_offset_below_epsilon_over_100=upper(native) < lower(log(eps/100)),
        continuation_angular_offset_below_epsilon_over_100=upper(axial_change) < lower(log(eps/100)),
        activation_width_below_one=upper(b['activation_width']) < 0,
        core_CE_over_sqrtX_below_expB=upper(laxis_amplitude-B) < 0,
        continuation_log_CE_C2_below_B=upper(lcontinuation-b['B']) < 0,
        shape_log_CE_C2_below_204B=upper(201*B+3-204*B) < 0,
        E_C2_below_one=upper(204*B-lc) < 0,
        pressure_length_factor_below_exp4B=upper(cp_length_factor-4*B) < 0,
        separation_before_entry=upper(b['xsep']) < -8,
        core_rows_C2_below_epsilon_over_4=upper(core_rows) < lower(log(eps/4)),
    )
    return dict(interval_digits=digits, checks=checks,
                log_bounds=dict(weighted_C2_evaluation_D2=ld2,
                    log_phi_star_C2=lstar, log_Phi_C2=llogphi,
                    continuation_log_CE_C2=lcontinuation,
                    axis_CE_over_sqrtX_C2=laxis_amplitude,
                    native_U_minus_Ustar_C2=native, continuation_U_change_C2=axial_change,
                    shape_CE_C2=204*B, Cp_C2=412*B-2*lc,
                    normalized_core_moments_C2=core_rows),
                fixed_actual_scale='Md64', raw_C_or_width_materialized=False,
                uses_exact_core_existence_argument=True)


def normalization_constants(digits=80):
    mp.iv.dps = digits
    exp = mp.iv.exp
    # Each row <= A_i*xsep^(1/5) + B_i*C^-2*exp(412B).
    coefficients = [
        [18*exp(6), p(0)],
        [p(60)/16*exp(p('9.6')), p(0)],
        [(p(10)/48+p(5)/8)*exp(p('9.6')), p(0)],
        [(p(6500)/256+p(5)/12)*exp(p('7.2')), p(0)],
        [p('2.5')*exp(p('1.2')), p(20)/256*exp(p('1.2'))],
    ]
    return dict(interval_digits=digits, row_coefficients=coefficients,
                inverse_f_C2=5, inverse_f_squared_C2=20,
                checks=dict(all_coefficients_below_one_million=all(
                    upper(v) < 1000000 for row in coefficients for v in row),
                    old_C1_factor_12_not_a_C2_bound=20 > 12))


def matching_transfer(digits=80):
    mp.iv.dps = digits
    exp = mp.iv.exp
    eps, gb = p('1e-16'), p('3e-18')
    core_budget = eps/4  # Certified by core_transfer, not inferred from C1 data.
    rows = [core_budget+exp(-2)*gb,
            core_budget+p(5)/8*exp(p('-3.2'))*gb,
            core_budget,
            core_budget+p(20)/256*exp(p('-.8'))*gb*gb,
            core_budget]
    zeta = p(20)/256*exp(p('1.2'))
    # beta^2 <= 64 beta, exact bump mass .08, all supports disjoint.
    qu = 64*p('.08')*(exp(p('.34'))+exp(p('.74')))
    qe = 32*p('.08')*sum((exp(p('1.2')*p(v)) for v in ('.19', '.54', '.89')), p(0))
    qc = 32*p('.08')*sum((exp(p('.2')*p(v)) for v in ('.19', '.54', '.89')), p(0))
    quadratic = zeta*qu+qe
    target = p('.286')*eps
    radius = 2000*target
    Z2 = 2*target  # Norm has the factorial denominator 2!.
    raw_curvature = previous.curvature_bound(Z2, digits)
    u_inverse = matching.inverse_bounds(['1', '1.6'], '.3', '.4')
    e_inverse = matching.inverse_bounds(['1.6', '1.2', '.2'], '.15', '.35')
    checks = dict(
        five_entry_rows_below_p255_epsilon=all(upper(v) < lower(p('.255')*eps) for v in rows),
        axial_offset_below_entry_epsilon=upper(gb) < lower(eps),
        target_after_axial_restoration=upper(p('.255')*eps+gb+gb*gb) < lower(target),
        zeta_C2_below_one=upper(zeta) < 1,
        quadratic_C2_below_100=max(upper(quadratic), upper(qc)) < 100,
        inverse_below_1000=max(u_inverse+e_inverse) < 1000,
        exact_C2_contraction=upper(8*1000**2*100*target) < 1,
        new_C2_ball_inside_old_C1_ball=upper(radius) < lower(p('6e-13')),
        raw_curvature_below_5p721e_minus14=upper(raw_curvature) < lower(p('5.721e-14')),
        implicit_curvature_improves_C2_norm_bound=upper(raw_curvature) < lower(2*radius),
    )
    return dict(interval_digits=digits, checks=checks,
                entry_row_C2_bounds=rows, g_C2_bound=gb, zeta_C2_bound=zeta,
                target_C2_bound=target, target_raw_second_derivative_bound=Z2,
                U_inverse_upper=u_inverse, E_inverse_upper=e_inverse,
                quadratic_C2_upper=quadratic, pressure_quadratic_C2_upper=qc,
                contraction_upper=8*1000**2*100*target,
                coefficient_C2_bound=radius, coefficient_raw_second_from_norm=2*radius,
                coefficient_raw_second_from_implicit_equation=raw_curvature,
                same_exact_root_as_previous_C1_construction=True)
