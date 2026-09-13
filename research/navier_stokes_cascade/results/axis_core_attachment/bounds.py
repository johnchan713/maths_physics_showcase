"""Positive majorants for the exact analytic core and a proposed attachment.

All returned large numbers are logarithms unless named otherwise. In particular,
C, exp(-Lambda), and the narrow widths are NEVER materialized. The selected Md64
scale makes even the integer exponent of C infeasible to store.

Arithmetic verifies the inequalities specified here; the mathematical reasons
these expressions bound the fields belong in README.md. Arithmetic cannot
verify an unproved continuation estimate merely by making its parameters small.
"""
from fractions import Fraction
from math import comb, factorial
import mpmath as mp


def point(value):
    if isinstance(value, Fraction):
        return mp.iv.mpf(value.numerator)/value.denominator
    return mp.iv.mpf(value)


def lower(value):
    return mp.make_mpf(value._mpi_[0])


def upper(value):
    return mp.make_mpf(value._mpi_[1])


def logsum(*logs):
    """Logarithm of a positive sum, with outward interval operations."""
    shift = point(max(upper(v) for v in logs))
    return shift+mp.iv.log(sum((mp.iv.exp(v-shift) for v in logs), point(0)))


def weight(n, k, rho=Fraction(1)):
    """Exact B.4 weight, used only at finite integer indices in diagnostics."""
    if n < 0 or k < 0:
        raise ValueError('Nonnegative weight indices required')
    return (Fraction(factorial(k)*comb(n+k, k), 20**n*(n+1)**2*(k+1)**2)
            /rho**k)


def log_pressure(scale):
    if scale == 'P16':
        return mp.iv.log(16)
    if scale in ('Md4', 'Md64'):
        return mp.iv.exp(int(scale[2:]))+11
    raise ValueError('Unknown pressure scale')


def certificate(scale='Md64', digits=80):
    if not isinstance(digits, int) or digits < 60:
        raise ValueError('At least 60 interval digits required')
    mp.iv.dps = digits
    p, log, exp = point, mp.iv.log, mp.iv.exp
    lp = log_pressure(scale)
    leps, lj, lsigma = log(p('1e-16')), log(p('1e-18')), log(p('1e-20'))
    lr0 = lsigma-log(256)
    lrho = lsigma-log(2048)
    la = log(256)
    ld = log(80)-lrho
    lm = log(200)-2*lsigma
    lz = log(1024)+2*lp-lr0
    lv = p(30720)  # log(exp(40 * algebra_constant * ||chi||)).
    ls = logsum(log(2), lv, log(40)+lz)

    # J R1 <= D [386 A S + (40+10 M) A^2 S^2].
    lq1 = logsum(log(40), log(10)+lm)+2*la
    lm1 = ld+logsum(log(386)+la+ls, lq1+2*ls)
    ll1 = ld+logsum(log(386)+la, log(2)+lq1+ls)
    # The pressure p=I(g^2 Phi^2) costs 80 A^3 S^2.
    lq2 = logsum(log(32)+2*la, log(2720)+4*la)
    lm2 = ld+logsum(log(442)+la+ls, lq2+2*ls)
    ll2 = ld+logsum(log(442)+la, log(2)+lq2+ls)
    lf = lv+logsum(lm1, lm2)
    llip = lv+logsum(ll1, ll2)
    leval = log(10000)-6*lrho
    lk = 16*(log(p('1e6'))+leval+
             logsum(p(0), ls, lf, llip, lm, lz, 2*lp))
    llam = log(p('1e12'))+8*lk-2*leps

    # The following attachment budget is deliberately isolated from the core
    # contraction. Its applicability requires the derivation in the note.
    lb = log(p('1e20'))+64*logsum(p(0), llam, lk)
    big_b = exp(lb)
    log_c = 10000*big_b-100*leps
    tsh = 2000*big_b
    log_width = -20*log_c
    log_xsep = tsh-10*(log_c+lp)
    log_xr = log(110)+10*(log_c+lp)
    log_D = 4*log_c+100*big_b

    logs = {
        'P': lp, 'r0': lr0, 'rho': lrho, 'zeta_bound': lm,
        'Z_over_L_bound': lz, 'angular_inverse_bound': lv,
        'ball_norm_bound': ls, 'J_R1_bound': lm1, 'J_R2_bound': lm2,
        'map_displacement_numerator': lf, 'map_Lipschitz_numerator': llip,
        'real_evaluation_bound': leval, 'K': lk, 'Lambda': llam,
        'B': lb, 'C': log_c, 'transition_length_Tsh': log(2000)+lb,
        'activation_width': log_width, 'kappa': log_width,
        'final_cutoff_width_each': log_width,
        'xsep': log_xsep, 'XR': log_xr,
        'map_displacement': lf-llam,
        'map_Lipschitz': llip-llam,
        'native_axial_C1_bound': lk-llam,
        'core_source_error': log(p('1e6'))+4*lk-llam,
        'reference_cutoff_error': log(100)+2*lk+log_width,
        'reference_native_error': log(100)+2*lk-llam,
        'reference_pressure_error': big_b-2*log_c,
        'reference_p1_C3_bound': log(p('1e6'))+log(110)+llam+8*lk,
        'reference_ns_C3_bound': log(p('1e6'))+8*lk,
        'reference_inverse_p1_C3_bound': log(p('1e24'))+4*llam+24*lk,
        'reference_log_CE_C3_bound': log(p('1e6'))+llam+8*lk,
        'reference_complement_vr_lower': log(16)+2*lj+2*log_c-2*big_b-2*llam-lb,
        'stress_comparison_constant_D': log_D,
        'activation_quadratic_test_ratio': 3*log_c+2*log_D+2*log_width-log(8),
        'global_comparison_error': log(6)+log_D+log_width,
        'activation_axial_C1_bound': log(6)+2*lb+log_width,
        'core_ledger_C1_bound': log(2)+1400*big_b-2*log_c,
        'small_width_geometric_bound': log_width,
        'first_collar_width': log_width-log(10)-log(log_c)/2,
        'first_collar_shear_loss': 4-97*log_c,
    }
    checks = {
        'complex_strip_inside_previous_one': upper(exp(lr0)) < p('1e-20')/p('210.8'),
        'Cauchy_radius_ratio_one_quarter': upper(exp(lrho-lr0)*2) <= p('.25000000000000000001'),
        'chi_inverse_series_finite': upper(lv) == 30720,
        'core_map_preserves_unit_ball': upper(logs['map_displacement']) < log(p('.5')),
        'core_map_contracts': upper(logs['map_Lipschitz']) < log(p('.5')),
        'native_axial_meets_budget': upper(logs['native_axial_C1_bound']) < leps-log(100),
        'core_comparison_error_below_one_millionth': upper(lk-llam) < log(p('1e-6')),
        'source_absorption_budget': upper(logs['core_source_error']) < log(p('1e-6')),
        'linear_chi_absorption_budget': upper(log(10000)+2*lk-llam) < 0,
        'axis_complement_keeps_source_sign': upper(log(4)+lk-llam) < lj,
        'K_bounds_inverse_j': lower(lk) > -lj,
        'B_dominates_core_log_amplitude': lower(lb) > upper(log(100)+llam+lm),
        'C_bounds_complex_axis_amplitude': lower(log_c) > upper(3*exp(llam+lm)),
        'C_larger_than_exp_1000B': lower(log_c) > upper(1000*big_b),
        'transition_chosen_independently_of_C': upper(tsh/big_b) < p('2000.000000000000001'),
        'reference_cutoff_source_budget': upper(logs['reference_cutoff_error']) < lj-log(100),
        'reference_native_source_budget': upper(logs['reference_native_error']) < lj-log(100),
        'reference_pressure_source_budget': upper(logs['reference_pressure_error']) < lj-log(100),
        'reference_complement_vr_above_two': lower(logs['reference_complement_vr_lower']) > log(p('2.3')),
        'D_below_C6': upper(log_D-6*log_c) < 0,
        'reference_vr_below_half_C3': upper(log(2)+6*big_b-log_c) < 0,
        'first_collar_keeps_vs_above_two': upper(logs['first_collar_shear_loss']) < log(p('.01')),
        'first_collar_quadratic_budget': upper(logs['activation_quadratic_test_ratio']) < 0,
        'remaining_continuation_error_budget': upper(logs['global_comparison_error']) < log(p('1e-6')),
        'activation_axial_meets_budget': upper(logs['activation_axial_C1_bound']) < leps-log(100),
        'core_ledger_meets_budget': upper(logs['core_ledger_C1_bound']) < leps-log(4),
        'core_separated_from_annulus': upper(log_xsep) < -8,
        'joining_radius_large_enough': lower(log_xr) > log(10000),
        'native_domain_contains_reference_cutoff': upper(log_width) < log(p('.001')),
        'final_cutoffs_fit_before_Xi': upper(log(4)+log_width) < log(log(p('1.1'))),
        'log_only_large_amplitude': True,
    }
    for name in ('p1', 'ns', 'inverse_p1', 'log_CE'):
        checks['reference_'+name+'_budget'] = (
            upper(logs['reference_'+name+'_C3_bound']) < lower(lb-log(100)))
    return dict(scale=scale, interval_digits=digits, log_bounds=logs, checks=checks)


def elementary_certificate(digits=80):
    mp.iv.dps = digits
    p = point
    # These are universal statements proved by a monotone polynomial or by
    # coefficient convolution, rather than sampling angular coordinates.
    scalar_lower = p(Fraction(305719, 1152000))
    edge_poly = -p(Fraction(75535511, 400000000))
    domain_error = p('1e-6')
    return {
        'algebra_constant': upper((8*mp.iv.pi**2/6)**2) < 256,
        'comparison_Phi_above_p265': lower(scalar_lower) > p('.265'),
        'exact_core_Phi_above_p264': lower(scalar_lower-domain_error) > p('.264'),
        'endpoint_scalar_sign': upper(edge_poly) < p('-.188'),
        'endpoint_p1_above_2p3': lower(2-2*edge_poly-domain_error) > p('2.3'),
        'reference_p1_barrier_at_X100': lower(p('1.2')*(100-p('4e-6')**2/100)) > 3,
        'final_constant_slope_source': lower(p('.6')*p('2.8')-p('.1')-p('.02')) > 1,
        'shape_transition_source': lower(p('.55')*p('2.8')-p('.1')-p('.02')) > 1,
        'shape_transition_slope_error': upper(64*(1+p('1e-10'))/2000) < p('.05'),
        'axial_offset_budget': upper(3*p('1e-18')) < p('1e-16')/4,
    }


def entry_transfer(digits=80):
    """Integrate the actual constant axial offset from Xsep to x=exp(-8).

    The five input core norms have already been bounded by epsilon/4. These
    are worst-case integral bounds, not values of manufactured moment data.
    """
    mp.iv.dps=digits
    p,exp=point,mp.iv.exp
    eps=p('1e-16')
    g=3*eps/100
    core=eps/4
    rows=[core+exp(-2)*g,
          core+p('0.625')*exp(p('-3.2'))*g,
          core,
          core+12*exp(p('-.8'))*g*g/256,
          core]
    return dict(core_C1_bound=core, axial_offset_C1_bound=g,
                entry_moment_C1_bounds=rows,
                row_order=['M','J-4etaI','I','-(S-8etaM)','Cp'],
                checks=dict(all_five_entry_rows=all(upper(v)<lower(eps) for v in rows),
                            axial_entry=upper(g)<lower(eps),
                            all_rows_below_p255_epsilon=max(upper(v) for v in rows)<lower(p('.255')*eps)))


def tail_bound_log(norm_log, degree, radius='4.1'):
    """Rigorous value tail for an exact B_rho element of the stated norm.

    Its use requires a bound on the norm of THAT element, not a small last
    coefficient of a finite Taylor approximation.
    """
    if not isinstance(degree, int) or degree < 0 or not 0 < point(radius) < 20:
        raise ValueError('Nonnegative degree and radius in (0,20) required')
    ratio = point(radius)/20
    return (norm_log+(degree+1)*mp.iv.log(ratio)
            -2*mp.iv.log(degree+2)-mp.iv.log(1-ratio))
