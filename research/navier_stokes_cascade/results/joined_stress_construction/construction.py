"""Actual support geometry, derivative transfer, and a REJECTED promotion.

Large quantities are represented by logarithms or symbolic expressions. The
candidate frequency is not certified by the arithmetic in this module.
"""
from fractions import Fraction
import importlib.util
from pathlib import Path
import sys

import mpmath as mp

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent


def module_at(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


core = module_at('joined_stress_core_bounds', RESULTS/'axis_core_attachment'/'bounds.py')
p, lower, upper = core.point, core.lower, core.upper


def geometry(digits=80):
    old = core.certificate('Md64', digits)
    log, exp = mp.iv.log, mp.iv.exp
    b = old['log_bounds']
    lc, lp = b['C'], b['P']
    td = lp-1
    tw = 240*td
    log_tc = b['first_collar_width']
    log_xa = log(4)-b['Lambda']
    log_rhop = b['XR']+td+tw-23
    # Do not form log(Xa) + tc/2: tc itself cannot be materialized at Md64,
    # and finite-precision addition would erase this protected separation.
    checks = dict(
        inherited_parameter_checks=all(old['checks'].values()),
        left_offsets_strictly_ordered=Fraction(0) < Fraction(1, 4)
            < Fraction(1, 2) < Fraction(1),
        compact_starts_above_inverse_C=lower(log_xa+lc) > 0,
        modulation_stops_before_repair=lower(tw-26) > 0,
        reserved_repair_fits=Fraction(0) < Fraction(45, 100)
            < Fraction(255, 100) < Fraction(5),
        repair_ends_before_heat_patch=Fraction(255, 100) < Fraction(5),
        compact_ends_below_C12=upper(log_rhop+5-12*lc) < 0,
        lambda_inverse_below_C=upper(4*td-lc) < 0,
        h_inverse_below_C=upper(8*td-lc) < 0,
        C_at_least_2pow16=lower(lc) > upper(16*log(2)),
        activation_step_lower_C401=lower(lc-log(4)) > 0,
        activation_ratio_below_Cminus28=upper(
            b['activation_quadratic_test_ratio']+28*lc) < 0,
        activation_quadratic_gap_above_Cminus1024=lower(222*lc) > 0,
        first_collar_vs_above_2p29=upper(b['first_collar_shear_loss']) < log(p('.001')),
        reference_error_below_one_millionth=upper(b['global_comparison_error']) < log(p('1e-6')),
        strengthened_shape_barrier=lower(100-p('.65')*p('2.5')) > 0,
        first_outer_ramp_floor=lower(exp(p('-.6'))*(p('1.1')
            -p('.01')/p('.98')*(exp(1)-1))) > p('.5'),
    )
    return dict(
        interval_digits=digits,
        locations={
            'protected_core_end': 'Xa * exp(tc/4)',
            'modulation_left': 'Xa * exp(tc/2)',
            'known_inner_collar_end': 'Xa * exp(tc)',
            'modulation_right': 'XR * exp(Td+3)',
            'repair_left_rhop': 'XR * exp(Td+Tw-23)',
            'repair_support': 'rhop * exp((.45,2.55))',
            'compact_right': 'rhop * exp(5)',
            'Xa': '4/Lambda', 'tc': 'C^(-20)/(10*sqrt(log(C)))',
            'Tw': '240*Td',
        },
        log_scales=dict(C=lc, Xa=log_xa, tc=log_tc,
                        left_offset=log_tc-log(2),
                        protected_offset=log_tc-log(4)),
        checks=checks,
        endpoint_offsets_materialized=False,
        field_derivative_envelope_verified=False,
    )


def derivative_constants(digits=80):
    mp.iv.dps = digits
    exp = mp.iv.exp
    quarter, three_quarters = p('.25'), p('.75')
    g = 2/quarter**3+2/three_quarters**3
    gp = 6/quarter**4+6/three_quarters**4
    gpp = 24/quarter**5+24/three_quarters**5
    middle = (g**3+3*g*gp+gpp)/4
    edge = 32*4**9*exp(p(16)/9-16)
    # zeta = exp(6/5) * (1+eta^2)^2 / P^2, P >= 16.
    zeta = exp(p('1.2'))*p(20)/256  # factorial-weighted C2 norm
    beta, q, radius = p(1000), p(100), p('6e-13')
    jac_inv = beta/(1-2*beta*q*radius)
    curvature_offset = jac_inv*7*q*radius**2
    return dict(
        interval_digits=digits,
        step_third_middle_upper=middle, step_third_edge_upper=edge,
        joining_zeta_C2_upper=zeta,
        small_root_Jacobian_inverse_upper=jac_inv,
        curvature_bound_at_zero_target_curvature=curvature_offset,
        checks=dict(step_third_below_1e6=max(upper(middle), upper(edge)) < 1000000,
                    joining_zeta_C2_below_one=upper(zeta) < 1,
                    Jacobian_inverse_below_1001=upper(jac_inv) < 1001,
                    curvature_offset_below_3e_minus19=upper(curvature_offset) < p('3e-19')),
        actual_target_curvature_supplied=False,
    )


def curvature_bound(target_second, digits=80):
    """Conditional supremum bound. It does not infer Z2 from a C1 bound."""
    mp.iv.dps = digits
    z = p(target_second)
    if not mp.isfinite(lower(z)) or not mp.isfinite(upper(z)) or lower(z) < 0:
        raise ValueError('A finite nonnegative actual target curvature bound is required')
    r = p('6e-13')
    return p(1000)/(1-2*1000*100*r)*(z+700*r*r)


def c1_counterexamples():
    """Analytic controls, not data from the joined profile.

    g_k=eps*sin(k*eta)/(2k), k>=2. The displayed C1 quantity is a uniform
    upper bound; the C2 derivative peak is attained at eta=pi/(2k).
    """
    eps = Fraction(1, 10**16)
    return [dict(k=k, value_sup=eps/(2*k), derivative_sup=eps/2,
                 C1_upper=eps*(k+1)/(2*k), second_derivative_sup=eps*k/2)
            for k in (100, 10**8, 10**16, 10**24)]


def finite_integer_comparison():
    """Exact inequalities used ONLY in the proposed cap/delta comparison.

    For A>=2^16, Z<2^28 A^10, log(64 A^(3/2))<=3 A,
    T<=exp(2 A^16), J<=exp(3 A^16). The tests below dominate coefficients
    at the left endpoint; every monomial ratio decreases as A increases.
    They do not establish a bound on the actual fields or pressure errors.
    """
    a = 2**16
    return dict(
        cap_Z_coefficient=16*2049**2 < 2**28,
        cap_log_polynomial=(2**28+3) < a**6,
        delta_quadratic_budget=6+8 < a**4,
        delta_below_gamma_half=2 < a**19,
        G_dominates_delta_inverse=2 < a**12,
        G_dominates_shear_excursion=8 < a**16,
    )


CANDIDATE = {
    'envelope_proposal': 'A=C^4096 (UNPROVED for the required actual jets)',
    'd0': '1/(16*A)',
    'Z': '16*(1+2048*A^5)^2',
    'mu_cap': '64*A^(3/2)*exp(Z)',
    'delta_loop': 'exp(-A^20)',
    'G': 'exp(A^32)',
    'H': 'exp(G^4)',
    'epsilon_proposal': 'G^(-64)',
    'Cstate_proposal': 'H^16',
    'D_proposal_including_inverse_lambda': 'H^16',
    'Ccorr_proposal': 'A^128',
    'coefficient_tolerance_proposal': 'A^(-128)',
    'N_proposal': '1+floor(H^32)',
    'certified_for_actual_joined_profile': False,
}

OPEN_OBLIGATIONS = {
    'actual_target_C2_transfer': 'Bound the actual five incoming moment functions and their second angular derivatives from the exact core and shape transition, including normalization derivatives.',
    'actual_compact_jet_envelope': 'Derive a uniform mixed-derivative bound and positive denominator bounds for the actual fields, cumulative moments, and ps throughout the selected compact interval.',
    'actual_modulation_error_constants': 'Propagate fixed-phase loop derivatives through the modulated fields and the original pressure formula to certify Cstate and the transformed discrepancy D, retaining 1/lambda.',
    'actual_repair_error_constants': 'Bound the correction-to-state map by Ccorr on a positive-field neighborhood and justify its coefficient tolerance.',
}


def promotion_review():
    """No caller-supplied boolean can turn an unproved estimate into evidence."""
    return dict(accepted=False, missing=list(OPEN_OBLIGATIONS),
                reason='The symbolic frequency is a proposal; its actual input and error bounds have not been proved.',
                actual_joined_profile_frequency_selected=False,
                full_admissible_stress_realized=False,
                full_PDE_corrections_verified=False,
                smooth_force_verified=False, blowup_verified=False)
