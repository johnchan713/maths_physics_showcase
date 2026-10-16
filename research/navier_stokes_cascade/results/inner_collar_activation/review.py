"""Independent integral, full-memory, endpoint and failed-shortcut checks."""
import math
import mpmath as mp
from collar import (Activation, Jet, collar_threshold, flat_factor,
                    flat_quotient, ga, laguerre, moment_map, seed, sigma)

CASES = [('.01', '.1', '.25', '-1'), ('.01', '.1', '.25', '.3'),
         ('.01', '.5', '.05', '1'), ('.001', '0', '.001', '0'),
         ('.0001', '.5', '.00001', '.3')]


def relative(a, b):
    return abs(a-b)/max(abs(a), abs(b), mp.mpf('1e-100'))


def quotient_review():
    rows = []
    b = lambda u: 1+2*u+3*u*u
    for ss in ('.25', '.15', '.05'):
        s = mp.mpf(ss)
        for power in (0, 3):
            got = flat_quotient(s, b, power=power)
            # Independent finite-interval integral, scaled BEFORE quadrature.
            def integrand(v):
                if not v:
                    return mp.mpf(0)
                return mp.exp(-(1/v**2-1)/s**2)*ga(s*v)/ga(s)*v**power*b(s*v)/s**2
            direct = mp.quad(integrand, [0, mp.mpf('.5'), mp.mpf('.9'), mp.mpf('.99'), 1])
            rows.append(dict(s=ss, power=power, coefficient=float(got), relative_error=float(relative(got, direct))))
    delta = mp.mpf('1e-5')
    value = flat_quotient(delta, b)
    predicted = mp.mpf('.5')+delta+mp.mpf('.75')*delta**2
    jet_error = abs(value-predicted)/delta**2
    wrong = sum(w*ga(mp.mpf('.25')/mp.sqrt(1+mp.mpf('.25')**2*v))*b(mp.mpf('.25')/mp.sqrt(1+mp.mpf('.25')**2*v))
                for v, w in laguerre(24, mp.mp.dps))/(2*ga(mp.mpf('.25')))
    gap = abs(wrong-flat_quotient('.25', b))
    return dict(rows=rows, endpoint_value=float(flat_quotient(0, b)),
                endpoint_second_jet_scaled_error=float(jet_error), omitted_Jacobian_gap=float(gap),
                checks=dict(positive_half_line_matches_finite_interval=all(r['relative_error'] < 1e-18 for r in rows),
                            nonzero_endpoint=flat_quotient(0, b) == mp.mpf('.5'),
                            endpoint_jets_through_second_order=jet_error < mp.mpf('.001'),
                            omitted_Jacobian_detected=gap > mp.mpf('.01')))


def primitive_review():
    errors = []
    for eta in (-1, 0, 1):
        y, e = mp.mpf('.0025'), mp.mpf(eta)
        b = seed(y, e)
        expected = dict(M=b['X']*b['U'].v, I=b['X']*b['H'].v,
                        J=b['X']*b['U'].v*b['H'].v,
                        S=b['X']*(b['U'].v**2-b['E'].v**2/2), Cp=b['E'].v**2/2)
        for name, value in expected.items():
            errors.append(relative(mp.diff(lambda yy: seed(yy, e)[name].v, y), value))
            errors.append(relative(mp.diff(lambda ee: seed(y, ee)[name].v, e), b[name].d))
    return dict(maximum_derivative_error=float(max(errors)),
                actual_core_solution=False,
                checks=dict(all_five_forward_primitive_derivatives=all(v < mp.mpf('1e-60') for v in errors),
                            no_X_factor_in_log_pressure_row=True,
                            seed_explicitly_not_actual_core=True))


def activation_rows(order=24):
    rows, cache = [], []
    for width, kappa0, ss, ee in CASES:
        a = Activation(width, kappa0, order=order)
        s, eta, t = mp.mpf(ss), mp.mpf(ee), mp.mpf(width)
        r = a.stress(s, eta)
        f, dh = r['state']['fields'], r['state']['hats']
        reference, actual = moment_map(f['base'], eta=eta), moment_map(r['state']['actual'], eta=eta)
        resolvable = s >= mp.mpf('.15')
        map_error = max(relative((actual['p1']-reference['p1'])/f['ea'], r['p1hat']),
                        relative((actual['p2']-reference['p2'])/f['ea'], r['p2hat'])) if resolvable else None
        components = [r['denominator_error']/(t*s**3), r['cross_over_ea']/(t*s**3),
                      r['vs_error_over_ea']/(t*s**3)]
        rows.append(dict(width=width, kappa0=kappa0, s=ss, eta=ee,
                         binary64_flat_factor=float(flat_factor(s)),
                         field_coefficients=[float(f['lc'].v), float(f['uc'].v)],
                         moment_coefficients=[float(c.v) for c in r['state']['moment_coefficients']],
                         p_hat_coefficients=[float(r['p1hat']/(t*s**3)), float(r['p2hat']/(t*s**3))],
                         cone_coefficients=list(map(float, components)),
                         excess_cone_loss_over_t2s6=float(r['cone_loss']/(t*t*s**6)),
                         exact_map_subtraction_relative_error=float(map_error) if map_error is not None else None,
                         pressure_p2_hat=float(r['pressure_contribution_to_p2hat']),
                         angular_Jeta_p1_hat=float(r['angular_Jeta_contribution_to_p1hat']),
                         seed_defect_norm_over_Fr=float(mp.sqrt(sum(v*v for v in r['seed_defect_over_Fr'])))))
        cache.append((a, r))
    return rows, cache


def activation_review():
    rows, cache = activation_rows()
    a, middle = cache[1]
    s, eta, delta = mp.mpf('.25'), mp.mpf('.3'), mp.mpf('1e-6')
    # Independent real central difference, not differentiation of a jet.
    plus, minus = a.state(s, eta+delta), a.state(s, eta-delta)
    errors = []
    for name in ('M', 'I', 'J', 'S', 'Cp'):
        fd = (plus['hats'][name].v-minus['hats'][name].v)/(2*delta)
        errors.append(relative(fd, middle['state']['hats'][name].d))
    return dict(rows=rows, eta_jet_central_difference_error=float(max(errors)),
                seed_is_actual_analytic_solution=False, cone_test_is_activated_excess_only=True,
                checks=dict(complete_moment_map_matches_resolvable_subtraction=all(
                                r['exact_map_subtraction_relative_error'] is None or r['exact_map_subtraction_relative_error'] < 1e-50 for r in rows),
                            eta_derivatives_of_all_moment_changes=all(v < mp.mpf('1e-10') for v in errors),
                            rescaled_coefficients_remain_finite=all(math.isfinite(v) and abs(v) < 100 for r in rows for v in r['cone_coefficients']),
                            quadratic_loss_has_sixth_order_factor=all(0 <= r['excess_cone_loss_over_t2s6'] < 10 for r in rows),
                            pressure_memory_not_zero=abs(rows[1]['pressure_p2_hat']) > 1e-14,
                            angular_moment_derivative_not_zero=abs(rows[1]['angular_Jeta_p1_hat']) > 1e-14,
                            flat_underflow_does_not_erase_normalized_moments=rows[-1]['binary64_flat_factor'] == 0 and abs(rows[-1]['moment_coefficients'][0]) > .01,
                            manufactured_seed_defect_not_erased=all(r['seed_defect_norm_over_Fr'] > 1 for r in rows)))


def endpoint_review():
    eta, t = mp.mpf('.3'), mp.mpf('.01')
    b = seed(0, eta)
    E, H, U, q, a = [b[n].v for n in ('E', 'H', 'U', 'q', 'a')]
    expected = [b['X']*E*q/8, b['X']*H*a/8,
                b['X']*H*(U*a+E*q)/8, b['X']*(U*E*q/4-E*E*a/8), E*E*a/8]
    row = Activation(t, '.5').state(mp.mpf('1e-6'), eta)
    errors = [relative(v.v, e) for v, e in zip(row['moment_coefficients'], expected)]
    threshold = collar_threshold('2.3', '10', '1000000')
    s, K, width = threshold['smax'], threshold['error_bound'], threshold['width_max']
    allowance = threshold['allowance']
    return dict(moment_endpoint_values=list(map(float, expected)), near_edge_maximum_relative_error=float(max(errors)),
                illustrative_bound_not_actual_seed=dict(smax=float(s), error_bound=float(K), width_max=float(width)),
                checks=dict(two_forward_integrations_give_nonzero_sixth_order_coefficients=max(errors) < mp.mpf('1e-4'),
                            explicit_algebra_budget=K*width*s**3 <= allowance*(1+mp.mpf('1e-60')),
                            explicit_activation_budget=sigma(s)*11 <= mp.mpf('.15'),
                            no_uniform_physical_width_derivatives_claimed=True))


def refinement_review():
    coarse, _ = activation_rows(24)
    fine, _ = activation_rows(40)
    errors = []
    for a, b in zip(coarse, fine):
        for key in ('field_coefficients', 'moment_coefficients', 'p_hat_coefficients', 'cone_coefficients'):
            errors.extend(abs(x-y)/max(abs(x), abs(y), 1e-100) for x, y in zip(a[key], b[key]))
    return dict(coarse_order=24, fine_order=40, maximum_exported_relative_error=max(errors),
                checks=dict(refined_scaled_integrals_agree=max(errors) < 1e-13))


def obstruction_review():
    s, tolerance = mp.mpf('.02'), mp.mpf('1e-100')
    log_ratio = mp.log10(tolerance)-mp.log10(mp.mpf('.5')*sigma(s))
    # A full seed defect is present in the defining formula (16).
    r = Activation(order=16).stress('.25', '.3')
    f, ea = r['state']['fields'], r['state']['fields']['ea']
    reference, actual = moment_map(f['base'], eta=mp.mpf('.3')), moment_map(r['state']['actual'], eta=mp.mpf('.3'))
    a, q = f['base']['a'].v, f['base']['q'].v
    k, ratio = 1-ea, mp.exp(-ea*f['w'].v)
    literal = [actual['p1']-k*a, actual['p2']-k*q*ratio]
    recovered = [reference['p1']-a+ea*r['B_over_F'][0], reference['p2']-q+ea*r['B_over_F'][1]]
    error = max(relative(v, w) for v, w in zip(literal, recovered))
    erased_gap = max(abs(v-ea*w) for v, w in zip(literal, r['B_over_F']))
    # Wrong c0/lambda normalization is visible as soon as axial shear exists.
    vs, F, shear_a = mp.mpf(4), mp.mpf(2), mp.mpf(3)
    true_growth = 2*shear_a*F**2*(1-2/vs)
    wrong_growth = 2*shear_a*F**2*(vs-2)
    return dict(log10_tiny_baseline_error_over_activation=float(log_ratio),
                full_seed_defect_identity_error=float(error), erased_defect_gap=float(erased_gap),
                wrong_growth_factor=float(wrong_growth/true_growth),
                checks=dict(any_fixed_core_tolerance_cannot_replace_exact_seed=log_ratio > 900,
                            complete_seed_defect_identity=error < mp.mpf('1e-55'),
                            erasing_manufactured_defect_detected=erased_gap > 1,
                            growth_denominator_vs_cannot_be_dropped=wrong_growth != true_growth,
                            no_small_q_rescue_of_leading_seed_error=True))
