"""Finite diagnostics; the exact all-angle inequalities are proved in README."""
from fractions import Fraction
import math
import mpmath as mp
from tail import HeatTail, axial_coefficients, heat_derivative, positive_coefficient, step_jets


def relative(a, b):
    return abs(a-b)/max(abs(a), abs(b)) if a or b else mp.mpf(0)


def algebra_review():
    errors = []
    for h in (Fraction(1, 200), Fraction(9, 1000), Fraction(1, 1000000)):
        for eta in map(Fraction, (-1, Fraction(-3, 5), 0, Fraction(1, 2), 1)):
            D, L, d = Fraction(1, 2)-h, 1-2*h*eta*eta, 1-eta*eta
            alpha = 4*eta*eta/L**2
            beta = -2*(d-2*D*eta*eta)/L**2-8*h*eta*eta*d/L**3
            errors.append(alpha*(1-h)-beta == 2/L**2+8*h*eta*eta*d/L**3)
    # Exact coefficient identity before multiplying by denominators:
    # [4(1-h)-4D] eta^2 + 2(1-eta^2) = 2.
    h = Fraction(7, 1000)
    constant, eta_squared = 2, 4*(1-h)-4*(Fraction(1, 2)-h)-2
    return dict(exact_test_count=len(errors), polynomial_coefficients=[str(constant), str(eta_squared)],
                checks=dict(coefficient_cancellation=all(errors) and eta_squared == 0,
                            beta_at_eta_zero_is_minus_two=axial_coefficients('.005', 0)[2] == -2))


def uniform_bounds_review():
    mp.iv.dps = 70
    p, exp = mp.iv.mpf, mp.iv.exp
    Lmin, hmax = p('.98'), p('.01')
    tH_over_h = (2*(1+hmax)/100)/(1-2*hmax*(1+hmax)/100)
    slope_over_h = p('.032')/(1-p('.001')*hmax)
    # C=alpha*(1-h+tH)-beta; 1-h+tH<=1 and |d-2D eta^2|<=1.
    Cmax = 6/Lmin**2+8*hmax/Lmin**3
    tail_step = (128+p(128)/27)*exp(p(16)/9-16)
    amplitude_log_ratio = -p('7.5')*64+p('10.3')
    cone_loss = 32*hmax**9
    lo = lambda v: mp.make_mpf(v._mpi_[0])
    hi = lambda v: mp.make_mpf(v._mpi_[1])
    up = lambda v: math.nextafter(float(hi(v)), math.inf)
    return dict(heat_log_slope_over_h_upper=up(tH_over_h), cutoff_log_slope_over_h_upper=up(slope_over_h),
                positive_integral_coefficient_upper=up(Cmax),
                log_e_star_over_h14_upper=up(amplitude_log_ratio), cone_loss_upper=up(cone_loss),
                interval_digits=70, actual_h_or_Xtail_numerically_evaluated=False,
                checks=dict(heat_slope_below_h_quarter=hi(tH_over_h) < p('.25'),
                            cutoff_slope_below_h_quarter=hi(slope_over_h) < p('.25'),
                            positive_integral_coefficient_below_seven=hi(Cmax) < 7,
                            step_tail_derivative_below_64=hi(tail_step) < 64,
                            actual_schedule_amplitude_bound=hi(amplitude_log_ratio) < 0,
                            actual_schedule_outer_cone_gap_above_1p99=hi(cone_loss) < p('.01')))


def heat_review():
    h, rows, max_error = mp.mpf('.005'), [], mp.mpf(0)
    for Z in map(mp.mpf, ('0', '.001', '.02')):
        exact = heat_derivative(h, Z)
        integral = mp.quad(lambda v: mp.exp(-v)*v**h*(1+Z*v)**(-h), [0, 1, mp.inf])/mp.gamma(1+h)
        Hp, Hpp = heat_derivative(h, Z, 1), heat_derivative(h, Z, 2)
        ode = Z*Z*Hpp+(1+2*(1+h)*Z)*Hp+h*(1+h)*exact
        gap = relative(exact, integral)
        max_error = max(max_error, gap, abs(ode))
        rows.append(dict(Z=float(Z), H=float(exact), independent_integral_relative_error=float(gap),
                         heat_ODE_error=float(abs(ode))))
    return dict(rows=rows, exact_integral_not_raw_Taylor_series=True,
                checks=dict(independent_heat_integral_agrees=max_error < mp.mpf('1e-50'),
                            positive_heat_values=all(0 < r['H'] <= 1 for r in rows)))


def stress_review(refined=False):
    tail = HeatTail(order=32 if refined else 24, panels=10 if refined else 8)
    rows, signed_error, positive = [], mp.mpf(0), True
    for y, eta in (('.5', '-1'), ('1.4', '0'), ('2.2', '.6'), ('2.2', '1')):
        r = tail.direct(y, eta)
        gap = relative(r['T1signed'], r['T1positive'])
        signed_error = max(signed_error, gap)
        cone0 = 2-(r['a']-2)*r['leading_ratio']**2
        updated = []
        for c in map(mp.mpf, ('0', '0.01', '1', '1000000')):
            theta = r['T0theta']+c*r['T1positive']
            ratio = r['T0z']/theta
            updated.append(dict(coefficient=float(c), ratio=float(ratio),
                                directional_loss=float((r['a']-2)*ratio**2)))
        positive &= (r['T1positive'] > 0 and r['correction_ratio'] <= r['uniform_ratio_bound']
                     and 2+tail.h < r['a'] <= 2+2*tail.h and abs(r['leading_ratio']) <= 2*r['K'])
        positive &= all(row['directional_loss'] <= float(2-cone0)*(1+1e-14) for row in updated)
        rows.append(dict(y=y, eta=eta, T0theta=float(r['T0theta']), T0z=float(r['T0z']),
                         T1signed=float(r['T1signed']), T1positive=float(r['T1positive']),
                         correction_ratio=float(r['correction_ratio']), updated=updated))
    return dict(rows=rows, maximum_signed_positive_relative_error=float(signed_error),
                parameters_are_diagnostic_not_actual_giant_schedule=True,
                checks=dict(signed_and_positive_stress_agree=signed_error < mp.mpf('1e-12'),
                            positive_correction_and_improved_cone=bool(positive)))


def physical_source_review():
    tail, rows, source_worst, stress_worst = HeatTail(order=32, panels=10), [], mp.mpf(0), mp.mpf(0)
    for q0, eta in (('1', '-1'), ('0.4', '0'), ('0.2', '.6'), ('0.7', '1')):
        q0, eta = mp.mpf(q0), mp.mpf(eta)
        y, h, D = mp.mpf('2.2'), tail.h, tail.D
        s = tail.state(y, eta)
        tau, z0, radius = q0*(1-eta*eta), eta*q0**D, mp.sqrt(q0)*s['r']
        # Kphysical is independent of z at fixed radius and tau. Only the
        # terminal cutoff changes as q(z,t) changes.
        Kphysical = q0**(-tail.A)*s['K']
        def implicit_q(z):
            return mp.findroot(lambda q: q-z*z*q**(2*h)-tau, q0)
        def cutoff_at_z(z):
            q = implicit_q(z)
            yy = mp.log(radius*radius/(2*q*tail.xtail))
            return 1-tail.rho*step_jets((3-yy)/2)[0]
        observed = Kphysical*mp.diff(cutoff_at_z, z0, 2)
        L, alpha, beta = axial_coefficients(h, eta)
        expected = q0**(-tail.A-2*D)*s['K']*(alpha*s['fpp']+beta*s['fp'])
        normalizer = q0**(-tail.A-2*D)*s['K']*tail.rho
        gap = abs(observed-expected)/normalizer
        # The cylindrical stress identity is (d_r+2/r)T1=-r_theta,1.
        # A symmetric radial difference checks the integral's cylindrical
        # divergence independently of the source formula. Its finite step
        # and quadrature error are diagnostics, not continuum error bars.
        step_size = mp.mpf('1e-8')
        T = tail.direct(y, eta)['T1positive']
        derivative = (tail.direct(y+step_size, eta)['T1positive']
                      -tail.direct(y-step_size, eta)['T1positive'])/(2*step_size)
        divergence = 2/s['r']*(derivative+T)
        desired = s['K']*(alpha*s['fpp']+beta*s['fp'])
        stress_gap = abs(divergence-desired)/(s['K']*tail.rho)
        source_worst, stress_worst = max(source_worst, gap), max(stress_worst, stress_gap)
        omitted_beta = abs(beta*s['fp'])/tail.rho
        rows.append(dict(q=float(q0), eta=float(eta), physical_source_normalized_error=float(gap),
                         radial_divergence_normalized_error=float(stress_gap),
                         omitted_beta_source_gap=float(omitted_beta)))
    return dict(rows=rows, keeps_nonunit_q_normalization=True,
                checks=dict(independent_physical_axial_derivative=source_worst < mp.mpf('1e-40'),
                            independent_radial_difference_stress_sign=stress_worst < mp.mpf('1e-9'),
                            omitted_beta_detected_at_equator=rows[1]['omitted_beta_source_gap'] > 0))


def edge_review(refined=False):
    tail = HeatTail(order=48 if refined else 32, panels=8)
    rows, last_error, bounded = [], mp.mpf(0), True
    for eta in map(mp.mpf, ('-1', '0', '.5', '1')):
        endpoint = tail.endpoint(eta)
        for delta_string in ('.3', '.03', '.001', '1e-7'):
            delta = mp.mpf(delta_string)
            r = tail.factored(delta, eta)
            error = abs(r['correction_ratio']-endpoint['correction_ratio'])/max(1, abs(endpoint['correction_ratio']))
            if delta_string == '1e-7':
                last_error = max(last_error, error)
            bounded &= 0 <= r['correction_ratio'] <= 2*tail.xb/mp.mpf('.98')**2+7
            rows.append(dict(eta=str(eta), delta=delta_string, b0=float(r['b0']), b1=float(r['b1']),
                             bz=float(r['bz']), correction_ratio=float(r['correction_ratio']),
                             expected_endpoint_ratio=float(endpoint['correction_ratio']),
                             binary64_flat_factor=r['float_flat_factor']))
    # The direct signed integral and flat-coordinate integral use different
    # variables and quadrature rules. Compare where raw values are resolvable.
    direct = tail.direct(mp.mpf('2.7'), mp.mpf('.5'))
    factored = tail.factored(mp.mpf('.3'), mp.mpf('.5'))
    scale = mp.exp(-4/mp.mpf('.3')**2)*mp.mpf('.3')**-3
    crosscheck = max(relative(direct['T0theta'], scale*factored['b0']),
                     relative(direct['T1signed'], scale*factored['b1']))
    return dict(rows=rows, smallest_delta='1e-7', maximum_last_relative_endpoint_error=float(last_error),
                direct_factored_relative_error=float(crosscheck),
                checks=dict(smooth_endpoint_ratio_matches=last_error < mp.mpf('2e-7'),
                            correction_ratio_uniform_to_flat_edge=bool(bounded),
                            direct_signed_and_factored_integrals_agree=crosscheck < mp.mpf('1e-9'),
                            positive_coefficient_survives_float_underflow=all(r['b0'] > 0 and r['binary64_flat_factor'] == 0
                                                                           for r in rows if r['delta'] == '1e-7')))


def failed_shortcuts_review():
    # A generic weighted signed perturbation can escape the positive cone,
    # even though every fixed q power is small away from the edge.
    q, h, delta = mp.mpf('.5'), mp.mpf('.005'), mp.mpf('1e-6')
    relative_signed = 1-q**(2*h)/delta
    # A negative Borel cutoff can also reverse the otherwise positive sum.
    wrong_cutoff_sum = 1-mp.mpf(2)
    return dict(generic_signed_perturbation_over_leading=float(relative_signed),
                negative_cutoff_counterexample=float(wrong_cutoff_sum),
                checks=dict(generic_weighted_smallness_does_not_preserve_cone=relative_signed < 0,
                            negative_cutoff_is_not_covered=wrong_cutoff_sum < 0,
                            no_global_or_complete_residual_claim=True))
