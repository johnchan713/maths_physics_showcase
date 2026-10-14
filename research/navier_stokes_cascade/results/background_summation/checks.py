"""Independent coefficient, conservative-moment, curl and summation controls."""
from fractions import Fraction
from itertools import product
import math
from types import SimpleNamespace
import sys
import mpmath as mp
import numpy as np

from hierarchy import (HERE, DIAGONAL, matrices, lower_sources, regular_v,
                       finite_residual_exponent)
from summation import (cutoff_log_scale, scalar_coefficients,
                       scalar_residual_coefficients, loss_budget)
from residual import (Poly, Jet, compose, coordinates, physical_point,
                      profile_jets, series_coefficients, axial_second)
sys.path.insert(0, str(HERE.parent/'first_correction_construction'))
from system import correction_vector, matrices as first_matrices
sys.path.insert(0, str(HERE.parent/'outer_pressure_pilot'))
from schedule import step
sys.path.insert(0, str(HERE))


def fixtures(count):
    """Compact polynomial F,U,Pi fixtures; no coefficient equations assumed."""
    x, e = Poly({(1, 0): 1}), Poly({(0, 1): 1})
    wall = Poly.cast(1)
    for _ in range(4):
        wall *= 1+(-1)*x
    profiles = []
    for n in range(count):
        F = wall*(1+e*(1/(n+2)))
        F = F*(1+e*e*(.03*(n+1)))*(1/(n+1))
        U = x*wall*(.3+e*(.4/(n+1))+e*e*(.02*(n+1)))
        Pi = wall*(.1/(n+1)+e*.03+e*e*.04)
        profiles.append(SimpleNamespace(F=F, U=U, Pi=Pi, average_U=U.average()))
    return profiles


def matrix_review():
    worst, first_error, nilpotence, axis_max, cases = 0., 0., 0., 0., 0
    for n in (1, 2, 3, 5):
        data = fixtures(n+1)
        for h, xi, eta, nu in product((.005, .009), (.1, .5, 1.1), (-.6, 0., .5), (.4, 1., 1.7)):
            a0, a1, f = matrices(data[:n], xi, eta, h, nu)
            w, wx, we = correction_vector(data[n], xi, eta)
            defect = wx+DIAGONAL*w/xi-a0@w-a1@we-f
            coef = series_coefficients(data, xi*xi, eta, h, nu)
            r, theta, axial = [coef[k][n] for k in ('radial', 'theta', 'axial')]
            expected = np.array([0., 0., 0., r/xi, -2*theta/nu,
                                 -2*(axial+eta*r/(1-2*h*eta*eta))/nu])
            worst = max(worst, float(np.max(np.abs(defect-expected)/np.maximum(1., abs(expected)))))
            _, right, _ = matrices(data[:n], xi/2, eta/2, h, nu)
            nilpotence = max(nilpotence, float(np.max(np.abs(a1@np.diag(np.arange(1., 7.))@right))))
            if n == 1:
                older = first_matrices(data[0], xi, eta, h, nu)
                first_error = max(first_error, *(float(np.max(abs(a-b))) for a, b in zip((a0, a1, f), older)))
            cases += 1
        at_axis = matrices(data[:n], 0., .4, .005)
        if not all(np.isfinite(value).all() for value in at_axis):
            raise ValueError('A regular general-order axis source became singular')
        axis_max = max(axis_max, *(float(np.max(abs(value))) for value in at_axis))
    controls = []
    for n in (2, 3, 5):
        data = fixtures(n+1)
        source = lower_sources(data[:n], .21, .43, .005)
        wrong = lower_sources(data[:n], .21, .43, .005, omit_radial_axial_diffusion=True)
        controls.append(dict(order=n, omitted_radial_axial_pressure_gap=abs(source['pressure']-wrong['pressure']),
                             expected_half_axial_source=abs(source['radial_axial_term_over_X'])/2))
    return dict(manufactured=True, sample_count=cases, tested_orders=[1, 2, 3, 5],
                maximum_relative_system_error=worst, first_order_regression_error=first_error,
                sparse_block_error=nilpotence, largest_finite_axis_entry=axis_max,
                controls=controls, checks=dict(independent_complete_residual_system=worst < 3e-11,
                                              first_order_system_preserved=first_error < 3e-12,
                                              sparse_block_at_every_order=nilpotence == 0.,
                                              shifted_radial_viscosity_detected=all(x['omitted_radial_axial_pressure_gap'] > 1e-5 for x in controls)))


def endpoint(poly):
    values = {}
    for (_, j), value in poly.c.items():
        values[(0, j)] = values.get((0, j), 0)+value
    return Poly(values)


def moment_review():
    """Conservative identities include lower viscous moments when nonzero."""
    x = Poly({(1, 0): 1})
    nodes, weights = np.polynomial.legendre.leggauss(40)
    rows = []
    for n in (1, 2, 3, 4):
        data = fixtures(n+1)
        products = [(i, n-i) for i in range(n+1)]
        m1 = endpoint(data[n].U.integral())
        m2 = endpoint((2*x*data[n].F).integral())
        m4 = endpoint((2*x*sum((data[i].U*data[j].F for i, j in products), Poly.cast(0))).integral())
        m5 = endpoint((sum((data[i].U*data[j].U for i, j in products), Poly.cast(0))+data[n].Pi).integral())
        lower1 = endpoint(data[n-1].U.integral())
        lower2 = endpoint((2*x*data[n-1].F).integral())
        for h, eta, nu in product((.005, .009), (-.7, .1, .6), (.4, 1.7)):
            A, D, lam = .5+h, .5-h, 2*h*n
            L, d = 1-2*h*eta*eta, 1-eta*eta
            vals = [p(0., eta) for p in (m1, m2, m4, m5)]
            angular = [p.derivative(1)(0., eta) for p in (m1, m2, m4, m5)]
            ktheta, ltheta, kz, lz = 1.5-A+lam, 1.5-2*A+lam, 1-A+lam, 1-2*A+lam
            expected = np.array([(-ktheta*vals[1]+D*eta*angular[1]+2*ltheta*eta*vals[2]+d*angular[2])/L,
                                 (-kz*vals[0]+D*eta*angular[0]+2*lz*eta*vals[3]+d*angular[3])/L])
            e = Jet.variable(eta, 1)
            viscous = nu*np.array([axial_second(ktheta-2*h, Jet.cast(lower2(0., e)), 0., eta, h),
                                   axial_second(kz-2*h, Jet.cast(lower1(0., e)), 0., eta, h)])
            expected -= viscous
            integrated = np.zeros(2)
            for node, weight in zip(nodes, weights):
                X = (node+1)/2
                coeff = series_coefficients(data, X, eta, h, nu)
                integrated += weight/2*np.array([2*X*coeff['theta'][n], coeff['axial'][n]])
            rows.append(dict(order=n, h=h, eta=eta, nu=nu,
                             identity_error=float(np.max(abs(integrated-expected))),
                             lower_viscous_moment_size=float(np.max(abs(viscous)))))
    return dict(manufactured=True, sample_count=len(rows), rows=rows,
                checks=dict(all_general_conservative_identities=all(x['identity_error'] < 3e-11 for x in rows),
                            nonzero_lower_viscous_moments_exercised=max(x['lower_viscous_moment_size'] for x in rows) > 1e-3))


def cutoff_derivatives(sigma):
    """The actual smooth step, with exact zero/one plateau handling."""
    with mp.workdps(60):
        s = mp.mpf(sigma)
        if s <= mp.mpf('.5'):
            return [1., 0., 0., 0.]
        if s >= 1:
            return [0., 0., 0., 0.]
        function = lambda v: 1-step(2*v-1)
        return [float(mp.diff(function, s, k)) for k in range(4)]


def cutoff_review():
    worst, curl_error, omitted, sign_gap, count = 0., 0., 0., 0., 0
    profile = fixtures(2)[1]
    for n, h, eta, X, sigma_value in product((1, 2, 4), (.005, .009), (-.6, .5), (0., .4), (.2, .65, .85, 1.2)):
        q_value, angle = .12, .4
        point = physical_point(q_value, eta, X, angle, h)
        x, y, z, t, q, xx, ee = coordinates(point, h)
        A, D, lam = .5+h, .5-h, 2*h*n
        scale = sigma_value/q_value
        sigma = scale*q
        c0, c1, c2, c3 = cutoff_derivatives(sigma.v)
        chi = compose(sigma, c0, c1, c2)
        chi_prime = compose(sigma, c1, c2, c3)
        average = profile.average_U(xx, ee)
        average_eta = profile.average_U.derivative(1)(xx, ee)
        v = (2*ee*profile.U(xx, ee)-2*(D+lam)*ee*average-(1-ee*ee)*average_eta)/(1-2*h*ee*ee)
        radial = q**(lam-1)*v/2
        swirl = q**(-A-.5+lam)*profile.F(xx, ee)
        direct = [chi*(radial*x-swirl*y), chi*(radial*y+swirl*x), chi*q**(-A+lam)*profile.U(xx, ee)]
        extra = -ee*sigma*chi_prime*q**(lam-1)*average/(1-2*h*ee*ee)
        corrected = [direct[0]+extra*x, direct[1]+extra*y, direct[2]]
        wrong_sign = [direct[0]-extra*x, direct[1]-extra*y, direct[2]]
        divscale = max(1., *(abs(u.g[i]) for i, u in enumerate(corrected)))
        worst = max(worst, abs(sum(u.g[i] for i, u in enumerate(corrected)))/divscale)
        omitted = max(omitted, abs(sum(u.g[i] for i, u in enumerate(direct)))/divscale)
        sign_gap = max(sign_gap, abs(sum(u.g[i] for i, u in enumerate(wrong_sign)))/divscale)
        potential = [-chi*q**(-A+lam)*average*y/2, chi*q**(-A+lam)*average*x/2, Jet(0.)]
        curl = [potential[2].g[1]-potential[1].g[2], potential[0].g[2]-potential[2].g[0], potential[1].g[0]-potential[0].g[1]]
        B = [-chi*swirl*y, chi*swirl*x, Jet(0.)]
        for i in range(3):
            curl_error = max(curl_error, abs(curl[i]+B[i].v-corrected[i].v))
        curl_gradient = np.array([potential[2].H[1]-potential[1].H[2],
                                  potential[0].H[2]-potential[2].H[0],
                                  potential[1].H[0]-potential[0].H[1]])
        curl_error = max(curl_error, float(np.max(abs(curl_gradient+np.array([v.g for v in B])-np.array([v.g for v in corrected])))))
        count += 1
    return dict(manufactured=True, sample_count=count, includes_axis_and_both_plateaus=True,
                maximum_relative_divergence=worst, maximum_curl_value_or_gradient_error=curl_error,
                missing_commutator_divergence=omitted, wrong_commutator_sign_divergence=sign_gap,
                checks=dict(exact_curl_cutoff_is_divergence_free=worst < 3e-11,
                            separate_potential_curl_agrees=curl_error < 3e-11,
                            direct_velocity_cutoff_rejected=omitted > 1e-4,
                            wrong_commutator_sign_rejected=sign_gap > 1e-4))


def scale_review():
    records = []
    for digits in (80, 110):
        with mp.workdps(digits):
            rows, previous = [], mp.mpf(0)
            for index, gain, constraints in ((1, '.005', [(10, 0), (20, 3)]),
                                             (3, '1e-20', [('1e6', 10000)]),
                                             (10, '2.5', [(3000, 7)])):
                s = cutoff_log_scale(index, gain, constraints, previous, 2)
                gaps = []
                for sample in (s, s+1, 2*s, s+10000):
                    for logC, power in constraints:
                        gap = -index*mp.log(2)-(mp.mpf(logC)+power*mp.log1p(sample)-mp.mpf(gain)*sample)
                        gaps.append(gap)
                rows.append(dict(index=index, half_gain=gain, log_scale=mp.nstr(s, 45),
                                 smallest_log_margin=mp.nstr(min(gaps), 30),
                                 scale_spacing=mp.nstr(s-previous, 30),
                                 checks_pass=bool(min(gaps) >= 0 and s-previous >= mp.log(2))))
                previous = s
            records.append(dict(digits=digits, rows=rows))
    return dict(records=records, actual_NS_coefficient_bounds_or_scales_evaluated=False,
                checks=dict(logarithmic_recipe_dominates_every_sample=all(x['checks_pass'] for r in records for x in r['rows'])))


def scalar_summation_review():
    with mp.workdps(80):
        coefficients = scalar_coefficients(32)
        scales, previous = [mp.mpf(0)], mp.mpf(0)
        for n in range(1, 33):
            logC = mp.log(coefficients[n])+100*mp.log(10)+3*mp.log(n+3)
            previous = cutoff_log_scale(n, mp.mpf(n)/2, [(logC, 0)], previous)
            scales.append(previous)
        rows = []
        for active in (4, 8, 12, 16):
            logq = -scales[active]-mp.mpf('.3')
            cuts, derivatives = [mp.mpf(0)], [mp.mpf(0)]
            for n in range(1, active+1):
                sigma = mp.exp(scales[n]+logq)
                if sigma <= mp.mpf('.5'):
                    chi, derivative = mp.mpf(1), mp.mpf(0)
                else:
                    function = lambda s: 1-step(2*s-1)
                    chi, derivative = function(sigma), sigma*mp.diff(function, sigma)
                cuts.append(chi); derivatives.append(derivative)
            grouped = scalar_residual_coefficients(coefficients[:active+1], cuts, derivatives)
            omitted = scalar_residual_coefficients(coefficients[:active+1], cuts, [0]*(active+1))
            value = mp.fsum(c*mp.exp(k*logq) for k, c in enumerate(grouped) if c)
            order = mp.log(abs(value))/logq
            gap = max(abs(a-b) for a, b in zip(grouped, omitted))
            rows.append(dict(active_index=active, log_q=mp.nstr(logq, 30),
                             effective_scalar_residual_order=mp.nstr(order, 30),
                             lower_coefficients_exactly_zero=all(c == 0 for c in grouped[:active]),
                             missing_cutoff_derivative_coefficient_gap=mp.nstr(gap, 30),
                             next_cutoff_is_zero=bool(scales[active+1]+logq >= 0)))
        polynomial_cancellations = []
        for N in (1, 2, 5, 10):
            grouped = scalar_residual_coefficients(coefficients[:N+1], [0]+[1]*N, [0]*(N+1))
            polynomial_cancellations.append(all(value == 0 for value in grouped[:N+1]))
    return dict(scalar_ODE_only=True, differential_polynomial='q^2*f_prime-f+q+f^2',
                formal_series_is_divergent_at_every_positive_q=True,
                factorial_lower_bound_verified_through=32, rows=rows,
                checks=dict(formal_recurrence_cancels_all_retained_orders=all(polynomial_cancellations),
                            factorial_growth=all(coefficients[n] >= math.factorial(n-1) for n in range(1, 33)),
                            locally_finite_cutoff_sum=all(x['next_cutoff_is_zero'] for x in rows),
                            retained_order_cancellations_survive=all(x['lower_coefficients_exactly_zero'] for x in rows),
                            improving_scalar_residual_order=all(mp.mpf(x['effective_scalar_residual_order']) > x['active_index']-1 for x in rows),
                            cutoff_derivative_omission_detected=all(mp.mpf(x['missing_cutoff_derivative_coefficient_gap']) > 1 for x in rows)))


def exponent_review():
    rows = []
    for h, m, target in product((Fraction(1, 200), Fraction(1, 1000)), (0, 1, 3), (0, 2, 5)):
        row = loss_budget(h, m, target)
        row['direct_finite_exponent'] = finite_residual_exponent(h, row['comparison_order'], m)
        rows.append(row)
    return dict(rows=rows,
                checks=dict(derivative_losses_do_not_grow_with_truncation=all(x['input_tuple_loss'] == 1+2*x['h']+x['derivative_order'] for x in rows),
                            finite_residual_and_tail_margins=all(x['check_finite_margin'] and x['check_tail_margin'] for x in rows),
                            independent_finite_exponents_agree=all(x['direct_finite_exponent'] == x['finite_residual_exponent'] for x in rows)))
