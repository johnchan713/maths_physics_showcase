"""Separate continuum bookkeeping from finite, omitted-term diagnostics."""
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
import sys
import mpmath as mp
import numpy as np

from extension import (HERE, Repair, normalized_target, solve_patch,
                       divergence_identity, radial_flux, geometry_ledger)
from fixtures import MemoryFixture, Rule, direct_patch_changes
sys.path.insert(0, str(HERE.parent/'physical_residual_budget'))
from residual import (Poly, Jet, series_coefficients, profile_jets,
                      time_operator, axial_operator, axial_second)
from review import fixtures as polynomial_fixtures
sys.path.insert(0, str(HERE))


def kinematic_review():
    rows = []
    for h, eta, X in ((Fraction(1, 200), Fraction(2, 3), Fraction(5, 7)),
                      (Fraction(1, 125), Fraction(-3, 4), Fraction(17, 2))):
        error = divergence_identity(Fraction(3, 5), Fraction(-2, 7), Fraction(4, 9), X, eta, h)
        bad = divergence_identity(Fraction(3, 5), Fraction(-2, 7), Fraction(4, 9), X, eta, h, omit_shift=True)
        memory_flux = radial_flux(0, eta, 1, X, eta, h, order=0)
        rows.append(dict(h=str(h), eta=str(eta), exact_divergence=str(error),
                         wrong_shift_defect=str(bad), leading_constant_memory_flux=str(memory_flux)))
    return dict(rows=rows, checks=dict(exact_positive_order_divergence=all(x['exact_divergence'] == '0' for x in rows),
                                      omitted_lambda_detected=all(Fraction(x['wrong_shift_defect']) != 0 for x in rows),
                                      zero_U_not_zero_V=all(Fraction(x['leading_constant_memory_flux']) == -1 for x in rows)))


def global_matching_review():
    rows = []
    with mp.workdps(80):
        # Refined integration is separate from the matrix's own quadrature.
        rule = Rule(48, 8)
        direct_rule = Rule(48, 12)
        for lam in ('.001', '.01'):
            repair = Repair(lam, order=48, panels=8)
            for radius in ('2', '17'):
                fixture = MemoryFixture(radius)
                for angle in ('-.6', '0', '.6'):
                    eta = mp.mpf(angle)
                    amplitude = mp.sqrt(2)/(1+eta*eta)*fixture.radius**(-mp.mpf('.5')-repair.lam)
                    incoming = fixture.incoming(eta, rule)
                    c = solve_patch(incoming, fixture.radius, amplitude, repair)
                    changes = direct_patch_changes(c, repair, fixture.radius, amplitude, direct_rule)
                    closed = incoming+changes
                    normalized = normalized_target(closed, fixture.radius, amplitude, repair.lam)
                    wrong_input = fixture.incoming(eta, rule, omit_omega=True)
                    wrong_c = solve_patch(wrong_input, fixture.radius, amplitude, repair)
                    wrong_closed = incoming+direct_patch_changes(wrong_c, repair, fixture.radius, amplitude, rule)
                    on_mean = fixture.radius*mp.exp(8)
                    tail = fixture.tail_pressure(on_mean, eta, rule)
                    patch_X = fixture.radius*mp.exp(2)
                    omega_error = abs(fixture.omega(patch_X, eta)+1/(2*patch_X))
                    rows.append(dict(lambda_value=lam, radius=radius, eta=angle,
                                     incoming_ordinary=[mp.nstr(x, 45) for x in incoming],
                                     coefficients=[mp.nstr(x, 45) for x in c],
                                     maximum_normalized_closure_error=mp.nstr(max(map(abs, normalized)), 30),
                                     omitted_omega_pressure_moment=mp.nstr(abs(wrong_closed[2]), 30),
                                     omitted_omega_fifth_moment=mp.nstr(abs(wrong_closed[4]), 30),
                                     pressure_on_reserved_mean_patch=mp.nstr(tail, 30),
                                     exact_patch_omega_error=mp.nstr(omega_error, 30)))
    return dict(manufactured=True, fixture_smoothness='C4 memory; C5 cutoff profiles; C-infinity repair bumps',
                actual_huge_frequency_profile_evaluated=False,
                lower_order_source_retained_through_outer_pulse=True,
                separate_matrix_and_integrand_quadrature=True,
                continuum_existence_not_inferred_from_quadrature=True,
                rows=rows, checks=dict(five_global_moments_close=all(mp.mpf(x['maximum_normalized_closure_error']) < mp.mpf('1e-35') for x in rows),
                                       missing_pressure_source_detected_off_midplane=all(mp.mpf(x['omitted_omega_pressure_moment']) > mp.mpf('.01') for x in rows if x['eta'] != '0'),
                                       midplane_total_pressure_source_cancels=all(mp.mpf(x['omitted_omega_pressure_moment']) < mp.mpf('1e-35') for x in rows if x['eta'] == '0'),
                                       missing_fifth_source_detected=all(mp.mpf(x['omitted_omega_fifth_moment']) > mp.mpf('.1') for x in rows),
                                       pressure_need_not_vanish_on_mean_patch=all(abs(mp.mpf(x['pressure_on_reserved_mean_patch'])) > mp.mpf('1e-8') for x in rows),
                                       nonzero_memory_patch_source=all(mp.mpf(x['exact_patch_omega_error']) < mp.mpf('1e-65') for x in rows)))


def memory_source_review():
    """Compare the scalar memory source with the general residual machinery."""
    fixture = MemoryFixture()
    start, end = float(fixture.pulse_start), float(fixture.support_end)
    coefficients = {5: 126, 6: -420, 7: 540, 8: -315, 9: 70}
    def B(X, derivative=0):
        value = X.v if isinstance(X, Jet) else X
        if value <= 0 or value >= end:
            return 0.
        if 1 <= value <= start:
            return 1. if derivative == 0 else 0.
        width = 1. if value < 1 else end-start
        t = X if value < 1 else (X-start)/width
        polynomial = sum(v*float(mp.factorial(n)/mp.factorial(n-derivative))*t**(n-derivative)
                         for n, v in coefficients.items() if n >= derivative)/width**derivative
        return polynomial if value < 1 else ((1-polynomial) if derivative == 0 else -polynomial)
    class Average:
        def __call__(self, X, eta):
            return eta*B(X)/X
        def derivative(self, axis):
            if axis != 1:
                raise ValueError('Only the angular derivative is needed here')
            return lambda X, eta: B(X)/X
    profile = SimpleNamespace(F=lambda X, eta: 0., U=lambda X, eta: eta*B(X, 1),
                              Pi=lambda X, eta: 0., average_U=Average())
    errors = []
    for X in (.2, .65, float(fixture.radius*mp.exp(2)), (start+end)/2, 2*end):
        for eta in (-.6, 0., .6):
            actual = series_coefficients([profile], X, eta, float(fixture.h))['omega'][0]
            expected = float(fixture.omega(mp.mpf(X), mp.mpf(eta)))
            errors.append(abs(actual-expected)/max(1., abs(expected)))
    return dict(manufactured=True, sample_count=len(errors), maximum_relative_source_error=max(errors),
                checks=dict(independent_radial_source=max(errors) < 1e-10))


def end_polynomial(poly):
    result = {}
    for (_, j), value in poly.c.items():
        result[(0, j)] = result.get((0, j), 0)+value
    return Poly(result)


def moment_identity_review():
    """Direct integrated residual vs conservative moment identities.

    These are different polynomial fixtures from the global patch example.
    Their leading weighted angular and axial moments are zero. Their first
    moments are deliberately NOT repaired, to test every nonzero RHS term.
    """
    x, e = Poly({(1, 0): 1}), Poly({(0, 1): 1})
    factor = Poly.cast(1)+(-1)*x
    wall = factor*factor*factor*factor
    F0 = wall*(1+(-3.5)*x)*(1+e*(1/3))
    F1 = x*wall*(1+(-.2)*e)
    U1 = x*wall*(1+e+e*e*(1/7))
    P0 = (F0*F0).integral()
    P0 = P0+(-1)*end_polynomial(P0)
    P1 = (2*F0*F1).integral()
    P1 = P1+(-1)*end_polynomial(P1)
    zero = Poly.cast(0)
    base = SimpleNamespace(F=F0, U=zero, Pi=P0, average_U=zero)
    correction = SimpleNamespace(F=F1, U=U1, Pi=P1, average_U=U1.average())
    moments = [end_polynomial(p.integral()) for p in
               (U1, 2*x*F1, 2*F0*F1, 2*x*U1*F0, (-2)*x*F0*F1)]
    nodes, weights = np.polynomial.legendre.leggauss(48)
    rows = []
    for h in (.005, .009):
        A, D, lam = .5+h, .5-h, 2*h
        for eta in (-.7, 0., .6):
            L, d = 1-2*h*eta*eta, 1-eta*eta
            m1, m2, _, m4, m5 = [p(0., eta) for p in moments]
            me = [p.derivative(1)(0., eta) for p in moments]
            expected = np.array([(-(1+h)*m2+D*eta*me[1]+eta*m4+d*me[3])/L,
                                 (-A*m1+D*eta*me[0]+d*me[4])/L])
            for nu in (.4, 1., 1.7):
                integrated = np.zeros(2)
                for node, weight in zip(nodes, weights):
                    X = (node+1)/2
                    coef = series_coefficients([base, correction], X, eta, h, nu)
                    integrated += weight/2*np.array([2*X*coef['theta'][1], coef['axial'][1]])
                rows.append(dict(h=h, eta=eta, nu=nu,
                                 maximum_absolute_identity_error=float(np.max(np.abs(integrated-expected))),
                                 nonzero_moment_prediction=list(map(float, expected))))
    return dict(manufactured=True, sample_count=len(rows), rows=rows,
                leading_weighted_angular_moment=float(end_polynomial((2*x*F0).integral())(0., .3)),
                checks=dict(direct_residual_integrals_match=all(row['maximum_absolute_identity_error'] < 1e-10 for row in rows),
                            leading_lower_angular_moment_zero=abs(end_polynomial((2*x*F0).integral())(0., .3)) < 1e-13,
                            nonzero_RHS_is_exercised=all(max(map(abs, row['nonzero_moment_prediction'])) > 1e-4 for row in rows)))


def remainder_review():
    """Retain the quadratic and coefficient axial-diffusion remainder terms."""
    base, first = polynomial_fixtures()
    worst, quadratic_gap, diffusion_gap, count = 0., 0., 0., 0
    for h in (.005, .009):
        A, lam = .5+h, 2*h
        b, c = -A-.5, -A
        for X, eta in ((.2, -.6), (.7, 0.), (1.4, .5)):
            p0, p1 = [profile_jets(p, X, eta, h, lam*n) for n, p in enumerate((base, first))]
            V0, V1, U0, U1, F1 = p0['V'], p1['V'], p0['U'], p1['U'], p1['F']
            for nu in (.4, 1., 1.7):
                omega1 = (time_operator(lam, V1, X, eta, h)-2*nu*X*V1.H[0, 0]
                          +V0.v*(V1.g[0]-V1.v/(2*X))+V1.v*(V0.g[0]-V0.v/(2*X))
                          +U0.v*axial_operator(lam, V1, X, eta, h)+U1.v*axial_operator(0., V0, X, eta, h))
                omega2 = V1.v*(V1.g[0]-V1.v/(2*X))+U1.v*axial_operator(lam, V1, X, eta, h)
                theta_quadratic = V1.v*(F1.g[0]+F1.v/X)+U1.v*axial_operator(b+lam, F1, X, eta, h)
                axial_quadratic = V1.v*U1.g[0]+U1.v*axial_operator(c+lam, U1, X, eta, h)
                theta_viscous = nu*axial_second(b+lam, F1, X, eta, h)
                axial_viscous = nu*axial_second(c+lam, U1, X, eta, h)
                expected = np.array([-2*X*F1.v**2+omega1-nu*axial_second(0., V0, X, eta, h),
                                     omega2-nu*axial_second(lam, V1, X, eta, h),
                                     theta_quadratic-theta_viscous, axial_quadratic-axial_viscous])
                coeff = series_coefficients([base, first], X, eta, h, nu)
                actual = np.array([coeff['radial'][2], coeff['radial'][3], coeff['theta'][2], coeff['axial'][2]])
                worst = max(worst, float(np.max(np.abs(actual-expected))))
                quadratic_gap = max(quadratic_gap, abs(theta_quadratic), abs(axial_quadratic))
                diffusion_gap = max(diffusion_gap, abs(theta_viscous), abs(axial_viscous))
                count += 1
    return dict(manufactured=True, first_two_coefficients_not_asserted_zero_for_these_fixtures=True,
                sample_count=count, maximum_absolute_identity_error=worst,
                missing_quadratic_gap=quadratic_gap, missing_coefficient_axial_diffusion_gap=diffusion_gap,
                checks=dict(exact_finite_remainder=worst < 2e-12,
                            quadratic_terms_cannot_be_dropped=quadratic_gap > 1e-4,
                            correction_axial_diffusion_cannot_be_dropped=diffusion_gap > 1e-4))
