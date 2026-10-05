"""Manufactured crosschecks and exact counterexamples for the residual ledger."""
from fractions import Fraction
import math

import numpy as np

from residual import (HERE, Jet, Poly, ManufacturedProfile, axial_operator,
                      axial_second, compose, cylindrical, cutoff_velocity,
                      evaluate_series, field, minimum_order_for_decay, momentum,
                      physical_point, physical_series, profile_jets, radial_profile,
                      series_coefficients, stress_pair,
                      symmetric_cross_stress_divergence)
from profiles import finite_difference_residual

BASE = dict(F=[[0, 0, 1.], [1, 0, .12], [0, 1, .08], [2, 0, .03], [1, 1, -.02]],
            U=[[0, 0, .1], [0, 1, .7], [1, 0, .25], [1, 1, .12], [0, 2, -.05], [2, 0, .03]],
            axis_pressure=[[0, 0, -1.], [0, 1, .2], [0, 2, .1]])
CORRECTION = dict(F=[[0, 0, .07], [1, 1, -.04], [0, 2, .025], [2, 0, .015]],
                  U=[[0, 0, -.09], [0, 1, .11], [1, 0, .06], [1, 2, -.02], [0, 3, .013]],
                  Pi=[[0, 0, .04], [0, 1, .02], [1, 1, .08], [2, 0, .03]])


def fixtures():
    base = ManufacturedProfile(BASE)
    correction = ManufacturedProfile({**CORRECTION, 'axis_pressure': []})
    correction.Pi = Poly({(i, j): v for i, j, v in CORRECTION['Pi']})
    return base, correction


def to_cartesian(vector, theta):
    ct, st = math.cos(theta), math.sin(theta)
    return np.array([ct*vector[0]-st*vector[1], st*vector[0]+ct*vector[1], vector[2]])


def physical_review():
    base, correction = fixtures()
    maxima = dict(base=0., finite_correction=0., divergence=0., old_ansatz=0.,
                  leading_stress=0., full_symmetric_stress=0.)
    cases = 0
    for h in (.005, .2):
        for nu in (.4, 1., 1.7):
            for q in (.08, .3, .8):
                for eta in (-.6, 0., .5):
                    for X in (.2, .7, 1.4):
                        for theta in (0., .4):
                            point = physical_point(q, eta, X, theta, h)
                            for profiles, key in (([base], 'base'), ([base, correction], 'finite_correction')):
                                u, p = physical_series(profiles, point, h)
                                actual, scales = momentum(u, p, nu)
                                coefficients = series_coefficients(profiles, X, eta, h, nu)
                                expected = to_cartesian(evaluate_series(coefficients, q, X, h), theta)
                                maxima[key] = max(maxima[key], float(np.max(np.abs(actual-expected)/scales)))
                                divscale = max(1., *(abs(u[i].g[i]) for i in range(3)))
                                maxima['divergence'] = max(maxima['divergence'], abs(sum(u[i].g[i] for i in range(3)))/divscale)
                            if nu == 1.:
                                old_u, old_p, _ = field(point, h, base)
                                new_u, new_p = physical_series([base], point, h)
                                maxima['old_ansatz'] = max(maxima['old_ansatz'], max(abs(a.v-b.v) for a, b in zip(old_u, new_u)), abs(old_p.v-new_p.v))
                            x, e = Jet.variable(X, 0), Jet.variable(eta, 1)
                            Ttheta, Tz = stress_pair(base, x, e, h, nu)
                            R = math.sqrt(2*X)
                            coefficient = series_coefficients([base], X, eta, h, nu)
                            leading = np.array([R*coefficient['theta'][0], coefficient['axial'][0]])
                            stress = -np.array([R*Ttheta.g[0]+2*Ttheta.v/R, R*Tz.g[0]+Tz.v/R])
                            maxima['leading_stress'] = max(maxima['leading_stress'], float(np.max(np.abs(leading-stress)/np.maximum(1., np.abs(leading)))))
                            if nu == 1.:
                                A = .5+h
                                predicted = np.array([q**(-1.5)*axial_operator(-A-.5, Tz, X, eta, h),
                                                      -q**(-A-1)*leading[0], -q**(-A-1)*leading[1]])
                                actual = cylindrical(symmetric_cross_stress_divergence(base, point, h), theta)
                                maxima['full_symmetric_stress'] = max(maxima['full_symmetric_stress'], float(np.max(np.abs(actual-predicted)/np.maximum(1., np.abs(predicted)))))
                            cases += 1
    return dict(manufactured=True, physical_sample_count=cases, maxima=maxima,
                checks={key: value < 2e-11 for key, value in maxima.items()})


def omission_review():
    base, correction = fixtures()
    h, q, eta, X, theta = .005, .3, .5, .7, .4
    point = physical_point(q, eta, X, theta, h)
    coefficients = series_coefficients([base], X, eta, h)
    full = evaluate_series(coefficients, q, X, h)
    without_axial = {key: list(row) for key, row in coefficients.items()}
    without_axial['theta'][1] = without_axial['axial'][1] = without_axial['radial'][2] = 0.
    missing_axial = evaluate_series(without_axial, q, X, h)
    jets = profile_jets(base, X, eta, h, 0.)
    shifted = axial_second(-.5-h, jets['U'], X, eta, h)
    unshifted = axial_second(-.5-h, jets['U'], X, eta, h, shift_power=False)
    u_bad, _ = physical_series([base, correction], point, h, omit_lambda=True)
    bad_divergence = abs(sum(u_bad[i].g[i] for i in range(3)))
    x, e = Jet.variable(X, 0), Jet.variable(eta, 1)
    _, Tz = stress_pair(base, x, e, h)
    radial_stress = q**(-1.5)*axial_operator(-1-h, Tz, X, eta, h)
    missing_curvature = math.sqrt(2*X)*q**(-1.5-h)*jets['F'].v/(2*X)
    wrong_sign = 2*abs(math.sqrt(2*X)*q**(-1.5-h)*coefficients['theta'][0])
    gaps = dict(omitted_axial_diffusion=float(np.max(np.abs(full-missing_axial)/np.maximum(1., np.abs(full)))),
                wrong_second_axial_power=abs(shifted-unshifted),
                missing_positive_order_lambda=bad_divergence,
                missing_symmetric_rz_axial_divergence=abs(radial_stress),
                omitted_swirl_curvature=missing_curvature,
                reversed_leading_stress_sign=wrong_sign)
    scalar_u, scalar_p = physical_series([base], point, h)
    exact, scale = momentum(scalar_u, scalar_p)
    errors = [float(np.max(np.abs(finite_difference_residual(point, h, base, step)-exact)/scale))
              for step in (.02, .01, .005)]
    return dict(manufactured=True, negative_control_gaps=gaps,
                finite_difference_errors=errors,
                checks={**{key: value > 1e-4 for key, value in gaps.items()},
                        'scalar_only_stencil_accuracy': errors[-1] < 1e-6,
                        'scalar_only_stencil_refinement': errors[0] > 4*errors[-1]})


class AngularSine:
    """U=m^-2 sin(m eta); all derivatives are evaluated without phase offsets."""
    def __init__(self, m, order=0):
        self.m, self.order = m, order

    def derivative(self, axis):
        return AngularSine(self.m, self.order+1) if axis == 1 else Poly({})

    def __call__(self, X, eta):
        value = eta.v if isinstance(eta, Jet) else eta
        trig = (math.sin(self.m*value), math.cos(self.m*value),
                -math.sin(self.m*value), -math.cos(self.m*value))
        def at(order):
            return self.m**(order-2)*trig[order % 4]
        if isinstance(eta, Jet):
            return compose(eta, at(self.order), at(self.order+1), at(self.order+2))
        return at(self.order)


def derivative_gap_review():
    rows = []
    for m in (8, 32, 128, 512):
        class Profile:
            U = average_U = AngularSine(m)
        x, eta = Jet.variable(.7, 0), Jet.variable(0., 1)
        V = radial_profile(Profile(), x, eta, .005)
        actual = axial_second(0., V, .7, 0., .005)
        expected = Fraction(7, 10)*(m+Fraction(6, m))
        rows.append(dict(m=m, C2_component_bounds=[str(Fraction(1, m*m)), str(Fraction(1, m)), '1'],
                         third_angular_amplitude=m,
                         radial_axial_diffusion_coefficient=str(expected),
                         formula_error=abs(actual-float(expected))))
    # A value-small oscillation has a second logarithmic derivative growing with N.
    phase = [dict(N=n, value_amplitude=str(Fraction(1, n)),
                  first_y_amplitude=1, second_y_amplitude=n) for n in (8, 64, 512)]
    return dict(manufactured=True, C2_counterexample=rows, fixed_phase_counterexample=phase,
                checks=dict(C2_is_uniform=all(all(Fraction(v) <= 1 for v in row['C2_component_bounds']) for row in rows),
                            third_derivative_grows=rows[-1]['third_angular_amplitude'] > 50*rows[0]['third_angular_amplitude'],
                            radial_viscosity_formula=all(row['formula_error'] < 1e-12 for row in rows),
                            radial_viscosity_unbounded_in_C2_family=Fraction(rows[-1]['radial_axial_diffusion_coefficient']) > 50*Fraction(rows[0]['radial_axial_diffusion_coefficient']),
                            larger_N_not_a_C2_smallness_argument=phase[-1]['second_y_amplitude'] > phase[0]['second_y_amplitude']))


def cutoff_review():
    base, _ = fixtures()
    point = physical_point(.3, .5, .7, .4, .005)
    naive = cutoff_velocity(base, point, .005, 2., False)
    preserved = cutoff_velocity(base, point, .005, 2., True)
    bad = abs(sum(naive[i].g[i] for i in range(3)))
    good = abs(sum(preserved[i].g[i] for i in range(3)))
    return dict(manufactured=True, diagnostic_cutoff_is_compact=False,
                naive_divergence=bad, streamfunction_divergence=good,
                checks=dict(naive_product_rejected=bad > 1e-4,
                            potential_product_preserves_divergence=good < 1e-12))


def exponent_review():
    h = Fraction(1, 200)
    rows = [dict(time_derivatives=m, minimum_sufficient_correction_order=minimum_order_for_decay(h, time=m))
            for m in (0, 1, 2)]
    return dict(manufactured_parameter_example=True, h=str(h), rows=rows,
                checks=dict(radial_residual_still_singular_at_order_one=2*h-Fraction(3, 2) < 0,
                            tangential_axial_diffusion_still_singular=2*h-(Fraction(3, 2)+h) < 0,
                            every_extra_time_derivative_needs_more_orders=rows[1]['minimum_sufficient_correction_order']-rows[0]['minimum_sufficient_correction_order'] == 100,
                            finite_order_not_all_derivatives=rows[-1]['minimum_sufficient_correction_order'] > rows[0]['minimum_sufficient_correction_order']))
