"""Manufactured full wave/mean potentials with the ACTUAL implicit q cutoff.

These exact Taylor-jet algebra diagnostics do not evaluate the global flow.
Their phase/amplitudes are the inherited compact test fixtures, not its bands.
"""
from pathlib import Path
import sys
import math
import numpy as np

PARENT = Path(__file__).resolve().parent.parent
sys.path.append(str(PARENT/'physical_wave_interface'))
from jets import Jet, INDICES, arrays, cross, variables
import waves

POINTS = ((1.12, -.4, -.19, .36), (1.3, .27, .2, .47), (1.66, -.9, .31, .64))
H = .005


def relative(a, b):
    return float(np.linalg.norm(a-b)/(1+np.linalg.norm(a)+np.linalg.norm(b)))


def q_value(z, time, h=H):
    """Monotone bisection of q-(1-t)-z^2*q^(2h)=0 on its positive branch."""
    tau = 1-time
    if tau <= 0 or not 0 < h < .01:
        raise ValueError('Diagnostic uses tau>0 and inherited h')
    lower = max(tau, abs(z)**(1/(.5-h)))
    upper = 2*tau+(2*z*z)**(1/(1-2*h))
    for _ in range(90):
        middle = (lower+upper)/2
        if middle-tau-z*z*middle**(2*h) > 0:
            upper = middle
        else:
            lower = middle
    return (lower+upper)/2


def q_jet(z, time, h=H):
    q = Jet(q_value(z.value, time.value, h))
    for _ in range(4):
        q = q-(q-(1-time)-z*z*q**(2*h))/(1-2*h*z*z*q**(2*h-1))
    return q


def chi(s):
    if s.value <= .5:
        return Jet(1)
    if s.value >= 1:
        return Jet(0)
    left, right = (-1/(s-.5)).exp(), (-1/(1-s)).exp()
    return right/(left+right)


def basis_theta(coords, cartesian):
    if cartesian:
        x, y = coords[:2]
        r = (x*x+y*y)**.5
        return [-y/r, x/r, Jet(0)]
    return [Jet(0), Jet(1), Jet(0)]


def representatives(coords, frequency=19, cartesian=False):
    if not waves.in_support(coords, cartesian):
        return [Jet(0) for _ in range(3)], [Jet(0) for _ in range(3)], Jet(0), Jet(0)
    data = [waves.carrier_data(coords, sign, cartesian=cartesian) for sign in (1, -1)]
    wave = [sum(-d['coefficient'][i]*(k*d['phase']).sin()/k
                for d, k in zip(data, (frequency, 2*frequency))) for i in range(3)]
    _, _, z, time, *rest = waves.geometry(coords, cartesian)
    envelope = data[0]['envelope']
    psi = .23*envelope*(1+.17*z.sin()+.11*time.cos())
    b = .19*envelope*(1+.13*z+.07*time)
    etheta = basis_theta(coords, cartesian)
    potential = [a+psi*e for a, e in zip(wave, etheta)]
    direct = [b*e for e in etheta]
    pressure = waves.wave_pressure(coords, frequency, cartesian)+.31*envelope*(1+.2*z)
    return potential, direct, pressure, psi


def cut_field(coords, scale, frequency=19, cartesian=False):
    A, B, p, psi = representatives(coords, frequency, cartesian)
    r = None if cartesian else coords[0]
    q = q_jet(coords[2], coords[3])
    cutoff = chi(scale*q)
    velocity = [a+cutoff*b for a, b in zip(waves.curl([cutoff*a for a in A], r, cartesian), B)]
    uncut = [a+b for a, b in zip(waves.curl(A, r, cartesian), B)]
    grad_cutoff = [cutoff.diff(j) for j in range(3)]
    if not cartesian:
        grad_cutoff[1] = grad_cutoff[1]/r
    commutator = cross(grad_cutoff, A)
    return dict(A=A, B=B, p=p, psi=psi, q=q, chi=cutoff,
                velocity=velocity, pressure=cutoff*p, uncut=uncut,
                grad_chi=grad_cutoff, commutator=commutator)


def cutoff_geometry_review():
    """Weighted derivatives agree as q shrinks at fixed eta, including mixed jets."""
    worst, rows, first = 0., [], None
    for q0 in (.2, .02, .002):
        eta, sigma = .43, .73
        z0, time0 = q0**(.5-H)*eta, 1-q0*(1-eta*eta)
        coords = variables((1.3, .27, z0, time0))
        q = q_jet(coords[2], coords[3])
        defect = q-(1-coords[3])-coords[2]**2*q**(2*H)
        cutoff = chi(sigma*q/q0)
        weighted = []
        for multi, value in zip(INDICES, cutoff.c):
            if multi[0] or multi[1]:
                continue
            a, b = multi[2:]
            derivative = value*math.factorial(a)*math.factorial(b)
            weighted.append(derivative*q0**((.5-H)*a+b))
        weighted = np.array(weighted)
        if first is None:
            first = weighted
        worst = max(worst, relative(first, weighted))
        rows.append(dict(q=q0, relative_q_value_error=abs(q.value/q0-1),
                         weighted_cutoff_jet_difference=relative(first, weighted),
                         implicit_equation_jet_error=float(np.linalg.norm(defect.c)/(1+np.linalg.norm(q.c)))))
    return dict(actual_implicit_q=True, derivatives_through_order=3, rows=rows,
                checks=dict(implicit_q_equation_in_all_jet_coefficients=max(r['implicit_equation_jet_error'] for r in rows) < 2e-12,
                            q_value_matches_similarity_inputs=max(r['relative_q_value_error'] for r in rows) < 2e-12,
                            scale_independent_weighted_cutoff_jets=worst < 2e-10))


def curl_review():
    error = divergence = cart_error = mean_error = 0.
    naive_gap = wrong_sign = commutator_size = 0.
    for point in POINTS:
        r, theta, z, time = point
        scale = .73/q_value(z, time)
        cc = variables(point)
        xx = variables((r*np.cos(theta), r*np.sin(theta), z, time))
        rotation = np.array([[np.cos(theta), -np.sin(theta), 0.],
                             [np.sin(theta), np.cos(theta), 0.], [0., 0., 1.]])
        for coords, cartesian in ((cc, False), (xx, True)):
            data = cut_field(coords, scale, cartesian=cartesian)
            exact = arrays(data['velocity'])
            expanded = arrays([data['chi']*v+c for v, c in zip(data['uncut'], data['commutator'])])
            error = max(error, *(relative(a, b) for a, b in zip(exact, expanded)))
            div = waves.divergence(exact, None if cartesian else r, cartesian)
            divergence = max(divergence, abs(div)/(1+np.linalg.norm(exact[1])))
            naive = arrays([data['chi']*v for v in data['uncut']])
            naive_gap = max(naive_gap, abs(waves.divergence(naive, None if cartesian else r, cartesian)))
            wrong_sign = max(wrong_sign, np.linalg.norm(exact[0]-arrays([
                data['chi']*v-c for v, c in zip(data['uncut'], data['commutator'])])[0]))
            commutator_size = max(commutator_size, np.linalg.norm(arrays(data['commutator'])[0]))
            if not cartesian:
                cyl = exact
                mean = [Jet(0), data['psi'], Jet(0)]
                full_mean = waves.curl([data['chi']*v for v in mean], cc[0])
                raw_mean = waves.curl(mean, cc[0])
                observed = arrays([a-data['chi']*b for a, b in zip(full_mean, raw_mean)])[0]
                expected = np.array([-data['chi'].diff(2).value*data['psi'].value, 0., 0.])
                mean_error = max(mean_error, relative(observed, expected))
            else:
                cart_error = max(cart_error, relative(exact[0], rotation@cyl[0]))
    plateau = []
    coords = variables(POINTS[1])
    for sigma in (.25, 1.25):
        data = cut_field(coords, sigma/q_value(POINTS[1][2], POINTS[1][3]))
        desired = data['uncut'] if sigma < .5 else [Jet(0) for _ in range(3)]
        plateau.append(max(relative(a, b) for a, b in zip(arrays(data['velocity']), arrays(desired))))
    support = []
    for point in ((.5, 0., .2, .47), (2.2, 0., .2, .47), (0., 0., .2, .47), (1.3, 0., 1.3, .47)):
        data = cut_field(variables(point), 1, cartesian=True)
        support.append(max(float(np.max(abs(a))) for a in arrays(data['velocity'])))
    return dict(manufactured_potentials=True, actual_q_cutoff=True,
                full_product_rule_error=error, divergence_error=divergence,
                cylindrical_cartesian_velocity_error=cart_error,
                mean_radial_commutator_error=mean_error,
                missing_curl_commutator_divergence=naive_gap,
                reversed_commutator_velocity_gap=wrong_sign,
                curl_commutator_size=float(commutator_size), plateau_errors=plateau,
                outside_support_errors=support,
                checks=dict(full_curl_product_rule_through_two_velocity_derivatives=error < 2e-10,
                            complete_cut_velocity_is_divergence_free=divergence < 2e-12,
                            cartesian_and_cylindrical_values_agree=cart_error < 2e-11,
                            mean_cutoff_adds_negative_axial_derivative_radial_velocity=mean_error < 2e-12,
                            naive_velocity_cutoff_breaks_divergence=naive_gap > 1e-4,
                            wrong_commutator_sign_is_detected=wrong_sign > 1e-4,
                            curl_commutator_is_nonzero=commutator_size > 1e-4,
                            both_plateaus_are_exact=max(plateau) < 1e-14,
                            core_axis_and_exterior_support_preserved=max(support) == 0.))


def residual_review():
    error = cart_error = time_gap = visc_gap = pressure_gap = 0.
    terms_max = {k: 0. for k in ('time', 'viscosity', 'pressure', 'old_new', 'new_old', 'new_square')}
    for point in POINTS:
        r, theta, z, time = point
        cc = variables(point)
        xx = variables((r*np.cos(theta), r*np.sin(theta), z, time))
        rotation = np.array([[np.cos(theta), -np.sin(theta), 0.],
                             [np.sin(theta), np.cos(theta), 0.], [0., 0., 1.]])
        for coords, cartesian in ((cc, False), (xx, True)):
            radius = None if cartesian else r
            first = cut_field(coords, .73/q_value(z, time), 11, cartesian)
            second = cut_field(coords, .81/q_value(z, time), 29, cartesian)
            e = [a+.37*b for a, b in zip(first['velocity'], second['velocity'])]
            pi = first['pressure']+.37*second['pressure']
            u, p = waves.background(coords, cartesian)
            U, E = arrays(u), arrays(e)
            combined = arrays([a+b for a, b in zip(u, e)])
            grad_pi = waves.pressure_gradient(pi, radius, cartesian)
            direct = waves.momentum(combined, waves.pressure_gradient(p+pi, radius, cartesian), radius, cartesian)-waves.momentum(U, waves.pressure_gradient(p, radius, cartesian), radius, cartesian)
            terms = dict(time=E[1][:, 3], viscosity=-waves.laplacian(E, radius, cartesian),
                         pressure=grad_pi, old_new=waves.advection(U, E, radius, cartesian),
                         new_old=waves.advection(E, U, radius, cartesian),
                         new_square=waves.advection(E, E, radius, cartesian))
            error = max(error, relative(direct, sum(terms.values())))
            for name, value in terms.items():
                terms_max[name] = max(terms_max[name], float(np.linalg.norm(value)))
            if not cartesian:
                cyl = direct
            else:
                cart_error = max(cart_error, relative(direct, rotation@cyl))
            ct, A = first['chi'].diff(3), first['A']
            grad_ct = [ct.diff(i) for i in range(3)]
            if not cartesian:
                grad_ct[1] = grad_ct[1]/coords[0]
            cutoff_time = arrays([ct*v+c for v, c in zip(first['uncut'], cross(grad_ct, A))])[0]
            time_gap = max(time_gap, float(np.linalg.norm(cutoff_time)))
            visc_gap = max(visc_gap, float(np.linalg.norm(waves.laplacian(arrays(first['velocity']), radius, cartesian)-first['chi'].value*waves.laplacian(arrays(first['uncut']), radius, cartesian))))
            pressure_gap = max(pressure_gap, float(np.linalg.norm(arrays([first['p']*v for v in first['grad_chi']])[0])))
    return dict(manufactured=True, two_distinct_cutoff_scales=True,
                full_residual_difference_error=error, cartesian_cylindrical_residual_error=cart_error,
                retained_term_sizes=terms_max, missing_cutoff_time_term_gap=time_gap,
                missing_cutoff_viscosity_commutator_gap=visc_gap,
                missing_pressure_cutoff_gradient_gap=pressure_gap,
                checks=dict(exact_complete_residual_difference=error < 2e-10,
                            independent_cartesian_residual_agrees=cart_error < 2e-10,
                            every_linear_and_nonlinear_term_is_nonzero=min(terms_max.values()) > 1e-5,
                            cutoff_time_terms_cannot_be_dropped=time_gap > 1e-4,
                            cutoff_viscosity_commutators_cannot_be_dropped=visc_gap > 1e-4,
                            pressure_cutoff_gradient_cannot_be_dropped=pressure_gap > 1e-4))


def divergent_series_review():
    # C_j=exp(j^2) grows faster than exp(s*h*j/10) for ANY fixed s.
    s, h = 1000., H
    j0 = 2
    while j0 <= h*s/10:
        j0 += 1
    log_growth = [j*j-h*j*s/10 for j in (j0, 2*j0, 20*j0)]
    return dict(synthetic_scalar_coefficients=True, coefficients='C_j=exp(j^2)',
                log_term_values=log_growth,
                all_fixed_q_argument='j^2-(h*s/10)*j tends to +infinity for each fixed s=-log q',
                checks=dict(arbitrarily_high_orders_do_not_imply_raw_series_convergence=log_growth[-1] > log_growth[0] > 0,
                            no_substitution_of_j_of_q_into_uncontrolled_constants=True,
                            flat_stage_errors_are_compared_not_independently_summed=True))
