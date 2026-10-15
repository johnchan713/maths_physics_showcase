"""Independent local identities for complete curl velocities and interactions.

The compact field is manufactured, not a simulation of the actual global
band construction. Coordinate identities test algebra; the conditional
all-scale estimates are proved separately in README.md.
"""
from pathlib import Path
import sys
import numpy as np

PARENT = Path(__file__).resolve().parent.parent
for folder in ('physical_wave_interface', 'compact_mean_corrections'):
    if str(PARENT/folder) not in sys.path:
        sys.path.append(str(PARENT/folder))
from jets import Jet, arrays, variables
from mean_system import mean_correction, mean_increment
import waves


def relative(a, b):
    return float(np.linalg.norm(a-b)/(1+np.linalg.norm(a)+np.linalg.norm(b)))


def coefficient(coords, sign, frequency, kind='primary', cartesian=False):
    """a=t+i*curl(n cross t/|n|^2)/frequency for a complex curl harmonic."""
    data = waves.carrier_data(coords, sign, kind, cartesian)
    t, r = arrays(data['amplitude']), arrays(data['correction'])
    return tuple(a+1j*b/frequency for a, b in zip(t, r)), data


def complete_harmonic(data, phase, frequency):
    a, grad, _ = data
    phi_gradient = np.array([phase.diff(i).value for i in range(4)])
    exponential = np.exp(1j*frequency*phase.value)
    return a*exponential, (grad+1j*frequency*a[:, None]*phi_gradient)*exponential


def expanded_transport(a, b, normal, radius, target_frequency):
    """Coefficient of (a e^ikmPhi . grad)(b e^ikm'Phi) in cylindrical form."""
    av, bv = a[0], b[0]
    frame = np.array([-av[1]*bv[1]/radius, av[1]*bv[0]/radius, 0j])
    return (av[0]*b[1][:, 0]+av[2]*b[1][:, 2]+frame
            +1j*target_frequency*np.dot(av, normal)*bv)


def longitudinal_and_transport_review():
    longitudinal = cartesian_error = incompressibility = frame_gap = 0.
    rows, principal_divergence_gap = [], 0.
    for point in ((1.12, -.4, -.19, .36), (1.3, .27, .2, .47), (1.66, -.9, .31, .64)):
        radius, theta, z, time = point
        cc = variables(point)
        xx = variables((radius*np.cos(theta), radius*np.sin(theta), z, time))
        rotation = np.array([[np.cos(theta), -np.sin(theta), 0.],
                             [np.sin(theta), np.cos(theta), 0.], [0., 0., 1.]])
        for m, other in ((1, 2), (2, -1), (1, -1)):
            k = 19
            a, ad = coefficient(cc, 1, k*m)
            b, bd = coefficient(cc, 1, k*other, 'increment')
            normal = np.array([n.value for n in ad['normal']])
            div_a = a[1][0, 0]+a[0][0]/radius+a[1][2, 2]
            longitudinal = max(longitudinal, relative(np.dot(a[0], normal), 1j*div_a/(k*m)))
            principal = arrays(ad['amplitude'])
            principal_divergence_gap = max(principal_divergence_gap,
                                           abs(principal[1][0, 0]+principal[0][0]/radius+principal[1][2, 2]))
            expanded = (expanded_transport(a, b, normal, radius, k*other)
                        +expanded_transport(b, a, normal, radius, k*m))
            ax, axd = coefficient(xx, 1, k*m, cartesian=True)
            bx, bxd = coefficient(xx, 1, k*other, 'increment', cartesian=True)
            axv, axg = complete_harmonic(ax, axd['phase'], k*m)
            bxv, bxg = complete_harmonic(bx, bxd['phase'], k*other)
            direct_cartesian = bxg[:, :3]@axv+axg[:, :3]@bxv
            physical_expanded = rotation@expanded*np.exp(1j*k*(m+other)*ad['phase'].value)
            cartesian_error = max(cartesian_error, relative(direct_cartesian, physical_expanded))
            incompressibility = max(incompressibility, abs(np.trace(axg[:, :3]))/(1+k*np.linalg.norm(axv)))
            frame = np.array([-2*a[0][1]*b[0][1]/radius,
                               (a[0][1]*b[0][0]+b[0][1]*a[0][0])/radius, 0j])
            frame_gap = max(frame_gap, np.linalg.norm(frame))
    # Complete same-label amplitudes have bounded nonlinear coefficients as
    # k grows; an explicitly longitudinal, non-curl perturbation loses that.
    cc = variables((1.3, .27, .2, .47))
    for k in (16, 64, 256, 1024):
        a, ad = coefficient(cc, 1, k)
        b, _ = coefficient(cc, 1, 2*k, 'increment')
        n = np.array([x.value for x in ad['normal']])
        exact = expanded_transport(a, b, n, 1.3, 2*k)+expanded_transport(b, a, n, 1.3, k)
        bad = (a[0]+.1*n, a[1], a[2])
        wrong = expanded_transport(bad, b, n, 1.3, 2*k)+expanded_transport(b, bad, n, 1.3, k)
        rows.append(dict(k=k, full_curl_transport_coefficient=float(np.linalg.norm(exact)),
                         longitudinal_perturbation_coefficient=float(np.linalg.norm(wrong))))
    return dict(maximum_exact_longitudinal_identity_relative_error=longitudinal,
                maximum_cartesian_cylindrical_transport_relative_error=cartesian_error,
                maximum_scaled_complex_curl_divergence=float(incompressibility),
                principal_only_divergence_gap=float(principal_divergence_gap),
                omitted_cylindrical_connection_gap=float(frame_gap), frequency_rows=rows,
                checks=dict(full_curl_longitudinal_identity=longitudinal < 2e-10,
                            all_harmonic_signs_and_zero_output_checked=cartesian_error < 2e-10,
                            direct_cartesian_transport_agrees=cartesian_error < 2e-10,
                            exact_complex_curl_incompressibility=incompressibility < 2e-10,
                            dropping_curl_breaks_divergence=principal_divergence_gap > 1e-3,
                            cylindrical_connection_is_present=frame_gap > 1e-3,
                            same_label_transport_has_no_half_power_loss=rows[-1]['full_curl_transport_coefficient'] < 3*rows[0]['full_curl_transport_coefficient'],
                            without_longitudinal_control_fast_loss_returns=rows[-1]['longitudinal_perturbation_coefficient'] > 20*rows[0]['longitudinal_perturbation_coefficient']))


def full_increment_review():
    worst = nonlinear = 0.
    omitted_self = omitted_pressure = omitted_viscosity = omitted_old_new = 0.
    for point in ((1.12, -.4, -.19, .36), (1.3, .27, .2, .47), (1.66, -.9, .31, .64)):
        radius, theta, z, time = point
        coords = variables(point)
        base, pbase = waves.background(coords)
        old_mean = mean_correction(coords)
        delta_mean = mean_increment(coords)
        old_wave = waves.field(coords, 13)
        delta_wave = waves.field(coords, 13, 'increment')
        old = arrays([a+b+c for a, b, c in zip(base, old_mean, old_wave)])
        delta = arrays([a+b for a, b in zip(delta_mean, delta_wave)])
        new = tuple(a+b for a, b in zip(old, delta))
        pressure = pbase+waves.wave_pressure(coords, 13)
        dp = .17*waves.geometry(coords)[8]+.09*coords[0]*coords[2]
        pgrad, dpgrad = waves.pressure_gradient(pressure, radius), waves.pressure_gradient(dp, radius)
        actual = waves.momentum(new, pgrad+dpgrad, radius)-waves.momentum(old, pgrad, radius)
        old_new = waves.advection(old, delta, radius)+waves.advection(delta, old, radius)
        square = waves.advection(delta, delta, radius)
        viscous = -waves.laplacian(delta, radius)
        expanded = delta[1][:, 3]+viscous+dpgrad+old_new+square
        worst = max(worst, relative(actual, expanded))
        omitted_self = max(omitted_self, np.linalg.norm(square))
        omitted_pressure = max(omitted_pressure, np.linalg.norm(dpgrad))
        omitted_viscosity = max(omitted_viscosity, np.linalg.norm(viscous))
        omitted_old_new = max(omitted_old_new, np.linalg.norm(old_new))
        # This second check is the exact nonlinear increment, in the same
        # cylindrical coordinates. The separate transport review above uses
        # independently differentiated Cartesian curl fields.
        nonlinear = max(nonlinear, relative(waves.advection(new, new, radius)-waves.advection(old, old, radius), old_new+square))
    return dict(maximum_same_field_full_residual_increment_relative_error=worst,
                maximum_exact_nonlinear_increment_relative_error=nonlinear,
                omitted_new_velocity_square_gap=float(omitted_self),
                omitted_pressure_gap=float(omitted_pressure), omitted_viscosity_gap=float(omitted_viscosity),
                omitted_old_new_cross_gap=float(omitted_old_new),
                checks=dict(full_actual_velocity_increment_identity=worst < 2e-10,
                            nonlinear_increment_has_all_three_products=nonlinear < 2e-10,
                            new_velocity_square_is_not_discarded=omitted_self > 1e-3,
                            pressure_is_not_discarded=omitted_pressure > 1e-3,
                            viscosity_is_not_discarded=omitted_viscosity > 1e-3,
                            old_new_cross_terms_are_not_discarded=omitted_old_new > 1e-3))


def covariance_increment_review():
    """Signed primary inverse with complete old/new curls and all remainders."""
    identity = prescribed = 0.
    old_remainder_gap = curl_remainder_gap = square_gap = 0.
    def covariance(a):
        return np.mean(a[:, 0, None]*a[:, 1:], axis=0)
    def bilinear(a, b):
        return np.mean(a[:, 0, None]*b[:, 1:]+b[:, 0, None]*a[:, 1:], axis=0)
    for point in ((1.12, -.4, -.19, .36), (1.3, .27, .2, .47), (1.66, -.9, .31, .64)):
        coords = variables(point)
        for k in (11, 29):
            primary = waves.harmonic_samples(point, k, principal_only=True)[0][0]
            exact = waves.harmonic_samples(point, k)[0][0]
            transverse = waves.harmonic_samples(point, k, mode='increment', principal_only=True)[0][0]
            increment = waves.harmonic_samples(point, k, mode='increment')[0][0]
            old = exact+.17*increment
            old_remainder, curl_remainder = old-primary, increment-transverse
            geometry = waves.geometry(coords)
            sn, sk = .65+.1*coords[2].sin(), 1.1+.1*coords[3].cos()
            target = np.array([(geometry[8]**2*(sn*geometry[6][i]+sk*geometry[5][i])).value for i in range(2)])
            leading = bilinear(primary, transverse)
            prescribed = max(prescribed, relative(leading, target))
            a, b, c = bilinear(old_remainder, transverse), bilinear(old, curl_remainder), covariance(increment)
            actual = covariance(old+increment)-covariance(old)
            expected = target+a+b+c
            identity = max(identity, relative(actual, expected))
            old_remainder_gap = max(old_remainder_gap, np.linalg.norm(a))
            curl_remainder_gap = max(curl_remainder_gap, np.linalg.norm(b))
            square_gap = max(square_gap, np.linalg.norm(c))
    return dict(maximum_prescribed_signed_cross_covariance_relative_error=prescribed,
                maximum_complete_covariance_increment_relative_error=identity,
                old_remainder_cross_gap=float(old_remainder_gap), new_curl_cross_gap=float(curl_remainder_gap),
                signed_increment_square_gap=float(square_gap),
                checks=dict(signed_primary_linear_inverse_matches_target=prescribed < 2e-10,
                            complete_covariance_change_has_all_remainders=identity < 2e-10,
                            old_wave_remainder_is_retained=old_remainder_gap > 1e-4,
                            new_curl_remainder_is_retained=curl_remainder_gap > 1e-4,
                            signed_increment_square_is_retained=square_gap > 1e-3))
