"""Independent coordinate checks and deliberate failures for the wave interface."""
import math
import numpy as np
from scipy.integrate import quad_vec
from jets import Jet, arrays, dot, variables
from waves import (AC, NU, advection, background, carrier_data, covariance,
                   divergence, field, geometry, harmonic_samples, laplacian, momentum,
                   pressure_gradient, tensor_divergence, tensor_values, wave_pressure)

POINTS = [(1.04, -.63, -.34, .31), (1.3, .27, .2, .47),
          (1.61, 1.17, -.12, .62), (1.79, -2.1, .38, .71)]
FREQUENCIES = [8, 16, 32, 64]


def relative(a, b):
    return float(np.linalg.norm(a-b)/(1+np.linalg.norm(a)+np.linalg.norm(b)))


def add_data(a, b):
    return tuple(x+y for x, y in zip(a, b))


def coordinate_review():
    worst_curl = worst_residual = worst_div = worst_normal = 0.
    periodic = 0.
    omitted_curl = omitted_connection = omitted_advection = omitted_pressure = 0.
    omitted_radial_curl = 0.
    for p in POINTS:
        r, theta, z, time = p
        coords = variables(p)
        xcoords = variables((r*math.cos(theta), r*math.sin(theta), z, time))
        rot = np.array([[math.cos(theta), -math.sin(theta), 0.],
                        [math.sin(theta), math.cos(theta), 0.], [0., 0., 1.]])
        for k in (8, 23):
            wcyl = field(coords, k)
            wpot = field(coords, k, via_potential=True)
            wcart = field(xcoords, k, cartesian=True, via_potential=True)
            cart_direct = field(xcoords, k, cartesian=True)
            # Compare cylindrical and Cartesian derivatives by evaluating
            # the cylindrical construction in the Cartesian chart as well.
            for aa, bb in ((arrays(wcyl), arrays(wpot)), (arrays(wcart), arrays(cart_direct))):
                for a, b in zip(aa, bb):
                    worst_curl = max(worst_curl, relative(a, b))
            worst_curl = max(worst_curl, relative(arrays(wcart)[0], rot@arrays(wcyl)[0]))
            wdata, wxdata = arrays(wcyl), arrays(wcart)
            worst_div = max(worst_div, abs(divergence(wdata, r))/(1+k*np.linalg.norm(wdata[0])),
                            abs(divergence(wxdata, cartesian=True))/(1+k*np.linalg.norm(wxdata[0])))
            vp, pressure = background(coords)
            vx, px = background(xcoords, True)
            ucyl = arrays([a+b for a, b in zip(vp, wcyl)])
            ucart = arrays([a+b for a, b in zip(vx, wcart)])
            pg = pressure_gradient(pressure+wave_pressure(coords, k), r)
            pxg = pressure_gradient(px+wave_pressure(xcoords, k, True), cartesian=True)
            rcyl = momentum(ucyl, pg, r)
            rcart = momentum(ucart, pxg, cartesian=True)
            worst_residual = max(worst_residual, relative(rot@rcyl, rcart))
            shifted = variables((r, theta+2*math.pi, z, time))
            for a, b in zip(arrays(wcyl), arrays(field(shifted, k))):
                periodic = max(periodic, relative(a, b))
            principal = [Jet(0)]*3
            missing_radial_curl = Jet(0)
            for sign, freq in ((1, k), (-1, 2*k)):
                d = carrier_data(coords, sign)
                worst_normal = max(worst_normal, abs(dot(d['normal'], d['amplitude']).value))
                principal = [a+b*(freq*d['phase']).cos() for a, b in zip(principal, d['amplitude'])]
                missing_radial_curl = missing_radial_curl+d['coefficient'][1]*(freq*d['phase']).sin()/(r*freq)
            omitted_curl = max(omitted_curl, abs(divergence(arrays(principal), r)))
            bad = [wcyl[0], wcyl[1], wcyl[2]+missing_radial_curl]
            omitted_radial_curl = max(omitted_radial_curl, abs(divergence(arrays(bad), r)))
            scalar_lap = (ucyl[2][:, 0, 0]+ucyl[1][:, 0]/r
                          +ucyl[2][:, 1, 1]/r**2+ucyl[2][:, 2, 2])
            omitted_connection = max(omitted_connection, float(np.linalg.norm(laplacian(ucyl, r)-scalar_lap)))
            naive_adv = (ucyl[0][0]*ucyl[1][:, 0]+ucyl[0][1]*ucyl[1][:, 1]/r
                         +ucyl[0][2]*ucyl[1][:, 2])
            omitted_advection = max(omitted_advection, float(np.linalg.norm(advection(ucyl, ucyl, r)-naive_adv)))
            omitted_pressure = max(omitted_pressure, float(np.linalg.norm(pg)))
    outside = [(0.7, .1, .2, .5), (2.1, .2, .1, .4), (1.4, .2, 1.4, .5),
               (1.4, .2, .1, 1.1), (.85, .2, .1, .5), (0., .2, .1, .5)]
    support_error = max(float(np.linalg.norm(a)) for p in outside for a in arrays(field(variables(p), 11)))
    support_error = max(support_error, max(float(np.linalg.norm(a)) for a in
                       arrays(field(variables((0., 0., .1, .5)), 11, cartesian=True))))
    return dict(samples=len(POINTS)*2,
                maximum_curl_value_gradient_hessian_relative_error=worst_curl,
                maximum_cartesian_cylindrical_full_residual_relative_error=worst_residual,
                maximum_scaled_divergence_error=float(worst_div), maximum_normal_amplitude_error=worst_normal,
                maximum_periodic_derivative_relative_error=periodic, outside_support_jet_error=support_error,
                principal_only_divergence_gap=omitted_curl,
                radial_curl_connection_omission_gap=omitted_radial_curl,
                scalar_laplacian_omission_gap=omitted_connection,
                advection_connection_omission_gap=omitted_advection, pressure_omission_gap=omitted_pressure,
                checks=dict(exact_curl_through_required_derivatives=worst_curl < 2e-10,
                            full_physical_residual_independent_coordinates=worst_residual < 2e-10,
                            exact_divergence=worst_div < 2e-10,
                            transverse_principal_amplitude=worst_normal < 2e-11,
                            global_angular_periodicity=periodic < 2e-10,
                            smooth_compact_support_fixture=support_error == 0.,
                            curl_remainder_cannot_be_discarded=omitted_curl > 1e-3,
                            radial_curl_connection_cannot_be_discarded=omitted_radial_curl > 1e-3,
                            vector_laplacian_connections_cannot_be_discarded=omitted_connection > 1e-3,
                            advection_connections_cannot_be_discarded=omitted_advection > 1e-3,
                            pressure_cannot_be_discarded=omitted_pressure > 1e-3))


def covariance_review():
    average_error = quadratic_error = mean_error = signed_error = 0.
    smallest_psd = math.inf
    cross_mean_gap = nonlinear_remainder_gap = 0.
    for p in POINTS:
        coords = variables(p)
        reference = None
        for k in FREQUENCIES:
            w, _ = harmonic_samples(p, k)
            q, q0 = covariance(coords, k)
            qv, q0v = tensor_values(q), tensor_values(q0)
            empirical = w[0].T@w[0]/len(w[0])
            average_error = max(average_error, relative(qv, empirical))
            scaled = (qv-q0v)*k*k
            if reference is None:
                reference = scaled
            quadratic_error = max(quadratic_error, relative(scaled, reference))
            smallest_psd = min(smallest_psd, float(np.linalg.eigvalsh(scaled).min()))
            mean_error = max(mean_error, relative(np.mean(advection(w, w, p[0]), axis=0), tensor_divergence(q, p[0])))
        k = 16
        w, _ = harmonic_samples(p, k)
        wp, _ = harmonic_samples(p, k, principal_only=True)
        v, _ = harmonic_samples(p, k, mode='increment')
        lp, _ = harmonic_samples(p, k, mode='increment', principal_only=True)
        e, rem = w[0]-wp[0], v[0]-lp[0]
        def c(a):
            return np.mean(a[:, 0, None]*a[:, 1:], axis=0)
        def b(a, z):
            return np.mean(a[:, 0, None]*z[:, 1:]+z[:, 0, None]*a[:, 1:], axis=0)
        geo = geometry(coords)
        z, time, f, kt, nt = geo[2], geo[3], geo[8], geo[5], geo[6]
        target = np.array([(f*f*((.65+.1*z.sin())*nt[j]
                                 +(1.1+.1*time.cos())*kt[j])).value for j in range(2)])
        signed_error = max(signed_error, relative(b(wp[0], lp[0]), target),
                           relative(c(w[0]+v[0])-c(w[0]), target+b(e, lp[0])+b(w[0], rem)+c(v[0])))
        nonlinear_remainder_gap = max(nonlinear_remainder_gap, float(np.linalg.norm(c(v[0]))))
        same, _ = harmonic_samples(p, k, principal_only=True, same_frequency=True)
        q0v = tensor_values(covariance(coords, k)[1])
        cross_mean_gap = max(cross_mean_gap, float(np.linalg.norm(same[0].T@same[0]/len(same[0])-q0v)))
    return dict(maximum_angular_covariance_relative_error=average_error,
                maximum_k_squared_covariance_remainder_relative_error=quadratic_error,
                minimum_k_squared_remainder_eigenvalue=smallest_psd,
                maximum_full_tensor_divergence_relative_error=mean_error,
                maximum_signed_full_increment_relative_error=signed_error,
                unsupported_same_frequency_mean_cross_gap=cross_mean_gap,
                omitted_signed_quadratic_increment_gap=nonlinear_remainder_gap,
                checks=dict(full_six_entry_covariance=average_error < 2e-10,
                            real_curl_covariance_error_is_quadratic=quadratic_error < 2e-10,
                            curl_covariance_error_positive_semidefinite=smallest_psd > -2e-10,
                            full_tensor_mean_divergence=mean_error < 2e-10,
                            signed_inverse_and_all_covariance_increment_terms=signed_error < 2e-10,
                            angular_mean_cross_terms_require_a_reason=cross_mean_gap > 1e-3,
                            signed_linear_inverse_does_not_remove_quadratic_remainder=nonlinear_remainder_gap > 1e-3))


def residual_review():
    p = POINTS[1]
    coords = variables(p)
    base, bp = background(coords)
    bd = arrays(base)
    bg = pressure_gradient(bp, p[0])
    rows = []
    decomposition_error = mean_error = 0.
    for k in FREQUENCIES:
        w, pg = harmonic_samples(p, k)
        linear = w[1][:, :, 3]+advection(bd, w, p[0])+advection(w, bd, p[0])-NU*laplacian(w, p[0])+pg
        nonlinear = advection(w, w, p[0])
        total = momentum(add_data(bd, w), bg+pg, p[0])-momentum(bd, bg, p[0])
        decomposition_error = max(decomposition_error, relative(total, linear+nonlinear))
        q, q0 = covariance(coords, k)
        mean_error = max(mean_error, relative(np.mean(total, axis=0), tensor_divergence(q, p[0])))
        rms = lambda a: float(np.sqrt(np.mean(np.sum(a*a, axis=-1))))
        rows.append(dict(k=k, full_increment_rms=rms(total), linear_increment_rms=rms(linear),
                         viscous_rms=rms(NU*laplacian(w, p[0])), nonlinear_rms=rms(nonlinear),
                         full_increment_rms_divided_by_k_squared=rms(total)/k**2,
                         mean_increment=np.mean(total, axis=0).tolist(),
                         principal_cross_covariance=[q0[0][j].value for j in (1, 2)]))
    limiting_square = 0.
    for sign, multiple in ((1, 1), (-1, 2)):
        data = carrier_data(coords, sign)
        n2 = dot(data['normal'], data['normal']).value
        a2 = dot(data['amplitude'], data['amplitude']).value
        limiting_square += .5*NU**2*multiple**4*n2**2*a2
    limit = math.sqrt(limiting_square)
    return dict(rows=rows, maximum_full_residual_decomposition_relative_error=decomposition_error,
                maximum_full_mean_residual_relative_error=mean_error,
                exact_leading_viscous_rms_limit_divided_by_k_squared=limit,
                checks=dict(full_residual_includes_all_wave_background_and_wave_wave_terms=decomposition_error < 2e-10,
                            full_mean_residual_is_full_tensor_divergence=mean_error < 2e-10,
                            kinematic_covariance_matching_does_not_reduce_residual=rows[-1]['full_increment_rms'] > 30*rows[0]['full_increment_rms'],
                            uncancelled_viscosity_has_k_squared_scale=all(.1 < r['full_increment_rms_divided_by_k_squared'] < 100 for r in rows),
                            finite_frequency_values_approach_derived_limit=abs(rows[-1]['full_increment_rms_divided_by_k_squared']-limit) < .02*limit),
                actual_background_or_pulse_evaluated=False,
                full_actual_residual_improvement_verified=False)


def mean_defect(q, r):
    """What survives after only Q_rtheta=S_theta, Q_rz=S_z are matched."""
    return [q[0][0].diff(0)+(q[0][0]-q[1][1])/r+q[0][2].diff(2),
            q[1][2].diff(2), q[2][2].diff(2)]


def mean_pressure_review():
    r, theta, z, time = POINTS[1]
    def radial_integrand(radius):
        coords = variables((radius, 0., z, time))
        q0 = covariance(coords, 16)[1]
        dr = mean_defect(q0, coords[0])[0]
        return np.array([dr.value, dr.diff(2).value])
    primitive, error = quad_vec(radial_integrand, r, 1.95, epsabs=2e-11, epsrel=2e-11)
    inner, inner_error = quad_vec(radial_integrand, .8, 1.95, epsabs=2e-11, epsrel=2e-11)
    q0 = covariance(variables((r, theta, z, time)), 16)[1]
    defect = np.array([x.value for x in mean_defect(q0, Jet.variable(r, 0))])
    after_pressure = np.array([0., defect[1], defect[2]+primitive[1]])
    return dict(principal_mean_defect_before_pressure=defect.tolist(),
                radial_pressure_primitive=float(primitive[0]), pressure_z=float(primitive[1]),
                mean_defect_after_radial_pressure=after_pressure.tolist(),
                inner_pressure_constant=float(inner[0]), inner_pressure_z=float(inner[1]),
                quadrature_error_estimates=[float(error), float(inner_error)],
                checks=dict(tangential_flux_matching_leaves_full_mean_terms=np.linalg.norm(defect) > 1e-3,
                            radial_pressure_leaves_axial_defect=abs(after_pressure[2]) > 1e-3,
                            pressure_primitive_need_not_vanish_inside_wave_support=abs(inner[0]) > 1e-3,
                            theta_z_flux_cannot_be_discarded=abs(after_pressure[1]) > 1e-3),
                actual_mean_corrector_constructed=False)
