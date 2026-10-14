"""Checks of the conditional estimate's algebra; all fluid data are fixtures."""
from fractions import Fraction
import math
import mpmath as mp
import numpy as np
from scipy.integrate import quad_vec, solve_ivp
from bounds import (barrier_certificate, cone_certificate, derivative_exponents,
                    endpoint_flat_threshold, eventual_band_threshold, exponent_ledger,
                    phase_majorant)
from charts import (H, KAPPA, REP, U, Jet, background, background_jet, band, frame,
                    label, local_coefficients, nearest_nonzero_integer, qstar,
                    rational, variables)


def fixture_first_band():
    n = eventual_band_threshold(rational(H), 3*10**6)
    return 400*((n+399)//400)


def point_for(data):
    return [mp.mpf(v)+mp.mpf(a)/data['scale']**3 for v, a in zip(REP, ('.3', '-.2', '.1'))]


def phase_identity_review():
    errors, omitted_fast_time, omitted_axial_phase = [], [], []
    eps, h = .03, .005
    p, pz, x0 = .37, -.43, .21
    for values in ((1.4, .12, .62, 3.), (1.31, -.16, .7, 11.), (1.6, .23, .55, 17.)):
        r, z, time, v = variables(values)
        b, f, g = background_jet(r, z, time, eps, h)
        transport = p*f+pz*g
        phi = pz*z/eps+x0*r-v*transport
        normal = np.array([phi.diff(0).value, p/r.value, eps*phi.diff(1).value])
        explicit = np.array([x0-v.value*transport.diff(0).value, p/r.value,
                             pz-eps*v.value*transport.diff(1).value])
        actual_eik = (-eps*phi.diff(2)+phi.diff(3)+b*phi.diff(0)+p*f+eps*g*phi.diff(1)).value
        explicit_eik = (eps*v*(transport.diff(2)-g*transport.diff(1))+b*phi.diff(0)).value
        errors.extend([float(np.linalg.norm(normal-explicit)), abs(actual_eik-explicit_eik)])
        omitted_fast_time.append(abs(actual_eik-(-eps*phi.diff(2)+b*phi.diff(0)+p*f+eps*g*phi.diff(1)).value))
        omitted_axial_phase.append(abs(normal[2]-pz))
    # Profile/chart normalization is independently checked by differentiating
    # the actual implicit coordinate map, before its simplified formula.
    normalization_error = 0.
    for values in ((1.4, .12, .62, 0.), (1.7, -.22, .71, 0.)):
        r, z, time, _ = variables(values)
        from charts import qstar_jet
        q = qstar_jet(z, time, h)
        a, d = .5+h, .5-h
        rp, eta = r/q**.5, z/q**d
        fp = (1+.1*eta+.03*eta*eta)/rp**3
        up = -1.5*(1+.1*eta+.03*eta*eta)/rp**2
        chart_f, chart_g = q**(-a-.5)*fp, q**(-a)*up
        simplified = background_jet(r, z, time, 0., h)[1:]
        for x, y in zip((chart_f, chart_g), simplified):
            normalization_error = max(normalization_error, float(np.max(np.abs(x.c-y.c))))
    return dict(maximum_phase_gradient_and_transport_identity_error=max(errors),
                maximum_chart_normalization_third_jet_error=normalization_error,
                omitted_fast_time_gap=max(omitted_fast_time), omitted_axial_phase_gap=max(omitted_axial_phase),
                checks=dict(actual_phase_formula_chain_rule=max(errors) < 2e-11,
                            implicit_profile_to_chart_derivatives=normalization_error < 2e-10,
                            fast_time_term_required=max(omitted_fast_time) > 1e-3,
                            axial_phase_derivative_required=max(omitted_axial_phase) > 1e-3))


def band_frame_review():
    rows = []
    max_frame_identity = max_basis_derivative = max_rounding = mp.mpf(0)
    omitted_basis = mp.mpf(0)
    first = fixture_first_band()
    for ell in (first, first+4000, first+12000):
        data = band(ell)
        for sign in (1, -1):
            lab = label(data, sign)
            coefficients = local_coefficients(data, lab, point_for(data))
            rounding = abs(data['k']*(lab['p']-lab['ptilde']))
            max_rounding = max(max_rounding, rounding)
            normal_error = frame_error = damping_error = normal_prime = mp.mpf(0)
            for fraction in (0, mp.mpf('.25'), mp.mpf('.5'), mp.mpf('.75'), 1):
                v = fraction*data['length']
                fr = frame(data, lab, coefficients, v)
                ident = fr['projected']*fr['basis']-fr['basis_prime']-fr['basis']*(mp.diag([fr['growth'], -fr['growth']])+fr['error'])
                max_frame_identity = max(max_frame_identity, mp.norm(ident), mp.norm(fr['left']*fr['basis']-mp.eye(2)),
                                         mp.norm(fr['normal'].T*fr['basis']))
                normal_error = max(normal_error, data['scale']*mp.norm(fr['normal']-fr['reference_normal']))
                frame_error = max(frame_error, data['scale']*max(abs(x) for x in fr['error']))
                damping_error = max(damping_error, data['scale']*abs(fr['damping_error']))
                normal_prime = max(normal_prime, data['scale']*mp.norm(fr['normal_prime']))
                derivative = mp.matrix(3, 2)
                for i in range(3):
                    for j in range(2):
                        derivative[i, j] = mp.diff(lambda vv: frame(data, lab, coefficients, vv)['basis'][i, j], v)
                max_basis_derivative = max(max_basis_derivative, mp.norm(derivative-fr['basis_prime']))
                omitted_basis = max(omitted_basis, data['scale']*mp.norm(fr['left']*fr['basis_prime']))
            rows.append(dict(ell=ell, sign=sign, S_squared_smallness_majorant=phase_majorant(ell, data['h']),
                             epsilon_times_k_squared=data['epsilon']*data['k']**2,
                             ci_times_S=data['ci']*data['scale'], length_over_S=data['length']/data['scale'],
                             radial_evaluation_scale_ratio=data['mi']/(data['epsilon']**(-rational(KAPPA))*data['scale']**(-data['rho'])),
                             scaled_normal_error=normal_error, scaled_frame_entry_error=frame_error,
                             scaled_damping_error=damping_error, scaled_normal_derivative=normal_prime,
                             angular_rounding_error_in_integer_units=rounding,
                             angular_frequency_nonzero=lab['angular_integer'] != 0))
    zero_guard = nearest_nonzero_integer(mp.mpf(0))
    half_bound_counterexample = abs(nearest_nonzero_integer(mp.mpf('.1'))-mp.mpf('.1'))
    return dict(rows=rows, maximum_exact_frame_identity_error=max_frame_identity,
                maximum_analytic_basis_prime_error=max_basis_derivative,
                maximum_angular_rounding_error=max_rounding,
                omitted_basis_prime_scaled_gap=omitted_basis,
                nearest_nonzero_integer_at_zero=zero_guard,
                false_half_integer_rounding_bound_gap=half_bound_counterexample-mp.mpf('.5'),
                checks=dict(full_moving_frame_identity=max_frame_identity < mp.mpf('1e-85'),
                            analytic_moving_basis_derivative=max_basis_derivative < mp.mpf('1e-85'),
                            correct_rounding_error_is_at_most_one=max_rounding <= 1,
                            angular_periodicity_nonzero=all(r['angular_frequency_nonzero'] for r in rows),
                            accepted_fixture_bands_are_small=all(r['S_squared_smallness_majorant'] < mp.mpf('1e-6') for r in rows),
                            phase_and_frame_fixture_have_inverse_S_error=all(max(r['scaled_normal_error'], r['scaled_frame_entry_error'], r['scaled_damping_error'], r['scaled_normal_derivative']) < 10000 for r in rows),
                            cover_and_carrier_scales=all(0 < r['ci_times_S'] <= 1 and 1 <= r['epsilon_times_k_squared'] <= 4 for r in rows),
                            moving_basis_derivative_required=omitted_basis > mp.mpf('.01'),
                            nearest_nonzero_guard_required=zero_guard != 0,
                            universal_half_integer_rounding_bound_rejected=half_bound_counterexample > mp.mpf('.5')),
                actual_global_background_numerically_evaluated=False)


def pulse_covariance_review(tighter=False):
    """Actual phase formulas on a huge fixture band, solved without P underflow.

    Ratio/log-amplitude evolution removes the enormous stable eigenvalue
    from the growing amplitude. Radau handles the remaining initial layer.
    Gaussian coordinates resolve the covariance peak instead of uniform v.
    """
    data = band(fixture_first_band())
    length, scale = float(data['length']), data['scale']
    rows, columns = [], []
    maximum_ratio = maximum_log_amplitude = 0.
    for sign in (1, -1):
        lab = label(data, sign)
        coeff = local_coefficients(data, lab, point_for(data))
        def matrices(u):
            fr = frame(data, lab, coeff, mp.mpf(float(u))*data['length'])
            return np.array(fr['error'].tolist(), dtype=float), float(fr['growth']), float(fr['damping_error'])
        def rhs(u, y):
            e, growth, delta = matrices(u)
            ratio = y[0]
            return length*np.array([e[1, 0]+(-2*growth+e[1, 1]-e[0, 0])*ratio-e[0, 1]*ratio**2,
                                    -delta+e[0, 0]+e[0, 1]*ratio])
        def jac(u, y):
            e, growth, delta = matrices(u)
            return length*np.array([[-2*growth+e[1, 1]-e[0, 0]-2*e[0, 1]*y[0], 0.], [e[0, 1], 0.]])
        sol = solve_ivp(rhs, (0., 1.), [0., 0.], method='Radau', jac=jac, dense_output=True,
                        rtol=2e-10 if tighter else 2e-9, atol=[2e-16 if tighter else 2e-15, 2e-12 if tighter else 2e-11],
                        max_step=.025)
        if not sol.success:
            raise ValueError('Normalized actual-phase fixture ODE failed')
        mu = lab['lambda0']/mp.sqrt(1+(mp.mpf('1.5')*U)**2)
        certificate = barrier_certificate(mu, 10000, 10000, scale, data['length']/scale)
        values = sol.sol(np.linspace(0., 1., 301))
        maximum_ratio = max(maximum_ratio, float(np.max(np.abs(values[0]))))
        maximum_log_amplitude = max(maximum_log_amplitude, float(np.max(np.abs(values[1]))))
        root_length = mp.sqrt(data['length'])
        def integrand(w):
            v = data['length']/2+mp.mpf(float(w))*root_length
            u = float(v/data['length'])
            ratio, loga = sol.sol(u)
            fr = frame(data, lab, coeff, v)
            y = np.array(fr['basis'].tolist(), dtype=float)@np.array([math.exp(loga), math.exp(loga)*ratio])
            rho = U*(mp.mpf('.5')+v/data['length'])
            logp = lab['lambda0']*data['length']/U*(mp.asinh(rho)-mp.asinh(U)
                    -((rho+rho**3/3)-(U+U**3/3))/(1+U*U)**mp.mpf('1.5'))
            return .5*math.exp(2*float(logp))*np.r_[y[0]**2, y[0]*y[1:]]
        integrated, error = quad_vec(integrand, -20., 20., epsabs=2e-11, epsrel=2e-10)
        hsqrt, columnsqrt = integrated[0], integrated[1:]
        h = hsqrt/math.sqrt(length)
        column = columnsqrt/math.sqrt(length)
        columns.append(column)
        ac = -float(lab['c0']*mp.sqrt(1+U*U))
        ref_direction = -ac*np.array(list(lab['nt']), dtype=float)-sign*float(U)*np.array(list(lab['kt']), dtype=float)
        rows.append(dict(sign=sign, maximum_ratio=float(np.max(np.abs(values[0]))),
                         log_amplitude_range=[float(np.min(values[1])), float(np.max(values[1]))],
                         supplied_invariant_ratio_bound=certificate['rho'], log_amplitude_bound=certificate['log_growing_amplitude_bound'],
                         h_times_sqrt_length=float(hsqrt), covariance_column=column.tolist(),
                         normalized_direction_error=float(np.linalg.norm(column/h-ref_direction)),
                         quadrature_error_estimate=float(error), successful_stiff_solve=sol.success))
    matrix = np.column_stack(columns)
    ac = -float(lab['c0']*mp.sqrt(1+U*U))
    nt, kt = np.array(list(lab['nt']), dtype=float), np.array(list(lab['kt']), dtype=float)
    weights = np.linalg.solve(matrix, -ac*nt+.2*float(U)*kt)
    signed = np.linalg.solve(matrix, .7*nt+1.1*kt)
    return dict(ell=data['ell'], length=length, rows=rows,
                positive_squared_amplitudes=weights.tolist(), signed_inverse=signed.tolist(),
                largest_ratio=maximum_ratio, largest_abs_log_amplitude=maximum_log_amplitude,
                checks=dict(ratio_fits_supplied_entry_bound_interval=all(r['maximum_ratio'] < float(r['supplied_invariant_ratio_bound']) for r in rows),
                            growing_amplitude_keeps_positive_lower_bound=all(max(abs(x) for x in r['log_amplitude_range']) < float(r['log_amplitude_bound']) for r in rows),
                            gaussian_peak_resolved=all(.01 < r['h_times_sqrt_length'] < 10 for r in rows),
                            actual_phase_fixture_covariance_directions=all(r['normalized_direction_error'] < 1e-5 for r in rows),
                            integrated_primary_weights_positive=min(weights) > 0,
                            integrated_signed_inverse_has_both_signs=min(signed) < 0 < max(signed)),
                covariance_quadrature_is_truncated_diagnostic=True,
                diagnostic_entry_bound_is_supplied_not_interval_certified=True,
                actual_global_background_numerically_evaluated=False)


def threshold_and_barrier_review():
    thresholds = []
    for h in (mp.mpf(1)/200, mp.mpf(1)/1000):
        for constant in (1, 10**6, 10**20):
            first = eventual_band_threshold(h, 3*constant)
            values = [constant*phase_majorant(ell, h) for ell in (first, first+1, 2*first, 10*first)]
            thresholds.append(dict(h=h, supplied_constant=constant, first_band=first,
                                   sampled_all_later_majorants=values))
    barriers = []
    for scale in (40, 100, 1000):
        barriers.append(barrier_certificate('.6', 2, '.4', scale, '.3'))
    return dict(thresholds=thresholds, barriers=barriers,
                low_band_majorant=phase_majorant(100, mp.mpf(1)/200),
                checks=dict(all_later_band_inequality=all(max(r['sampled_all_later_majorants']) < 1 for r in thresholds),
                            low_bands_not_silently_accepted=phase_majorant(100, mp.mpf(1)/200) > 1,
                            riccati_barriers_point_inward=all(r['checks']['invariant_interval_strict'] for r in barriers),
                            strictly_positive_radial_lower_bound=all(r['lower_radial_factor'] > 0 for r in barriers)),
                actual_profile_constants_or_actual_qstar_evaluated=False)


def cone_and_edge_review():
    ac, u, margin = mp.mpf('1.9'), mp.mpf('1.5'), mp.mpf('.2')
    error = margin/(8*(1/ac+1/u))
    certificate = cone_certificate(ac, u, margin, error)
    reference = np.array([[-float(ac), -float(ac)], [-float(u), float(u)]])
    minimum = math.inf
    max_inverse_error = 0.
    # The certificate proves the continuum ball. The grid only tests signs
    # and matrix conventions, including both edges of the target cone.
    for target_k in (-.8, 0., .8):
        target = np.array([-float(ac), target_k*float(u)])
        z0 = np.linalg.solve(reference, target)
        for angle1 in np.linspace(0., 2*np.pi, 9)[:-1]:
            for angle2 in np.linspace(0., 2*np.pi, 9)[:-1]:
                perturbation = float(error)*np.array([[np.cos(angle1), np.cos(angle2)], [np.sin(angle1), np.sin(angle2)]])
                actual = np.linalg.solve(reference+perturbation, target)
                minimum = min(minimum, float(min(actual)))
                max_inverse_error = max(max_inverse_error, float(np.max(np.abs(actual-z0))))
    log_edge = []
    derivative_error = mp.mpf(0)
    # y=S^(1/2)*zeta*v_positive is a scalar amplitude fixture. Differentiate
    # sqrt(y) and compare with its first two explicit chain-rule formulas.
    for distance in (mp.mpf(1)/8, mp.mpf(1)/16, mp.mpf(1)/32, mp.mpf(1)/64):
        def zeta(x):
            return mp.exp(-1/x**2-2/(1-x)**2)
        def y(x):
            return 10*zeta(x)*(2+mp.sin(x))
        root = mp.sqrt(y(distance))
        yp, ypp = mp.diff(y, distance), mp.diff(y, distance, 2)
        exact1 = yp/(2*root)
        exact2 = ypp/(2*root)-yp**2/(4*root**3)
        calculated1 = mp.diff(lambda x: mp.sqrt(y(x)), distance)
        calculated2 = mp.diff(lambda x: mp.sqrt(y(x)), distance, 2)
        derivative_error = max(derivative_error, abs(exact1-calculated1)/(root*distance**-3),
                               abs(exact2-calculated2)/(root*distance**-6))
        log_edge.append(dict(distance=distance, log_sqrt_zeta_times_delta_to_minus_12=
                             -mp.mpf('.5')/distance**2-1/(1-distance)**2-12*mp.log(distance)))
    bad_reference = np.array([[-float(ac), -float(ac)], [-float(u), float(u)]])
    # Outside the strict cone, the required squared weight is negative.
    bad_weights = np.linalg.solve(bad_reference, np.array([-float(ac), 1.2*float(u)]))
    return dict(cone_certificate=certificate, cone_grid_cases=3*8*8,
                minimum_positive_normalized_weight=minimum, maximum_observed_weight_perturbation=max_inverse_error,
                maximum_scaled_square_root_derivative_error=derivative_error,
                flat_edge_log_values=log_edge, outside_cone_negative_weight=float(min(bad_weights)),
                checks=dict(continuum_cone_certificate=all(certificate['checks'].values()),
                            finite_cone_matrix_crosscheck=minimum > float(margin)/4,
                            finite_neumann_error_crosscheck=max_inverse_error < float(certificate['relative_inverse_error']),
                            square_root_weight_derivatives_keep_half_weight=derivative_error < mp.mpf('1e-85'),
                            edge_weight_absorbs_every_fixed_inverse_power=all(log_edge[i+1]['log_sqrt_zeta_times_delta_to_minus_12'] < log_edge[i]['log_sqrt_zeta_times_delta_to_minus_12'] for i in range(len(log_edge)-1)),
                            outside_cone_not_promoted_to_positive=min(bad_weights) < 0))


def derivative_cutoff_review():
    homogeneous = derivative_exponents(8)
    forced = derivative_exponents(8, source_powers=[0]*9)
    cutoffs = []
    for target in (0, 1, 10, 100):
        ell0 = endpoint_flat_threshold(10**20, '.001', 40, 70, target)
        logs = []
        for ell in (ell0, ell0+1, 2*ell0):
            logs.append(mp.log(10**20)+70*mp.log(ell)+(40+target*mp.log(2))*ell-mp.mpf('.001')*ell**2)
        cutoffs.append(dict(target=target, first_band=ell0, log_error_over_target=logs))
    exponents = exponent_ledger()
    return dict(homogeneous_slow_derivative_powers=homogeneous,
                unit_source_slow_derivative_powers=forced,
                localized_primary_powers=[Fraction(3*j, 1)+Fraction(1, 4) for j in range(9)],
                curl_correction_powers=[Fraction(3*j, 1)+Fraction(13, 4) for j in range(9)],
                averaged_real_curl_covariance_powers=[3*j+6 for j in range(9)],
                cutoff_flatness=cutoffs, exponents=exponents,
                checks=dict(homogeneous_differentiation_has_finite_polynomial_loss=homogeneous == [2*j for j in range(9)],
                            duhamel_adds_one_length_power=forced == [2*j+1 for j in range(9)],
                            all_fixed_endpoint_losses_absorbed=all(max(r['log_error_over_target']) < 0 for r in cutoffs),
                            conditional_real_covariance_has_positive_relative_gain=exponents['relative_physical_stress_gain'] > 0,
                            covariance_gain_is_not_absolute_residual_decay=exponents['physical_covariance_remainder_power'] < 0 and not exponents['full_residual_gain']))


def support_review():
    bg, tg = mp.sqrt(2)-1, 4+mp.sqrt(2)
    vr, vt = mp.matrix([1, -bg]), mp.matrix([bg, 1])
    dual = vt.T/(1+bg*bg)
    projector = mp.eye(2)-vt*dual
    matrix = mp.matrix([[3, 1], [1, 5]])
    max_error = omitted_cover = mp.mpf(0)
    data = band(fixture_first_band())
    ci, r0 = data['ci'], mp.mpf('.04')
    for level_difference in (0, 1, 3):
        current_v, earlier_v = mp.mpf('.63')*data['length'], mp.mpf('.27')*data['length']
        transverse = mp.mpf('.01')
        center = mp.matrix([mp.mpf('.2'), mp.mpf('.3')])
        eta = ci*current_v-r0
        current = matrix**(-level_difference)*(center+transverse*vr+eta*vt)
        earlier = current+tg**(-level_difference)*ci*(earlier_v-current_v)*vt
        band_point = matrix**level_difference*earlier-center
        max_error = max(max_error, abs((dual*band_point)[0]-(ci*earlier_v-r0)),
                        abs((vr.T*band_point)[0]/(1+bg*bg)-transverse),
                        mp.norm(projector*projector-projector), mp.norm(projector*vt), mp.norm(projector*vr-vr))
        bad = current+ci*(earlier_v-current_v)*vt
        if level_difference:
            omitted_cover = max(omitted_cover, abs((dual*(matrix**level_difference*bad-center))[0]-(ci*earlier_v-r0)))
    # Two common-torus points can be identical on the finer band torus.
    # A source descending only to the common torus can distinguish them.
    shift = matrix**-1*mp.matrix([1, 0])
    source_preimage_gap = abs(mp.cos(0)-mp.cos(2*mp.pi*shift[0]))
    return dict(maximum_affine_common_path_identity_error=max_error,
                omitted_cover_factor_gap=omitted_cover, source_on_distinct_preimages_gap=source_preimage_gap,
                checks=dict(pulse_path_fixes_slow_and_transverse_coordinates=max_error < mp.mpf('1e-85'),
                            common_torus_cover_factor_required=omitted_cover > mp.mpf('.01'),
                            allowed_sources_must_not_be_averaged_over_preimages=source_preimage_gap > mp.mpf('.01')),
                global_color_centers_or_all_rectangles_materialized=False)


def linear_remainder_review():
    """Differentiate the full real harmonic, then compare with (9.1)+(9.2).

    A local exponential in the transverse coordinate exercises the fast
    radial operator as Dr f=f_R+M(R)*lambda_xi*f, including Dr^2.
    The phase is independent of that transverse coordinate.
    """
    eps, h, k, p, pz, x0 = .03, .005, 6, 1/3, -.43, .21
    maximum = omitted_angular = omitted_mixed = omitted_pressure = 0.
    for point in ((1.4, .12, .62, 3.), (1.31, -.16, .7, 11.), (1.6, .23, .55, 17.)):
        for multiplier in (.0, .7):
            r, z, time, v = variables(point)
            b, f, g = background_jet(r, z, time, eps, h)
            swirl = r*f
            phase = pz*z/eps+x0*r-v*(p*f+pz*g)
            n = [phase.diff(0), p/r, eps*phase.diff(1)]
            amplitudes = [(r+.1*z)*(.3*v).sin(), (1+.2*time)*(.2*v).cos(), .4*r*z+.1*v]
            pressure_amp = .17*r*z+.03*time*v
            transverse_rates, pressure_rate = (.3, -.2, .4), .15
            radial_multiplier = multiplier*r**.3
            def dr(a, rate):
                return a.diff(0)+radial_multiplier*rate*a
            def drr(a, rate):
                return (a.diff(0).diff(0)+2*radial_multiplier*rate*a.diff(0)
                        +(radial_multiplier.diff(0)*rate+radial_multiplier**2*rate**2)*a)
            cosine, sine = (k*phase).cos(), (k*phase).sin()
            velocity = [a*cosine for a in amplitudes]
            pressure = -pressure_amp*sine
            theta = [-k*p*a*sine for a in amplitudes]
            scalar_lap = [drr(w, rate)+dr(w, rate)/r+eps**2*w.diff(1).diff(1)
                          -(k*p)**2*w/r**2 for w, rate in zip(velocity, transverse_rates)]
            vector_lap = [scalar_lap[0]-velocity[0]/r**2-2*theta[1]/r**2,
                          scalar_lap[1]-velocity[1]/r**2+2*theta[0]/r**2, scalar_lap[2]]
            adv = [-eps*w.diff(2)+w.diff(3)+b*dr(w, rate)+f*th+g*eps*w.diff(1)
                   for w, rate, th in zip(velocity, transverse_rates, theta)]
            direct = [adv[0]+velocity[0]*b.diff(0)+velocity[2]*eps*b.diff(1)-2*f*velocity[1]
                       -eps*vector_lap[0]+dr(pressure, pressure_rate),
                      adv[1]+(swirl.diff(0)+f)*velocity[0]+velocity[2]*eps*swirl.diff(1)+b*velocity[1]/r
                       -eps*vector_lap[1]-k*p*pressure_amp*cosine/r,
                      adv[2]+velocity[0]*g.diff(0)+velocity[2]*eps*g.diff(1)
                       -eps*vector_lap[2]+eps*pressure.diff(1)]
            aa = amplitudes
            kmat_a = [-2*f*aa[1], (2*f+r*f.diff(0))*aa[0], g.diff(0)*aa[0]]
            n2 = sum(x*x for x in n)
            principal = [a.diff(3)+ka+eps*k*k*n2*a-k*ni*pressure_amp
                         for a, ka, ni in zip(aa, kmat_a, n)]
            defect = (-eps*phase.diff(2)+phase.diff(3)+b*phase.diff(0)+f*p+g*eps*phase.diff(1))
            base_derivative = [aa[0]*b.diff(0)+aa[2]*eps*b.diff(1),
                               aa[2]*eps*swirl.diff(1)+b*aa[1]/r, aa[2]*eps*g.diff(1)]
            amplitude_lap = [drr(a, rate)+dr(a, rate)/r+eps**2*a.diff(1).diff(1)
                             -(a/r**2 if i < 2 else 0) for i, (a, rate) in enumerate(zip(aa, transverse_rates))]
            real = [-eps*a.diff(2)+b*dr(a, rate)+g*eps*a.diff(1)+bd-eps*lap
                    for a, rate, bd, lap in zip(aa, transverse_rates, base_derivative, amplitude_lap)]
            phase_div = n[0].diff(0)+n[0]/r+eps*n[2].diff(1)
            angular_connection = [-2*n[1]*aa[1]/r, 2*n[1]*aa[0]/r, Jet(0)]
            pressure_slow = [dr(pressure_amp, pressure_rate), Jet(0), eps*pressure_amp.diff(1)]
            imaginary, mixed = [], []
            for a, rate, ac, pg in zip(aa, transverse_rates, angular_connection, pressure_slow):
                mix = 2*(n[0]*dr(a, rate)+n[2]*eps*a.diff(1))+phase_div*a
                mixed.append(mix)
                imaginary.append(k*defect*a+pg-eps*k*(mix+ac))
            expanded = np.array([((pr+rr)*cosine-ii*sine).value for pr, rr, ii in zip(principal, real, imaginary)])
            actual = np.array([a.value for a in direct])
            maximum = max(maximum, float(np.linalg.norm(actual-expanded)/(1+np.linalg.norm(actual)+np.linalg.norm(expanded))))
            omitted_angular = max(omitted_angular, float(np.linalg.norm([(eps*k*x*sine).value for x in angular_connection])))
            omitted_mixed = max(omitted_mixed, float(np.linalg.norm([(eps*k*x*sine).value for x in mixed])))
            omitted_pressure = max(omitted_pressure, float(np.linalg.norm([(x*sine).value for x in pressure_slow])))
    return dict(cases=6, maximum_full_linear_harmonic_identity_relative_error=maximum,
                angular_viscous_connection_omission_gap=omitted_angular,
                mixed_phase_amplitude_omission_gap=omitted_mixed, pressure_amplitude_gradient_omission_gap=omitted_pressure,
                checks=dict(full_principal_plus_linear_remainder_identity=maximum < 2e-10,
                            angular_viscous_connection_retained=omitted_angular > 1e-3,
                            mixed_phase_amplitude_terms_retained=omitted_mixed > 1e-3,
                            pressure_amplitude_gradient_retained=omitted_pressure > 1e-3),
                nonlinear_or_mean_correction_terms_checked=False)


def partition_scaling_review():
    h = rational(H)
    a = mp.mpf('.5')+h
    q = mp.power(2, -100)
    target = mp.matrix([mp.mpf('.7'), mp.mpf('-.3')])
    expected = q**(-a-mp.mpf('.5'))*target
    total, omitted_square, omitted_epsilon = mp.matrix([0, 0]), mp.matrix([0, 0]), mp.matrix([0, 0])
    for level, square_weight in ((100, mp.mpf(2)/5), (101, mp.mpf(3)/5)):
        scale = mp.power(2, -level)
        eps = scale**h
        local_target = (scale/q)**(a+mp.mpf('.5'))*target
        total += scale**(-2*a)*square_weight*eps*local_target
        omitted_square += scale**(-2*a)*eps*local_target
        omitted_epsilon += scale**(-2*a)*square_weight*local_target
    relative = mp.norm(total-expected)/mp.norm(expected)
    return dict(physical_leading_covariance_relative_error=relative,
                omitted_squared_partition_relative_gap=mp.norm(omitted_square-expected)/mp.norm(expected),
                omitted_epsilon_relative_gap=mp.norm(omitted_epsilon-expected)/mp.norm(expected),
                checks=dict(squared_partition_recovers_physical_target=relative < mp.mpf('1e-85'),
                            squared_slow_cutoff_cannot_be_omitted=mp.norm(omitted_square-expected)/mp.norm(expected) > mp.mpf('.1'),
                            wave_amplitude_normalization_cannot_be_omitted=mp.norm(omitted_epsilon-expected)/mp.norm(expected) > mp.mpf('.1')))
