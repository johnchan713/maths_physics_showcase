"""Independent finite diagnostics; none samples the actual giant-N profile."""
from math import factorial
from fractions import Fraction
from pathlib import Path
import sys
import mpmath as mp

from higher_bounds import HERE, module_at
from higher_jets import (regular_variance_factor, signed_F, root_jet,
                         implicit_phase_jet, pulled_back_third,
                         fixed_frequency_derivative)

sys.path.insert(0, str(HERE.parent/'stress_realization_audit'))
from repair import Repair
from schedule import step
sys.path.insert(0, str(HERE))
profiles = module_at('higher_jet_manufactured_profiles', HERE.parent/'paper_profile_audit'/'profiles.py')


def small(value):
    return mp.nstr(value, 16)


def variance_review():
    rows = []
    with mp.workdps(90):
        for z in map(mp.mpf, ('-2', '-.7', '.4', '1.5')):
            direct = (mp.besseli(0, 2*z)/mp.besseli(0, z)**2-1)/z**2
            rows.append(dict(z=str(z), error=small(abs(regular_variance_factor(z)-direct))))
        d0 = mp.mpf('.2')
        origin = mp.diff(lambda mu: signed_F(mu, mp.mpf('.7'), d0), 0)
        origin_error = abs(origin-d0/mp.sqrt(2))
        wrong_denominator = mp.diff(lambda mu: signed_F(mu, mp.mpf('.7'), d0)**2, 0)
        second = mp.diff(regular_variance_factor, 0, 2)
    return dict(manufactured=True, samples=rows, signed_F_mu_at_zero=small(origin),
                variance_mu_at_zero=small(wrong_denominator), R_second_at_zero=small(second),
                checks=dict(positive_series_matches_independent_Bessel_formula=all(mp.mpf(row['error']) < mp.mpf('1e-70') for row in rows),
                            regular_origin_value=regular_variance_factor(0) == mp.mpf('.5'),
                            signed_root_derivative=origin_error < mp.mpf('1e-60'),
                            direct_variance_division_is_singular=wrong_denominator == 0,
                            even_factor_second_derivative=abs(second+mp.mpf(3)/16) < mp.mpf('1e-60')))


def root_review():
    rows = []
    with mp.workdps(80):
        for constant in ('.1', '0'):
            d0 = mp.mpf('.2')
            target = lambda x: mp.mpf(constant)+mp.mpf('.03')*x+mp.mpf('.007')*x*x-mp.mpf('.002')*x**3
            pressure = lambda x: mp.mpf('.4')+mp.mpf('.2')*x+mp.mpf('.04')*x*x
            jet = root_jet(target, pressure, d0, mp.mpf(0))
            solve = lambda x: mp.findroot(lambda mu: signed_F(mu, pressure(x), d0)-target(x), target(x)*mp.sqrt(2)/d0)
            error = max(abs(jet[k]-mp.diff(solve, 0, k)) for k in range(4))
            frozen = root_jet(target, lambda x: pressure(0), d0, mp.mpf(0))
            rows.append(dict(target_constant=constant, maximum_jet_error=small(error),
                             frozen_pressure_third_gap=small(abs(frozen[3]-jet[3]))))
    return dict(manufactured=True, samples=rows,
                checks=dict(third_implicit_root_matches_independent_solved_function=all(mp.mpf(row['maximum_jet_error']) < mp.mpf('1e-55') for row in rows),
                            moving_pressure_terms_detected=mp.mpf(rows[0]['frozen_pressure_third_gap']) > mp.mpf('1e-6'),
                            zero_variance_boundary_included=rows[1]['target_constant'] == '0'))


def phase_review():
    with mp.workdps(80):
        a = lambda x: mp.mpf('.8')+mp.mpf('.1')*x+mp.mpf('.04')*x*x
        v = lambda x: mp.mpf('2.1')+mp.mpf('.07')*x+mp.mpf('.03')*x*x+mp.mpf('.01')*x**3
        E = lambda x: 1+mp.mpf('.1')*x+mp.mpf('.02')*x*x
        k = lambda x: mp.sqrt(2*(v(x)/a(x)-1))
        forward = lambda x, theta: theta/(2*mp.pi)-a(x)*k(x)**2*mp.sin(2*theta)/(8*mp.pi*v(x))
        A = lambda x, theta: -a(x)**2*k(x)**2*mp.sin(2*theta)/(16*mp.pi*v(x))
        B = lambda x, theta: E(x)*a(x)*k(x)*mp.cos(theta)/(4*mp.pi)
        rows = []
        for phase in map(mp.mpf, ('.17', '.29', '.43')):
            theta = mp.findroot(lambda t: forward(0, t)-phase, 2*mp.pi*phase)
            jet = implicit_phase_jet(forward, 0, theta)
            solve = lambda x: mp.findroot(lambda t: forward(x, t)-phase, theta)
            inverse_error = max(abs(jet[j]-mp.diff(solve, 0, j)) for j in range(4))
            for name, function in (('A', A), ('B', B)):
                expected = mp.diff(lambda x: function(x, solve(x)), 0, 3)
                actual = pulled_back_third(function, 0, jet)
                frozen = mp.diff(lambda x: function(x, theta), 0, 3)
                rows.append(dict(phase=str(phase), primitive=name, inverse_error=small(inverse_error),
                                 composition_error=small(abs(actual-expected)),
                                 frozen_phase_gap=small(abs(frozen-expected))))
    return dict(manufactured_zero_pressure=True, samples=rows,
                checks=dict(third_inverse_phase_derivative=all(mp.mpf(row['inverse_error']) < mp.mpf('1e-55') for row in rows),
                            pulled_back_third_chain_rule=all(mp.mpf(row['composition_error']) < mp.mpf('1e-55') for row in rows),
                            freezing_the_phase_is_detected=all(mp.mpf(row['frozen_phase_gap']) > mp.mpf('1e-8') for row in rows)))


def modulation_review():
    with mp.workdps(80):
        primitive = lambda y, e, s: (mp.mpf('.2')+y/32+e/16+y*e/64+e*e/32)*mp.sin(2*mp.pi*s)+(1+e**3)*mp.cos(4*mp.pi*s)/100
        logE = lambda y, e: y/8+e/4+y*e/32+e*e/64
        y, eta = mp.mpf('.13'), mp.mpf('.27')
        chain_error, exponential_error, cases = mp.mpf(0), mp.mpf(0), 0
        for N in (8, 32, 128):
            beta = lambda yy, ee: primitive(yy, ee, N*yy)/N
            field = lambda yy, ee: mp.exp(logE(yy, ee)+beta(yy, ee))
            for r, s in ((0, 3), (1, 2), (2, 0), (1, 1), (0, 2)):
                actual = fixed_frequency_derivative(primitive, y, eta, N, r, s)
                chain_error = max(chain_error, abs(actual-mp.diff(beta, (y, eta), (r, s))))
                cases += 1
            L = lambda yy, ee: logE(yy, ee)+beta(yy, ee)
            L1, L2, L3 = [mp.diff(lambda ee: L(y, ee), eta, j) for j in (1, 2, 3)]
            predicted = field(y, eta)*(L3+3*L1*L2+L1**3)
            exponential_error = max(exponential_error, abs(predicted-mp.diff(lambda ee: field(y, ee), eta, 3)))
        N = 32
        radial_full = fixed_frequency_derivative(primitive, y, eta, N, 2, 0)
        without_fast_second = radial_full-N*mp.diff(primitive, (y, eta, N*y), (0, 0, 2))
        missing_radial = abs(radial_full-without_fast_second)
        L1 = mp.diff(lambda ee: logE(y, ee)+primitive(y, ee, N*y)/N, eta)
        missing_exponential = abs(field(y, eta)*L1**3)
    return dict(manufactured=True, chain_sample_count=cases, chain_error=small(chain_error),
                exponential_error=small(exponential_error),
                missing_N_phase_second_gap=small(missing_radial), missing_exponential_cube_gap=small(missing_exponential),
                checks=dict(all_mixed_and_angular_chain_terms=chain_error < mp.mpf('1e-55'),
                            exponential_third_product_rule=exponential_error < mp.mpf('1e-55'),
                            positive_N_radial_loss_detected=missing_radial > 1,
                            omitted_exponential_cube_detected=missing_exponential > mp.mpf('1e-4')))


def fourth_implicit_review():
    with mp.workdps(90):
        c = lambda t: mp.sin(t)+t*t/5+mp.mpf('.01')
        q = lambda t: mp.exp(t)/7
        eta = mp.mpf('.37')
        cc, qq = [[mp.diff(f, eta, j) for j in range(5)] for f in (c, q)]
        z = lambda t: 2*c(t)+q(t)*c(t)**2
        terms = [8*qq[0]*cc[1]*cc[3], 6*qq[0]*cc[2]**2,
                 8*qq[1]*cc[0]*cc[3], 24*qq[1]*cc[1]*cc[2],
                 12*qq[2]*cc[0]*cc[2], 12*qq[2]*cc[1]**2,
                 8*qq[3]*cc[0]*cc[1], qq[4]*cc[0]**2]
        error = abs((2+2*qq[0]*cc[0])*cc[4]-mp.diff(z, eta, 4)+mp.fsum(terms))
    return dict(identity_error=small(error), omitted_term_gaps=[small(abs(v)) for v in terms],
                checks=dict(fourth_implicit_identity=error < mp.mpf('1e-70'),
                            all_eight_product_terms_detected=all(abs(v) > mp.mpf('1e-5') for v in terms)))


def repair_C3_review():
    with mp.workdps(80):
        repair, amplitude, eta = Repair('.001'), mp.mpf('1e-16'), mp.mpf('.23')
        e = lambda x: mp.exp(x/3+x*x/7)
        raw = lambda x: [mp.sin(x), mp.cos(2*x), mp.exp(x), x**3, mp.log(2+x)]
        weights = lambda x: [1/e(x), 1/e(x)**2, 1/e(x), 1/e(x)**2, 1/e(x)**2]
        target = lambda x: mp.matrix([amplitude*a*b for a, b in zip(raw(x), weights(x))])
        solve = lambda x: repair.solve(target(x))
        jets = [mp.matrix([mp.diff(lambda x: solve(x)[i], eta, j) for i in range(5)]) for j in range(4)]
        target3 = mp.matrix([mp.diff(lambda x: target(x)[i], eta, 3) for i in range(5)])
        residual = repair.B*jets[3]+2*repair.quadratic(jets[0], jets[3])+6*repair.quadratic(jets[1], jets[2])-target3
        error = max(abs(v) for v in residual)/amplitude
        nonlinear_gap = max(abs(v) for v in 6*repair.quadratic(jets[1], jets[2]))/amplitude**2
        frozen = mp.matrix([amplitude*weights(eta)[i]*mp.diff(lambda x: raw(x)[i], eta, 3) for i in range(5)])
        normalization_gap = max(abs(v) for v in target3-frozen)/amplitude
        sample_target_norm = max(mp.fsum(abs(mp.diff(lambda x: target(x)[i], eta, j))/factorial(j) for j in range(4)) for i in range(5))
        sample_coefficient_norm = max(mp.fsum(abs(jets[j][i])/factorial(j) for j in range(4)) for i in range(5))
    return dict(manufactured=True, quadrature_is_not_a_continuum_certificate=True,
                normalized_third_identity_error=small(error), nonlinear_gap_scaled_by_target_squared=small(nonlinear_gap),
                frozen_normalization_gap_scaled_by_target=small(normalization_gap),
                sampled_C3_coefficient_to_target_ratio=small(sample_coefficient_norm/sample_target_norm),
                checks=dict(exact_nonlinear_third_jet=error < mp.mpf('1e-50'),
                            quadratic_cross_term_retained=nonlinear_gap > mp.mpf('1e-5'),
                            normalization_derivatives_retained=normalization_gap > mp.mpf('.01'),
                            sample_inside_inherited_inverse_cap=sample_coefficient_norm < 2000*sample_target_norm))


def moment_review():
    data = dict(F=[[0, 0, 1], [1, 0, .125], [0, 4, .0625], [1, 5, .03125]],
                U=[[0, 0, .25], [0, 1, .5], [0, 4, .03125], [1, 3, .015625]], axis_pressure=[])
    profile = profiles.ManufacturedProfile(data)
    f_exact = {(i, j): Fraction(str(c)) for i, j, c in data['F']}
    u_exact = {(i, j): Fraction(str(c)) for i, j, c in data['U']}
    def product(a, b):
        result = {}
        for (i, j), c in a.items():
            for (k, l), d in b.items():
                key = (i+k, j+l)
                result[key] = result.get(key, Fraction())+c*d
        return result
    def plus(a, b):
        result = dict(a)
        for key, c in b.items():
            result[key] = result.get(key, Fraction())+c
        return result
    times_x = lambda a, factor: {(i+1, j): factor*c for (i, j), c in a.items()}
    exact_densities = dict(M=u_exact, I=times_x(f_exact, 2),
                           J=times_x(product(u_exact, f_exact), 2),
                           S=plus(product(u_exact, u_exact), times_x(product(f_exact, f_exact), -1)),
                           Pi=product(f_exact, f_exact))
    rows = []
    with mp.workdps(70):
        X, eta = mp.mpf('.7'), mp.mpf('.3')
        f, u = profile.F, profile.U
        densities = dict(M=lambda x, e: u(x, e), I=lambda x, e: 2*x*f(x, e),
                         J=lambda x, e: 2*x*u(x, e)*f(x, e),
                         S=lambda x, e: u(x, e)**2-x*f(x, e)**2,
                         Pi=lambda x, e: f(x, e)**2)
        for name, density in densities.items():
            polynomial = getattr(profile, name)
            for _ in range(4):
                polynomial = polynomial.derivative(1)
            actual = polynomial(X, eta)
            exact = mp.fsum(mp.mpf(c.numerator)/c.denominator/(i+1)*X**(i+1)
                           *factorial(j)/factorial(j-4)*eta**(j-4)
                           for (i, j), c in exact_densities[name].items() if j >= 4)
            expected = mp.quad(lambda x: mp.diff(lambda e: density(x, e), eta, 4), [0, X])
            rows.append(dict(moment=name, fourth_derivative=small(expected), error=small(abs(exact-expected)),
                             legacy_binary64_polynomial_error=small(abs(actual-expected))))
    return dict(manufactured=True, samples=rows,
                axis_pressure_integrand_uses_F_squared=True,
                checks=dict(five_original_fourth_moment_identities=all(mp.mpf(row['error']) < mp.mpf('1e-55') for row in rows),
                            legacy_float_polynomial_agrees_to_roundoff=all(mp.mpf(row['legacy_binary64_polynomial_error']) < mp.mpf('2e-14') for row in rows),
                            all_moments_have_nonzero_fourth_jet=all(abs(mp.mpf(row['fourth_derivative'])) > mp.mpf('.001') for row in rows)))


def derivative_gap_review():
    rows = [dict(frequency=m, factorial_C3_upper=str(mp.mpf(1)/m**3+mp.mpf(1)/m**2+mp.mpf(1)/(2*m)+mp.mpf(1)/6),
                 raw_fourth_peak=m, stress_third_peak=m) for m in (8, 32, 128, 512)]
    return dict(exact_family='f_m=m^(-3)*sin(m*eta)', rows=rows,
                missing_C4_normalization_contribution=1,
                checks=dict(C3_uniform=all(mp.mpf(row['factorial_C3_upper']) < 1 for row in rows),
                            input_C4_and_stress_C3_unbounded=rows[-1]['stress_third_peak'] == 64*rows[0]['stress_third_peak'],
                            f_inverse_square_norm_must_be_25=24+1 == 25))


def reference_review():
    """Original cumulative p1 formula versus its reference source integral."""
    with mp.workdps(80):
        h, D = mp.mpf('.005'), mp.mpf('.495')
        a = lambda e: 1+e**4/16
        b = lambda e: mp.mpf(1)/8+e**5/32
        c = lambda e: mp.mpf(1)/4+e/2+e**4/32
        dd = lambda e: e**3/64
        R = lambda x, e: a(e)+b(e)*x
        F = lambda x, e: mp.exp(e/5)*R(x, e)
        U = lambda x, e: c(e)+dd(e)*x
        M = lambda x, e: c(e)*x+dd(e)*x*x/2
        I = lambda x, e: mp.exp(e/5)*(a(e)*x*x+mp.mpf(2)/3*b(e)*x**3)
        J = lambda x, e: mp.exp(e/5)*(a(e)*c(e)*x*x+mp.mpf(2)/3*(a(e)*dd(e)+b(e)*c(e))*x**3+b(e)*dd(e)*x**4/2)
        L = lambda e: 1-2*h*e*e
        W = lambda x, e: 1-(2*D*e*M(x, e)+(1-e*e)*mp.diff(lambda ee: M(x, ee), e))/x
        def p1(x, e):
            numerator = ((1-h)*I(x, e)-D*e*mp.diff(lambda ee: I(x, ee), e)
                         -(1-e*e)*mp.diff(lambda ee: J(x, ee), e)+2*(h-D)*e*J(x, e))
            return x/L(e)*(-W(x, e)+numerator/(2*x*x*F(x, e)))
        lr = lambda x, e: 1+x*b(e)/R(x, e)
        def source(x, e, freeze_axis_factor=False):
            log_eta = mp.diff(lambda ee: mp.log(R(x, ee)), e)+(0 if freeze_axis_factor else mp.mpf(1)/5)
            return -W(x, e)*lr(x, e)-h*(1-2*e*U(x, e))-(D*e+(1-e*e)*U(x, e))*log_eta
        xa, X, eta = mp.mpf('.2'), mp.mpf('.9'), mp.mpf('.3')
        ode = lambda e: X*mp.diff(lambda x: p1(x, e), X)-X*source(X, e)/L(e)+lr(X, e)*p1(X, e)
        ode_error = max(abs(mp.diff(ode, eta, j)) for j in range(5))
        def integral(e, wrong_factor=False, frozen=False):
            factor = (lambda x: R(x, e)) if wrong_factor else (lambda x: x*R(x, e))
            # For x*R the forcing in y becomes R*x*Sq dx; for the
            # deliberately wrong R integrating factor the x is absent.
            return (factor(xa)*p1(xa, e)+mp.quad(lambda x: factor(x)*source(x, e, frozen)/L(e), [xa, X]))/factor(X)
        errors = [abs(mp.diff(lambda e: integral(e)-p1(X, e), eta, j)) for j in range(5)]
        wrong_factor = abs(integral(eta, wrong_factor=True)-p1(X, eta))
        frozen_axis = abs(integral(eta, frozen=True)-p1(X, eta))
    return dict(manufactured=True, maximum_original_source_jet_error=small(ode_error),
                maximum_integral_jet_error=small(max(errors)),
                missing_X_integrating_factor_gap=small(wrong_factor),
                frozen_axis_gradient_gap=small(frozen_axis),
                checks=dict(original_source_identity_through_C4=ode_error < mp.mpf('1e-55'),
                            exact_reference_integral_through_C4=max(errors) < mp.mpf('1e-55'),
                            correct_integrating_factor_required=wrong_factor > mp.mpf('1e-4'),
                            axis_gradient_remains_in_source=frozen_axis > mp.mpf('1e-4')))
