"""Moving-normal reference pulse; no actual band-uniform estimate asserted.

Divide out the scalar Gaussian envelope before integration. This avoids
accepting an exponentially tiny zero solution at the endpoints.
"""
import math
import numpy as np
from scipy.integrate import quad_vec, solve_ivp

F, GABS, U, K_CARRIER = 1., 3., 1.5, 16
N = np.array([-1/math.sqrt(2), 1/math.sqrt(2)])
KT = np.array([-N[1], N[0]])
LAMBDA0 = math.sqrt(-2*F*N[0]*(2*F*N[0]+GABS))
C0 = LAMBDA0/(2*F*N[0])
AC = -C0*math.sqrt(1+U*U)
BS2 = LAMBDA0/(1+U*U)**1.5  # epsilon*k^2 = nu = 1
BS = math.sqrt(BS2)
SHEAR = GABS*N
KMAT = np.array([[0., -2*F, 0.], [2*F+SHEAR[0], 0., 0.], [SHEAR[1], 0., 0.]])


def coefficients(v, length, sign):
    s = sign*U*(.5+v/length)
    normal = BS*np.r_[s, KT]
    normal_prime = BS*np.array([sign*U/length, 0., 0.])
    n2 = normal@normal
    projected = -KMAT+np.outer(normal, normal@KMAT-normal_prime)/n2
    growth = LAMBDA0/math.sqrt(1+s*s)
    damping = BS2*(1+s*s)
    return normal, normal_prime, projected, growth, damping


def log_envelope(v, length):
    rho = U*(.5+v/length)
    return LAMBDA0*length/U*(math.asinh(rho)-math.asinh(U)
        -((rho+rho**3/3)-(U+U**3/3))/(1+U*U)**1.5)


def gaussian_constants():
    low, high = U/2, 3*U/2
    common = 2/(1+U*U)**1.5
    return (U*LAMBDA0*low*((1+high*high)**-1.5+common),
            U*LAMBDA0*high*((1+low*low)**-1.5+common))


def step(x):
    if x <= 0:
        return 0., 0.
    if x >= 1:
        return 1., 0.
    a, b = math.exp(-1/x), math.exp(-1/(1-x))
    return a/(a+b), (a*b/(a+b)**2)*(1/x**2+1/(1-x)**2)


def cutoff(v, length):
    a, da = step((v/length-.1)/.15)
    b, db = step((.9-v/length)/.15)
    return a*b, (da*b-a*db)/(.15*length)


def solve_pulse(length, sign, tighter=False):
    s = sign*U/2
    initial = np.r_[1., -s*KT+C0*math.sqrt(1+s*s)*N]

    def rhs(v, y):
        _, _, projected, growth, _ = coefficients(v, length, sign)
        return projected@y-growth*y

    solution = solve_ivp(rhs, (0., length), initial, method='DOP853', dense_output=True,
                         rtol=2e-12 if tighter else 2e-11,
                         atol=2e-13 if tighter else 2e-12, max_step=length/80)
    if not solution.success:
        raise ValueError('Reference pulse integration failed')
    return solution, rhs


def pulse_review(lengths=(40, 100, 240), tighter=False):
    rows = []
    max_constraint = max_principal = max_energy = max_gaussian = 0.
    min_moving_gap = min_pressure_gap = min_viscosity_gap = math.inf
    for length in lengths:
        columns, errors, weights, quadrature_errors = [], [], [], []
        cutoff_bound = 0.
        for sign in (1, -1):
            sol, rhs = solve_pulse(length, sign, tighter)

            def integrand(v):
                y = sol.sol(v)
                psi, _ = cutoff(v, length)
                factor = .5/length*psi**2*math.exp(2*log_envelope(v, length))
                return factor*np.r_[y[0]**2, y[0]*y[1:]]

            integrated, integration_error = quad_vec(integrand, 0., length,
                                                     epsabs=2e-13, epsrel=2e-11)
            h, column = integrated[0], integrated[1:]
            quadrature_errors.append(float(integration_error))
            weights.append(float(h))
            columns.append(column)
            errors.append(float(np.linalg.norm(column/h+AC*N+sign*U*KT)))
            for v in np.linspace(0., length, 121):
                normal, normal_prime, projected, growth, damping = coefficients(v, length, sign)
                y = sol.sol(v)
                yp = rhs(v, y)
                n2 = normal@normal
                pi_over_p = 1j*(normal@KMAT@y-normal_prime@y)/(K_CARRIER*n2)
                principal = yp+growth*y+KMAT@y+1j*K_CARRIER*normal*pi_over_p
                max_constraint = max(max_constraint, abs(normal@y))
                max_principal = max(max_principal, float(np.linalg.norm(principal)))
                energy = y@yp+SHEAR@(y[0]*y[1:])+growth*(y@y)
                max_energy = max(max_energy, abs(energy))
                c, bigc = gaussian_constants()
                offset = (v-length/2)**2/(2*length)
                logp = log_envelope(v, length)
                max_gaussian = max(max_gaussian, -bigc*offset-logp, logp+c*offset, logp)
                _, dpsi = cutoff(v, length)
                cutoff_bound = max(cutoff_bound, abs(dpsi)*math.exp(logp)*np.linalg.norm(y))
                if .2*length < v < .8*length:
                    bad_a = projected+np.outer(normal, normal_prime)/n2
                    moving_gap = abs(normal_prime@y+normal@((bad_a-growth*np.eye(3))@y))
                    bad_pi = (normal@KMAT@y-normal_prime@y)/(K_CARRIER*n2)
                    pressure_gap = np.linalg.norm(yp+growth*y+KMAT@y+1j*K_CARRIER*normal*bad_pi)
                    min_moving_gap = min(min_moving_gap, moving_gap)
                    min_pressure_gap = min(min_pressure_gap, float(pressure_gap))
                    min_viscosity_gap = min(min_viscosity_gap, float(damping*np.linalg.norm(y)))
        matrix = np.column_stack(columns)
        target = -AC*N+.2*U*KT
        squared_amplitudes = np.linalg.solve(matrix, target)
        signed_target = .7*N+1.1*KT
        signed_solution = np.linalg.solve(matrix, signed_target)
        base = np.sqrt(squared_amplitudes)
        delta = signed_solution/(2*base)
        inverse_error = np.linalg.norm(matrix@(2*base*delta)-signed_target)
        rows.append(dict(length=length, weights=weights,
                         sqrt_length_times_weights=(np.sqrt(length)*np.array(weights)).tolist(),
                         direction_errors=errors, covariance_columns=matrix.tolist(),
                         positive_squared_amplitudes=squared_amplitudes.tolist(),
                         signed_linear_solution=signed_solution.tolist(),
                         signed_inverse_error=float(inverse_error),
                         sampled_cutoff_principal_residual_max=float(cutoff_bound),
                         covariance_quadrature_error_estimates=quadrature_errors))
    checks = dict(moving_normal_constraint=max_constraint < 2e-9,
                  complex_pressure_principal_equation=max_principal < 2e-9,
                  shear_viscosity_energy_identity=max_energy < 2e-9,
                  continuum_gaussian_bound_sample_crosscheck=max_gaussian < 2e-12,
                  positive_covariance_inverse_in_reference_fixture=all(min(r['positive_squared_amplitudes']) > 0 for r in rows),
                  signed_linear_inverse_retained=all(min(r['signed_linear_solution']) < 0 < max(r['signed_linear_solution']) for r in rows),
                  covariance_directions_approach_reference=max(rows[-1]['direction_errors']) < max(rows[0]['direction_errors']),
                  weight_has_reference_sqrt_length_scale=all(.1 < w < 10 for r in rows for w in r['sqrt_length_times_weights']),
                  moving_normal_omission_detected=min_moving_gap > 1e-5,
                  pressure_i_omission_detected=min_pressure_gap > 1e-3,
                  viscous_damping_omission_detected=min_viscosity_gap > 1e-3)
    return dict(rows=rows, maximum_normal_constraint_error=float(max_constraint),
                maximum_normalized_principal_equation_error=float(max_principal),
                maximum_normalized_energy_error=float(max_energy),
                maximum_gaussian_inequality_defect=float(max_gaussian),
                minimum_moving_normal_omission_gap=float(min_moving_gap),
                minimum_pressure_i_omission_gap=float(min_pressure_gap),
                minimum_viscosity_omission_gap=float(min_viscosity_gap),
                gaussian_constants=list(gaussian_constants()), checks=checks,
                actual_band_coefficients_evaluated=False,
                uniform_actual_pulse_or_parameter_derivative_bound_verified=False)
