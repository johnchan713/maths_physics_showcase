"""Independent residual comparisons, sparse-word and first-moment controls."""
import itertools
from pathlib import Path
import sys
import mpmath as mp
import numpy as np

from system import DIAGONAL, axis_slopes, matrices, correction_vector
from series import build
from majorant import allowed_word, log_term, log_sum_bound
from review import fixtures
from residual import series_coefficients

sys.path.insert(0, str(Path(__file__).resolve().parent.parent/'stress_realization_audit'))
from repair import Repair, precondition_target
from schedule import step_prime
sys.path.insert(0, str(Path(__file__).resolve().parent))


def matrix_review():
    base, correction = fixtures()
    worst, nilpotence, cases = 0., 0., 0
    for h, xi, eta, nu in itertools.product((.005, .2), (.1, .5, 1.2), (-.7, 0., .6), (.4, 1., 1.7)):
        a0, a1, f = matrices(base, xi, eta, h, nu)
        w, wx, we = correction_vector(correction, xi, eta)
        defect = wx+DIAGONAL*w/xi-a0@w-a1@we-f
        coeff = series_coefficients([base, correction], xi*xi, eta, h, nu)
        r, t, z = [coeff[k][1] for k in ('radial', 'theta', 'axial')]
        L = 1-2*h*eta*eta
        expected = np.array([0., 0., 0., r/xi, -2*t/nu, -2*(z+eta*r/L)/nu])
        worst = max(worst, float(np.max(np.abs(defect-expected)/np.maximum(1., np.abs(expected)))))
        diagonal = np.diag([.1, .3, .7, 1., 2., 3.])
        right = matrices(base, xi/2, eta/2, h, nu)[1]
        nilpotence = max(nilpotence, float(np.max(np.abs(a1@diagonal@right))))
        cases += 1
    xi, eta, h, nu = .7, .6, .005, 1.
    a0, a1, f = matrices(base, xi, eta, h, nu)
    w, wx, we = correction_vector(correction, xi, eta)
    bad = a0.copy(); bad[5, 3] = 0
    pressure_gap = float(np.max(np.abs(a0@w-bad@w)))
    # Changing h also changes the leading operators; instead remove the exact
    # positive-order lambda terms from the independent kinematic residual.
    from residual import profile_jets
    actual_v = profile_jets(correction, xi*xi, eta, h, 2*h)['V'].v
    frozen_v = profile_jets(correction, xi*xi, eta, h, 0.)['V'].v
    derivative_gap = float(np.max(np.abs(a1@we)))
    parity_error = 0.
    J = np.diag([1., 1., 1., 1., -1., -1.])
    # The evaluators take X=xi^2; odd coefficients explicitly contain xi.
    # Their formula extends to negative xi even though the public validation
    # deliberately restricts radial coordinates to nonnegative values.
    odd_entries = ((3, 0), (4, 4), (5, 5))
    a0negative = a0.copy()
    for i, j in odd_entries:
        a0negative[i, j] *= -1
    fnegative = f.copy(); fnegative[:4] *= -1
    parity_error = max(float(np.max(np.abs(a0negative+J@a0@J))),
                       float(np.max(np.abs(a1+J@a1@J))),
                       float(np.max(np.abs(fnegative+J@f))))
    return dict(manufactured=True, sample_count=cases, maximum_relative_error=worst,
                nilpotence_error=nilpotence, parity_error=parity_error,
                negative_controls=dict(pressure_coupling=pressure_gap,
                                       positive_order_lambda=abs(actual_v-frozen_v),
                                       missing_angular_derivatives=derivative_gap),
                checks=dict(independent_full_residual=worst < 2e-12,
                            derivative_block_square_zero=nilpotence == 0.,
                            radial_parity=parity_error == 0.,
                            all_missing_term_controls_detected=min(pressure_gap, abs(actual_v-frozen_v), derivative_gap) > 1e-5))


def sparse_review():
    rows = []
    for k in range(11):
        words = [word for word in itertools.product((0, 1), repeat=k) if allowed_word(word)]
        rows.append(dict(length=k, allowed_count=len(words), maximum_derivatives=max(map(sum, words))))
    samples = []
    with mp.workdps(80):
        for C, a, delta in ((1, '.01', '.1'), (10, '.02', '.03'), (100, '.1', '.2')):
            bound = log_sum_bound(C, a, delta)
            total = mp.fsum(mp.exp(log_term(k, C, a, delta)) for k in range(401))
            samples.append(dict(C=str(C), a=a, strip_loss=delta,
                                partial_log_sum=mp.nstr(mp.log(total), 25),
                                log_sum_bound=mp.nstr(bound, 25), below_bound=mp.log(total) <= bound))
        roots = [dict(k=k, log_kth_root=mp.nstr(log_term(k, 10, '.02', '.03')/k, 25))
                 for k in (20, 100, 500, 2000)]
    return dict(exact_word_enumeration=rows, finite_sum_checks=samples, term_roots=roots,
                checks=dict(sparse_derivative_count=all(row['maximum_derivatives'] == (row['length']+1)//2 for row in rows),
                            sum_majorant_controls_samples=all(row['below_bound'] for row in samples),
                            kth_root_falls=mp.mpf(roots[-1]['log_kth_root']) < mp.mpf(roots[0]['log_kth_root'])-1))


def recurrence_review():
    base, _ = fixtures()
    rows, slope_error = [], 0.
    with mp.workdps(60):
        for degree in (3, 5, 7):
            correction, budget = build(base, .005, 1., degree)
            errors = []
            for X, eta in ((.01, .2), (.006, -.3), (.004, .4)):
                coeff = series_coefficients([base, correction], X, eta, .005)
                errors.append(max(abs(coeff[key][1]) for key in ('radial', 'theta', 'axial')))
                slopes = axis_slopes(base, eta, .005)
                for name, value in slopes.items():
                    actual = getattr(correction, name.split('1_')[0]).derivative(0)(0., eta)
                    slope_error = max(slope_error, abs(actual-value))
            rows.append(dict(degree=degree, errors=list(map(float, errors)), budget=budget))
    refinement = max(rows[0]['errors'])/max(rows[-1]['errors'])
    return dict(manufactured=True, rows=rows, axis_slope_error=float(slope_error),
                passed_refinement_lower_bound=1000000,
                checks=dict(axis_slopes_match_independent_operators=slope_error < 2e-12,
                            final_coefficient_residual_small=max(rows[-1]['errors']) < 2e-11,
                            radial_refinement=refinement > 1e6,
                            actual_large_parameter_profile_not_evaluated=all(not row['budget']['actual_large_parameter_profile_evaluated'] for row in rows)))


def first_moment_review():
    rows = []
    with mp.workdps(80):
        target = mp.matrix(list(map(mp.mpf, ('.3', '-.4', '.2', '.1', '-.2'))))
        for lam in ('.0001', '.001', '.01', '1e-50'):
            repair = Repair(lam)
            c = repair.inverse*target
            m, i, j, s, cp = [mp.mpf(0) for _ in range(5)]
            for index in range(5):
                for z, w in repair.rule.nodes:
                    y = repair.centers[index]+repair.width*(z-mp.mpf('.5'))
                    bump = step_prime(z)/repair.width
                    du, de = (c[index]*bump, 0) if index < 2 else (0, c[index]*bump)
                    e0 = mp.exp((-mp.mpf('.5')-repair.lam)*y)
                    weight = repair.width*w
                    m += weight*mp.exp(y)*du
                    i += weight*mp.exp(mp.mpf('1.5')*y)*de
                    j += weight*mp.exp(mp.mpf('1.5')*y)*e0*du
                    s -= weight*mp.exp(y)*e0*de
                    cp += weight*e0*de
            ordinary = mp.matrix([m, i, j, s, cp])
            normalized = precondition_target(ordinary, repair.lam)
            err = max(abs(v) for v in normalized-target)
            algebra = max(abs(v) for v in repair.B*c-target)
            nonlinear_gap = max(abs(v) for v in repair.value(c)-target)
            inverse_norm = max(mp.fsum(abs(repair.inverse[k, l]) for l in range(5)) for k in range(5))
            rows.append(dict(lambda_value=lam, linear_matrix_error=mp.nstr(algebra, 15),
                             independent_integrand_error=mp.nstr(err, 15),
                             wrong_nonlinear_map_gap=mp.nstr(nonlinear_gap, 15),
                             numerical_inverse_norm=mp.nstr(inverse_norm, 15),
                             coefficients_finite=all(mp.isfinite(v) for v in c)))
    return dict(manufactured_power_law_patch=True, positive_order_patch='I_3',
                quadrature_is_not_a_continuum_certificate=True,
                actual_annular_support_hypotheses_fully_verified=False,
                finite_target_need_not_be_small=True,
                quadratic_leading_repair_must_not_be_reused=True, rows=rows,
                checks=dict(linear_matrix_inversion=all(mp.mpf(row['linear_matrix_error']) < mp.mpf('1e-65') for row in rows),
                            independent_linear_moments=all(mp.mpf(row['independent_integrand_error']) < mp.mpf('1e-25') for row in rows),
                            inverse_below_inherited_analytic_cap=all(mp.mpf(row['numerical_inverse_norm']) < 1000 for row in rows),
                            leading_nonlinear_map_rejected=all(mp.mpf(row['wrong_nonlinear_map_gap']) > mp.mpf('.01') for row in rows),
                            arbitrary_target_has_finite_solution=all(row['coefficients_finite'] for row in rows)))


def implicit_review():
    # Independent raw third differentiation of B*c+q(eta)*c^2=z(eta).
    with mp.workdps(80):
        c = lambda t: mp.sin(t)+t*t/5+mp.mpf('.01')
        q = lambda t: mp.exp(t)/7
        B, eta = mp.mpf(2), mp.mpf('.37')
        z = lambda t: B*c(t)+q(t)*c(t)**2
        cc = [mp.diff(c, eta, k) for k in range(4)]
        qq = [mp.diff(q, eta, k) for k in range(4)]
        L = B+2*qq[0]*cc[0]
        terms = [6*qq[0]*cc[1]*cc[2], 6*qq[1]*cc[0]*cc[2],
                 6*qq[1]*cc[1]**2, 6*qq[2]*cc[0]*cc[1], qq[3]*cc[0]**2]
        exact = mp.diff(z, eta, 3)-mp.fsum(terms)
        error = abs(L*cc[3]-exact)
        gaps = [abs(value) for value in terms]
    return dict(identity_error=mp.nstr(error, 15), omitted_term_gaps=[mp.nstr(v, 15) for v in gaps],
                checks=dict(raw_third_product_rule=error < mp.mpf('1e-65'),
                            each_missing_term_detected=min(gaps) > mp.mpf('1e-4')))
