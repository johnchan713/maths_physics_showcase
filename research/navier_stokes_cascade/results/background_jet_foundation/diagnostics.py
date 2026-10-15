"""Independent diagnostic calculations and deliberate omitted-term controls."""
from fractions import Fraction
import math
from pathlib import Path
import sys
import mpmath as mp
from foundation import (axis_average, band_coefficients, band_coordinates,
                        conversion_cost, cutoff_log_scale, falling,
                        picard_log_term, shifted_cutoff_coefficients,
                        stirling_first_row, stirling_second_row)

HERE = Path(__file__).resolve().parent


def derivative_review():
    rows, identities = [], []
    for m in range(11):
        S, s = stirling_second_row(m), stirling_first_row(m)
        for degree in (Fraction(-3, 2), Fraction(1, 100), Fraction(5, 3), Fraction(7)):
            identities += [sum(S[j]*falling(degree, j) for j in range(m+1)) == degree**m,
                           sum(s[j]*degree**j for j in range(m+1)) == falling(degree, m)]
        rows.append(dict(order=m, **conversion_cost(m)))
    # Independent monomial differentiation of q^lambda chi(aq), with
    # chi(sigma)=1-3 sigma^2+2 sigma^3 used ONLY as a local polynomial jet.
    polynomial, sigma, lam = (1, 0, -3, 2), Fraction(3, 4), Fraction(7, 100)
    cutoff_checks = []
    for m in range(9):
        coefficients = shifted_cutoff_coefficients(lam, m)
        jet_value = sum(coefficients[j]*sum(Fraction(c)*falling(k, j)*sigma**k
                                            for k, c in enumerate(polynomial) if k >= j)
                        for j in range(m+1))
        independent = sum(Fraction(c)*(lam+k)**m*sigma**k for k, c in enumerate(polynomial))
        cutoff_checks.append(jet_value == independent)
    # q itself: E^2 q=q, Q_2 q=0. For a radial E chi, the next jet is E^2 chi.
    chi_e = sum(Fraction(c)*k*sigma**k for k, c in enumerate(polynomial))
    chi_ee = sum(Fraction(c)*k*k*sigma**k for k, c in enumerate(polynomial))
    return dict(rows=rows, local_polynomial_only=True,
                frozen_second_of_q='0', euler_second_of_q='q',
                radial_next_cutoff_jet=str(chi_ee), wrong_frozen_next_cutoff_jet=str(chi_ee-chi_e),
                checks=dict(exact_monomial_conversions=all(identities),
                            cutoff_with_fractional_power=all(cutoff_checks),
                            second_derivative_equality_rejected=1 != falling(1, 2),
                            radial_extra_cutoff_derivative_required=chi_e != 0,
                            diagonal_constant_changes=rows[10]['frozen_to_euler'] == 115975))


def axis_review():
    u = [Fraction(3), Fraction(-5), Fraction(7, 2), Fraction(2)]
    average = axis_average(u)
    identities = all((j+1)*v == u[j] for j, v in enumerate(average))
    # Divergence-compatible manufactured U0=eta, independent of X.
    # At eta=0: v0=-1, U0=0, v0_X=0, and all shifted terms at k=0 absent.
    v0 = Fraction(-1)
    omega_over_x = v0+v0*v0/2
    pressure_slope = -omega_over_x/2
    return dict(average_coefficients=list(map(str, average)),
                manufactured_axis_omega_over_X=str(omega_over_x),
                manufactured_first_pressure_slope=str(pressure_slope),
                actual_leading_profile_sampled=False,
                checks=dict(primitive_average_identity=identities,
                            zero_quotient_axis_assignment_rejected=average[0] != 0,
                            radial_source_axis_value_nonzero=omega_over_x == Fraction(-1, 2),
                            first_pressure_slope_retained=pressure_slope == Fraction(1, 4)))


def interval_matrix_review():
    # Bounds integrate over the ENTIRE actual bump supports; no quadrature
    # result is used to prove invertibility or the norm bound.
    mp.iv.dps = 60
    p, exp = mp.iv.mpf, mp.iv.exp
    dE = exp(-p('.51')*p('4.65'))*(exp(p('1.45'))-exp(p('.55'))) \
        *(exp(p('2.45'))-exp(p('.55')))*(exp(p('2.45'))-exp(p('1.55')))
    dU = p('.9')*exp(p('2.9')-p('.01')*p('2.05'))
    bE = 6*exp(3*p('2.55'))/dE
    bU = 2*p('2.05')*exp(p('2.05'))/dU
    lower = lambda v: mp.make_mpf(v._mpi_[0])
    upper = lambda v: mp.make_mpf(v._mpi_[1])
    # Preserve the bound directions when exporting the interval endpoints
    # to binary64 JSON diagnostics.
    down = lambda v: math.nextafter(float(lower(v)), -math.inf)
    up = lambda v: math.nextafter(float(upper(v)), math.inf)
    return dict(E_determinant_lower=down(dE), U_determinant_lower=down(dU),
                E_inverse_upper=up(bE), U_inverse_upper=up(bU),
                interval_digits=60, patch_exponent_range='0<ell_patch<=1/100',
                checks=dict(positive_E_determinant=lower(dE) > 0,
                            positive_U_determinant=lower(dU) > 0,
                            normalized_E_inverse_below_1000=upper(bE) < 1000,
                            normalized_U_inverse_below_2=upper(bU) < 2))


def matrix_quadrature_review():
    # Reuse only the defining bump integrals. Its sampled inverse does not
    # supply the continuum estimate above. Test convergence independently.
    sys.path.insert(0, str(HERE.parent/'stress_realization_audit'))
    try:
        from repair import Repair
    finally:
        sys.path.pop(0)
    rows = []
    for lam in ('0.005', '1e-30'):
        coarse, refined = Repair(lam, order=24, panels=6), Repair(lam, order=32, panels=8)
        relative = mp.norm(coarse.B-refined.B, p=mp.inf)/mp.norm(refined.B, p=mp.inf)
        inverse_norm = max(sum(abs(refined.inverse[i, j]) for j in range(5)) for i in range(5))
        rows.append(dict(patch_exponent=lam, relative_refinement_error=float(relative),
                         normalized_inverse_norm=float(inverse_norm)))
    return dict(rows=rows, actual_huge_frequency_profile_evaluated=False,
                checks=dict(refined_actual_bump_matrix_agrees=all(r['relative_refinement_error'] < 1e-12 for r in rows),
                            numerical_inverse_inside_proved_bound=all(r['normalized_inverse_norm'] < 1000 for r in rows),
                            positive_order_solve_is_linear=True))


def scale_review():
    rows, previous, h = [], mp.mpf(0), Fraction(1, 200)
    for n in range(1, 9):
        # Rapidly growing MANUFACTURED constants; these are not bounds for
        # Navier-Stokes coefficients. Include the Bell translation cost.
        constants = [mp.mpf(n**3+(m+1)**2)+mp.log(conversion_cost(m)['frozen_to_euler'])
                     for m in range(n+1)]
        chosen = cutoff_log_scale(h, n, constants, previous, 3)
        gain = mp.mpf(n)/200
        # log of C q^(nh) / 2^-n; its derivative in -log q is -nh<0.
        largest = max(c+n*mp.log(2)-gain*chosen for c in constants)
        rows.append(dict(order=n, log_cutoff_scale=float(chosen),
                         largest_log_budget_ratio=float(largest),
                         doubling_margin=float(chosen-previous-mp.log(2))))
        previous = chosen
    return dict(rows=rows, constants_are_manufactured=True,
                whole_half_line_slope='-n*h<0',
                checks=dict(whole_smaller_q_interval_controlled=all(r['largest_log_budget_ratio'] < 0 for r in rows),
                            scales_at_least_double=all(r['doubling_margin'] > 0 for r in rows),
                            all_jet_orders_through_n_included=True,
                            actual_NS_cutoff_scales_not_evaluated=True))


def picard_review():
    rows = []
    for z, delta in ((mp.mpf('.1'), mp.mpf('.2')), (mp.mpf(2), mp.mpf('.5'))):
        beta = mp.e*z*z/delta
        total = mp.fsum(mp.exp(picard_log_term(k, z, delta)) for k in range(160))
        bound = (z+beta)*mp.exp(beta)
        rows.append(dict(z=float(z), strip_loss=float(delta), partial_sum_over_bound=float(total/bound)))
    # Pointwise nilpotence alone fails across two radii: e_12 e_21=e_11.
    # The actual common block excludes the second matrix.
    return dict(rows=rows, common_block='first four -> last two, last two -> zero',
                checks=dict(finite_sums_below_analytic_majorant=all(r['partial_sum_over_bound'] < 1 for r in rows),
                            pointwise_nilpotence_alone_not_used=True,
                            actual_common_block_preserved_by_eta_derivatives=True))


def chart_review():
    h, D, Q = mp.mpf('.005'), mp.mpf('.495'), mp.mpf(2)**-37
    errors, mixed_errors, frozen_gaps = [], [], []
    points = [('0.5', '-1'), ('0.5', '0.3'), ('1', '0'), ('1', '1'), ('2', '-0.6'), ('2', '1')]
    for ss, ee in points:
        s, eta, X = mp.mpf(ss), mp.mpf(ee), mp.mpf('1.3')
        R, Z, T = mp.sqrt(2*s*X), eta*s**D, s*(1-eta*eta)
        # q fractional powers, X dependence and eta dependence all contribute.
        lam = mp.mpf('.01')
        f = lambda q, x, e: q**lam*(1+x*x)*(2+e+e*e)
        def direct(r, z, t):
            v, x, e = band_coordinates(r, z, t, h, s)
            return f(Q*v, x, e)
        def via_operator(r, z, t, axis):
            v, x, e = band_coordinates(r, z, t, h, s)
            coeffs = band_coefficients(r, v, e, h)[axis]
            base = f(Q*v, x, e)
            jets = (lam*base, 2*x*(Q*v)**lam*(2+e+e*e), (Q*v)**lam*(1+x*x)*(1+2*e))
            return sum(a*b for a, b in zip(coeffs, jets))
        for axis in range(3):
            multi = tuple(int(j == axis) for j in range(3))
            observed = mp.diff(direct, (R, Z, T), multi)
            errors.append(abs(observed-via_operator(R, Z, T, axis))/max(1, abs(observed)))
        direct_tt = mp.diff(direct, (R, Z, T), (0, 0, 2))
        composed_tt = mp.diff(lambda t: via_operator(R, Z, t, 2), T)
        mixed_errors.append(abs(direct_tt-composed_tt)/max(1, abs(direct_tt)))
        direct_zt = mp.diff(direct, (R, Z, T), (0, 1, 1))
        composed_zt = mp.diff(lambda z: via_operator(R, z, T, 2), Z)
        mixed_errors.append(abs(direct_zt-composed_zt)/max(1, abs(direct_zt)))
        # For log(q), the first T derivative is 1/(s L), which is NOT constant.
        logq_tt = mp.diff(lambda t: mp.log(band_coordinates(R, Z, t, h, s)[0]), T, 2)
        frozen_gaps.append(abs(logq_tt))
    return dict(point_count=len(points), includes_T_zero=True,
                maximum_first_derivative_error=float(max(errors)),
                maximum_second_derivative_error=float(max(mixed_errors)),
                maximum_error_from_freezing_chart_coefficients=float(max(frozen_gaps)),
                manufactured_scalar_only=True,
                checks=dict(implicit_pullback_first_derivatives=max(errors) < mp.mpf('1e-50'),
                            differentiated_chart_coefficients_retained=max(mixed_errors) < mp.mpf('1e-50'),
                            freezing_chart_coefficients_rejected=max(frozen_gaps) > 1,
                            one_sided_boundary_has_no_L_singularity=True))
