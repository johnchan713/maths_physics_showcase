"""Derivative controls independent of the analytic supremum estimates.

Manufactured functions test formulas. They are never passed off as values of
the actual core or its five incoming moments.
"""
from fractions import Fraction
import mpmath as mp
from jets import Jet2, physical_to_rows


def normalization_review(digits=80):
    with mp.workdps(digits):
        pressure, radius = mp.mpf(16), mp.mpf(10000)

        def defects(e):
            return [(1+e+e*e)/32, (2-e+e**3)/40, (1+2*e+e*e-e**3)/27,
                    (3-e+2*e*e)/19, (1+e*e)/23]

        def scalar_rows(e):
            # Assemble scalars independently; mp.diff differentiates everything.
            M, I, J, S, Cp = defects(e)
            xc, f = mp.exp(-6), 1/(1+e*e)
            angular = mp.sqrt(2)*radius**mp.mpf('1.5')*pressure*f*xc**mp.mpf('1.6')
            return [M/(radius*xc), (J-4*e*I)/angular, I/angular,
                    -(S-8*e*M)/(radius*pressure**2*f*f*xc**mp.mpf('1.2')),
                    Cp/(pressure**2*f*f*xc**mp.mpf('.2'))]

        error = mp.mpf(0)
        for e in map(mp.mpf, ['-1', '-.37', '0', '.5', '1']):
            inputs = [Jet2(defects(e)[j], mp.diff(lambda t: defects(t)[j], e),
                           mp.diff(lambda t: defects(t)[j], e, 2)) for j in range(5)]
            rows = physical_to_rows(inputs, e, pressure, radius)
            for j, row in enumerate(rows):
                expected = [scalar_rows(e)[j], mp.diff(lambda t: scalar_rows(t)[j], e),
                            mp.diff(lambda t: scalar_rows(t)[j], e, 2)]
                error = max(error, *(abs(a-b) for a, b in zip((row.v, row.d, row.dd), expected)))
        # At eta=0, constant Cp has zero input derivatives. Its NORMALIZED
        # second derivative is nevertheless nonzero because f^-2 varies.
        inputs = [Jet2(mp.mpf(0)) for _ in range(4)]+[Jet2(mp.mpf(1))]
        correct = physical_to_rows(inputs, mp.mpf(0), pressure, radius)[4]
        frozen = physical_to_rows(inputs, mp.mpf(0), pressure, radius, True)[4]
        e = Jet2(Fraction(1), Fraction(1))
        inverse_f_squared = (1+e*e)*(1+e*e)
        actual_C2 = abs(inverse_f_squared.v)+abs(inverse_f_squared.d)+abs(inverse_f_squared.dd)/2
        return dict(manufactured_inputs=True, maximum_jet_oracle_error=error,
                    constant_pressure_row_raw_second=correct.dd,
                    frozen_normalization_raw_second=frozen.dd,
                    inverse_f_squared_exact_C2=actual_C2,
                    checks=dict(all_five_rows_keep_second_derivatives=error < mp.mpf('1e-65'),
                                frozen_normalization_failure_detected=correct.dd > mp.mpf('.01') and frozen.dd == 0,
                                old_factor_12_rejected=actual_C2 == 20 and actual_C2 > 12))


def cancellation_review():
    e = Jet2(Fraction(3, 10), Fraction(1))
    g = Fraction(1, 10**16)*(1+e+e*e)
    U = 4*e+g
    actual = U*U-16*e*e-8*e*(U-4*e)
    expected = g*g
    missing_cross = U*U-16*e*e
    return dict(manufactured_inputs=True, actual=actual, expected=expected,
                checks=dict(energy_row_cancellation_through_C2=actual == expected,
                            omitted_cross_term_detected=missing_cross.dd != expected.dd))


def exponential_review(digits=80):
    with mp.workdps(digits):
        eta, cutoff = mp.mpf('.3'), mp.mpf('.37')
        e = Jet2(eta, 1)
        logCE = mp.mpf('.3')+(1-cutoff)*(e+e*e)-cutoff*(1+e*e).log()
        value = logCE.exp()

        def scalar(t):
            return mp.exp(mp.mpf('.3')+(1-cutoff)*(t+t*t)-cutoff*mp.log(1+t*t))

        expected = mp.diff(scalar, eta, 2)
        omitted_square = value.v*logCE.dd
        return dict(manufactured_shape_input=True, correct_error=abs(value.dd-expected),
                    missing_square_error=abs(omitted_square-expected),
                    checks=dict(exponentiation_keeps_second_derivative=abs(value.dd-expected) < mp.mpf('1e-65'),
                                missing_first_derivative_square_detected=abs(omitted_square-expected) > mp.mpf('.1')))


def axis_integrability_review():
    # Toy E=sqrt(2X)*A(eta), so E^2/(2X)=A^2, including angular jets.
    e = Jet2(Fraction(1, 3), Fraction(1))
    amplitude = 1+e*e
    X = Fraction(1, 7)
    density = (2*X*amplitude*amplitude)/(2*X)
    # A uniform E<=1 estimate alone gives integral_{exp(-k)}^1 dX/(2X)=k/2.
    bad_bounds = [Fraction(k, 2) for k in (10, 100, 1000)]
    return dict(manufactured_axis=True, uniform_only_logarithmic_bounds=bad_bounds,
                checks=dict(axis_factor_cancels_in_all_three_jets=density == amplitude*amplitude,
                            dropping_axis_factor_has_no_finite_integral_bound=bad_bounds[0] < bad_bounds[1] < bad_bounds[2]))
