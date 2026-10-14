"""Positive-order matching with the lower-order pressure source retained.

These finite-dimensional helpers implement the analytic construction in the
README. Numerical callers supply manufactured data, not the huge-N profile.
"""
from fractions import Fraction
from pathlib import Path
import sys
import mpmath as mp

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
sys.path.insert(0, str(HERE.parent / 'stress_realization_audit'))
from repair import Repair
sys.path.insert(0, str(HERE))


def normalized_target(moments, radius, amplitude, lam):
    """Map ordinary (m1,m2,m3,m4,m5) to the inherited B row order.

    amplitude is E0 at the LEFT endpoint of I3, not a constant independent
    of eta. The divided difference still costs lam**-1 for arbitrary data.
    """
    if len(moments) != 5:
        raise ValueError('Exactly five ordinary discrepancies are required')
    radius, amplitude, lam = map(mp.mpf, (radius, amplitude, lam))
    if not all(mp.isfinite(v) for v in (radius, amplitude, lam)) or min(radius, amplitude, lam) <= 0:
        raise ValueError('Positive finite patch radius, amplitude and slope required')
    values = list(map(mp.mpf, moments))
    if not all(mp.isfinite(v) for v in values):
        raise ValueError('Moment discrepancies must be finite')
    angular_scale = mp.sqrt(2)*radius**mp.mpf('1.5')
    m, i, cp, j, s = (values[0]/radius, values[1]/angular_scale,
                      values[2]/amplitude, values[3]/(angular_scale*amplitude),
                      values[4]/(radius*amplitude))
    return mp.matrix([m, (m-j)/lam, i, -s, cp])


def solve_patch(moments, radius, amplitude, repair):
    """The FIRST coefficient uses B*c=-z, never B*c+Q(c,c)=-z."""
    return -repair.inverse*normalized_target(moments, radius, amplitude, repair.lam)


def radial_flux(U, M, M_eta, X, eta, h, order=1):
    """Reconstruct V_n from a streamfunction moment, retaining D+2*n*h."""
    D, lam = Fraction(1, 2)-h, 2*order*h
    return (2*eta*X*U-2*(D+lam)*eta*M-(1-eta*eta)*M_eta)/(1-2*h*eta*eta)


def divergence_identity(U, U_X, U_eta, X, eta, h, order=1, omit_shift=False):
    """Differentiate V using M_X=U, then add the independent Z operator."""
    A, D, lam = Fraction(1, 2)+h, Fraction(1, 2)-h, 2*order*h
    L, d = 1-2*h*eta*eta, 1-eta*eta
    shift = 0 if omit_shift else lam
    V_X = (2*eta*U+2*eta*X*U_X-2*(D+shift)*eta*U-d*U_eta)/L
    Z_U = (2*(-A+lam)*eta*U-2*eta*X*U_X+d*U_eta)/L
    return V_X+Z_U


def geometry_ledger():
    """Exact offsets, not floating-point additions to log(Xa)."""
    offsets = [Fraction(1, 16), Fraction(1, 8), Fraction(3, 16), Fraction(1, 4), Fraction(1, 2)]
    return dict(inner_offsets_in_units_tc=offsets,
                offset_names=['stress_free_through', 'keep', 'cut', 'analytic_endpoint', 'leading_modulation_start'],
                positive_patch_stage_offsets=[-14, -9], mean_patch_stage_offsets=[-8, -3],
                positive_patch_radius='rho3=XR*exp(Td+Tw-12)',
                actual_bump_log_supports=[['.45', '.55'], ['.95', '1.05'],
                                         ['1.45', '1.55'], ['1.95', '2.05'], ['2.45', '2.55']],
                global_pressure_support_edge='Xv (after all reserved patches and the axial pulse)',
                checks=dict(positive_and_strict_inner_order=all(0 < a < b for a, b in zip(offsets, offsets[1:])),
                            extension_cutoff_outside_stress_free_core=offsets[0] > 0,
                            cutoff_before_protected_endpoint=offsets[2] < offsets[3],
                            protected_endpoint_before_modulation=offsets[3] < offsets[4],
                            patch_separation=-9 < -8))


def construction_ledger():
    return dict(dependencies=['inherited leading analytic-eta collar', 'inherited leading exact moments and heat exterior',
                              'inherited continuum inverse bound for B', 'first six-variable system and sparse Picard proof'],
                collar_operator_bound='C_col=max(1,sup_Sr ||A0||,sup_Sr ||A1||,sup_Sr ||f1||)<infinity',
                collar_bound_is_C16=False, collar_bound_uniform_in_large_parameters=False,
                collar_endpoint='a^2=Xa*exp(tc/4)',
                sum_bound='(z+beta)*exp(beta), z=2*C_col*a, beta=e*z^2/Delta',
                same_inner_coefficient_by_analytic_uniqueness=True,
                cutoff_fields=['F1', 'U1'], integrated_fields=['M1', 'V1', 'Pi1'],
                zero_axis_pressure_retained=True,
                pressure_equation='Pi1_X=2*F0*F1-Omega0/(2*X)',
                leading_patch_U0_zero=True, leading_patch_M0_zero_required=False,
                leading_patch_Omega0_zero_required=False,
                velocities_vanish_on_Imean=True, pressure_vanishes_on_Imean_not_asserted=True,
                coefficient_bounds='finite for every fixed profile derivative; not numerically materialized',
                stress_pair_is_full_wave_covariance=False,
                background_stress_support='[Xa*exp(tc/16),Xb]',
                checks=dict(same_frequency=True, no_small_target_assumption=True,
                            no_zero_memory_shortcut=True, pressure_not_cut_off_directly=True))


def residual_gain_ledger(h=Fraction(1, 200), derivative_order=0):
    """One relative gain for the augmented BACKGROUND, not a smooth force."""
    h = Fraction(h)
    if not 0 < h < Fraction(1, 100) or type(derivative_order) is not int or derivative_order < 0:
        raise ValueError('Use the construction range and a nonnegative derivative order')
    return dict(relative_gain=2*h, radial_physical_exponent=-Fraction(3, 2)+2*h,
                tangential_physical_exponent=-Fraction(3, 2)+3*h,
                sufficient_all_component_exponent=-Fraction(3, 2)+2*h-derivative_order,
                constants='finite C_m on each fixed compact profile range; not evaluated',
                cancels_first_order_by_global_pressure_and_retained_stress=True,
                error_is_full_wave_or_localized_residual=False,
                absolute_error_proved_to_decay=False)
