"""Support, Fourier and deliberately failed interaction controls.

These are manufactured examples. They do not materialize or independently
certify the actual global auxiliary rectangles and coloring.
"""
import math
import numpy as np
from scipy.integrate import quad


def bump(x):
    return math.exp(-1/(1-x*x)) if abs(x) < 1 else 0.


def bump_derivative(x):
    return -2*x*math.exp(-1/(1-x*x))/(1-x*x)**2 if abs(x) < 1 else 0.


def support_and_flat_review():
    width = .18
    mass = width*quad(bump, -1, 1, epsabs=1e-14, epsrel=1e-13)[0]
    mean = mass*mass
    chi = lambda x, y: bump((x-.5)/width)*bump((y-.5)/width)
    n1, n2 = math.sqrt(2)-1, 1.
    source = lambda x, y: (n1*bump_derivative((x-.5)/width)*bump((y-.5)/width)
                           +n2*bump((x-.5)/width)*bump_derivative((y-.5)/width))/width
    outside = (chi(.1, .1)-mean)
    source_outside = source(.1, .1)
    # A finite Fourier coefficient check of N^-1(N chi)=chi-mean.
    fourier_error = 0.
    for k1, k2 in ((1, 0), (0, 1), (1, -2), (3, 1)):
        one = []
        for k in (k1, k2):
            one.append(width*quad(lambda s: bump(s)*math.cos(2*math.pi*k*(.5+width*s)), -1, 1,
                                  epsabs=1e-13, epsrel=1e-12)[0]
                       -1j*width*quad(lambda s: bump(s)*math.sin(2*math.pi*k*(.5+width*s)), -1, 1,
                                     epsabs=1e-13, epsrel=1e-12)[0])
        coefficient = one[0]*one[1]
        n_coefficient = 2j*math.pi*(n1*k1+n2*k2)*coefficient
        recovered = n_coefficient/(2j*math.pi*(n1*k1+n2*k2))
        fourier_error = max(fourier_error, abs(recovered-coefficient))
    product_outside = (chi(.1, .1)-mean)*chi(.1, .1)
    disjoint = max(abs(bump((x-.2)/.1)*bump((x-.8)/.1)) for x in np.linspace(0, 1, 201))
    # A zero-mean correction can contribute a nonzero mean in a later
    # quadratic product. Even a flat factor must stay in the actual field.
    flat_rows = []
    y = np.arange(65)/65
    for ell in (1, 3, 8):
        eta = math.exp(-ell*ell)
        correction = eta*np.cos(2*math.pi*y)
        other = np.cos(2*math.pi*y)
        flat_rows.append(dict(ell=ell, flat_factor=eta, mean_of_old_correction=float(np.mean(correction)),
                              mean_after_product=float(np.mean(correction*other)), expected_mean=eta/2))
    # Finite angular sets grow by addition; means do not introduce phases.
    modes = {-1, 1}
    harmonic_rows = []
    for stage in range(7):
        residual = {a+b for a in modes|{0} for b in modes|{0}}-{0}
        harmonic_rows.append(dict(stage=stage, maximum_residual_harmonic=max(map(abs, residual)),
                                  count=len(residual)))
        modes |= residual|{-1, 1}
    return dict(auxiliary_bump_mean=mean, temporal_inverse_outside_original_rectangle=outside,
                temporal_source_outside_original_rectangle=source_outside,
                maximum_finite_fourier_inverse_error=fourier_error,
                mean_times_wave_outside_wave_support=product_outside,
                separated_closed_support_product=disjoint,
                flat_product_rows=flat_rows, finite_harmonic_rows=harmonic_rows,
                checks=dict(temporal_inverse_can_fill_auxiliary_torus=abs(outside) > .001 and source_outside == 0.,
                            temporal_inverse_fourier_identity=fourier_error < 1e-12,
                            mean_wave_product_keeps_wave_support=product_outside == 0.,
                            separated_support_products_vanish=disjoint == 0.,
                            zero_mean_flat_velocity_is_retained_in_later_products=all(abs(r['mean_after_product']-r['expected_mean']) < 1e-12*r['flat_factor'] for r in flat_rows),
                            angular_harmonics_remain_finite_and_band_independent=all(r['maximum_residual_harmonic'] <= 2**(r['stage']+1) for r in harmonic_rows)))


def sharp_loss_controls():
    """Cartesian examples isolating losses; neither is the actual band flow."""
    rows = []
    for k in (4, 16, 64):
        # Distinct overlapping phases: w1=cos(kx)e_y, w2=cos(ky)e_x.
        # Each is exactly divergence-free, but cross transport is order k.
        nodes = np.arange(8*k+1)*2*math.pi/(8*k+1)
        sine, cosine = np.sin(k*nodes), np.cos(k*nodes)
        squared = k*k*(np.outer(cosine*cosine, sine*sine)+np.outer(sine*sine, cosine*cosine))
        rms = float(np.sqrt(np.mean(squared)))
        # A constant axial mean advects w=epsilon^alpha cos(kz)e_x.
        # epsilon=k^-2, so the coefficient has order alpha+mu-1/2.
        eps, alpha, mu = k**-2, .7, .9
        mean_wave = k*eps**(alpha+mu)
        rows.append(dict(k=k, overlapping_distinct_label_cross_rms=rms,
                         cross_rms_over_k=rms/k,
                         mean_wave_coefficient=mean_wave,
                         mean_wave_exact_scale_ratio=mean_wave/eps**(alpha+mu-.5),
                         omitted_phase_derivative_scale_ratio=mean_wave/eps**(alpha+mu)))
    return dict(rows=rows, checks=dict(
        closed_support_separation_is_needed=rows[-1]['overlapping_distinct_label_cross_rms'] > 15*rows[0]['overlapping_distinct_label_cross_rms'],
        same_label_incompressibility_cannot_be_used_for_different_phases=all(abs(r['cross_rms_over_k']-1/math.sqrt(2)) < 1e-12 for r in rows),
        mean_wave_half_power_loss_is_sharp=all(abs(r['mean_wave_exact_scale_ratio']-1) < 1e-12 for r in rows),
        dropping_mean_advection_of_phase_is_invalid=rows[-1]['omitted_phase_derivative_scale_ratio'] > 60))
