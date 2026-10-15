"""Exact order bookkeeping for the CONDITIONAL finite correction cycle.

These computations check the inequalities after the analytic operator and
interaction estimates in README.md. Rational arithmetic does not prove
those PDE estimates or the inherited leading profile and support geometry.
"""
from dataclasses import dataclass
from fractions import Fraction as F

KAPPA = F(1, 100000)
MIN_B, STEP, RHO, MU = F(7, 10), F(1, 10), F(17, 25), F(9, 10)


@dataclass(frozen=True)
class Affine:
    """slope*B+constant, on the WHOLE half-line B>=0.7."""
    slope: F
    constant: F

    def at(self, b):
        return self.slope*F(b)+self.constant

    def margin_above(self, target):
        slope = self.slope-target.slope
        if slope < 0:
            raise ValueError('An order inequality fails at arbitrarily large stages')
        return slope*MIN_B+self.constant-target.constant


def line(slope, constant):
    return Affine(F(slope), F(constant))


def initialization(kappa=KAPPA):
    k = F(kappa)
    if not 0 < k < F(1, 100):
        raise ValueError('Small positive graph loss required')
    h0 = 1-k
    # The sharper real-primary covariance bound from the saved checkpoint
    # yields an auxiliary-averaged mean residual of order 2-3*kappa.
    orders = dict(primary_linear_wave=1-3*k, primary_wave_nonlinearity=1-k,
                  raw_full_mean=1-k, primary_averaged_mean=2-3*k,
                  temporal_mean_increment=h0, induced_radial_increment=h0+1,
                  updated_wave_residual=1-3*k,
                  updated_tangential_mean=2-3*k,
                  updated_compatibility_defects=h0+MU-2*k,
                  total_mean_velocity_and_pressure=h0,
                  total_radial_mean_velocity=h0+1,
                  exact_primary_wave_remainder=1-k)
    margins = dict(wave=orders['updated_wave_residual']-MIN_B,
                   tangential_mean=orders['updated_tangential_mean']-(MIN_B+F(1, 2)),
                   defects=orders['updated_compatibility_defects']-(MIN_B+F(1, 2)),
                   total_mean=h0-MU, total_radial_mean=h0+1-(MU+1),
                   primary_remainder=1-k-RHO)
    if min(margins.values()) <= 0:
        raise ValueError('Initialization does not meet the declared state bounds')
    return dict(orders=orders, margins=margins, B0=MIN_B, C0=MIN_B+F(1, 2),
                theorem_is_conditional=True, actual_global_fields_evaluated=False)


def cycle_certificate(kappa=KAPPA, requested_gain=STEP):
    k, gain = F(kappa), F(requested_gain)
    if k <= 0 or gain <= 0:
        raise ValueError('Positive graph loss and requested gain required')
    wave_target, mean_target = line(1, gain), line(1, F(1, 2)+gain)
    step1_wave = dict(linear=line(1, F(1, 2)-3*k),
                      old_wave_cross=line(1, F(1, 2)-k),
                      new_wave_square=line(2, -k),
                      total_mean_cross=line(1, MU-F(1, 2)))
    step2_wave = dict(linear=line(1, F(1, 2)-4*k),
                      total_mean_cross=line(1, MU-F(1, 2)-k),
                      old_wave_cross=line(1, F(1, 2)-2*k),
                      signed_square=line(2, -3*k))
    # Covariance remainders AFTER radial divergence. Here C=B+1/2.
    averaged_remainders = dict(old_primary_remainder=line(1, RHO-2*k),
                               new_curl_cross=line(1, 1-3*k),
                               signed_square=line(2, -3*k),
                               retained_moment_bump=line(1, F(3, 2)-k),
                               axial_pressure_and_flux=line(1, F(3, 2)-2*k))
    averaged_target = line(1, F(1, 2)+F(17, 100))
    h = line(1, F(1, 2)-2*k)
    final_mean_rows = dict(improved_auxiliary_average=averaged_target,
                           temporal_and_slow_mean_errors=line(1, F(3, 2)-4*k))
    final_defect = line(1, F(1, 2)+MU-4*k)
    groups = dict(step1_wave=step1_wave, step2_wave=step2_wave,
                  averaged_covariance_remainders=averaged_remainders,
                  final_mean_rows=final_mean_rows)
    margins = {}
    for name, rows in groups.items():
        target = (averaged_target if name == 'averaged_covariance_remainders'
                  else mean_target if name == 'final_mean_rows' else wave_target)
        for term, value in rows.items():
            margins[name+':'+term] = value.margin_above(target)
    margins['mean_updates:wave_cross'] = h.margin_above(wave_target)
    margins['final_defects'] = final_defect.margin_above(mean_target)
    margins['preserve_total_wave_remainder'] = MIN_B-k-RHO
    margins['preserve_total_mean'] = h.at(MIN_B)-MU
    margins['five_map_input_at_least_0.9'] = h.at(MIN_B)-MU
    if min(margins.values()) < 0:
        raise ValueError('Requested all-stage cycle gain or a required input estimate fails')
    return dict(kappa=k, requested_gain=gain, minimum_B=MIN_B,
                parameter_domain='every B>=7/10, not a finite list of stages',
                term_groups=groups, temporal_and_moment_increment_order=h,
                wave_target=wave_target, mean_and_defect_target=mean_target,
                final_defect_order=final_defect, all_half_line_margins=margins,
                minimum_margin=min(margins.values()),
                fixed_background_primary_and_inverse_maps=True,
                theorem_is_conditional=True,
                unconditionally_verified_blowup=False)


def physical_loss(h=F(1, 200), derivative_order=0):
    """A conservative stage-INDEPENDENT loss for fixed Cartesian jets.

    Phase/graph derivatives cost at most q^-2 per physical derivative on
    the annular support. One extra q^-1 absorbs each fixed logarithmic
    polynomial. Its constant may depend on stage and derivative order.
    """
    h = F(h)
    if not 0 < h < F(1, 100) or type(derivative_order) is not int or derivative_order < 0:
        raise ValueError('Construction h and a nonnegative integer jet order required')
    a = F(1, 2)+h
    return 2*a+F(3, 2)+2*derivative_order


def physical_residual_power(stage, derivative_order=0, h=F(1, 200)):
    if type(stage) is not int or stage < 0:
        raise ValueError('A nonnegative finite stage is required')
    return F(h)*(MIN_B+STEP*stage)-physical_loss(h, derivative_order)
