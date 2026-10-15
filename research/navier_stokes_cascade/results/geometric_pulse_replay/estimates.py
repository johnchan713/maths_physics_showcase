"""Radius-explicit continuum certificates; numerical fixtures live separately."""
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parent.parent


def frozen_module(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


charts = frozen_module('replay_frozen_charts', 'actual_band_pulse_bounds/charts.py')
old_bounds = frozen_module('replay_frozen_bounds', 'actual_band_pulse_bounds/bounds.py')
geometry = frozen_module('replay_frozen_geometry', 'auxiliary_geometry_audit/geometry.py')
RADIUS = geometry.radius_bound()


def real(x):
    return mp.mpf(x.numerator) / x.denominator if isinstance(x, F) else mp.mpf(x)


def radius(x):
    x = real(x)
    if not mp.isfinite(x) or not 0 < x <= real(RADIUS):
        raise ValueError('A positive radius within the constructive global bound is required')
    return x


def length_interval(ell, r=RADIUS):
    if type(ell) is not int or ell < 16:
        raise ValueError('An integer band at least 16 is required')
    r = radius(r)
    return 2*r*ell**2, 2*r*(4+mp.sqrt(2))*ell**2


def all_later_threshold(h, r, minimum_length, preliminary_band=16):
    """Every later band has L>=minimum_length and S^2*(eps+eps^2+1/k)<=1.

    The exponential majorant is an all-half-line proof, not a band scan.
    Profile value/gap conditions must already hold from preliminary_band.
    """
    r, h, minimum_length = radius(r), real(h), real(minimum_length)
    if not mp.isfinite(minimum_length) or minimum_length < 1:
        raise ValueError('A finite minimum length at least one is required')
    if type(preliminary_band) is not int or preliminary_band < 16:
        raise ValueError('A fixed preliminary band at least 16 is required')
    exponential = old_bounds.eventual_band_threshold(h, 3, 4)
    length = int(mp.ceil(mp.sqrt(minimum_length/(2*r))))+1
    return dict(first_band=max(preliminary_band, exponential, length),
                preliminary_band=preliminary_band, exponential_band=exponential,
                pulse_length_band=length, minimum_length=minimum_length,
                actual_profile_constants_numerically_certified=False)


def gaussian_constants(lam_min, lam_max, u):
    lo, hi, u = map(real, (lam_min, lam_max, u))
    if not all(mp.isfinite(x) for x in (lo, hi, u)) or not 0 < lo <= hi or u <= 0:
        raise ValueError('Finite positive spectral interval and tilt required')
    c = u*lo*(u/2)*((1+(3*u/2)**2)**(-mp.mpf('1.5'))
                      +2*(1+u*u)**(-mp.mpf('1.5')))
    upper = u*hi*(3*u/2)*((1+(u/2)**2)**(-mp.mpf('1.5'))
                            +2*(1+u*u)**(-mp.mpf('1.5')))
    return dict(lower=c, upper=upper)


def pulse_barrier(mu, ce, cd, length):
    """Entrywise |E|<=ce/L and |d-dref|<=cd/L, hence radius-free factors."""
    mu, ce, cd, length = map(real, (mu, ce, cd, length))
    if not all(mp.isfinite(x) for x in (mu, ce, cd, length)) or min(mu, ce, length) <= 0 or cd < 0:
        raise ValueError('Finite positive gap, entry bound and length required')
    if ce/length > mu/8:
        raise ValueError('Pulse too short for the Riccati invariant interval')
    ratio = ce/(mu*length)
    outward = -2*mu*ratio+ce/length*(1+2*ratio+ratio**2)
    log_bound = cd+ce*(1+ratio)
    return dict(ratio_bound=ratio, outward_derivative_upper=outward,
                log_amplitude_bound=log_bound,
                radial_lower=(1-ratio)*mp.exp(-log_bound),
                radial_upper=(1+ratio)*mp.exp(log_bound),
                log_forward_propagator_factor=2*ce+cd)


def minimum_pulse_length(mu, ce, direction_constant, ac, u, cone_margin):
    mu, ce, direction, ac, u, margin = map(real, (mu, ce, direction_constant, ac, u, cone_margin))
    if not all(mp.isfinite(x) for x in (mu, ce, direction, ac, u, margin)) or min(mu, ce, direction, ac, u, margin) <= 0 or margin > 1:
        raise ValueError('Positive finite gaps/constants and cone margin at most one required')
    cone = (direction*(1/ac+1/u)*(4+margin)/margin)**2
    return max(mp.mpf(25), 8*ce/mu, cone)


def covariance_interval(r, length, lower_gaussian, upper_gaussian, radial_lower, radial_upper):
    """Full Haar weights, including transverse area and angular factor 1/2.

    chi(xi)=beta(2*xi/r), supported in [-r/2,r/2] strictly inside
    the rectangle. Its squared integral lies between (r/2)*exp(-2/3)
    and r. psi=1 centrally.
    Bounds apply once L>=25 so |v-L/2|<=sqrt(L) is central.
    """
    r = radius(r)
    length, c, upper, xlo, xhi = map(real, (length, lower_gaussian, upper_gaussian, radial_lower, radial_upper))
    if not all(mp.isfinite(x) for x in (length, c, upper, xlo, xhi)) or length < 25 or not 0 < c <= upper or not 0 < xlo <= xhi:
        raise ValueError('Length at least 25 and positive ordered Gaussian/radial bounds required')
    determinant = 1+(mp.sqrt(2)-1)**2
    factor = r*r/mp.sqrt(length)
    lo = determinant*mp.exp(-mp.mpf(2)/3)*xlo*xlo*mp.exp(-upper)
    hi = determinant*xhi*xhi*mp.sqrt(mp.pi/c)
    return dict(lower=lo*factor, upper=hi*factor,
                lower_over_r_squared_length_inverse_half=lo,
                upper_over_r_squared_length_inverse_half=hi)


def derivative_degrees(order, source=None):
    """Differentiated Duhamel: one S per integral, S^i per coefficient jet."""
    return old_bounds.derivative_exponents(order, source_powers=source)


def flat_tail_threshold(r, c_gaussian, constant, physical_linear_loss, log_power, target):
    r, cg = radius(r), real(c_gaussian)
    if not mp.isfinite(cg) or cg <= 0:
        raise ValueError('A positive Gaussian rate is required')
    decay = cg*r/25  # P<=exp(-cg*L/50), L>=2*r*ell^2.
    threshold = old_bounds.endpoint_flat_threshold(constant, decay, physical_linear_loss, log_power, target)
    return dict(first_band=threshold, quadratic_decay=decay)


def radius_power_ledger():
    """Exact exponents in powers of r and S; no floating fitted scalings."""
    return dict(haar_covariance_r=F(3, 2), haar_covariance_S=F(-1, 2),
                covariance_inverse_r=F(-3, 2), covariance_inverse_S=F(1, 2),
                primary_root_r=F(-3, 4), primary_root_S=F(1, 4),
                squared_primary_root_times_haar_r=2*F(-3, 4)+F(3, 2),
                finite_radius_does_not_change_epsilon_orders=True)


def exact_scale_review():
    powers = radius_power_ledger()
    return dict(radius=RADIUS, powers=powers, maximum_length_over_scale=12*RADIUS,
                checks=dict(global_radius_is_positive=RADIUS>0,
                            length_upper_bound_is_less_than_S=12*RADIUS<1,
                            transverse_area_factor_is_retained=powers['haar_covariance_r']==F(3, 2),
                            primary_root_compensates_small_support=powers['squared_primary_root_times_haar_r']==0,
                            no_uniform_in_radius_covariance_inverse_claim=powers['covariance_inverse_r']<0,
                            fixed_radius_preserves_epsilon_orders=powers['finite_radius_does_not_change_epsilon_orders']))


def threshold_review():
    rows = []
    for r in (RADIUS, RADIUS/4):
        row = all_later_threshold(F(1, 200), r, 1000)
        ell = row['first_band']
        row.update(radius=r, lower_length=length_interval(ell, r)[0],
                   last_majorant=old_bounds.phase_majorant(ell, real(F(1, 200))))
        rows.append(row)
    low = length_interval(16)[0]
    return dict(rows=rows, low_band_lower_length=low,
                checks=dict(all_later_lengths_follow_lower_bound=all(row['lower_length']>=1000 for row in rows),
                            exponential_half_line_condition=all(row['last_majorant']<1 for row in rows),
                            low_geometry_bands_are_not_analytic_pulse_bands=low<1,
                            smaller_radius_requires_a_later_length_threshold=rows[1]['pulse_length_band']>rows[0]['pulse_length_band']))


def derivative_and_tail_review():
    homogeneous=derivative_degrees(6)
    source=derivative_degrees(6,[0]*7)
    growing=derivative_degrees(6,[3*j for j in range(7)])
    rows=[]
    with mp.workdps(70):
        for r in (RADIUS,RADIUS/2):
            row=flat_tail_threshold(r,mp.mpf('.1'),mp.mpf('1e80'),9,12,5)
            a=row['quadratic_decay']
            first=row['first_band']
            logs=[mp.log(mp.mpf('1e80'))+12*mp.log(ell)+(9+5*mp.log(2))*ell-a*ell**2
                  for ell in (first,first+1,2*first)]
            rows.append(dict(radius=r,first_band=first,maximum_log_tail_over_q_fifth=float(max(logs))))
    return dict(homogeneous_degrees=homogeneous,constant_source_degrees=source,
                growing_source_degrees=growing,tail_rows=rows,
                checks=dict(every_fixed_homogeneous_jet_has_finite_degree=homogeneous==[2*j for j in range(7)],
                            source_integral_cost_is_retained=source==[2*j+1 for j in range(7)],
                            arbitrary_finite_source_jet_growth_is_allowed=growing==[3*j+1 for j in range(7)],
                            radius_dependent_Gaussian_beats_all_fixed_physical_losses=all(r['maximum_log_tail_over_q_fifth']<0 for r in rows),
                            smaller_radius_does_not_have_uniform_tail_threshold=rows[1]['first_band']>rows[0]['first_band'],
                            cutoff_derivatives_absorbed_by_polynomial_scale=F(12,25)<1))
