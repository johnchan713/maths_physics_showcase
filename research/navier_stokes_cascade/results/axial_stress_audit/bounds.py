"""Outward interval bounds for A.24 on the complete axial transition.

See README.md for the reduction from the exact finite-h moment identities.
No P*, h, exp(Td), or inner log-amplitude is evaluated at large Md.
"""
from fractions import Fraction
import mpmath as mp


def lo(value):
    return mp.make_mpf(value._mpi_[0])


def hi(value):
    return mp.make_mpf(value._mpi_[1])


def point(value):
    if isinstance(value, Fraction):
        return mp.iv.mpf(value.numerator) / value.denominator
    return mp.iv.mpf(value)


def step_at(value):
    """Interval enclosure at an exact rational point, using reflection."""
    value = Fraction(value)
    if value <= 0:
        return point(0)
    if value >= 1:
        return point(1)
    if value > Fraction(1, 2):
        return 1 - step_at(1 - value)
    x = point(value)
    return 1 / (1 + mp.iv.exp(1/x**2 - 1/(1-x)**2))


def derivative_upper(a, b):
    """Enclose sup sigma' over [a,b]; no sampled derivative is substituted."""
    if not 0 <= a < b <= 1:
        raise ValueError('Expected an ordered subinterval of [0,1]')
    if a >= Fraction(3, 4):
        return derivative_upper(1-b, 1-a)
    if b <= Fraction(1, 4):
        # This majorant is increasing for 0<t<=1/4 and vanishes at zero.
        x = point(b)
        return min(mp.mpf(8), hi((2/x**3+16)*mp.iv.exp(4-1/x**2)))
    if a == 0 or b == 1:
        return mp.mpf(8)
    # The rational factor is convex; sigma*(1-sigma) increases toward 1/2.
    closest = min(max(Fraction(1, 2), a), b)
    factors = [(2/point(x)**3+2/(1-point(x))**3) for x in (a, b)]
    rational_upper = max(hi(v) for v in factors)
    product = step_at(closest)*step_at(1-closest)
    return min(mp.mpf(8), hi(point(rational_upper)*product))


def axial_bound(md, boxes=1024, digits=50):
    """All y in [0,Td], all eta in [-1,1], for the stated reference prefix.

    Returns sufficient A.24 margins, not a full profile or finite-XR cone.
    The finite remainder delta=exp(-Td) is bounded by 2e-28 for every Md>=4.
    The bound holds uniformly for 0<h<=exp(-2*Td), with h<lambda<0.1.
    """
    if not isinstance(md, int) or md < 4:
        raise ValueError('Md must be an integer >=4')
    if not isinstance(boxes, int) or boxes < 4 or digits < 30:
        raise ValueError('Use at least four boxes and 30 interval digits')
    mp.iv.dps = digits
    delta = point('2e-28')
    q_loss = 217*delta
    bounds = []
    for i in range(boxes):
        a, b = Fraction(i, boxes), Fraction(i+1, boxes)
        slope = derivative_upper(a, b)
        r = mp.iv.exp(-md*point(a))
        k = 4*step_at(1-b)
        angular = max(hi(2*r), hi((2+4*r)/(1+2*k)))
        bound = 8*point(slope)*(point(angular)+370*delta)/(md*(1-q_loss))
        bounds.append(hi(bound))
    index = max(range(boxes), key=bounds.__getitem__)
    B = point(bounds[index])
    bs2 = 1024*delta
    first = 2-B
    second = 2-2*B-bs2/2
    return dict(Md=md, boxes=boxes, interval_digits=digits,
                maximum_bsw_absolute_upper=hi(B),
                first_A24_margin_lower=lo(first),
                second_A24_margin_lower=lo(second),
                Pc_over_ps1_lower=lo(1-B/2),
                relative_Q_loss_upper=hi(q_loss),
                bs_squared_upper=hi(bs2),
                worst_t_box=[str(Fraction(index, boxes)), str(Fraction(index+1, boxes))],
                status='axial-ratio-bound-passed' if lo(second)>0 and lo(first)>0
                       else 'bound-inconclusive')


def elementary_checks():
    """Interval checks of constants used in the analytic inequalities."""
    mp.iv.dps = 50
    e = mp.iv.exp(1)
    memory_floor = 45/(128*e)
    return {
        'delta_below_2e-28': hi(mp.iv.exp(-mp.iv.exp(4)-10)) < mp.mpf('2e-28'),
        'energy_memory_below_1p5': hi(1+point(Fraction(5,6))*mp.iv.exp(point('-.6'))) < mp.mpf('1.5'),
        'Q_loss_constant_below_217': hi(28/memory_floor) < 217,
        'geometric_constant_below_45': hi(44+64*point('4e-56')) < 45,
    }
