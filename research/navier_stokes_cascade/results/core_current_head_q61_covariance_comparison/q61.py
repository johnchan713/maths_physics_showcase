"""Value gates for the conditional same-label comparison."""
from fractions import Fraction


def branch(eps, floor_cost, volterra_cost, inverse_cost):
    if not 0 < eps <= 1 or min(floor_cost, volterra_cost, inverse_cost) < 0:
        raise ValueError("invalid scale or comparison cost")
    B = max(1, floor_cost, volterra_cost, inverse_cost)
    if 2*eps**2*B <= 1:
        return 'compare', Fraction(1)/(1-eps**2*volterra_cost)
    return 'direct', 2*B


def record():
    return {"correction_order_offset": 6, "leading_order_offset": 99998,
            "physical_enclosure_verified": False,
            "boundary": str(branch(Fraction(1,4),1,8,1)[1])}
