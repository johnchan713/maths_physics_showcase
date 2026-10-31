"""Finite constants for the inverse leading viscosity multiplier on T^2."""
from fractions import Fraction as Q

# Lattice separation and two sides of the origin; Parseval with normalized Haar.
SHELL_CONSTANT = 4 * 12**3 * 2 * 2**5
INVERSE_CONSTANT = 2 * SHELL_CONSTANT


def necessary_target(coefficient, graph_lower, floor_inverse):
    """Coefficient at fixed S: floor_inverse >= 1/c >= 1; no epsilon here."""
    coefficient, graph_lower, floor_inverse = map(Q,(coefficient, graph_lower, floor_inverse))
    if coefficient < 0 or graph_lower <= 0 or floor_inverse < 1:
        raise ValueError('Require a nonnegative jet bound and positive graph/floor bounds')
    return INVERSE_CONSTANT * coefficient * floor_inverse**2 / graph_lower**3


def target_order(viscous_order,kappa=Q(1,100000)):
    return Q(viscous_order)-1+3*kappa


def pell_modes(count):
    """(a+b sqrt(2))=(1+sqrt(2))^j; k1+k2=a, k2=b."""
    a,b=1,0
    for _ in range(count):
        a,b=a+2*b,a+b
        yield a-b,b


def record():
    k=Q(1,100000)
    return dict(auxiliary_derivatives=5, inverse_constant=INVERSE_CONSTANT,
                critical_target_order=str(target_order(3-6*k,k)),
                example_modes=[list(k) for k in pell_modes(5)], physical_target_lower_bound=False,
                complete_residual_or_blowup=False)
