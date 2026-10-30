"""Exact projected-label bounds; no evaluation of actual pulse seminorms."""
from fractions import Fraction as Q
from math import ceil, comb

BANDS = 3
TRANSVERSE_MESH = 3**2
SIGNS = 2
DECK = 14**3
S_RATIO = Q(25, 16)

def nonnegative(x):
    if isinstance(x, (float, bool)):
        raise ValueError('exact nonnegative input required')
    x = Q(x)
    if x < 0:
        raise ValueError('negative input')
    return x

def projected_count(S, radial_length):
    S = nonnegative(S); length = nonnegative(radial_length)
    if S < 1:
        raise ValueError('S >= 1 required')
    return BANDS*SIGNS*TRANSVERSE_MESH*DECK*(ceil(length*(S_RATIO*S)**3)+3)

def count_constant(radial_length):
    length = nonnegative(radial_length)
    return BANDS*SIGNS*TRANSVERSE_MESH*DECK*(length*S_RATIO**3+4)

def aggregate_polynomial(coefficients, radial_length):
    """C S^3 P(S): a safe sum bound after radial/auxiliary support fills."""
    cs = [nonnegative(c) for c in coefficients]
    if not cs:
        raise ValueError('nonempty polynomial required')
    return [Q(0)]*3+[count_constant(radial_length)*c for c in cs]

def product_jet_bound(left, right):
    """Order-m directional Leibniz majorant, retaining both cross orientations."""
    if len(left) != len(right) or not left:
        raise ValueError('equal nonempty jet lists required')
    a = [nonnegative(x) for x in left]; b = [nonnegative(x) for x in right]
    m = len(a)-1
    return sum(comb(m,k)*a[k]*b[m-k] for k in range(m+1))

def evaluate(coefficients, x):
    return sum(Q(c)*Q(x)**i for i,c in enumerate(coefficients))

# Rational polynomial fixtures below test the pressure moment mechanism.
# They are finite-regularity compact fixtures, not actual pulse profiles.
def integral(p, power=0):
    return sum(Q(c,Q(i+power+1)) for i,c in enumerate(p))

def pressure_fixture(g, rho):
    """On [0,1], p'=g-rho*P, p(0)=0, int rho=1."""
    g = list(map(Q,g)); rho = list(map(Q,rho))
    if integral(rho) != 1:
        raise ValueError('unit-mass rho required')
    P = integral(g)
    n = max(len(g),len(rho))
    f = [(g[i] if i<len(g) else 0)-P*(rho[i] if i<len(rho) else 0)
         for i in range(n)]
    p = [Q(0)]+[v/Q(i+1) for i,v in enumerate(f)]
    return P,p

def record():
    return {'projected_label_bound':'54*2744*(ceil(L_R*(25*S/16)^3)+3)',
        'polynomial_bound':'54*2744*(L_R*(25/16)^3+4)*S^3',
        'extra_log_degree_per_label_sum':6,
        'linear_maps_preserve_slow_projection':True,
        'radial_or_auxiliary_disjointness_of_means':False,
        'zero_Haar_A0_is_not_zero_A0':True,
        'same_base_pressure_memory_retained':['Pi1','Pi2'],
        'actual_normalized_coefficient_estimates_certified':False}
