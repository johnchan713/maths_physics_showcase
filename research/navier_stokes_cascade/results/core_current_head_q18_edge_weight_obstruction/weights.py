"""Exact radial edge-weight bookkeeping; not actual-field certification."""
from fractions import Fraction as Q

def weight(x):
    if isinstance(x,(float,bool)): raise ValueError('exact weight required')
    x=Q(x)
    if x<0: raise ValueError('nonnegative weight required')
    return x

def product(*xs): return sum((weight(x) for x in xs),Q(0))
def sum_weight(*xs):
    if not xs: raise ValueError('nonempty sum required')
    return min(map(weight,xs))
def embeds(source,target): return weight(source)>=weight(target)
def signed_quotient(target): return weight(target)-Q(1,2)

def radial_log_lower(n,p):
    """w=1 fixture: log(zeta^(-1/2) x^p) >= n^2/2-p*n at x=1/n.

    Uses log n <= n. The omitted outer-edge exponential is positive.
    """
    if type(n) is not int or n<1 or type(p) is not int or p<0:
        raise ValueError('integer n>=1, p>=0 required')
    return Q(n*n,2)-p*n

def record():
    return {'linear_wave_source': '1/2', 'wave_wave_source':'1',
        'mean_wave_source_with_full_weight_mean':'3/2',
        'sum_of_these_sources':'1/2', 'mean_covariance':'1',
        'full_weight_signed_target_after_root_division':'1/2',
        'half_weight_target_after_root_division':'0',
        'Q5_blanket_full_weight_source_promotion_justified':False,
        'actual_source_violates_full_weight_bound_proved':False,
        'typed_all_stage_induction_completed':False}
