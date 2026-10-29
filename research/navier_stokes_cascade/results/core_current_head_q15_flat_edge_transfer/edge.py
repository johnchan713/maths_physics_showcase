"""Fixed-width flat-edge transfer; constants remain symbolic at giant requests."""
from fractions import Fraction as Q
from math import factorial

def ceiling(x):
    x=Q(x);return -(-x.numerator//x.denominator)

def request(p,vanishing_order=0):
    if isinstance(p,(float,bool)) or isinstance(vanishing_order,(float,bool)):
        raise ValueError('exact nonnegative orders required')
    p,s=Q(p),Q(vanishing_order)
    if min(p,s)<0:raise ValueError('nonnegative orders required')
    n=ceiling((p+s)/2)
    return {'pole':str(p),'vanishing_order':str(s),'n':n,
            'power_margin':str(2*n-p-s),
            'constant':'factorial(n)*w^(-2*n)',
            'materialized_constant':False}

def fixture_constant(w,p,s=0):
    """Small rational diagnostic, never invoked for a root-stage pole archive."""
    if isinstance(w,(float,bool)) or not 0<Q(w)<=1:
        raise ValueError('exact 0<w<=1 required')
    n=request(p,s)['n']
    if n>1000:raise ValueError('large factorial must remain symbolic')
    return Q(factorial(n))/Q(w)**(2*n)

def normalized_majorant(w,p,s,delta):
    """Rational upper bound obtained from exp(u)>=u^n/n!, without exp evaluation."""
    delta=Q(delta)
    if not 0<delta<=1:raise ValueError('0<delta<=1 required')
    n=request(p,s)['n']
    # Integer exponents are used only by this diagnostic evaluator.
    power=2*n-Q(p)-Q(s)
    if power.denominator!=1:raise ValueError('integer exponent fixture required')
    return fixture_constant(w,p,s)*delta**power.numerator
