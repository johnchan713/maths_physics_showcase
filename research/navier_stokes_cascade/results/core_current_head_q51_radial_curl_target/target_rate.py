"""Conditional first-real target bound with the component-specific curl operator."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,HERE.parent/path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
curl=load('q51_curl','core_complete_curl_operator_budget/bounds.py')
profile=load('q51_profile','core_current_head_q36_profile_signed_target/profile_target.py')

def radial_operator(op,j):
    """U-chart coefficient norm of B_r / epsilon, acting on A."""
    curl.order(j)
    def raw(k):
        return Q(2*5**k,factorial(k))*(op.lift(k)+op.lift(k+1))
    return op.compose(raw,j)

def target_bound(target,op,n,e):
    """Reuse the Q36 polynomial only after verifying domination at every needed order.

    op must describe the SAME normal on the SAME native hull (physical premise).
    """
    if not isinstance(target,profile.signed.SignedTarget) or not isinstance(target.source,profile.Source):
        raise TypeError('Q36 original first-real target required')
    if not isinstance(op,curl.CurlBudget):raise TypeError('Canonical curl operator required')
    target.check(n)
    if e not in (1,2):raise ValueError('Selected flux only')
    if (op.Rmin,op.Rmax)!=(target.Rmin,target.source.profile.Rmax):
        raise ValueError('Inconsistent native hull')
    for j in range(n+4):
        upper=target.root.data('curl',j).c
        if any(upper.get(i,0)<v for i,v in radial_operator(op,j).c.items()):
            raise ValueError('Existing curl polynomial does not cover the stripped radial operator')
    return dict(order=n,pole=target.pole+3*n,epsilon_order=2,
                normalized_target=target.native(n,e))

def record():
    k=Q(1,100000);alpha=2-3*k
    return dict(selected_curl_order=str(3-k),selected_defect_order=str(min(3,3-k)),
                haar_row_and_target_order=2, isolated_viscosity_order=str(3-3*k),
                full_fast_source_order=str(min(alpha,2-2*k,2-k,2-k)),
                actual_field_enclosures_or_complete_residual=False)
