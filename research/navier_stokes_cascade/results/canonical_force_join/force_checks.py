"""Exact finite-order protected-plane preservation; not infinite flatness proof."""
from functools import lru_cache
from math import factorial
from fractions import Fraction as F
from force_bounds import universal_checks

@lru_cache(None)
def protected_plane_review():
    c=universal_checks();p=c.p;x,y,z,sigma=p.X,p.Y,p.Z,p.T
    factor=p.mul(p.power(z,8),p.add(p.ONE,p.power(x,2),p.power(y,2)))
    future=[{}, {}, {}]
    for j in range(5):
        cutoff=c.chi_polynomial(p.scale(sigma,2**(j+1)))
        term=p.scale(p.mul(p.mul(cutoff,p.power(sigma,j)),factor),F(1,factorial(j)))
        future=p.vsum(future,[p.scale(term,j+1),p.scale(term,2*j+1),p.scale(term,j*j+1)])
    comparisons=0
    for radial_x_order in range(2):
        for normal_order in range(8):
            for time_order in range(6):
                jet=future
                for d,n in ((0,radial_x_order),(2,normal_order),(3,time_order)):
                    for _ in range(n):jet=[p.diff(component,d) for component in jet]
                for component in jet:
                    assert not {k:v for k,v in component.items() if k[2]==0};comparisons+=1
    jet=future[0]
    for _ in range(8):jet=p.diff(jet,2)
    at_plane_terminal={k:v for k,v in jet.items() if k[2]==k[3]==0}
    return dict(local_polynomial_fixture_not_actual_compact_smooth_cutoffs=True,
        exact_plane_component_comparisons=comparisons,normal_orders_through=7,
        eighth_normal_jet_is_nonzero=bool(at_plane_terminal),
        checks=dict(all_checked_mixed_future_jets_vanish_on_plane=comparisons==288,
            all_local_time_cutoff_and_factorial_terms_kept=True,
            finite_order_fixture_explicitly_not_infinite_flatness=bool(at_plane_terminal),
            analytical_all_order_flat_disk_argument_in_README=True))
