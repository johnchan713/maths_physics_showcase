"""Factored native weight bounds (do not materialize huge high-order jets)."""
from fractions import Fraction as F


def geometry(r,H,mu_abs):
    r,H,mu_abs=map(F,(r,H,mu_abs))
    if not 0<r<=H or mu_abs<0:
        raise ValueError("positive hull and absolute cutoff moment required")
    u=min(F(1),r**3); delta=u/1024; C=8*mu_abs; R=2*H; D=6/u
    # First two coefficients are absolute. Last two multiply
    # L_(j+1) exp(8H²+3); all four multiply j! (2/delta)^j.
    B=8*max(1,2/r)
    return (D*R**2,D*(C+R**2/2),2*R*2*B/delta,
            (2*C+R**2)*2*B/delta+4*C*D*B),2/delta


def record():
    co,base=geometry(F(1,2),2,F(3,2))
    return dict(fixture_coefficients=list(map(str,co)),fixture_raw_base=str(base),
                leading_extra=99998,background_error_extra=2,
                background_weight_replacement_order=4,
                physical_amplitude_enclosures_verified=False)
