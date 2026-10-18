"""Exact local gauge/cutoff controls; not actual compact canonical fields."""
from functools import lru_cache
from fractions import Fraction as F
from realization import universal_controls

@lru_cache(None)
def review():
    p=universal_controls();X,Y,Z,T=p.X,p.Y,p.Z,p.T;one=p.ONE
    phi=p.add(p.mul(X,Z),p.mul(T,p.power(Y,2)))
    q=p.add(one,p.scale(T,-1),p.power(Z,2))
    chi=p.add(one,p.scale(p.power(q,2),-3),p.scale(p.power(q,3),2))
    chi2=p.add(one,p.scale(p.power(q,2),-5),p.scale(p.power(q,3),3))
    gradphi=p.grad(phi);artifact=p.curl(p.vmul(chi,gradphi))
    product=p.cross(p.grad(chi),gradphi)
    common=p.add(one,T);c1=p.add(common,p.power(T,2));c2=p.add(common,p.power(T,3))
    pressure_artifact=p.grad(p.add(p.mul(chi,p.add(c1,p.scale(common,-1))),
        p.mul(chi2,p.add(c2,p.scale(c1,-1)))))
    common_increment=p.add(common,p.scale(common,-1))
    # Partial telescoping p2-p0=(p1-p0)+(p2-p1), including middle pressure.
    p0=p.mul(X,Z);p1=p.add(p0,p.mul(Y,T));p2=p.add(p1,p.mul(p.power(Z,2),T))
    increment=p.add(p1,p.scale(p0,-1),p2,p.scale(p1,-1))
    telescopes=not p.add(increment,p.scale(p.add(p2,p.scale(p0,-1)),-1))
    h=F(1,200);eta=F(1,2);sigma=F(3,4);L=1-2*h*eta*eta
    transition_factor=2*eta*sigma/L
    return dict(local_exact_polynomial_fixture=True,actual_supported_gauge_sequence_evaluated=False,
        potential_gauge_artifact_components=sum(bool(v) for v in artifact),
        stage_pressure_gauge_artifact_components=sum(bool(v) for v in pressure_artifact),
        actual_q_gauge_counterexample=dict(h=str(h),eta=str(eta),sigma=str(sigma),
            pressure_gradient_factor_per_chi_prime=str(transition_factor),power=str(-(F(1,2)-h)),
            actual_canonical_fields=False),
        checks=dict(uncut_gradient_gauge_has_zero_curl=not any(p.curl(gradphi)),
            potential_gauge_after_cutoff_changes_velocity=any(artifact),
            exact_potential_gauge_artifact_is_grad_chi_cross_grad_phi=not any(p.vsum(artifact,p.vscale(product,-1))),
            potential_gauge_artifact_still_divergence_free=not p.divergence(artifact),
            stagewise_time_gauges_have_zero_uncut_spatial_gradient=not any(p.grad(c1)) and not any(p.grad(c2)),
            independent_stage_time_gauges_create_cut_pressure_force=any(pressure_artifact),
            common_time_gauge_cancels_from_all_pressure_increments=not common_increment and not any(p.grad(common)),
            actual_q_counterexample_points_inside_t_less_than_one=1-eta*eta>0,
            strict_doubling_isolates_one_cutoff_derivative=sigma/2<F(1,2) and 2*sigma>1,
            nonflat_pressure_force_factor_and_negative_power=transition_factor>0 and F(1,2)-h>0,
            full_middle_pressure_telescope_exact=telescopes,
            prematurely_discarding_middle_pressure_change_is_detected=bool(p.add(p2,p.scale(p1,-1),p.scale(increment,-1)))))
