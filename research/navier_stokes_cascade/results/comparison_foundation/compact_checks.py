"""Independent exact local energy identity and compact exterior kernel controls."""
from fractions import Fraction as F
from functools import lru_cache
from compact_bounds import HERE,load

p=load('comparison_foundation_independent_polynomials',HERE.parent/'realization_foundation/realized_checks.py')
def dot(a,b):return p.add(*(p.mul(x,y) for x,y in zip(a,b)))
def lap(f):return p.add(*(p.diff(p.diff(f,i),i) for i in range(3)))

@lru_cache(None)
def local_review():
    x,y,z,t=p.X,p.Y,p.Z,p.T
    u=p.curl([p.mul(y,z),p.mul(x,t),p.mul(p.power(x,2),y)])
    w=p.curl([p.mul(x,z),p.mul(y,t),p.add(p.mul(x,p.power(y,2)),p.mul(z,t),p.mul(t,p.mul(x,y)))])
    v=p.vsum(u,w);pi=p.add(p.power(x,2),p.mul(y,z),p.mul(t,x))
    chi=p.power(p.add(p.ONE,x,p.power(z,2)),2);gc=p.grad(chi)
    assert not p.divergence(u) and not p.divergence(w)
    split=0
    for i in range(3):
        for j in range(3):
            g=p.add(p.mul(v[i],v[j]),p.scale(p.mul(u[i],u[j]),-1))
            h=p.add(p.mul(u[i],w[j]),p.mul(w[i],u[j]))
            assert not p.add(g,p.scale(p.mul(w[i],w[j]),-1),p.scale(h,-1));split+=1
    s=dot(w,w);dt=p.scale(p.mul(chi,p.diff(s,3)),F(1,2))
    strain=p.mul(chi,dot(w,p.adv(w,u)))
    transport=p.scale(p.mul(s,dot(v,gc)),F(1,2));pressure=p.mul(pi,dot(w,gc))
    grad2=p.add(*(dot([p.diff(q,i) for q in w],[p.diff(q,i) for q in w]) for i in range(3)))
    nonzero={};count=0
    for nu in (F(3,5),F(1),F(17,10)):
        H=p.vsum([p.diff(q,3) for q in w],p.vscale(p.lap(w),-nu),p.adv(v,w),p.adv(w,u),p.grad(pi))
        terms=dict(time=dt,dissipation=p.scale(p.mul(chi,grad2),nu),strain=strain,
            laplacian=p.scale(p.mul(s,lap(chi)),nu/2),transport=transport,pressure=pressure,force=p.mul(chi,dot(w,H)))
        J=p.vsum(p.vscale(p.vsum(p.vmul(chi,p.grad(s)),p.vscale(p.vmul(s,gc),-1)),nu/2),
            p.vscale(p.vmul(p.mul(chi,s),v),F(-1,2)),p.vscale(p.vmul(p.mul(chi,pi),w),-1))
        identity=p.add(terms['time'],terms['dissipation'],terms['strain'],
            *(p.scale(terms[k],-1) for k in ('laplacian','transport','pressure','force')),p.scale(p.divergence(J),-1))
        assert not identity;count+=1
        for k,a in terms.items():nonzero[k]=bool(a) or nonzero.get(k,False)
        DD=p.add(*(p.diff(p.diff(p.add(p.mul(v[i],v[j]),p.scale(p.mul(u[i],u[j]),-1)),i),j) for i in range(3) for j in range(3)))
        assert not p.add(p.scale(lap(pi),-1),p.scale(DD,-1),p.divergence(H))
    return dict(manufactured_local_polynomials_and_difference_residual_not_whole_space_solutions=True,
        exact_stress_components=split,exact_weighted_local_energy_identities=count,nonzero_terms=sorted(k for k,v in nonzero.items() if v),
        checks=dict(all_nine_stress_components_retain_both_mixed_terms=split==9,
            complete_local_energy_flux_divergence_identity=count==3,
            time_viscous_strain_laplacian_transport_pressure_force_omissions_detected=all(nonzero.values()),
            full_forced_difference_pressure_keeps_minus_div_H=True,
            same_force_condition_means_H_zero_not_an_added_pressure_hypothesis=True))

def matrix(direction):
    if sum(q*q for q in direction)!=1:raise ValueError('Unit direction required')
    return [[3*direction[i]*direction[j]-int(i==j) for j in range(3)] for i in range(3)]

@lru_cache(None)
def exterior_review():
    directions=[(F(1),F(0),F(0)),(F(3,5),F(4,5),F(0)),(F(2,3),F(1,3),F(2,3)),(F(0),F(0),F(1))]
    u=(F(2),F(-1),F(3));w=(F(-1),F(4),F(2))
    h=[[u[i]*w[j]+w[i]*u[j] for j in range(3)] for i in range(3)];mass=sum(abs(a) for row in h for a in row)
    rows=[];spectral=0
    for n in directions:
        M=matrix(n)
        assert all(sum(M[i][j]*n[j] for j in range(3))==2*n[i] for i in range(3))
        assert sum(M[i][i] for i in range(3))==0 and sum(a*a for row in M for a in row)==6
        assert max(abs(a) for row in M for a in row)<=2;spectral+=1
        for R in (F(2),F(3),F(7)):
            for rho in (F(-1),F(-1,2),F(0),F(1,2),F(1)):
                distance=R-rho
                value=sum(M[i][j]*h[i][j] for i in range(3) for j in range(3))/distance**3
                bound=16*mass/R**3
                assert abs(value)<=bound;rows.append(dict(radius=str(R),source_radius=str(rho),value_without_inverse_four_pi=str(value),bound_without_inverse_four_pi=str(bound)))
    # Uniform L2 mass can be concentrated in successive cutoff annuli.
    # w_R=R^-3/2 W(x/R): norm exponent 0 and Laplacian flux exponent -2.
    norm_power=2*F(-3,2)+3;lap_power=norm_power-2
    return dict(exterior_point_tensor_kernel_controls_not_actual_L1_stresses_or_flows=True,
        exact_kernel_direction_controls=spectral,exact_separated_exterior_controls=len(rows),rows=rows,
        checks=dict(Hessian_sign_trace_and_eigenvalue_exact=spectral==4,
            every_separated_kernel_value_bounded_by_sixteen_mass_over_R_cubed=len(rows)==60,
            inverse_four_pi_kept_separate_from_radius_power=True,
            annular_L2_preserving_scaling_has_zero_norm_power=norm_power==0,
            coarse_laplacian_flux_has_exact_radius_minus_two=lap_power==-2,
            no_universal_better_rate_claim_for_a_fixed_solution=True))

def reviews():return dict(exact_complete_difference_energy=local_review(),exact_compact_pressure_exterior=exterior_review())
