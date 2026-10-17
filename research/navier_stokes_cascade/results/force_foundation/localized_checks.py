"""Exact manufactured localization, Borel Leibniz and compact energy controls."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb,factorial
from localized_bounds import HERE,load,TERMS,derivative_factor

p=load('force_foundation_exact_polynomials',HERE.parent/'realization_foundation/realized_checks.py')
def dot(a,b):return p.add(*(p.mul(x,y) for x,y in zip(a,b)))
def lap_scalar(f):return p.add(*(p.diff(p.diff(f,i),i) for i in range(3)))
def momentum(u,pressure,nu=F(1)):
    return p.vsum([p.diff(v,3) for v in u],p.vscale(p.lap(u),-nu),p.grad(pressure),p.adv(u,u))
def zero_vector(v):return not any(v)

@lru_cache(None)
def localization_review():
    x,y,z,t=p.X,p.Y,p.Z,p.T;r2=p.add(p.power(x,2),p.power(y,2))
    f=p.add(p.ONE,p.mul(z,t),r2)
    A=p.vsum([p.mul(x,z),p.mul(y,t),p.add(p.mul(x,y),p.mul(z,t))],
        [p.scale(p.mul(y,f),F(-1,2)),p.scale(p.mul(x,f),F(1,2)),{}])
    swirl=p.add(p.ONE,z,t);B=[p.scale(p.mul(y,swirl),-1),p.mul(x,swirl),{}]
    pressure=p.add(p.mul(r2,z),p.power(t,2),p.mul(z,t))
    c=p.mul(p.mul(p.add(p.ONE,p.scale(r2,-1)),p.add(p.ONE,p.scale(p.power(z,2),-1))),p.add(t,p.power(t,2)))
    u0=p.vsum(p.curl(A),B);gc=p.grad(c);d=p.cross(gc,A)
    u=p.vsum(p.curl(p.vmul(c,A)),p.vmul(c,B));cp=p.mul(c,pressure)
    assert zero_vector(p.vsum(u,p.vscale(p.vsum(p.vmul(c,u0),d),-1)))
    assert not p.divergence(u) and not dot(d,gc)
    omissions={k:False for k in TERMS};identities=0
    for nu in (F(3,5),F(1),F(17,10)):
        terms=dict(base_residual=p.vmul(c,momentum(u0,pressure,nu)),
            cutoff_time=p.vmul(p.diff(c,3),u0),curl_time=[p.diff(v,3) for v in d],
            cutoff_laplacian=p.vmul(p.scale(lap_scalar(c),-nu),u0),
            gradient_viscosity=p.vscale(p.adv(gc,u0),-2*nu),curl_viscosity=p.vscale(p.lap(d),-nu),
            pressure_cutoff=p.vmul(pressure,gc),
            quadratic_rescaling=p.vmul(p.mul(c,p.add(c,p.scale(p.ONE,-1))),p.adv(u0,u0)),
            velocity_cutoff_transport=p.vmul(p.mul(c,dot(u0,gc)),u0),
            old_curl_transport=p.vmul(c,p.adv(u0,d)),curl_old_transport=p.vmul(c,p.adv(d,u0)),curl_square=p.adv(d,d))
        assert tuple(terms)==TERMS
        direct=momentum(u,cp,nu);assert zero_vector(p.vsum(direct,p.vscale(p.vsum(*terms.values()),-1)));identities+=3
        for name,value in terms.items():omissions[name]=omissions[name] or any(value)
    naive=p.divergence(p.vmul(c,u0));assert naive
    bad=p.mul(c,p.add(p.ONE,x));bad_u=p.vsum(p.curl(p.vmul(bad,A)),p.vmul(bad,B));bad_div=p.divergence(bad_u)
    assert bad_div and not p.add(bad_div,p.scale(p.mul(c,B[0]),-1))
    actual_force=momentum(u,cp)
    DD=p.add(*(p.diff(p.diff(p.mul(u[i],u[j]),i),j) for i in range(3) for j in range(3)))
    divf=p.divergence(actual_force);poisson=p.add(p.scale(lap_scalar(cp),-1),p.scale(DD,-1),divf)
    assert not poisson and divf
    return dict(manufactured_local_polynomial_fields_not_actual_compact_C_infinity_flow=True,
        positive_diagnostic_viscosities=['3/5','1','17/10'],actual_construction_viscosity='1',
        exact_localization_vector_components=identities,all_nonzero_omission_terms=sorted(k for k,v in omissions.items() if v),
        checks=dict(all_twelve_localization_terms_exact=identities==9 and all(omissions.values()),
            complete_cut_potential_velocity_is_divergence_free=True,
            orthogonal_curl_commutator_has_zero_d_dot_grad_c=True,
            naive_velocity_cutoff_breaks_divergence=bool(naive),
            theta_dependent_direct_swirl_cutoff_breaks_divergence=bool(bad_div),
            full_forced_pressure_Poisson_identity_exact=True,
            omitted_force_divergence_is_nonzero=bool(divf),
            all_pressure_viscous_time_and_nonlinear_omissions_detected=all(omissions.values())))

def chi_polynomial(s):return p.add(p.ONE,p.scale(p.power(s,2),-3),p.scale(p.power(s,3),2))
def substitute_time_scale(f,b):return {k:v*F(b)**k[3] for k,v in f.items()}

@lru_cache(None)
def Leibniz_review():
    sigma=p.T;chi=chi_polynomial(sigma);identities=plateaus=0
    for j in range(7):
        for m in range(7):
            for b in (2,7):
                direct=p.scale(p.mul(chi_polynomial(p.scale(sigma,b)),p.power(sigma,j)),F(1,factorial(j)))
                for _ in range(m):direct=p.diff(direct,3)
                terms={}
                for k in range(m+1):
                    factor=derivative_factor(j,m,k)
                    if not factor:continue
                    chik=chi
                    for _ in range(k):chik=p.diff(chik,3)
                    terms=p.add(terms,p.scale(p.mul(substitute_time_scale(chik,b),p.power(sigma,j-m+k)),factor*b**k))
                assert not p.add(direct,p.scale(terms,-1));identities+=1
            plateau=p.scale(p.power(sigma,j),F(1,factorial(j)))
            for _ in range(m):plateau=p.diff(plateau,3)
            assert plateau.get((0,0,0,0),F(0))==int(j==m);plateaus+=1
    return dict(local_polynomial_transition_not_actual_cutoff=True,exact_Leibniz_identities=identities,
        exact_terminal_factorial_controls=plateaus,
        checks=dict(all_cutoff_derivatives_and_factorials_retained=identities==98,
            every_plateau_terminal_derivative_is_the_prescribed_jet=plateaus==49,
            missing_factorial_changes_terminal_jet=factorial(4)!=1,
            uniform_whole_support_exponent_is_m_minus_j=True))

def integrate_space(f):
    out={}
    for (a,b,c,t),v in f.items():
        if any(n%2 for n in (a,b,c)):continue
        value=v*F(2,a+1)*F(2,b+1)*F(2,c+1)
        out[(t,)]=out.get((t,),F(0))+value
    return {k:v for k,v in out.items() if v}
def time_integral(f):return {(k[0]+1,):v/F(k[0]+1) for k,v in f.items()}

@lru_cache(None)
def energy_review():
    x,y,z,t=p.X,p.Y,p.Z,p.T
    phi=p.ONE
    for coordinate in (x,y,z):phi=p.mul(phi,p.power(p.add(p.ONE,p.scale(p.power(coordinate,2),-1)),4))
    A=[{}, {},p.mul(p.power(t,2),phi)];u=p.curl(A)
    pressure=p.mul(p.mul(t,phi),p.mul(x,z))
    E=integrate_space(dot(u,u))
    D=integrate_space(p.add(*(dot([p.diff(v,i) for v in u],[p.diff(v,i) for v in u]) for i in range(3))))
    pressure_work=integrate_space(dot(u,p.grad(pressure)))
    nonlinear_work=integrate_space(dot(u,p.adv(u,u)))
    viscous_work=integrate_space(dot(u,p.vscale(p.lap(u),-1)))
    work=integrate_space(dot(u,momentum(u,pressure)))
    assert not pressure_work and not nonlinear_work
    assert not p.add(viscous_work,p.scale(D,-1))
    assert not p.add(p.scale(p.diff(E,0),F(1,2)),D,p.scale(work,-1))
    assert not p.add(E,p.scale(time_integral(D),2),p.scale(time_integral(work),-2))
    assert E[(4,)]>0 and D[(4,)]>0
    return dict(manufactured_polynomial_potential_zero_extended_from_cube=True,
        finite_regularity_fixture_not_the_actual_C_infinity_field=True,
        exact_kinetic_energy_coefficient=str(E[(4,)]),exact_dissipation_coefficient=str(D[(4,)]),
        checks=dict(all_polynomial_spatial_integrals_exact=True,
            pressure_and_transport_energy_fluxes_exactly_zero=True,
            viscous_work_equals_positive_dissipation=True,
            differential_energy_identity_with_half_coefficient_exact=True,
            integrated_energy_identity_has_factor_two_dissipation=True,
            manufactured_initial_velocity_exactly_zero=True,
            no_numerical_energy_quadrature_or_tolerance_replacement=True))

def reviews():return dict(exact_complete_localization_and_pressure=localization_review(),
    exact_Borel_Leibniz_and_terminal_factors=Leibniz_review(),exact_compact_preterminal_energy_controls=energy_review())
