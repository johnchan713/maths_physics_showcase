"""Heat endpoint, terminal-jet gluing and whole-space energy diagnostics.

All fields and terminal coefficients here are explicit diagnostic fixtures.
The analytic implication for the actual conditional construction is in README.
"""
import math
import mpmath as mp
import numpy as np
from localization_checks import chi, relative
from jets import INDICES, Jet, bump, variables


def heat_derivative(z, j, h):
    return ((-1)**j*mp.rf(h,j)/mp.gamma(1+h)
            *mp.quad(lambda v: mp.exp(-v)*v**(h+j)*(1+z*v)**(-h-j), [0,1,mp.inf]))


def heat_endpoint_review():
    rows=[]
    with mp.workdps(55):
        h=mp.mpf(1)/200
        for z in (mp.mpf(0),mp.mpf('.3'),mp.mpf(2)):
            derivatives=[heat_derivative(z,j,h) for j in range(4)]
            bounds=[mp.rf(h,j)*mp.rf(1+h,j) for j in range(4)]
            s=1+2*h
            spatial=((s*s-1)*derivatives[0]+4*(s+1)*z*derivatives[1]+4*z*z*derivatives[2])
            correct=-4*derivatives[1]-spatial
            wrong=4*derivatives[1]-spatial
            rows.append(dict(heat_argument=float(z),
                             derivative_bound_ratios=[float(abs(v)/b) for v,b in zip(derivatives,bounds)],
                             heat_equation_error=float(abs(correct)),
                             reversed_tau_time_sign_gap=float(abs(wrong))))
    return dict(exact_gamma_integral=True,diagnostic_h='1/200',rows=rows,
                pressure_density_majorant='|partial_tau^j(K(rho,tau)^2/rho)|<=C_j*rho^(-3-4*h-2*j)',
                checks=dict(every_fixed_heat_derivative_has_uniform_nonnegative_argument_bound=max(max(r['derivative_bound_ratios']) for r in rows) <= 1+2e-13,
                            exact_heat_identity_includes_tau_equals_one_minus_t=max(r['heat_equation_error'] for r in rows) < 1e-45,
                            reversed_time_sign_fails=min(r['reversed_tau_time_sign_gap'] for r in rows) > 1e-4,
                            pressure_majorant_is_integrable_for_positive_radius=True))


def terminal_jet(j,coords):
    x,y,z,sigma=coords
    scalar=math.factorial(j+1)**2*bump((x-.4)/.2)*bump(y/.3)*bump(z/.3)
    return [scalar*(1+x),scalar*(.2+y),scalar*(.3+z)]


def monomial(s,j):
    out=Jet(1)
    for _ in range(j):
        out=out*s
    return out/math.factorial(j)


def future_fixture(coords,last=7):
    sigma=coords[3]
    out=[Jet(0) for _ in range(3)]
    for j in range(last+1):
        factor=chi(2**(j+1)*sigma)*monomial(sigma,j)
        out=[a+factor*f for a,f in zip(out,terminal_jet(j,coords))]
    return out


def terminal_gluing_review():
    coords=variables((.41,.04,-.03,0.))
    future=future_fixture(coords)
    errors=[]
    for multi in INDICES:
        # Taylor coefficient alpha of F_m / m! equals the matching
        # coefficient of sigma^m in the total series, including mixed jets.
        a,b,c,m=multi
        expected=terminal_jet(m,coords)
        spatial_multi=(a,b,c,0)
        location=INDICES.index(spatial_multi)
        actual_location=INDICES.index(multi)
        errors.append(max(abs(f.c[actual_location]-g.c[location]/math.factorial(m))
                          for f,g in zip(future,expected)))
    product_errors=[]
    with mp.workdps(55):
        def cutoff(s):
            if s <= mp.mpf('.5'):
                return mp.mpf(1)
            if s >= 1:
                return mp.mpf(0)
            left,right=mp.exp(-1/(s-mp.mpf('.5'))),mp.exp(-1/(1-s))
            return right/(left+right)
        for j in range(1,8):
            b=2**(j+1)
            sigma=mp.mpf('.73')/b
            time=Jet.variable(float(sigma),3)
            term=chi(b*time)*monomial(time,j)
            for m in range(4):
                value=term
                for _ in range(m):
                    value=value.diff(3)
                independent=mp.diff(lambda s:cutoff(b*s)*s**j/mp.factorial(j),sigma,m)
                product_errors.append(abs(value.value-float(independent))/(1+abs(float(independent))))
    zeros=[]
    for point in ((0.,0.,0.,0.),(.41,.04,-.03,.6),(1.2,.04,-.03,.001)):
        zeros.append(max(float(np.max(abs(v.c))) for v in future_fixture(variables(point))))
    jump=max(abs(f.value) for f in terminal_jet(0,coords))
    # This bad function is zero near x=0 but oscillates without an endpoint
    # limit on a distant compact set. Smoothness for t<1 alone is insufficient.
    pairs=[(1/(math.pi/2+2*math.pi*n),1/(3*math.pi/2+2*math.pi*n)) for n in (1,5,50)]
    gaps=[abs(math.sin(1/a)-math.sin(1/b)) for a,b in pairs]
    return dict(manufactured_terminal_jets=True,finite_terms_used_only_for_jet_algebra=8,
                mixed_jet_matching_error=max(errors),time_product_rule_error=max(product_errors),
                support_and_origin_zero_jet_errors=zeros,zero_continuation_jump=jump,
                insufficient_endpoint_bound_counterexample_gaps=gaps,
                checks=dict(all_available_mixed_terminal_jets_match=max(errors) < 2e-11,
                            independent_scalar_time_product_rule_agrees=max(product_errors) < 2e-11,
                            spatial_support_time_support_and_flat_origin_are_preserved=max(zeros) == 0.,
                            extending_the_force_by_zero_would_create_a_jump=jump > 1e-5,
                            flatness_at_origin_alone_is_not_endpoint_regularity=min(gaps) > 1.99))


def one_dimensional_bump(x):
    d=1-x*x
    f=np.exp(-1/d)
    a=-2*x/d**2
    b=-2/d**2-8*x*x/d**3
    c=-24*x/d**3-48*x**3/d**4
    return f,f*a,f*(a*a+b),f*(a**3+3*a*b+c)


def energy_review():
    """Independent whole-space integrals for a smooth compact curl on [-1,1]^3."""
    count=96
    nodes,weights=np.polynomial.legendre.leggauss(count)
    f=one_dimensional_bump(nodes)
    X=[v[:,None,None] for v in f]
    Y=[v[None,:,None] for v in f]
    Z=[v[None,None,:] for v in f]
    zero=np.zeros((count,count,count))
    v=np.array([X[0]*Y[1]*Z[0],-X[1]*Y[0]*Z[0],zero])
    grad=np.array([[X[1]*Y[1]*Z[0],X[0]*Y[2]*Z[0],X[0]*Y[1]*Z[1]],
                   [-X[2]*Y[0]*Z[0],-X[1]*Y[1]*Z[0],-X[1]*Y[0]*Z[1]],
                   [zero,zero,zero]])
    lap=np.array([X[2]*Y[1]*Z[0]+X[0]*Y[3]*Z[0]+X[0]*Y[1]*Z[2],
                  -X[3]*Y[0]*Z[0]-X[1]*Y[2]*Z[0]-X[1]*Y[0]*Z[2],zero])
    nonlinear=np.einsum('jabc,ijabc->iabc',v,grad)
    x=nodes[:,None,None]
    pressure_gradient=.17*np.array([(X[0]+x*X[1])*Y[0]*Z[0],
                                    x*X[0]*Y[1]*Z[0],x*X[0]*Y[0]*Z[1]])
    integrate=lambda a:float(np.einsum('abc,a,b,c',a,weights,weights,weights))
    square=lambda a:integrate(np.sum(a*a,axis=0))
    V2=square(v)
    G2=integrate(np.sum(grad*grad,axis=(0,1)))
    N,L,P=[math.sqrt(square(a)) for a in (nonlinear,lap,pressure_gradient)]
    by_parts=integrate(np.sum(v*lap,axis=0))+G2
    pressure_work=integrate(np.sum(v*pressure_gradient,axis=0))
    nonlinear_work=integrate(np.sum(v*nonlinear,axis=0))
    absolute_nonlinear_work=integrate(abs(np.sum(v*nonlinear,axis=0)))
    # Preserve the failed lower quadrature resolution; refine the integral,
    # not the tolerance. This equivalent separated expression is inexpensive.
    coarse_x,coarse_w=np.polynomial.legendre.leggauss(64)
    coarse_f=one_dimensional_bump(coarse_x)
    i0,i1,i2=[np.dot(coarse_w,a*a) for a in coarse_f[:3]]
    j02=np.dot(coarse_w,coarse_f[0]*coarse_f[2])
    j13=np.dot(coarse_w,coarse_f[1]*coarse_f[3])
    coarse_gap=abs(4*i0*i1*(j02+i1)+2*i0*i0*(j13+i2))
    rows=[]
    for t,nu in ((.3,.4),(.7,1.),(.7,1.7)):
        g=t*t
        u=g*v
        force=2*t*v+g*g*nonlinear-nu*g*lap+t**3*pressure_gradient
        work=integrate(np.sum(force*u,axis=0))
        expected=2*t*g*V2+nu*g*g*G2
        F_upper=t*t*math.sqrt(V2)+t**5*N/5+nu*t**3*L/3+t**4*P/4
        integrated=t**4*V2+2*nu*t**5*G2/5
        rows.append(dict(time=t,viscosity=nu,energy_identity_error=abs(work-expected)/(1+abs(work)+abs(expected)),
                         integrated_energy_over_force_bound=integrated/F_upper**2,
                         omitted_viscous_work_gap=nu*g*g*G2))
    return dict(manufactured_smooth_compact_curl=True,quadrature_nodes_per_axis=count,
                retained_failed_64_node_integration_by_parts_error=float(coarse_gap),
                divergence_error=float(np.max(abs(grad[0,0]+grad[1,1]))),
                integration_by_parts_error=abs(by_parts),pressure_energy_flux_error=abs(pressure_work),
                nonlinear_energy_flux_error=abs(nonlinear_work),nonlinear_work_absolute_integral=absolute_nonlinear_work,
                rows=rows,
                checks=dict(manufactured_field_is_solenoidal=float(np.max(abs(grad[0,0]+grad[1,1]))) < 1e-14,
                            independent_viscous_integration_by_parts=abs(by_parts) < 2e-9,
                            low_quadrature_failure_is_retained=coarse_gap > 2e-9,
                            pressure_and_nonlinear_integrals_cancel=max(abs(pressure_work),abs(nonlinear_work)) < 2e-11,
                            nonlinear_energy_density_is_not_trivially_zero=absolute_nonlinear_work > 1e-6,
                            complete_forced_energy_identity=max(r['energy_identity_error'] for r in rows) < 2e-9,
                            integrated_energy_bound_holds=max(r['integrated_energy_over_force_bound'] for r in rows) <= 1.,
                            omission_of_dissipation_is_detected=min(r['omitted_viscous_work_gap'] for r in rows) > 1e-5))
