"""Full physical localization identity; compact manufactured wave/mean curls.

The fixture combines nonaxisymmetric inherited carriers with a smooth
Cartesian base. It is not the actual global field or auxiliary geometry.
"""
from pathlib import Path
import sys
import numpy as np

PARENT = Path(__file__).resolve().parent.parent
sys.path.append(str(PARENT/'full_potential_realization'))
from cutoff_checks import chi, representatives
from jets import Jet, arrays, cross, dot, variables
import waves

POINTS = ((1.12, -.4, -.19, .64), (1.3, .27, .25, .72), (1.44, -.9, .31, .8))


def relative(a, b):
    return float(np.linalg.norm(a-b)/(1+np.linalg.norm(a)+np.linalg.norm(b)))


def spatial_gradient(f, r=None, cartesian=False):
    grad = [f.diff(i) for i in range(3)]
    if not cartesian:
        grad[1] = grad[1]/r
    return grad


def scalar_laplacian(f, r=None, cartesian=False):
    if cartesian:
        return sum(f.diff(i).diff(i) for i in range(3))
    return f.diff(0).diff(0)+f.diff(0)/r+f.diff(1).diff(1)/r**2+f.diff(2).diff(2)


def physical_fixture(coords, cartesian=False):
    A, B, p, _ = representatives(coords, 19, cartesian)
    if cartesian:
        x, y, z, time = coords
        radius_squared = x*x+y*y
        smooth_theta = [-y, x, Jet(0)]
    else:
        r, theta, z, time = coords
        radius_squared = r*r
        smooth_theta = [Jet(0), r, Jet(0)]
    factor = (-.2*radius_squared-.3*z*z).exp()*(.8+.2*time)
    A = [a+.12*factor*e for a, e in zip(A, smooth_theta)]
    B = [b+.07*factor*e for b, e in zip(B, smooth_theta)]
    p = p+.13*factor*radius_squared+.17*z*factor
    return A, B, p


def physical_cutoff(coords, cartesian=False):
    if cartesian:
        r2 = coords[0]**2+coords[1]**2
    else:
        r2 = coords[0]**2
    z, time = coords[2:]
    return chi(r2/1.55**2)*chi(z*z/.35**2)*chi((1-time)/.45)


def momentum(data, pressure, radius, cartesian, nu):
    return (data[1][:, 3]+waves.advection(data, data, radius, cartesian)
            -nu*waves.laplacian(data, radius, cartesian)
            +waves.pressure_gradient(pressure, radius, cartesian))


def localization_data(coords, cartesian=False):
    radius = None if cartesian else coords[0]
    A, B, p = physical_fixture(coords, cartesian)
    c = physical_cutoff(coords, cartesian)
    u0 = [a+b for a, b in zip(waves.curl(A, radius, cartesian), B)]
    grad_c = spatial_gradient(c, radius, cartesian)
    d = cross(grad_c, A)
    localized = [a+c*b for a, b in zip(waves.curl([c*a for a in A], radius, cartesian), B)]
    return dict(A=A, B=B, p=p, c=c, grad_c=grad_c, d=d, u0=u0,
                localized=localized, localized_pressure=c*p)


def localization_review():
    worst = cart_error = product_error = div_error = naive_div = orthogonality = angular_failure = 0.
    sizes = {name: 0. for name in ('base_residual', 'cutoff_time', 'curl_time',
             'cutoff_laplacian', 'gradient_viscosity', 'curl_viscosity', 'pressure_cutoff',
             'quadratic_rescaling', 'velocity_cutoff_transport', 'old_curl_transport',
             'curl_old_transport', 'curl_square')}
    for point in POINTS:
        r, theta, z, time = point
        cylindrical = variables(point)
        cartesian_coords = variables((r*np.cos(theta), r*np.sin(theta), z, time))
        rotation = np.array([[np.cos(theta), -np.sin(theta), 0.],
                             [np.sin(theta), np.cos(theta), 0.], [0., 0., 1.]])
        for nu in (.4, 1., 1.7):
            for coords, cartesian in ((cylindrical, False), (cartesian_coords, True)):
                radius = None if cartesian else r
                data = localization_data(coords, cartesian)
                A, B, p, c, d, grad_c = [data[k] for k in ('A', 'B', 'p', 'c', 'd', 'grad_c')]
                U, D, cut = arrays(data['u0']), arrays(d), arrays(data['localized'])
                product_error = max(product_error, *(relative(a,b) for a,b in zip(cut, arrays([c*v+dd for v,dd in zip(data['u0'],d)]))))
                div_error = max(div_error, abs(waves.divergence(cut, radius, cartesian))/(1+np.linalg.norm(cut[1])))
                naive_div = max(naive_div, abs(waves.divergence(arrays([c*v for v in data['u0']]), radius, cartesian)))
                orthogonality = max(orthogonality, abs(dot(d,grad_c).value))
                # (grad c . grad)u uses the full vector/cylindrical connection.
                grad_c_data = arrays(grad_c)
                direct = momentum(cut, c*p, radius, cartesian, nu)
                gc = grad_c_data[0]
                terms = dict(base_residual=c.value*momentum(U,p,radius,cartesian,nu),
                    cutoff_time=c.diff(3).value*U[0], curl_time=D[1][:,3],
                    cutoff_laplacian=-nu*scalar_laplacian(c,radius,cartesian).value*U[0],
                    gradient_viscosity=-2*nu*waves.advection(grad_c_data,U,radius,cartesian),
                    curl_viscosity=-nu*waves.laplacian(D,radius,cartesian),
                    pressure_cutoff=p.value*gc,
                    quadratic_rescaling=c.value*(c.value-1)*waves.advection(U,U,radius,cartesian),
                    velocity_cutoff_transport=c.value*np.dot(U[0],gc)*U[0],
                    old_curl_transport=c.value*waves.advection(U,D,radius,cartesian),
                    curl_old_transport=c.value*waves.advection(D,U,radius,cartesian),
                    curl_square=waves.advection(D,D,radius,cartesian))
                worst = max(worst,relative(direct,sum(terms.values())))
                for name,value in terms.items():
                    sizes[name] = max(sizes[name],float(np.linalg.norm(value)))
                if not cartesian:
                    cyl = direct
                else:
                    cart_error = max(cart_error,relative(direct,rotation@cyl))
                    # Nonaxisymmetric localization of direct B would add divergence.
                    bad = c*(1+.1*coords[0])
                    bad_velocity = [a+bad*b for a,b in zip(waves.curl([bad*a for a in A],None,True),B)]
                    angular_failure = max(angular_failure,abs(waves.divergence(arrays(bad_velocity),None,True)))
    return dict(manufactured=True, positive_viscosities=[.4,1.,1.7],
                complete_localization_identity_error=worst,
                independent_cartesian_cylindrical_residual_error=cart_error,
                velocity_product_rule_error=product_error, localized_divergence_error=div_error,
                naive_velocity_cutoff_divergence=naive_div,
                d_dot_grad_c_error=orthogonality,
                nonaxisymmetric_direct_swirl_cutoff_divergence=angular_failure,
                retained_term_sizes=sizes,
                checks=dict(full_localization_identity_with_every_term=worst < 2e-10,
                            independent_cartesian_residual_agrees=cart_error < 2e-10,
                            all_velocity_jets_through_order_two_match_product_rule=product_error < 2e-10,
                            exact_localized_divergence=div_error < 2e-12,
                            each_retained_term_is_nonzero=min(sizes.values()) > 1e-7,
                            omitted_curl_cutoff_breaks_divergence=naive_div > 1e-4,
                            d_dot_grad_c_is_exactly_zero=orthogonality < 2e-12,
                            axisymmetry_is_required_for_direct_swirl=angular_failure > 1e-5))


def support_axis_review():
    rows=[]
    for point in ((0.,0.,.1,.72),(.1,-.1,.05,.8),(1.8,0.,.1,.72),
                  (.5,.2,.5,.72),(.5,.2,.1,.3),(.15,.1,.05,.9)):
        coords=variables(point)
        data=localization_data(coords,True)
        v=arrays(data['localized'])
        rows.append(dict(point=list(point),cutoff=data['c'].value,
                         velocity_norm=float(np.linalg.norm(v[0])),
                         velocity_jet_max=float(max(np.max(abs(a)) for a in v)),
                         divergence=abs(float(waves.divergence(v,None,True)))))
    return dict(rows=rows,
                checks=dict(nonzero_smooth_axis_velocity=rows[0]['velocity_norm'] > 1e-3,
                            Cartesian_axis_has_zero_divergence=rows[0]['divergence'] < 2e-12,
                            spatial_and_initial_zero_extensions_are_exact=all(rows[i]['velocity_jet_max']==0 for i in (2,3,4)),
                            inner_cutoff_plateau_preserves_actual_field=rows[-1]['cutoff']==1.))
