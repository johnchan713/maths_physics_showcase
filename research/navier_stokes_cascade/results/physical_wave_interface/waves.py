"""A compact physical fixture, exact curl carriers and the full NS residual.

This is a manufactured field on an annulus, not the manuscript's actual
band field. Two integer angular frequencies avoid mean cross products in
this fixture only. All oscillatory nonlinear interactions are retained.
"""
import math
import numpy as np
from jets import Jet, arrays, bump, cross, dot, variables

U, C0, BETA, OMEGA, NU = 1.5, -.75, -1., .37, 1.
AC = -C0*math.sqrt(1+U*U)


def in_support(coords, cartesian=False):
    radius = math.hypot(coords[0].value, coords[1].value) if cartesian else coords[0].value
    return .85 < radius < 1.95 and abs(coords[2].value) < 1.2 and .05 < coords[3].value < .95


def geometry(coords, cartesian=False):
    if cartesian:
        x, y, z, time = coords
        r = (x*x+y*y)**.5
        theta0 = math.atan2(y.value, x.value)
        along = math.cos(theta0)*x+math.sin(theta0)*y
        across = -math.sin(theta0)*x+math.cos(theta0)*y
        theta = theta0+(across/along).atan()
    else:
        r, theta, z, time = coords
    a = (1/r**2+BETA**2)**.5
    kt = [1/(r*a), BETA/a]
    nt = [kt[1], -kt[0]]
    primitive = (1+r*r)**.5+(r/(1+(1+r*r)**.5)).log()
    envelope = math.exp(3)*bump((r-1.4)/.55)*bump(z/1.2)*bump((time-.5)/.45)
    tau = .25+.1*z.sin()+.05*time.cos()
    return r, theta, z, time, a, kt, nt, primitive, envelope, tau


def to_cart(field, r, coords):
    x, y = coords[:2]
    return [(x*field[0]-y*field[1])/r, (y*field[0]+x*field[1])/r, field[2]]


def curl(field, r=None, cartesian=False):
    a, b, c = field
    if cartesian:
        return [c.diff(1)-b.diff(2), a.diff(2)-c.diff(0), b.diff(0)-a.diff(1)]
    return [c.diff(1)/r-b.diff(2), a.diff(2)-c.diff(0),
            b.diff(0)+b/r-a.diff(1)/r]


def carrier_data(coords, sign, mode='primary', cartesian=False):
    r, theta, z, time, at, kt, nt, primitive, f, tau = geometry(coords, cartesian)
    normal = [sign*U*at, 1/r, Jet(BETA)]
    direction = [Jet(1)]+[-sign*U*kt[j]-AC*nt[j] for j in range(2)]
    weight = (1-sign*tau)**.5
    if mode == 'increment':
        sn, sk = .65+.1*z.sin(), 1.1+.1*time.cos()
        signed_column = -sn/AC-sign*sk/U
        weight = signed_column/(2*weight)
    elif mode != 'primary':
        raise ValueError('Unknown amplitude mode')
    amplitude = [f*weight*d for d in direction]
    phi = theta+BETA*z+sign*U*primitive-OMEGA*time
    if cartesian:
        normal, amplitude = to_cart(normal, r, coords), to_cart(amplitude, r, coords)
    potential_coefficient = [c/dot(normal, normal) for c in cross(normal, amplitude)]
    correction = curl(potential_coefficient, r, cartesian)
    return dict(r=r, phase=phi, normal=normal, amplitude=amplitude,
                coefficient=potential_coefficient, correction=correction, envelope=f)


def carrier(coords, sign, frequency, mode='primary', cartesian=False, via_potential=False):
    data = carrier_data(coords, sign, mode, cartesian)
    phase = frequency*data['phase']
    if via_potential:
        potential = [-c*phase.sin()/frequency for c in data['coefficient']]
        return curl(potential, data['r'], cartesian)
    return [a*phase.cos()-b*phase.sin()/frequency
            for a, b in zip(data['amplitude'], data['correction'])]


def field(coords, k, mode='primary', cartesian=False, via_potential=False):
    if not in_support(coords, cartesian):
        return [Jet(0) for _ in range(3)]
    plus = carrier(coords, 1, k, mode, cartesian, via_potential)
    minus = carrier(coords, -1, 2*k, mode, cartesian, via_potential)
    return [a+b for a, b in zip(plus, minus)]


def wave_pressure(coords, k, cartesian=False):
    if not in_support(coords, cartesian):
        return Jet(0)
    data = carrier_data(coords, 1, cartesian=cartesian)
    return .13*data['envelope']*(k*data['phase']).sin()


def background(coords, cartesian=False):
    r, theta, z, time = geometry(coords, cartesian)[:4]
    alpha = .17
    omega = .43+.31*r*r/2+.19*z+.11*time
    velocity = [-alpha*r, r*omega, 2*alpha*z]
    if cartesian:
        velocity = to_cart(velocity, r, coords)
    pressure = .11*r*r+.07*r*z+.13*z*z+.03*time*r
    return velocity, pressure


def divergence(data, r=None, cartesian=False):
    val, grad, _ = data
    if cartesian:
        return grad[..., 0, 0]+grad[..., 1, 1]+grad[..., 2, 2]
    return grad[..., 0, 0]+val[..., 0]/r+grad[..., 1, 1]/r+grad[..., 2, 2]


def advection(left, right, r=None, cartesian=False):
    """(left dot grad)right, including every cylindrical frame derivative."""
    lv, rv, rg = left[0], right[0], right[1]
    if cartesian:
        return np.einsum('...j,...ij->...i', lv, rg[..., :3])
    directional = (lv[..., 0, None]*rg[..., :, 0]+lv[..., 1, None]*rg[..., :, 1]/r
                   +lv[..., 2, None]*rg[..., :, 2])
    frame = np.stack([-lv[..., 1]*rv[..., 1]/r, lv[..., 1]*rv[..., 0]/r,
                      np.zeros_like(lv[..., 0]+rv[..., 0])], axis=-1)
    return directional+frame


def laplacian(data, r=None, cartesian=False):
    val, grad, hess = data
    if cartesian:
        return sum(hess[..., :, j, j] for j in range(3))
    scalar = hess[..., :, 0, 0]+grad[..., :, 0]/r+hess[..., :, 1, 1]/r**2+hess[..., :, 2, 2]
    frame = np.stack([-val[..., 0]/r**2-2*grad[..., 1, 1]/r**2,
                      -val[..., 1]/r**2+2*grad[..., 0, 1]/r**2,
                      np.zeros_like(val[..., 0])], axis=-1)
    return scalar+frame


def pressure_gradient(pressure, r=None, cartesian=False):
    grad = arrays([pressure])[1][0, :3]
    return grad if cartesian else grad*np.array([1, 1/r, 1])


def momentum(data, pressure_grad, r=None, cartesian=False):
    return (data[1][..., :, 3]+advection(data, data, r, cartesian)
            -NU*laplacian(data, r, cartesian)+pressure_grad)


def covariance(coords, k, mode='primary'):
    q = [[Jet(0) for _ in range(3)] for _ in range(3)]
    q0 = [[Jet(0) for _ in range(3)] for _ in range(3)]
    if not in_support(coords):
        return q, q0
    for sign, freq in ((1, k), (-1, 2*k)):
        data = carrier_data(coords, sign, mode)
        a, b = data['amplitude'], data['correction']
        for i in range(3):
            for j in range(3):
                q0[i][j] = q0[i][j]+a[i]*a[j]/2
                q[i][j] = q[i][j]+(a[i]*a[j]+b[i]*b[j]/freq**2)/2
    return q, q0


def tensor_values(q):
    return np.array([[x.value for x in row] for row in q])


def tensor_divergence(q, r):
    r = float(r)
    return np.array([q[0][0].diff(0).value+(q[0][0].value-q[1][1].value)/r+q[0][2].diff(2).value,
                     q[0][1].diff(0).value+2*q[0][1].value/r+q[1][2].diff(2).value,
                     q[0][2].diff(0).value+q[0][2].value/r+q[2][2].diff(2).value])


def harmonic_samples(point, k, count=None, mode='primary', principal_only=False, same_frequency=False):
    """Trapezoidal angular samples from exact local harmonic derivative arrays."""
    coords = variables((point[0], 0., point[2], point[3]))
    count = count or 8*k+1
    angles = np.arange(count)*(2*np.pi/count)
    out = [np.zeros((count, 3)), np.zeros((count, 3, 4)), np.zeros((count, 3, 4, 4))]
    for sign, freq in ((1, k), (-1, k if same_frequency else 2*k)):
        data = carrier_data(coords, sign, mode)
        phase = freq*data['phase']
        a = data['amplitude']
        b = [Jet(0)]*3 if principal_only else data['correction']
        cosine = [x*phase.cos()-y*phase.sin()/freq for x, y in zip(a, b)]
        sine = [-x*phase.sin()-y*phase.cos()/freq for x, y in zip(a, b)]
        for i, (c, s) in enumerate(zip(arrays(cosine), arrays(sine))):
            shape = (count,)+(1,)*c.ndim
            out[i] += np.cos(freq*angles).reshape(shape)*c+np.sin(freq*angles).reshape(shape)*s
    # The chosen oscillatory pressure has the plus carrier's angular frequency.
    pdata = carrier_data(coords, 1)
    pc = .13*pdata['envelope']*(k*pdata['phase']).sin()
    ps = .13*pdata['envelope']*(k*pdata['phase']).cos()
    pgrad = (np.cos(k*angles)[:, None]*pressure_gradient(pc, point[0])
             +np.sin(k*angles)[:, None]*pressure_gradient(ps, point[0]))
    return tuple(out), pgrad
