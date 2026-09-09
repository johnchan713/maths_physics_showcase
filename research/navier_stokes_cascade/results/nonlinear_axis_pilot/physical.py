"""Full Cartesian momentum residual for the nonlinear pilot, at viscosity one.

Derivative arithmetic and a separate scalar-only finite-difference stencil
share the same finite profile. Neither checks the uncomputed global solution.
"""
import mpmath as mp
from inner import axis_series


class Jet:
    def __init__(self, value, gradient=None, hessian=None):
        self.v = mp.mpf(value)
        self.g = list(gradient) if gradient is not None else [mp.mpf(0)]*4
        self.H = [list(row) for row in hessian] if hessian is not None else [[mp.mpf(0)]*4 for _ in range(4)]

    @staticmethod
    def cast(value):
        return value if isinstance(value, Jet) else Jet(value)

    @staticmethod
    def variable(value, axis):
        result = Jet(value)
        result.g[axis] = mp.mpf(1)
        return result

    def __add__(self, other):
        other = Jet.cast(other)
        return Jet(self.v+other.v, [a+b for a,b in zip(self.g,other.g)],
                   [[self.H[i][j]+other.H[i][j] for j in range(4)] for i in range(4)])

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, [-a for a in self.g], [[-a for a in row] for row in self.H])

    def __sub__(self, other):
        return self+-Jet.cast(other)

    def __rsub__(self, other):
        return Jet.cast(other)+-self

    def __mul__(self, other):
        other = Jet.cast(other)
        return Jet(self.v*other.v, [self.g[i]*other.v+other.g[i]*self.v for i in range(4)],
                   [[self.H[i][j]*other.v+other.H[i][j]*self.v+self.g[i]*other.g[j]+other.g[i]*self.g[j]
                     for j in range(4)] for i in range(4)])

    __rmul__ = __mul__

    def __pow__(self, power):
        if power == 0:
            return Jet(1)
        if power == 1:
            return self
        return compose(self, self.v**power, power*self.v**(power-1), power*(power-1)*self.v**(power-2))

    def __truediv__(self, other):
        return self*Jet.cast(other)**-1

    def __rtruediv__(self, other):
        return Jet.cast(other)*self**-1


def compose(a, value, first, second):
    return Jet(value, [first*x for x in a.g],
               [[first*a.H[i][j]+second*a.g[i]*a.g[j] for j in range(4)] for i in range(4)])


def exponential(value):
    if isinstance(value, Jet):
        result = mp.exp(value.v)
        return compose(value, result, result, result)
    return mp.exp(value)


def polynomial(coefficients, Y, delta, eta_derivative=0, average=False):
    result = 0
    for n in reversed(range(len(coefficients))):
        angular = 0
        row = coefficients[n]
        for k in reversed(range(eta_derivative, len(row))):
            angular = angular*delta+row[k]*mp.factorial(k)/mp.factorial(k-eta_derivative)
        result = result*Y+angular/(n+1 if average else 1)
    return result


def q_value(z, time, p):
    if not (mp.isfinite(z) and mp.isfinite(time) and 0 <= time < 1):
        raise ValueError("Physical checks require 0 <= t < 1")
    tau = 1-time
    if not z:
        return tau
    guess = max(tau, abs(z)**(1/p.D))
    return mp.findroot(lambda q:q-z*z*q**(2*p.h)-tau, guess,
                       df=lambda q:1-2*p.h*z*z*q**(2*p.h-1), solver='newton')


def field(profile, point, differentiated=True):
    p = profile.parameters
    q = q_value(point[2],point[3],p)
    if differentiated:
        z = point[2]
        L = 1-2*p.h*z*z*q**(2*p.h-1)
        gradient = [mp.mpf(0),mp.mpf(0),2*z*q**(2*p.h)/L,-1/L]
        Gqq = -2*p.h*(2*p.h-1)*z*z*q**(2*p.h-2)
        Gqi = [0,0,-4*p.h*z*q**(2*p.h-1),0]
        H = [[-(Gqq*gradient[i]*gradient[j]+Gqi[i]*gradient[j]+Gqi[j]*gradient[i]
                +(-2*q**(2*p.h) if i == j == 2 else 0))/L for j in range(4)] for i in range(4)]
        q = Jet(q,gradient,H)
        x,y,z,time = [Jet.variable(v,i) for i,v in enumerate(point)]
    else:
        x,y,z,time = point
    eta = z*q**-p.D
    delta = eta-profile.eta0
    X = (x*x+y*y)/(2*q)
    Y = p.lam*X
    ph = polynomial(profile.phi,Y,delta)
    U = 4*eta+p.j+polynomial(profile.u,Y,delta)/p.lam
    avg = 4*eta+p.j+polynomial(profile.u,Y,delta,average=True)/p.lam
    avge = 4+polynomial(profile.u,Y,delta,eta_derivative=1,average=True)/p.lam
    zeta = axis_series(p,profile.eta0,profile.degree+3)['zeta']
    primitive = 0
    for k in reversed(range(len(zeta))):
        primitive = primitive*delta+zeta[k]/(k+1)
    logg = profile.log_g0+p.lam*delta*primitive
    F = exponential(logg)*ph
    v0 = (2*eta*U-2*p.D*eta*avg-(1-eta*eta)*avge)/(1-2*p.h*eta*eta)
    radial = v0/(2*q)
    swirl = q**(-p.A-mp.mpf('.5'))*F
    u = [radial*x-swirl*y, radial*y+swirl*x, q**-p.A*U]
    pressure = q**(-2*p.A)*(-p.pressure/(1+eta*eta)**2+polynomial(profile.pressure_increment,Y,delta)/p.lam)
    return u, pressure


def physical_point(profile,Y,q=mp.mpf('.2')):
    eta = profile.eta0
    if not abs(eta) < 1 or Y <= 0:
        raise ValueError("Choose an interior eta and positive radius")
    return [mp.sqrt(2*q*Y/profile.parameters.lam),mp.mpf(0),q**profile.parameters.D*eta,1-q*(1-eta*eta)]


def residual(profile,point):
    u,pressure = field(profile,point)
    terms = [[v.g[3],mp.fsum(u[j].v*v.g[j] for j in range(3)),
              -mp.fsum(v.H[j][j] for j in range(3)),pressure.g[i]] for i,v in enumerate(u)]
    total = [mp.fsum(row) for row in terms]
    scales = [max(abs(v) for v in row) for row in terms]
    if any(s == 0 for s in scales):
        raise ValueError("A component has no nonzero normalization")
    divergence_scale = max(mp.mpf(1),*(abs(u[i].g[i]) for i in range(3)))
    # At theta=0, Cartesian y is precisely the azimuthal direction.
    # Adding u_zz removes only axial diffusion from the tangential residual.
    leading = [total[i]+u[i].H[2][2] for i in (1,2)]
    return dict(total=total,terms=terms,scales=scales,
                normalized_full=[abs(v)/s for v,s in zip(total,scales)],
                normalized_leading=[abs(v)/scales[i] for v,i in zip(leading,(1,2))],
                divergence=abs(mp.fsum(u[i].g[i] for i in range(3)))/divergence_scale)


def cartesian_finite_difference(profile,point,relative_step):
    """Rejected control: mixed Cartesian samples absorb the extremely tiny swirl.

    Retained deliberately to reproduce the discovered numerical failure. This
    checker must not be used to validate this amplitude normalization.
    """
    p = profile.parameters
    if relative_step <= 0:
        raise ValueError("Positive stencil step required")
    q = q_value(point[2],point[3],p)
    eta = profile.eta0
    stiffness = max(mp.mpf(1),abs(p.lam*p.zeta(eta)),mp.sqrt(abs(p.lam*mp.diff(p.zeta,eta))))
    widths = [relative_step*w/stiffness for w in (mp.sqrt(q),mp.sqrt(q),q**p.D,q)]
    def values(where):
        u,pressure = field(profile,where,differentiated=False)
        return u+[pressure]
    center = values(point)
    first = [[mp.mpf(0)]*4 for _ in range(4)]
    laplacian = [mp.mpf(0)]*3
    for axis, step in enumerate(widths):
        samples = []
        for multiple in (-2,-1,1,2):
            shifted = list(point)
            shifted[axis] += multiple*step
            samples.append(values(shifted))
        a,b,c,d = samples
        for i in range(4):
            first[i][axis] = (a[i]-8*b[i]+8*c[i]-d[i])/(12*step)
        if axis < 3:
            for i in range(3):
                laplacian[i] += (-a[i]+16*b[i]-30*center[i]+16*c[i]-d[i])/(12*step*step)
    return [first[i][3]+mp.fsum(first[i][j]*center[j] for j in range(3))-laplacian[i]+first[3][i] for i in range(3)]


def finite_difference(profile,point,relative_step):
    """Component-separated cylindrical physical stencils, with all curvature.

    At theta=0 these components coincide with Cartesian x,y,z. Evaluating only
    on that ray avoids adding the tiny swirl to the much larger radial flow.
    This is NOT the rejected mixed-Cartesian stencil above.
    """
    p=profile.parameters
    if relative_step<=0 or point[1]!=0 or point[0]<=0:
        raise ValueError('Positive step and positive theta=0 ray required')
    q=q_value(point[2],point[3],p)
    eta=profile.eta0
    stiffness=max(mp.mpf(1),abs(p.lam*p.zeta(eta)),mp.sqrt(abs(p.lam*mp.diff(p.zeta,eta))))
    widths=[relative_step*w/stiffness for w in (mp.sqrt(q),q**p.D,q)]
    def values(where):
        u,pressure=field(profile,where,differentiated=False)
        return u+[pressure]
    center=values(point)
    first=[[mp.mpf(0)]*3 for _ in range(4)]
    second=[[mp.mpf(0)]*2 for _ in range(3)]
    for j,(axis,step) in enumerate(zip((0,2,3),widths)):
        rows=[]
        for multiple in (-2,-1,1,2):
            where=list(point);where[axis]+=multiple*step
            rows.append(values(where))
        a,b,c,d=rows
        for i in range(4):
            first[i][j]=(a[i]-8*b[i]+8*c[i]-d[i])/(12*step)
        if j<2:
            for i in range(3):
                second[i][j]=(-a[i]+16*b[i]-30*center[i]+16*c[i]-d[i])/(12*step*step)
    radius=point[0]
    result=[first[i][2]+center[0]*first[i][0]+center[2]*first[i][1]
            -second[i][0]-first[i][0]/radius-second[i][1] for i in range(3)]
    result[0]+=-center[1]**2/radius+center[0]/radius**2+first[3][0]
    result[1]+=center[0]*center[1]/radius+center[1]/radius**2
    result[2]+=first[3][1]
    return result
