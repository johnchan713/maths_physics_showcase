"""The unedited A.2 swirl schedule, A.21 pressure, and A.11 angular bumps.

All radial locations are logarithmic. Pressure masses are stored separately by
stage so that a tiny terminal contribution is never tested against a unit floor.
This is a finite numerical construction, not a certified matched outer flow.
"""
from dataclasses import dataclass
import mpmath as mp


def step(x):
    """A.5: a smooth step flat at both ends; reflect to avoid overflow."""
    if x <= 0:
        return mp.mpf(0)
    if x >= 1:
        return mp.mpf(1)
    if x > mp.mpf('.5'):
        return 1-step(1-x)
    return 1/(1+mp.exp(1/x**2-1/(1-x)**2))


def step_prime(x):
    if not 0 < x < 1:
        return mp.mpf(0)
    # Evaluating both complementary factors retains very small derivatives.
    return step(x)*step(1-x)*(2/x**3+2/(1-x)**3)


class Rule:
    """Positive composite Gauss-Legendre quadrature at the active precision."""
    def __init__(self, order=32, panels=4):
        if order not in (16, 24, 32, 48) or panels < 1:
            raise ValueError('Unsupported quadrature rule')
        self.order, self.panels = order, panels
        x, w = mp.gauss_quadrature(order, 'legendre')
        self.nodes = [((i+(1+x[j])/2)/panels, w[j]/(2*panels))
                      for i in range(panels) for j in range(order)]

    def integrate(self, f, a=0, b=1):
        a, b = mp.mpf(a), mp.mpf(b)
        return (b-a)*mp.fsum(w*f(a+(b-a)*x) for x,w in self.nodes)

    def step_integral(self, x):
        if x <= 0:
            return mp.mpf(0)
        if x >= 1:
            return x-mp.mpf('.5')
        if x > mp.mpf('.5'):
            return x-mp.mpf('.5')+self.step_integral(1-x)
        return self.integrate(step, 0, x)


def power_jet(eta, exponent, degree):
    """Taylor coefficients of (1+eta^2)^(-exponent), from its first-order ODE."""
    a = [(1+eta*eta)**(-exponent)]
    for k in range(degree):
        prev = a[k-1] if k else 0
        a.append(-(2*eta*(k+exponent)*a[k]+(k-1+2*exponent)*prev)
                 /((k+1)*(1+eta*eta)))
    return a


@dataclass
class Stage:
    name: str
    kind: str
    start: object
    length: object
    log_amplitude: object
    theta: object
    left_slope: object
    right_slope: object


class Schedule:
    def __init__(self, data, order=32, panels=4, moment_panels=8):
        self.rule = Rule(order, panels)
        self.moment_rule = Rule(order, moment_panels)
        self.md, self.lam, self.tf, self.co = [mp.mpf(data[k]) for k in ('Md','lambda','Tf','co')]
        if not (self.md > 0 and 0 < self.lam < mp.mpf('.1') and self.tf >= 200 and 0 < self.co <= mp.mpf('.001')):
            raise ValueError('Parameters outside the declared schedule domain')
        self.td = mp.exp(self.md)+10
        self.log_p = self.td+mp.mpf(data['pressure_log_offset'])
        self.h = mp.exp(-mp.mpf(data['h_exponent_factor'])*self.td)
        if not (self.log_p > self.td and 0 < self.h < min(self.lam,mp.exp(-self.td))):
            raise ValueError('Explicit A.6 inequalities not satisfied')
        self.tw, self.tu = 60*mp.log(1/self.lam), 30*mp.log(1/self.lam)
        if self.tw <= 25:
            raise ValueError('Reserved intervals do not fit')
        self.rho = self.co*self.h
        # A.13 simplifies to rho/(1-rho) times this positive integral.
        self.qp = self.rho/(1-self.rho)*self.rule.integrate(
            lambda v:mp.exp((1-self.h)*(1+2*v))*step_prime(v))
        self.q_initial = (self.lam-self.h)/(1-self.lam)
        q1 = self.propagate_q(self.q_initial,-self.lam,-1)
        self.hold = 4*mp.log(1/self.h)
        q2 = q1+(1-self.h)*self.hold
        self.q_before_wait = self.propagate_q(q2,-1,-self.h)
        self.wait = mp.log(self.q_before_wait/self.qp)/(1-self.h)
        if not self.wait > 0:
            raise ValueError('No positive terminal waiting interval')
        self.stages = []
        self._append('initial-slope-transition','ramp',1,1,mp.mpf('.6'),0)
        self._append('axial-transition','constant',self.td,1,0,0)
        self._append('intermediate-slope-transition','ramp',1,1,0,-self.lam)
        self._append('reserved-power-interval','constant',self.tw,1,-self.lam,-self.lam)
        self._append('axial-pulse-swirl','constant',13/self.lam,1,-self.lam,-self.lam)
        self._append('parameter-interpolation','interpolation',self.tf,1,-self.lam,-self.lam)
        self._append('angular-correction-interval','constant',self.tu,0,-self.lam,-self.lam)
        self._append('steep-decay-transition','ramp',1,0,-self.lam,-1)
        self._append('steep-decay-hold','constant',self.hold,0,-1,-1)
        self._append('tail-slope-transition','ramp',1,0,-1,-self.h)
        self._append('terminal-wait','constant',self.wait,0,-self.h,-self.h)
        self._append('terminal-cutoff','terminal',3,0,-self.h,-self.h)
        self._append('infinite-power-tail','constant',mp.inf,0,-self.h,-self.h)
        # Cache per-stage positive quadrature masses at active precision.
        self.components = self._pressure_components()

    def _append(self,name,kind,length,theta,left,right):
        if self.stages:
            previous = self.stages[-1]
            start = previous.start+previous.length
            loga = previous.log_amplitude+self.log_shape(previous,previous.length,0)
        else:
            start, loga = mp.mpf(0), mp.mpf(0)
        self.stages.append(Stage(name,kind,start,mp.mpf(length),loga,mp.mpf(theta),mp.mpf(left),mp.mpf(right)))

    def log_fo(self, y):
        return mp.log1p(-self.rho*step((3-y)/2))

    def log_shape(self,stage,y,eta):
        """Log(E/P*) minus the stage's stored log amplitude."""
        j = mp.log1p(eta*eta)
        if stage.kind == 'interpolation':
            theta = 1-step(y/self.tf)
            return -(mp.mpf('.5')+self.lam)*y-theta*j-(1-theta)*mp.log(2)
        if stage.kind == 'terminal':
            return -(mp.mpf('.5')+self.h)*y+self.log_fo(y)-self.log_fo(0)
        value = (stage.left_slope-mp.mpf('.5'))*y-stage.theta*j
        if stage.kind == 'ramp':
            value += (stage.right_slope-stage.left_slope)*self.rule.step_integral(y)
        return value

    def slope(self,stage,y,eta=0):
        if stage.kind == 'interpolation':
            return -self.lam-step_prime(y/self.tf)/self.tf*mp.log(2/(1+eta*eta))
        if stage.kind == 'terminal':
            return -self.h+self.rho*step_prime((y-1)/2)/(2*mp.exp(self.log_fo(y)))
        if stage.kind == 'ramp':
            return stage.left_slope+(stage.right_slope-stage.left_slope)*step(y)
        return stage.left_slope

    def propagate_q(self,initial,left,right):
        G = lambda t:(1+left)*t+(right-left)*self.rule.step_integral(t)
        integral = self.rule.integrate(lambda t:mp.exp(G(t))*(-left-(right-left)*step(t)-self.h))
        return mp.exp(-G(1))*(initial+integral)

    def _pressure_components(self):
        # A.7 contributes exactly integral exp(y/5) dy = 5.
        out = [dict(name='ideal-reference',log_scale=mp.mpf(0),nodes=[(mp.mpf(5),mp.mpf(1))])]
        for stage in self.stages:
            if stage.kind == 'constant':
                rate = 1-2*stage.left_slope
                mass = (1 if mp.isinf(stage.length) else -mp.expm1(-rate*stage.length))/rate
                nodes = [(mass,stage.theta)]
            else:
                nodes = []
                for x,w in self.rule.nodes:
                    y = stage.length*x
                    theta = 1-step(x) if stage.kind == 'interpolation' else stage.theta
                    nodes.append((stage.length*w*mp.exp(2*self.log_shape(stage,y,0)),theta))
            out.append(dict(name=stage.name,log_scale=2*stage.log_amplitude,nodes=nodes))
        return out

    def pressure_parts(self,eta,degree=0):
        parts = []
        for component in self.components:
            jets = [(weight,power_jet(eta,2*theta,degree)) for weight,theta in component['nodes']]
            # Small stage scales are kept separate from their O(1) integrals.
            normalized = [-mp.fsum(w*a[k] for w,a in jets)/2 for k in range(degree+1)]
            parts.append(dict(name=component['name'],log_scale=component['log_scale'],jet=normalized))
        return parts

    def pressure_jet(self,eta,degree,normalized=False):
        parts = self.pressure_parts(eta,degree)
        amplitude = mp.mpf(1) if normalized else mp.exp(2*self.log_p)
        return [amplitude*mp.fsum(mp.exp(v['log_scale'])*v['jet'][k] for v in parts)
                for k in range(degree+1)]

    def pressure(self,eta,normalized=False):
        return self.pressure_jet(eta,0,normalized)[0]

    def interpolation_moment_discrepancy(self,eta,rule=None):
        """r_I-1/(1-lambda) at interpolation end, with tiny differences retained."""
        beta = 1-self.lam
        def ramp_r(r,left,right):
            G = lambda t:(1+left)*t+(right-left)*self.rule.step_integral(t)
            return mp.exp(-G(1))*(r+self.rule.integrate(lambda t:mp.exp(G(t))))
        r = ramp_r(mp.mpf(5)/8,mp.mpf('.6'),mp.mpf(0))
        r = 1+(r-1)*mp.exp(-self.td)
        r = ramp_r(r,0,-self.lam)
        delta0 = (r-1/beta)*mp.exp(-beta*(self.tw+13/self.lam))
        a = mp.log(2/(1+eta*eta))
        driven = a/beta*(rule or self.moment_rule).integrate(lambda z:
            mp.exp(-beta*self.tf*(1-z)+a*step(1-z))*step_prime(z))
        return mp.exp(-beta*self.tf+a)*delta0+driven

    def angular_correction(self,eta):
        """Solve the actual A.11 discrepancy using two disjoint relative E bumps."""
        beta, rate, width = 1-self.lam, 1+2*self.lam, mp.mpf('.3')
        B, quadratic = mp.matrix(2,2), []
        for i,center in enumerate((mp.mpf(0),mp.mpf(2))):
            y = lambda z:center+width*(z-mp.mpf('.5'))
            B[0,i] = width*self.rule.integrate(lambda z:mp.exp(beta*y(z))*step_prime(z))
            B[1,i] = 2*width*self.rule.integrate(lambda z:mp.exp(-rate*y(z))*step_prime(z))
            quadratic.append(width*self.rule.integrate(lambda z:mp.exp(-rate*y(z))*step_prime(z)**2))
        discrepancy = -self.interpolation_moment_discrepancy(eta)*mp.exp(-beta*(self.tu-3))
        c = mp.matrix([0,0])
        for _ in range(8):
            q = mp.fsum(quadratic[i]*c[i]**2 for i in range(2))
            c = mp.lu_solve(B,mp.matrix([discrepancy,-q]))
        rows = [mp.fsum([B[0,i]*c[i] for i in range(2)]+[-discrepancy]),
                mp.fsum([B[1,i]*c[i] for i in range(2)]+[quadratic[i]*c[i]**2 for i in range(2)])]
        scales = [max(abs(discrepancy),*(abs(B[0,i]*c[i]) for i in range(2))),
                  max(abs(B[1,i]*c[i]) for i in range(2))]
        errors = [abs(v)/s if s else abs(v) for v,s in zip(rows,scales)]
        inverse = B**-1
        beta_bound = max(mp.fsum(abs(inverse[i,j]) for j in range(2)) for i in range(2))
        # A deliberately loose uniform bound: |delta_interp| <= 4/beta.
        d_bound = 4/beta*mp.exp(-beta*(self.tu-3))
        smallness_bound = 8*beta_bound**2*sum(quadratic)*d_bound
        return dict(coefficients=list(c),discrepancy=discrepancy,matrix=[list(B[i,:]) for i in range(2)],
                    quadratic=quadratic,relative_errors=errors,smallness_bound=smallness_bound,
                    negative_control_pressure_gap=abs(B[1,0]/B[0,0]),
                    maximum_relative_edit_bound=16*max(abs(v) for v in c))
