"""Finite Md=4 diagnostics using moments and a separate source-ODE assembly."""
import json
import math
from pathlib import Path
import sys
import mpmath as mp
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / 'outer_pressure_pilot'
sys.path.append(str(PREVIOUS))
from schedule import Schedule, step, step_prime  # noqa: E402


class Moments:
    def __init__(self, order=32):
        self.protocol = json.loads((PREVIOUS/'protocol.json').read_text())
        self.s = Schedule(self.protocol['parameters'], order)
        s, r = self.s, self.s.rule
        shape = lambda t: s.log_shape(s.stages[0], t, 0)
        self.rI1 = (mp.mpf(5)/8+r.integrate(lambda t: mp.exp(mp.mpf('1.5')*t+shape(t))))/mp.exp(mp.mpf('1.3'))
        self.energy1 = (mp.mpf(5)/6+r.integrate(lambda t: mp.exp(t+2*shape(t))))/mp.exp(mp.mpf('.6'))
        self.pressure1 = (5+r.integrate(lambda t: mp.exp(2*shape(t))))/2
        self.cache = {}

    def k(self, y):
        return 4*step(1-mp.log1p(y)/self.s.md)

    def kp(self, y):
        return -4*step_prime(mp.log1p(y)/self.s.md)/(self.s.md*(1+y))

    def averages(self, y):
        if y not in self.cache:
            end = min(y, mp.expm1(self.s.md))
            r = self.s.rule
            self.cache[y] = tuple(initial*mp.exp(-y)+r.integrate(
                lambda z: mp.exp(z-y)*self.k(z)**power, 0, end)
                for initial, power in ((4, 1), (16, 2)))
        return self.cache[y]

    def axial(self, y, eta):
        s, e = self.s, mp.mpf(eta)
        y = mp.mpf(y)
        if not 0 <= y <= s.td or not -1 <= e <= 1:
            raise ValueError('Point outside the axial interval')
        z, d, h = e*e, 1-e*e, s.h
        A, D, L = mp.mpf('.5')+h, mp.mpf('.5')-h, 1-2*h*z
        Jp = 2*e/(1+z)
        K, K2 = self.averages(y)
        memory = (self.rI1-1)*mp.exp(-y)
        rI, rK = 1+memory, K+4*memory
        k, kp = self.k(y), self.kp(y)
        E02 = mp.exp(2*s.log_p)/(1+z)**2
        E2 = E02*mp.exp(-mp.mpf('.4')-y)
        cp = E02*(self.pressure1+mp.exp(-mp.mpf('.4'))*(-mp.expm1(-y))/2)
        datum = s.pressure_jet(e, 1)
        pi, pie = datum[0]+cp, datum[1]-2*Jp*cp
        W, U = 1-L*K, e*k
        sx = z*K2-E2*(self.energy1+y)/2
        sxe = 2*e*K2+Jp*E2*(self.energy1+y)
        Q = -W+(1-h+D*e*Jp)*rI+(-d*(1-e*Jp)+2*(h-D)*z)*rK
        n = (-W*U+4*h*e*sx-d*sxe+4*A*e*pi-d*pie)/E2
        Q0 = (z*(1+2*d*K)+(1-self.rI1)*(3-6*z+8*z*z)*mp.exp(-y))/(1+z)
        dq = h*(-(1+3*z)/(1+z)+2*z*K+memory*(-(1+3*z)/(1+z)+16*z))
        bsw = 2*e*kp*n/Q
        return dict(y=y, eta=e, Q=Q, n=n, bsw=bsw, Pc_over_ps1=1-bsw/2,
                    Q_from_positive_identity=Q0+dq, Q_at_h_zero=Q0,
                    bs_squared=4*z*kp*kp/E2)

    def intermediate_start(self, eta):
        """At the start of the constant -lambda interval, with U=0, bs=0."""
        s, e = self.s, mp.mpf(eta)
        z, d, h = e*e, 1-e*e, s.h
        A, D, L = mp.mpf('.5')+h, mp.mpf('.5')-h, 1-2*h*z
        Jp = 2*e/(1+z)
        K, K2 = self.averages(s.td)
        memory = (self.rI1-1)*mp.exp(-s.td)
        r = s.rule
        g = lambda t: t-s.lam*r.step_integral(t)
        rI = mp.exp(-g(1))*(1+memory+r.integrate(lambda t: mp.exp(g(t))))
        rK = (K+4*memory)*mp.exp(-g(1))
        W = 1-L*K/mp.e
        Q = -W+(1-h+D*e*Jp)*rI+(-d*(1-e*Jp)+2*(h-D)*z)*rK
        energy = (self.energy1+s.td+r.integrate(lambda t: mp.exp(-2*s.lam*r.step_integral(t))))*mp.exp(s.lam)
        stage = s.stages[3]
        relative_log_E2 = 2*(stage.log_amplitude-mp.log1p(z))
        # Backward pressure avoids cancellation of two O(P*^2) quantities.
        parts = s.pressure_parts(e, 1)[4:]
        pi, pie = [mp.fsum(mp.exp(p['log_scale']-relative_log_E2)*p['jet'][j]
                          for p in parts) for j in range(2)]
        E2 = mp.exp(2*s.log_p+relative_log_E2)
        n = e*(4*h*z-2*d)*K2/(mp.e*E2)-(d*Jp+2*h*e)*energy+4*A*e*pi-d*pie
        return dict(eta=e, Q=Q, n=n, E_squared=E2, a=2+2*s.lam,
                    lambda_w_squared=s.lam*E2*(n/Q)**2,
                    necessary_condition='lambda*w^2 < 1',
                    status='rejected-for-every-positive-XR')


def independent_ode(moments, eta, ys, max_step=.01):
    """Integrate (4.9), without using the positive-Q or cumulative-N formula.

    The pressure datum and initial ideal moments are shared with the moment
    method. This is an equation-assembly crosscheck, not an independent paper.
    Only y<=10 is supported: double precision loses the pressure cancellation
    on long axial intervals.
    """
    if not ys or min(ys)<0 or max(ys)>10:
        raise ValueError('The double-precision ODE check requires 0<=y<=10')
    s, e = moments.s, float(eta)
    h, md = float(s.h), float(s.md)
    A, D, d, L = .5+h, .5-h, 1-e*e, 1-2*h*e*e
    Jp, E02 = 2*e/(1+e*e), float(mp.exp(2*s.log_p)/(1+mp.mpf(eta)**2)**2)
    datum = s.pressure_jet(mp.mpf(eta), 1)
    p0, pe0 = float(datum[0]/E02)+2.5, float(datum[1]/E02)-5*Jp
    W0, U0 = 1-4*L, 4*e
    sx0, sxe0 = 16*e*e-E02*5/12, 32*e+Jp*E02*5/6
    n0 = (-W0*U0+4*h*e*sx0-d*sxe0)/E02+4*A*e*p0-d*pe0
    Q0 = (-W0*.6-h*(1-8*e*e)+(D*e+d*U0)*Jp)/1.6

    def sigma(x):
        if x<=0: return 0.
        if x>=1: return 1.
        v = 1/x**2-1/(1-x)**2
        return math.exp(-v)/(1+math.exp(-v)) if v>=0 else 1/(1+math.exp(v))

    def rhs(axial):
        def evaluate(y, state):
            logE, p, pe, K, Q, n = state
            if axial:
                t = math.log1p(y)/md
                k = 4*sigma(1-t)
                sp = sigma(t)*sigma(1-t)*(2/t**3+2/(1-t)**3) if 0<t<1 else 0.
                kp, l = -4*sp/(md*(1+y)), 0.
            else:
                k, kp, l = 4., 0., .6*(1-sigma(y))
            U, W = e*k, 1-L*K
            Hc = D*e+d*U
            Sq = -W*l-h*(1-2*e*U)+Hc*Jp
            scale = math.exp(2*logE)
            Sn_over_E2 = (-W*e*kp-A*(1-2*e*U)*U-Hc*k)/(E02*scale)+(-d*pe+4*A*e*p)/scale+e
            return [l-.5, .5*scale, -Jp*scale, k-K,
                    Sq-(1+l)*Q, Sn_over_E2-2*l*n]
        return evaluate

    initial = [0., p0, pe0, 4., Q0, n0]
    ramp = solve_ivp(rhs(False), (0,1), initial, method='DOP853',
                     rtol=2e-13, atol=2e-15, max_step=max_step)
    if not ramp.success: raise RuntimeError(ramp.message)
    axial = solve_ivp(rhs(True), (0,max(ys)), ramp.y[:,-1], method='DOP853',
                      rtol=2e-13, atol=2e-15, max_step=max_step, dense_output=True)
    if not axial.success: raise RuntimeError(axial.message)
    return [dict(y=y, eta=e, Q=float(axial.sol(y)[4]), n=float(axial.sol(y)[5])) for y in ys]
