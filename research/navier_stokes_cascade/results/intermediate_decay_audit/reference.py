"""Computable Md=4 prefix: moments versus the source equations.

The future pressure is replaced by its continuing power law with a proved
absolute normalized error <=lambda^8. This is a diagnostic prefix, not the
complete matched outer profile. The all-angle certificate uses the exact
pressure envelope instead of this approximation.
"""
import math
from pathlib import Path
import sys
import mpmath as mp
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
OUTER = HERE.parent/'outer_pressure_pilot'
sys.path.append(str(OUTER))
from schedule import Rule, step  # noqa: E402


class Reference:
    def __init__(self, order=32, md=4):
        if md != 4:
            raise ValueError('The numerical prefix is restricted to Md=4; use bounds for Md=64')
        if mp.mp.dps < 260:
            raise ValueError('Use at least 260 digits to resolve the finite-h pole increment')
        self.md = mp.mpf(md)
        self.T = mp.exp(self.md)+10
        self.delta = mp.exp(-self.T)
        self.lam = mp.exp(-4*self.T)
        self.h = self.lam**2
        self.beta = 1-self.lam
        self.Tw = 240*self.T
        self.rule = Rule(order,4)
        r = self.rule
        shape = lambda t: mp.mpf('.1')*t-mp.mpf('.6')*r.step_integral(t)
        self.r1 = (mp.mpf(5)/8+r.integrate(lambda t: mp.exp(mp.mpf('1.5')*t+shape(t))))/mp.exp(mp.mpf('1.3'))
        self.F1 = (mp.mpf(5)/6+r.integrate(lambda t: mp.exp(t+2*shape(t))))/mp.exp(mp.mpf('.6'))
        k = lambda t: 4*step(1-mp.log1p(t)/self.md)
        end = mp.expm1(self.md)
        self.Ka, self.K2a = [initial*self.delta+r.integrate(
            lambda t: mp.exp(t-self.T)*k(t)**power,0,end)
            for initial,power in ((4,1),(16,2))]
        self.memory = (self.r1-1)*self.delta
        self.cache = {}

    def ramp_factors(self,x):
        """Integrated swirl, angular moments, energy moment and pressure deficit."""
        x = mp.mpf(x)
        if not 0 <= x <= 1:
            raise ValueError('Ramp coordinate must be in [0,1]')
        key = ('ramp',x)
        if key not in self.cache:
            r, lam = self.rule, self.lam
            S = lambda t: lam*r.step_integral(t)
            G = lambda t: t-S(t)
            rI = mp.exp(-G(x))*(1+self.memory+r.integrate(lambda t: mp.exp(G(t)),0,x))
            rK = (self.Ka+4*self.memory)*mp.exp(-G(x))
            F = mp.exp(2*S(x))*(self.F1+self.T+r.integrate(lambda t: mp.exp(-2*S(t)),0,x))
            # Z=(1-B)/(2lambda) avoids subtracting an O(lambda) deficit from 1.
            R = lambda t: t+2*S(t)
            Z = mp.exp(R(x)-R(1))/(1+2*lam)+r.integrate(
                lambda t: mp.exp(R(x)-R(t))*step(t),x,1)
            self.cache[key] = dict(x=x,S=S(x),rI=rI,rK=rK,F=F,Z=Z,
                                   K=self.Ka*mp.exp(-x),K2=self.K2a*mp.exp(-x))
        return self.cache[key]

    def power_factors(self,y):
        y = mp.mpf(y)
        if not 0 <= y <= self.Tw:
            raise ValueError('Power coordinate outside [0,Tw]')
        key = ('power',y)
        if key not in self.cache:
            b = self.ramp_factors(1)
            decay = mp.exp(-self.beta*y)
            self.cache[key] = dict(x=1+y,S=self.lam*(y+mp.mpf('.5')),
                rI=b['rI']*decay-mp.expm1(-self.beta*y)/self.beta,
                rK=b['rK']*decay,
                F=mp.exp(2*self.lam*y)*(b['F']-mp.expm1(-2*self.lam*y)/(2*self.lam)),
                Z=1/(1+2*self.lam),K=b['K']*mp.exp(-y),K2=b['K2']*mp.exp(-y))
        return self.cache[key]

    def state(self,stage,coordinate,eta,h_override=None):
        """Evaluate (4.16); v=N/(eta E^2) has its continuous value at eta=0."""
        if stage not in ('ramp','power'):
            raise ValueError('Unknown stage')
        f = self.ramp_factors(coordinate) if stage=='ramp' else self.power_factors(coordinate)
        e = mp.mpf(eta)
        if abs(e)>1: raise ValueError('Eta outside [-1,1]')
        z, d = e*e, 1-e*e
        h = self.h if h_override is None else mp.mpf(h_override)
        L, D = 1-2*h*z, mp.mpf('.5')-h
        ejp, c = 2*z/(1+z), 2*d/(1+z)
        W = 1-L*f['K']
        Q = -W+(1-h+D*ejp)*f['rI']+(-d*(1-ejp)+2*(h-D)*z)*f['rK']
        # P*^2 contributes 2(T+1), first ramp -.4, axial stage -T.
        E2 = mp.exp(self.T+mp.mpf('1.6')-f['x']-2*f['S'])/(1+z)**2
        B = 1-2*self.lam*f['Z']
        v = (4*h*z-2*d)*f['K2']/E2-(c+2*h)*f['F']-(1+2*h+c)*B
        danger = self.lam*E2*z*(v/Q)**2
        # Pressure-only uncertainty in v, propagated over both U=0 stages.
        v_tail_error = 4*self.lam**8*(1+f['x'])*mp.exp(2*f['S'])
        return dict(stage=stage,coordinate=mp.mpf(coordinate),eta=e,Q=Q,v=v,
                    E_squared=E2,lambda_w_squared=danger,v_tail_error_upper=v_tail_error)

    def source_power(self,y,eta):
        """Solve the scalar source equations on the power stage using expm1."""
        y, e = mp.mpf(y), mp.mpf(eta)
        b = self.state('ramp',1,e)
        f = self.ramp_factors(1)
        z, c = e*e, 2*(1-e*e)/(1+e*e)
        ceta = (1-2*self.h)*z/(1+z)
        decay = mp.exp(-self.beta*y)
        Q = b['Q']*decay+(self.lam-self.h+ceta)*(-mp.expm1(-self.beta*y))/self.beta
        Q -= (1-2*self.h*z)*f['K']*decay*(-mp.expm1(-self.lam*y))
        g = (2*self.lam-2*self.h-c)/(1+2*self.lam)
        increment = mp.expm1(2*self.lam*y)*(b['v']+g/(2*self.lam))
        return dict(Q=Q,v=b['v']+increment,v_increment=increment)

    def independent_ramp(self,eta,xs,max_step=.02):
        """Double-precision source IVP, with pressure deficit and Q scaled first."""
        e = mp.mpf(eta)
        lam,h = float(self.lam),float(self.h)
        z,c = float(e*e),float(2*(1-e*e)/(1+e*e))
        scale = e*e+self.delta+self.lam
        ceta = (1-2*self.h)*e*e/(1+e*e)
        lam_scaled,h_scaled,c_scaled = [float(v/scale) for v in (self.lam,self.h,ceta)]
        K0 = float(self.Ka)
        def sigma(x):
            if x<=0: return 0.
            if x>=1: return 1.
            v=1/x**2-1/(1-x)**2
            return math.exp(-v)/(1+math.exp(-v)) if v>=0 else 1/(1+math.exp(v))
        # Backward pressure-deficit equation, independently coded from its integral.
        pressure = solve_ivp(lambda x,Z:[(1+2*lam*sigma(x))*Z[0]-sigma(x)],
            (1,0),[float(1/(1+2*self.lam))],method='DOP853',rtol=2e-13,atol=2e-15,
            max_step=max_step,dense_output=True)
        initial=self.state('ramp',0,e)
        def rhs(x,state):
            alpha=lam*sigma(x)
            W=1-(1-2*h*z)*K0*math.exp(-x)
            g=-2*h-c+2*lam*float(pressure.sol(x)[0])*(1+2*h+c)
            return [-(1-alpha)*state[0]+lam_scaled*sigma(x)*W-h_scaled+c_scaled,
                    g+2*alpha*state[1]]
        flow=solve_ivp(rhs,(0,1),[float(initial['Q']/scale),float(initial['v'])],
            method='DOP853',rtol=2e-13,atol=2e-15,max_step=max_step,dense_output=True)
        if not pressure.success or not flow.success:
            raise RuntimeError('Independent source integration failed')
        return [dict(x=x,Q_over_scale=float(flow.sol(x)[0]),v=float(flow.sol(x)[1]),
                     Z=float(pressure.sol(x)[0]),scale=scale) for x in xs]
