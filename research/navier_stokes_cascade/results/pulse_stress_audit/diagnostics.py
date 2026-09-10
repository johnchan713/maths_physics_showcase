"""Independent finite-parameter identities and deliberately failed controls.

The ODE fixtures are manufactured profiles, not the selected reference
family. They make finite-h, angular, pressure and end-patch terms visible.
The cancellation fixture is an exactly solvable affine particular solution,
not an assertion that the complete pulse equals that particular solution.
"""
import math
from fractions import Fraction
import numpy as np
import mpmath as mp
from scipy.integrate import quad, solve_ivp
from scipy.special import expit


def step_jet(t):
    """The fixed flat step and two derivatives, for double-precision controls."""
    if t <= 0: return 0.0, 0.0, 0.0
    if t >= 1: return 1.0, 0.0, 0.0
    s = float(expit(1/(1-t)**2-1/t**2))
    g = 2/t**3+2/(1-t)**3
    gp = -6/t**4+6/(1-t)**4
    return s, s*(1-s)*g, s*(1-s)*((1-2*s)*g*g+gp)


def shape_jet(xi):
    """R0 and R0_xi; the primitive is exact outside the short initial ramp."""
    if xi <= 0 or xi >= 11: return 0.0, 0.0
    if xi < .02:
        phi = .02*quad(lambda t:step_jet(t)[0],0,xi/.02,epsabs=1e-14,epsrel=1e-13)[0]
    else: phi = xi-.01
    start = step_jet(xi/.02)[0]
    cut, cutp, _ = step_jet(xi-10)
    return phi*(1-cut), start*(1-cut)-phi*cutp


def moment_rhs(state, R, Re, lam):
    """Normalized derivatives of M,I,J,S and their eta derivatives."""
    m, me, r, re, j, je, s, se = state
    beta, alpha, gamma = .5-lam, 1-lam, .5-2*lam
    return [R-beta*m, Re-beta*me, 1-alpha*r, -alpha*re,
            R-gamma*j, Re-gamma*je, R*R-.5+2*lam*s, 2*R*Re+2*lam*se]


def moment_stress(state, E, R, eta, lam, h):
    """Full (4.16), with an infinite constant-slope pressure tail."""
    m, me, r, re, j, je, s, se = state
    D, Ap = .5-h, .5+h
    d, Jp = 1-eta*eta, 2*eta/(1+eta*eta)
    c = D*eta*Jp
    B = 2*D*eta*m+d*(me-Jp*m)
    W = 1-E*B
    pressure, pressure_eta = -1/(2*(1+2*lam)), Jp/(1+2*lam)
    q = -W+(1-h+c)*r-D*eta*re+E*(-d*je+(2*d*Jp+2*(h-D)*eta)*j)
    energy = 4*h*eta*s-d*(se-2*Jp*s)
    pressure_part = 4*Ap*eta*pressure-d*pressure_eta
    n = -W*R+(D+c)*m-D*eta*me+E*(energy+pressure_part)
    return q, n


def source_stress(state, E, R, Ry, Re, eta, lam, h):
    """Direct (4.9), without the integrated Q/N formulas."""
    m, me = state[:2]
    D, Ap = .5-h, .5+h
    d, Jp = 1-eta*eta, 2*eta/(1+eta*eta)
    W = 1-E*(2*D*eta*m+d*(me-Jp*m))
    Hc = D*eta+d*E*R
    pressure, pressure_eta = -1/(2*(1+2*lam)), Jp/(1+2*lam)
    Sq = lam*W-h*(1-2*eta*E*R)+Hc*Jp
    Sn_over_E = -W*(Ry-(.5+lam)*R)-Ap*(1-2*eta*E*R)*R
    Sn_over_E -= Hc*(Re-Jp*R)
    Sn_over_E += E*(-d*pressure_eta+4*Ap*eta*pressure+eta)
    return Sq, Sn_over_E


class SyntheticPulse:
    """Visible, signed end bumps and a non-even strength for identity checks."""
    def __init__(self, eta, lam=.02):
        self.eta, self.lam = float(eta), float(lam)
        if not -1 <= self.eta <= 1 or not 0 < self.lam < .1:
            raise ValueError('Expected eta in [-1,1] and lambda in (0,.1)')
        self.h = self.lam**2
        self.end = 13/self.lam
        self.centers = (self.end-3,self.end-1)

    def profile(self, y, omit_bumps=False):
        e, lam = self.eta, self.lam
        base, deriv = shape_jet(lam*y)
        amp, ampe = 1+.02*e+.01*e*e, .02+.02*e
        R, Ry, Re = amp*base, lam*amp*deriv, ampe*base
        if not omit_bumps:
            coeffs = ((1+.2*e)*1e-5,-(1-.1*e)*2e-5)
            for center, coeff, ce in zip(self.centers,coeffs,(2e-6,2e-6)):
                _, bump, bump_y = step_jet((y-center)/.3+.5)
                R += coeff*bump
                Ry += coeff*bump_y/.3
                Re += ce*bump
        E = .03/(1+e*e)*math.exp(-(.5+lam)*y)
        return E, R, Ry, Re

    def initial(self):
        e = self.eta
        state = [.1+.02*e*e,.04*e,.8+.01*e*e,.02*e,
                 .2-.03*e,-.03,-.2+.02*e*e,.04*e]
        E,R,_,_ = self.profile(0)
        return np.array(state+list(moment_stress(state,E,R,e,self.lam,self.h)))

    def rhs(self, y, state):
        E,R,Ry,Re = self.profile(y)
        moments = moment_rhs(state[:8],R,Re,self.lam)
        Sq,Sn = source_stress(state,E,R,Ry,Re,self.eta,self.lam,self.h)
        return moments+[Sq-(1-self.lam)*state[8], Sn-(.5-self.lam)*state[9]]

    def samples(self):
        values = [0,.25,.5,.75,1,.51/self.lam,1/self.lam,5/self.lam,
                  10/self.lam,10.3/self.lam,10.5/self.lam,10.7/self.lam,
                  11/self.lam,12/self.lam,self.end]
        values += [c+d for c in self.centers for d in (-.15,-.075,0,.075,.15)]
        return sorted(set(values))

    def integrate(self, bulk_step=2, patch_step=.04):
        y0, state, rows = 0.0, self.initial(), []
        for y1 in self.samples():
            if y1 > y0:
                near_patch = any(y0 < c+.3 and y1 > c-.3 for c in self.centers)
                result = solve_ivp(self.rhs,(y0,y1),state,method='DOP853',
                    rtol=2e-12,atol=2e-13,max_step=patch_step if near_patch else bulk_step)
                if not result.success: raise RuntimeError(result.message)
                state = result.y[:,-1]
            E,R,Ry,Re = self.profile(y1)
            qm,nm = moment_stress(state[:8],E,R,self.eta,self.lam,self.h)
            gap = max(abs(qm-state[8])/(1+abs(qm)),abs(nm-state[9])/(1+abs(nm)))
            full_bs = 2*(Ry-(.5+self.lam)*R)
            _,r_bad,ry_bad,_ = self.profile(y1,omit_bumps=True)
            drop_bs = 2*(ry_bad-(.5+self.lam)*r_bad)
            rows.append(dict(y=y1,eta=self.eta,Q_source=state[8],Q_moment=qm,
                N_over_E_source=state[9],N_over_E_moment=nm,
                normalized_gap=gap,bs=full_bs,omitted_bump_bs_gap=abs(full_bs-drop_bs)))
            y0 = y1
        return rows


def algebra_controls():
    """Differentiate (4.16) automatically, compare with independent (4.9)."""
    with mp.workdps(110):
        lam = mp.mpf('.02')
        h = lam**2
        state = list(map(mp.mpf,('.7','.2','.9','-.15','.3','.11','-.4','.17')))
        R,Ry,Re = map(mp.mpf,('.8','-.3','.13'))
        errors,omitted_h,omitted_eta,omitted_energy, rows = [],[],[],[],[]
        for eta in map(mp.mpf,('-1','-.3','0','.3','1')):
            E = mp.mpf('.03')/(1+eta**2)
            velocities = moment_rhs(state,R,Re,lam)
            moment_at = lambda t:moment_stress([v+t*w for v,w in zip(state,velocities)],
                E*mp.exp(-(.5+lam)*t),R+t*Ry,eta,lam,h)
            q,n = moment_at(0)
            Sq,Sn = source_stress(state,E,R,Ry,Re,eta,lam,h)
            residual_q = mp.diff(lambda t:moment_at(t)[0],0)+(1-lam)*q-Sq
            residual_n = mp.diff(lambda t:moment_at(t)[1],0)+(.5-lam)*n-Sn
            errors.extend((abs(residual_q),abs(residual_n)))
            q_bad,n_bad = moment_stress(state,E,R,eta,lam,mp.mpf(0))
            omitted_h.append(max(abs(q-q_bad),abs(n-n_bad)))
            omitted_eta.append(abs((.5-h)*eta*state[1]))
            d,Jp = 1-eta**2,2*eta/(1+eta**2)
            omitted_energy.append(abs(E*(4*h*eta*state[6]-d*(state[7]-2*Jp*state[6]))))
            rows.append(dict(eta=eta,Q_source_residual=residual_q,N_source_residual=residual_n))
        return dict(maximum_source_identity_error=max(errors),rows=rows,
            omitted_h_maximum_gap=max(omitted_h),omitted_eta_maximum_gap=max(omitted_eta),
            omitted_energy_maximum_gap=max(omitted_energy))


def affine_control(digits, angle_ratio='0'):
    """Exact affine particular solution exposes small-angle cancellation.

    eta=angle_ratio*sqrt(lambda); this is not the full closed reference pulse.
    At xi=.51 the startup transient is absent by definition of the particular
    solution. Positive finite remainders of the actual pulse are in bounds.py.
    """
    with mp.workdps(digits):
        T = mp.exp(4)+10
        lam = mp.exp(-4*T)
        h = lam**2
        eta = mp.mpf(angle_ratio)*mp.sqrt(lam)
        D,beta = mp.mpf('.5')-h,mp.mpf('.5')-lam
        c = (1-2*h)*eta**2/(1+eta**2)
        amp = mp.mpf('1.0100503')
        R = amp*mp.mpf('.5')
        q = (lam-h+c)/(1-lam)
        Cd = lam*(D+c)*(1-lam)/(beta**2*(lam-h+c))
        m = R/beta-lam*amp/beta**2
        direct = (-R+(D+c)*m)/q
        stable = (1-lam)/beta*R-Cd*amp
        missing_derivative = (1-lam)/beta*R
        return dict(digits=digits,angle_over_sqrt_lambda=angle_ratio,lambda_value=lam,
            Q=q,Cd=Cd,direct_w=direct,stable_w=stable,
            direct_relative_error=abs(direct-stable)/max(abs(stable),mp.mpf('1e-80')),
            omitted_derivative_absolute_error=abs(missing_derivative-stable),
            a_minus_two=2*lam,rounded_a_minus_two=2+2*lam-2,
            binary64_rounded_a_minus_two=2+2*float(lam)-2)


def finite_cone_controls():
    """Check the exact quadratic identity and a radius-too-small counterexample."""
    with mp.workdps(100):
        residuals = []
        for a,b,w in ((mp.mpf('2.02'),mp.mpf('-.5'),mp.mpf('-1.1')),
                      (mp.mpf('2.0001'),mp.mpf('-10'),mp.mpf('100')),
                      (mp.mpf('2.01'),mp.mpf('.001'),mp.mpf('-.02'))):
            c,j,v = 1-b*w/a,w+b/a,a+b*b/a
            second = 2*b*w+b*b/a+(a-2)*w*w
            residuals.append(abs(2*c*c-(v-2)*j*j-(1+b*b/(a*a))*(2-second)))
        # Positive A.24 margins alone are not the finite-radius cone.
        a,b,w,p = mp.mpf('2.01'),mp.mpf('-.5'),mp.mpf(-1),mp.mpf(1)
        v = a+b*b/a
        Pc = p*(1-b*w/a)
        return dict(identity_maximum_error=max(residuals),small_radius=dict(
            first_margin=a-b*w,second_margin=2-2*b*w-b*b/a-(a-2)*w*w,
            Pc=Pc,v=v,passes_ratio_test=a-b*w>0 and 2-2*b*w-b*b/a-(a-2)*w*w>0,
            passes_finite_cone=Pc>v))


def angular_equality_control():
    """Exact rational scaling avoids a false floating comparison at equality."""
    ratios = [Fraction(s) for s in ('0','0.1','1','3','1000')]
    exact = all(2*t<=1+t*t for t in ratios)
    with mp.workdps(120):
        lam = mp.mpf('1e-80')
        eta = mp.sqrt(lam)
        direct = eta/(lam+eta*eta)
        proposed_upper = 1/(2*mp.sqrt(lam))
        return dict(exact_rational_checks=exact,
            scaled_identity='1+t^2-2t=(t-1)^2>=0',
            equality_scale='1',digits=120,
            floating_comparison_passed=direct<=proposed_upper,
            positive_rounding_excess=direct-proposed_upper,
            relative_rounding_excess=(direct-proposed_upper)/proposed_upper)
