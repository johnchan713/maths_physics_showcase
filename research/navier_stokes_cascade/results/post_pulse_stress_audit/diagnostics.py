"""Actual angular patches and independently integrated terminal controls.

The moderate-lambda Md4 patch fixtures retain their historical axial failure.
They test equations and partial integrals; the continuum claim is bounds.py.
"""
from pathlib import Path
import sys
import mpmath as mp
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.append(str(HERE.parent/'outer_pressure_pilot'))
from schedule import Schedule, Rule, step, step_prime  # noqa: E402


def relative(a,b):
    return abs(a-b)/max(abs(a),abs(b)) if a or b else mp.mpf(0)


class PatchSchedule(Schedule):
    """Reuse the exact A.11 solver; cache its angle-independent input."""
    def _pressure_components(self):
        return None

    def pressure_parts(self,*args,**kwargs):
        raise RuntimeError('This diagnostic exposes only explicitly normalized patch pressure')

    def discrepancy_jet(self,eta,rule=None):
        beta = 1-self.lam
        if not hasattr(self,'_initial_delta'):
            # At eta=1, a=0; undo only the known interpolation factor.
            self._initial_delta = Schedule.interpolation_moment_discrepancy(self,mp.mpf(1))*mp.exp(beta*self.tf)
        a = mp.log(2/(1+eta*eta))
        ae = -2*eta/(1+eta*eta)
        r = rule or self.moment_rule
        weight = lambda z:mp.exp(-beta*self.tf*(1-z)+a*step(1-z))*step_prime(z)
        memory = mp.exp(-beta*self.tf+a)*self._initial_delta
        value = memory+a/beta*r.integrate(weight)
        derivative = ae*(memory+r.integrate(lambda z:weight(z)*(1+a*step(1-z)))/beta)
        return value,derivative

    def interpolation_moment_discrepancy(self,eta,rule=None):
        return self.discrepancy_jet(eta,rule)[0]


def fixture(order=32,lam='0.0001',panels=4,moment_panels=32):
    return PatchSchedule(dict(Md=4,**{'lambda':lam},Tf=1000,co='.001',
        pressure_log_offset=1,h_exponent_factor=2),order=order,panels=panels,moment_panels=moment_panels)


class Patch:
    """Partial backward pressure, including the actual quadratic corrections."""
    def __init__(self,schedule,eta):
        self.s,self.eta = schedule,mp.mpf(eta)
        result = schedule.angular_correction(self.eta)
        self.result,self.c = result,result['coefficients']
        self.scale = max(abs(c) for c in self.c)
        B = mp.matrix(result['matrix'])
        jac = B.copy()
        for j in range(2): jac[1,j] += 2*result['quadratic'][j]*self.c[j]
        factor = mp.exp(-(1-schedule.lam)*(schedule.tu-3))
        target_eta = -schedule.discrepancy_jet(self.eta)[1]*factor
        self.ce = list(mp.lu_solve(jac,mp.matrix([target_eta,0])))
        self.eta_scale = max([abs(target_eta)]+[abs(c) for c in self.ce])
        residual = jac*mp.matrix(self.ce)-mp.matrix([target_eta,0])
        self.derivative_error = max(abs(x) for x in residual)/(self.eta_scale or 1)

    def remaining(self,t,derivative=False):
        rate,width = 1+2*self.s.lam,mp.mpf('.3')
        terms = []
        for center,c,ce in zip((0,2),self.c,self.ce):
            low = max(mp.mpf(0),(mp.mpf(t)-center)/width+mp.mpf('.5'))
            if low>=1: continue
            def integrand(z):
                bump = step_prime(z)
                density = 2*ce*bump+2*c*ce*bump*bump if derivative else 2*c*bump+c*c*bump*bump
                return mp.exp(-rate*(center+width*(z-mp.mpf('.5'))))*density
            terms.append(-width*self.s.rule.integrate(integrand,low,1)/2)
        return mp.fsum(terms)

    def density(self,t):
        # Scale before conversion. A positive bound below accounts for any
        # quadratic term that binary64 rounds away in this diagnostic path.
        value = mp.mpf(0)
        for center,c in zip((0,2),self.c):
            bump = step_prime((mp.mpf(t)-center)/mp.mpf('.3')+mp.mpf('.5'))
            value += 2*(c/self.scale)*bump+(c*c/self.scale)*bump*bump
        return float(mp.exp(-(1+2*self.s.lam)*mp.mpf(t))*value/2)

    def independent(self,cadence):
        points = [-.2,-.15,-.075,0,.075,.15,1,1.85,1.925,2,2.075,2.15,2.2]
        state,t0,rows = 0.0,points[0],[]
        for t in points:
            if t>t0:
                solution = solve_ivp(lambda y,v:[self.density(y)],(t0,t),[state],
                    method='DOP853',rtol=2e-12,atol=2e-13,max_step=cadence)
                if not solution.success: raise RuntimeError(solution.message)
                state = solution.y[0,-1]
            actual = self.remaining(str(t))/self.scale
            rows.append(dict(t=t,backward_pressure_over_cscale=actual,
                forward_ODE_over_cscale=state,normalized_gap=abs(mp.mpf(state)-actual)/(1+abs(actual))))
            t0 = t
        return rows

    def record(self,cadences):
        rows = [dict(max_step=v,rows=self.independent(v)) for v in cadences]
        quadratic = sum(abs(c*c)/self.scale for c in self.c)*20*mp.exp(mp.mpf('.16'))
        # Integral beta^2 <= .3*64 <20. This bounds the full quadratic
        # pressure part, hence also the part omitted by binary64 arithmetic.
        return dict(eta=self.eta,coefficients=self.c,coefficient_eta=self.ce,
            coefficient_scale=self.scale,derivative_row_relative_error=self.derivative_error,
            moment_relative_errors=self.result['relative_errors'],
            complete_pressure_relative_residual=self.remaining('-.2')/self.scale,
            mid_patch_pressure_over_cscale=self.remaining(0)/self.scale,
            gap_pressure_over_cscale=self.remaining(1)/self.scale,
            post_patch_pressure=self.remaining('2.2'),
            mid_patch_pressure_eta=self.remaining(0,True),
            quadratic_normalized_positive_upper=quadratic,independent_runs=rows)


def terminal_q_over_h(y,h,rule,co=None):
    """Positive A.16, normalized before integration; zero endpoint is exact."""
    y,h = mp.mpf(y),mp.mpf(h)
    co = mp.mpf('.001') if co is None else mp.mpf(co)
    if y<0 or y>3: raise ValueError('Terminal y must be in [0,3]')
    if y==3: return mp.mpf(0)
    lower = max(mp.mpf(0),(y-1)/2)
    fo = 1-co*h*step((3-y)/2)
    return co/fo*rule.integrate(lambda z:mp.exp((1-h)*(1+2*z-y))*step_prime(z),lower,1)


def terminal_control(digits=320,order=48,panels=8):
    with mp.workdps(digits):
        T = mp.exp(4)+10
        h = mp.exp(-8*T)
        co,r = mp.mpf('.001'),Rule(order,panels)
        qh = terminal_q_over_h(0,h,r,co)
        qzero = terminal_q_over_h(0,0,r,co)
        direct = (qh-qzero)/h
        stable = co/(1-co*h)*r.integrate(lambda z:
            mp.exp(1+2*z)*(mp.expm1(-h*(1+2*z))/h+co)*step_prime(z))
        points = ('0','0.5','1','1.5','2','2.5','2.75','3')
        rows = [dict(y=y,Q_over_h=terminal_q_over_h(y,h,r,co)) for y in points]
        # Pure power tail: h*S/E^2=1/4, not zero even as h becomes small.
        s = 1/(4*h)
        pi = -1/(2*(1+2*h))
        energy_term = 4*h*s
        pressure_term = 4*(mp.mpf('.5')+h)*pi
        return dict(digits=digits,order=order,panels=panels,h=h,Qp_over_h=qh,
            terminal_rows=rows,direct_finite_h_difference_over_h=direct,
            stable_finite_h_difference_over_h=stable,
            direct_relative_gap=relative(direct,stable),
            binary64_Qp_over_h_difference=float(qh)-float(qzero),
            pure_power_energy_term=energy_term,pure_power_pressure_term=pressure_term,
            pure_power_N_residual=energy_term+pressure_term,
            omitted_h_energy_N_over_E2=pressure_term,
            tail_Q=0,tail_Pc=0,tail_cone_pass=False)


def independent_terminal(h,rule,cadence):
    """Solve (4.9) backwards from Q(3)=0, in the variable Q/h."""
    hf,co = float(h),.001
    def rhs(y,q):
        psi = float(step((mp.mpf(3)-y)/2))
        prime = float(step_prime((mp.mpf(y)-1)/2))/2
        fo = 1-co*hf*psi
        ell = -hf+co*hf*prime/fo
        return [-(1+ell)*q[0]-co*prime/fo]
    points = [3,2.75,2.5,2,1.5,1,.5,0]
    solution = solve_ivp(rhs,(3,0),[0],t_eval=points,method='DOP853',
        rtol=2e-12,atol=2e-14,max_step=cadence)
    if not solution.success: raise RuntimeError(solution.message)
    return [dict(y=y,Q_over_h_ODE=float(q),Q_over_h_formula=terminal_q_over_h(y,h,rule),
                 normalized_gap=abs(mp.mpf(q)-terminal_q_over_h(y,h,rule))/(1+abs(q)))
            for y,q in zip(points,solution.y[0])]


def exterior_control(s):
    """Separate source ODEs for both unit ramps; exact hold and log wait."""
    lam,h = s.lam,s.h
    def ramp(q,left,right):
        def rhs(t,v):
            ell = float(left+(right-left)*step(t))
            return [-ell-float(h)-(1+ell)*v[0]]
        result = solve_ivp(rhs,(0,1),[float(q)],method='DOP853',rtol=2e-12,atol=2e-14,max_step=.02)
        if not result.success: raise RuntimeError(result.message)
        return result.y[0,-1]
    q1 = s.propagate_q(s.q_initial,-lam,-1)
    q1_ode = ramp(s.q_initial,-lam,-1)
    q2 = q1+(1-h)*s.hold
    qr_ode = ramp(q2,-1,-h)
    wait_identity = mp.log(s.q_before_wait)-(1-h)*s.wait-mp.log(s.qp)
    return dict(q_initial=s.q_initial,q_after_first_ramp=q1,q_after_first_ramp_ODE=q1_ode,
        q_before_wait=s.q_before_wait,q_before_wait_ODE=qr_ode,
        maximum_ODE_relative_gap=max(relative(mp.mpf(q1_ode),q1),relative(mp.mpf(qr_ode),s.q_before_wait)),
        wait=s.wait,wait_log_residual=wait_identity,Qp=s.qp,
        E_hold_factor=mp.exp(-mp.mpf('1.5')*s.hold),expected_E_hold_factor=h**6,
        XE2_hold_factor=mp.exp(-2*s.hold),expected_XE2_hold_factor=h**8)
