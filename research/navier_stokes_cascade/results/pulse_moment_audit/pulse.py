"""Scaled numerical diagnostics for the exact reference moment equations.

The modest-lambda fixtures evaluate both axial corrections. At the much
smaller selected family, the theorem uses bounds; the numerical diagnostic
retains the amplitude increment and explicitly bounds omitted bump energy.
Neither diagnostic asserts the pulse stress cone or a complete PDE solution.
"""
from pathlib import Path
import math
import sys
import mpmath as mp
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
OUTER = HERE.parent/'outer_pressure_pilot'
sys.path.append(str(OUTER))
from schedule import Rule, Schedule, step, step_prime  # noqa: E402


def relative(a, b):
    return abs(a-b)/max(abs(a),abs(b)) if a or b else mp.mpf(0)


def pulse_constant(rule):
    """K_b using the exact middle integral and positive end integrals."""
    width, offset = mp.mpf('.02'), mp.mpf('.01')
    start = rule.integrate(lambda s: width**3*mp.exp(-2*width*s)*rule.step_integral(s)**2)
    F = lambda x: -mp.exp(-2*x)*((x-offset)**2/2+(x-offset)/2+mp.mpf('.25'))
    middle = F(10)-F(width)
    cutoff = rule.integrate(lambda s:mp.exp(-2*(10+s))*(10+s-offset)**2*step(1-s)**2)
    return start+middle+cutoff


def incoming_integrals(md, rule):
    """Positive incoming U/U^2 integrals; coarse panels miss relative digits."""
    T, end = mp.exp(md)+10, mp.expm1(md)
    k = lambda t:4*step(1-mp.log1p(t)/md)
    return [initial*mp.exp(-T)+rule.integrate(
        lambda t:mp.exp(t-T)*k(t)**power,0,end)
        for initial,power in ((4,1),(16,2))]


class Prefix:
    """Accumulate M, J and S from the actual A.7 prefix to pulse start.

    m0=M/(Xp Ep), j0=J/(Xp Hp Ep), s0=S/(Xp Ep^2).
    Here Ep=eb/(1+eta^2). The earlier Md=4 cone failure is not changed.
    """
    def __init__(self, lam=None, order=32, md=4, incoming_panels=16):
        if md != 4 or mp.mp.dps < 80:
            raise ValueError('Numerical prefix requires Md=4 and at least 80 digits')
        self.md = md
        self.T = mp.exp(md)+10
        self.lam = mp.exp(-4*self.T) if lam is None else mp.mpf(lam)
        if not 0 < self.lam <= mp.mpf('.001'):
            raise ValueError('Diagnostic lambda must lie in (0,.001]')
        self.Tw = 60*mp.log(1/self.lam)
        self.rule = r = Rule(order,4)
        first = lambda t:mp.mpf('.1')*t-mp.mpf('.6')*r.step_integral(t)
        self.r1 = (mp.mpf(5)/8+r.integrate(lambda t:mp.exp(mp.mpf('1.5')*t+first(t))))/mp.exp(mp.mpf('1.3'))
        self.F1 = (mp.mpf(5)/6+r.integrate(lambda t:mp.exp(t+2*first(t))))/mp.exp(mp.mpf('.6'))
        delta = mp.exp(-self.T)
        self.Ka,self.K2a = incoming_integrals(md,Rule(order,incoming_panels))
        self.Ja = self.Ka+4*(self.r1-1)*delta
        growth = self.lam*(self.Tw+mp.mpf('.5'))
        self.log_eb = self.T/2+mp.mpf('.3')-self.Tw/2-growth
        self.mcoef = self.Ka*mp.exp(-1-self.Tw-self.log_eb)
        self.jcoef = self.Ja*mp.exp(-1-self.Tw+growth-self.log_eb)
        self.scoef = self.K2a*mp.exp(-1-self.Tw-2*self.log_eb)
        ramp = r.integrate(lambda t:mp.exp(-2*self.lam*r.step_integral(t)))
        self.F = mp.exp(2*growth)*(self.F1+self.T+ramp)+mp.expm1(2*self.lam*self.Tw)/(2*self.lam)

    def moments(self, eta):
        e = mp.mpf(eta)
        if abs(e)>1:
            raise ValueError('Eta outside [-1,1]')
        inverse_f = 1+e*e
        return dict(m0=e*inverse_f*self.mcoef,j0=e*inverse_f*self.jcoef,
                    s0=e*e*inverse_f**2*self.scoef-self.F/2)


class _EnergySchedule(Schedule):
    """Reuse the prescribed stages without constructing unused pressure jets."""
    def __init__(self, prefix, data, **kwargs):
        self._prefix = prefix
        super().__init__(data, **kwargs)

    def _pressure_components(self):
        return None

    def pressure_parts(self, eta, degree=0):
        raise RuntimeError('This private schedule computes energy, not pressure')

    def interpolation_moment_discrepancy(self, eta, rule=None):
        """The same A.11 input, reusing the already integrated first ramp."""
        beta = 1-self.lam
        if not hasattr(self,'_pulse_delta0'):
            r = 1+(self._prefix.r1-1)*mp.exp(-self.td)
            G = lambda t:t-self.lam*self.rule.step_integral(t)
            r = mp.exp(-1+self.lam/2)*(r+self.rule.integrate(lambda t:mp.exp(G(t))))
            self._pulse_delta0 = (r-1/beta)*mp.exp(-beta*(self.tw+13/self.lam))
        a = mp.log(2/(1+eta*eta))
        driven = a/beta*(rule or self.moment_rule).integrate(lambda z:
            mp.exp(-beta*self.tf*(1-z)+a*step(1-z))*step_prime(z))
        return mp.exp(-beta*self.tf+a)*self._pulse_delta0+driven


class Continuation:
    """Full A.2 swirl after the pulse, including actual angular edits.

    Every stage's positive XE^2 integral is stored separately. The angular
    increments are separate signed terms, so tiny edits are not labelled zero
    just because adding them to the large unedited total loses their digits.
    """
    def __init__(self, prefix, order=32, family=False):
        self.prefix = prefix
        hfactor = 8 if family else 2
        data = dict(Md=4, **{'lambda':mp.nstr(prefix.lam,mp.mp.dps)}, Tf=1000,
                    co='.001',pressure_log_offset=1,h_exponent_factor=hfactor)
        self.schedule = _EnergySchedule(prefix,data,order=order,panels=4,moment_panels=16)
        self.cache = {}
        self._uniform_masses = {}

    def energy_mass(self, eta):
        # E and its angular edits are even, independently of the axial pulse.
        eta = abs(mp.mpf(eta))
        if eta in self.cache:
            return self.cache[eta]
        s, r = self.schedule, self.schedule.rule
        pulse = s.stages[4]
        end_y = pulse.start+pulse.length
        end_loga = pulse.log_amplitude+s.log_shape(pulse,pulse.length,0)
        j = mp.log1p(eta*eta)
        parts = []
        for stage in s.stages[5:]:
            logscale = stage.start-end_y+2*(stage.log_amplitude-end_loga)
            if stage.kind=='interpolation':
                value = mp.exp(logscale+2*j)*r.integrate(
                    lambda y:mp.exp(y+2*s.log_shape(stage,y,eta)),0,stage.length)
            else:
                # Every subsequent unedited stage is independent of eta.
                # Only division by the pulse-end f(eta)^2 changes its mass.
                if stage.name not in self._uniform_masses:
                    if stage.kind=='constant':
                        rate = -2*stage.left_slope
                        mass = (1 if mp.isinf(stage.length) else -mp.expm1(-rate*stage.length))/rate
                    else:
                        mass = r.integrate(lambda y:mp.exp(y+2*s.log_shape(stage,y,0)),0,stage.length)
                    self._uniform_masses[stage.name] = mp.exp(logscale)*mass
                value = self._uniform_masses[stage.name]*(1+eta*eta)**2
            parts.append(dict(stage=stage.name,mass=value))
        edit = s.angular_correction(eta)
        uniform = s.stages[6]
        increments = []
        for coefficient, center in zip(edit['coefficients'],(s.tu-3,s.tu-1)):
            width = mp.mpf('.3')
            base = uniform.start-end_y+2*(uniform.log_amplitude-end_loga)+2*j-2*s.lam*center
            for power, factor in ((1,2*coefficient),(2,coefficient**2)):
                local = width*r.integrate(lambda z:mp.exp(-2*s.lam*width*(z-mp.mpf('.5')))*step_prime(z)**power)
                increments.append(mp.exp(base)*factor*local)
        total = mp.fsum([p['mass'] for p in parts]+increments)
        result = dict(total=total,parts=parts,angular_energy_increments=increments,
                      angular_coefficients=edit['coefficients'],
                      angular_relative_errors=edit['relative_errors'])
        self.cache[eta] = result
        return result


def log_main_moment(lam, slope, rule):
    """Log of integral exp(s*(y-Y)) R0(lambda*y) dy, Y=13/lambda-3.

    The flat cutoff produces a narrow saddle. Locate it and integrate in
    its natural width; a uniform xi grid would miss the main mass. The
    middle piece is analytic. The remaining [0,.02] contribution is enclosed
    positively and returned as a relative upper bound, never called zero.
    """
    lam, s = mp.mpf(lam), mp.mpf(slope)
    if not mp.mpf('1e-5')<=lam<=mp.mpf('.001'):
        raise ValueError('Actual numerical end corrections use lambda in [1e-5,.001]')
    k = s/lam
    seed = (2/k)**(mp.mpf(1)/3)
    derivative = lambda d:-k+(1-step(d))*(2/d**3+2/(1-d)**3)-1/(mp.mpf('10.99')-d)
    a, b = seed/2, 2*seed
    for _ in range(mp.mp.prec+10):
        mid = (a+b)/2
        if derivative(mid)>0: a=mid
        else: b=mid
    center = (a+b)/2
    phase = lambda d:-k*d+mp.log(mp.mpf('10.99')-d)+mp.log(step(d))
    maximum = phase(center)
    width = center**2/mp.sqrt(6)
    knots = [mp.mpf(0)]
    knots += [center+i*width for i in range(-32,33,4) if 0<center+i*width<1]
    knots += [mp.mpf(1)]
    integral = mp.fsum(rule.integrate(lambda d:mp.exp(phase(d)-maximum),a,b)
                      for a,b in zip(knots,knots[1:]))
    # Exact integral on [.02,10], with exp(-k*(11-xi)) factored first.
    middle = mp.exp(-k)*(mp.mpf('9.99')/k-1/k**2)
    middle -= mp.exp(-mp.mpf('10.98')*k)*(mp.mpf('.01')/k-1/k**2)
    start_upper = mp.mpf('.02')**2/2*mp.exp(-mp.mpf('10.98')*k)
    body = integral+mp.exp(-maximum)*middle
    log_value = -2*s/lam+3*s-mp.log(lam)+maximum+mp.log(body)
    return dict(log_value=log_value,cutoff_saddle=center,cutoff_width=width,
                prefix_relative_error_upper=start_upper*mp.exp(-maximum)/body)


def stable_two_row(lam, b, rhs):
    """Solve [b_i, exp(2*s_i)*b_i]c=rhs with the small determinant retained."""
    s1 = mp.mpf('.5')-lam
    denominator = mp.exp(2*s1)*mp.expm1(-2*lam)
    if not 0<lam<mp.mpf('.1') or not all(v>0 for v in b):
        raise ValueError('Expected positive lambda and positive bump masses')
    r1, r2 = rhs[0]/b[0], rhs[1]/b[1]
    c2 = (r2-r1)/denominator
    return [r1-mp.exp(2*s1)*c2,c2]


class MomentSolver:
    """Actual axial end corrections for a computable moderate-lambda fixture."""
    def __init__(self, prefix, order=32):
        self.prefix, self.lam = prefix, prefix.lam
        self.rule = r = Rule(order,4)
        self.slopes = [mp.mpf('.5')-self.lam,mp.mpf('.5')-2*self.lam]
        self.main = [log_main_moment(self.lam,s,r) for s in self.slopes]
        self.log_scale = max(v['log_value'] for v in self.main)
        self.main_scaled = [mp.exp(v['log_value']-self.log_scale) for v in self.main]
        self.b = [mp.mpf('.3')*r.integrate(lambda z:mp.exp(s*mp.mpf('.3')*(z-mp.mpf('.5')))*step_prime(z)) for s in self.slopes]
        self.matrix = [[b,b*mp.exp(2*s)] for b,s in zip(self.b,self.slopes)]
        self.a_scaled = stable_two_row(self.lam,self.b,[-v for v in self.main_scaled])
        self.energy_weights = [mp.exp(-26+6*self.lam-2*self.lam*center)*mp.mpf('.3')*r.integrate(
            lambda z:mp.exp(-2*self.lam*mp.mpf('.3')*(z-mp.mpf('.5')))*step_prime(z)**2) for center in (0,2)]
        self.K = pulse_constant(r)
        self.C = -mp.expm1(-26)/4

    def affine(self, eta):
        m = self.prefix.moments(eta)
        Y = 13/self.lam-3
        pre = [v*mp.exp(-s*Y-self.log_scale) for v,s in zip((m['m0'],m['j0']),self.slopes)]
        b = stable_two_row(self.lam,self.b,[-v for v in pre])
        return dict(a=self.a_scaled,b=b,pre_scaled=pre,moments=m)

    def root(self, eta, continuation):
        affine = self.affine(eta)
        a,b = affine['a'],affine['b']
        tiny = self.lam*mp.exp(2*self.log_scale)
        dk = tiny*mp.fsum(w*x*x for w,x in zip(self.energy_weights,a))
        linear = 2*tiny*mp.fsum(w*x*y for w,x,y in zip(self.energy_weights,a,b))
        dc = tiny*mp.fsum(w*y*y for w,y in zip(self.energy_weights,b))
        post = continuation.energy_mass(eta)
        prefix_part = self.lam*affine['moments']['s0']
        post_part = -self.lam*mp.exp(-26)*post['total']/2
        constant = mp.fsum([-self.C,prefix_part,post_part,dc])
        quadratic = self.K+dk
        amplitude = -2*constant/(linear+mp.sqrt(linear**2-4*quadratic*constant))
        scaled = [x*amplitude+y for x,y in zip(a,b)]
        rows = [mp.fsum([affine['pre_scaled'][i],amplitude*self.main_scaled[i]]+
                       [self.matrix[i][j]*scaled[j] for j in range(2)]) for i in range(2)]
        scales = [max(abs(affine['pre_scaled'][i]),abs(amplitude*self.main_scaled[i]),
                      *(abs(self.matrix[i][j]*scaled[j]) for j in range(2))) for i in range(2)]
        quadratic_residual = mp.fsum([self.K*amplitude**2,-self.C,prefix_part,post_part,
                                     tiny*mp.fsum(w*c*c for w,c in zip(self.energy_weights,scaled))])
        # Setting only c1 from M is a deliberately incomplete correction.
        incomplete = -(affine['pre_scaled'][0]+amplitude*self.main_scaled[0])/self.b[0]
        bad_J = mp.fsum([affine['pre_scaled'][1],amplitude*self.main_scaled[1],self.b[1]*incomplete])
        bad_scale = max(abs(amplitude*self.main_scaled[1]),abs(self.b[1]*incomplete))
        return dict(eta=mp.mpf(eta),amplitude=amplitude,leading_amplitude=mp.sqrt(self.C/self.K),
            coefficient_log_scale=self.log_scale,scaled_coefficients=scaled,
            affine_a_scaled=a,affine_b_scaled=b,linear_relative_errors=[abs(x)/z for x,z in zip(rows,scales)],
            quadratic_relative_residual=abs(quadratic_residual)/self.C,
            S_parts=dict(prefix=prefix_part,post_swirl=post_part,bump_quadratic_increment=dk,
                         bump_linear_coefficient=linear,bump_constant_increment=dc),
            only_M_correction_J_relative_error=abs(bad_J)/bad_scale,
            omitted_end_corrections_relative_error=mp.mpf(1),
            leading_root_S_error=prefix_part+post_part,
            post_mass=post)

    def source_check(self, eta, amplitude, max_step=.01):
        """Independent DOP853 integration of both normalized moment ODEs.

        Only the final correction patch is integrated. The incoming moments
        share the main-pulse quadrature above. The bump formula and ODE
        assembly here are independent of the weighted linear-system solve.
        """
        affine = self.affine(eta)
        coeff = [x*amplitude+y for x,y in zip(affine['a'],affine['b'])]
        scale = max(abs(v) for v in coeff)
        coeff = [float(v/scale) for v in coeff]
        slopes = [float(s) for s in self.slopes]
        start = -.15
        initial = [float(mp.exp(-s*start)*(pre+amplitude*main)/scale)
                   for s,pre,main in zip(self.slopes,affine['pre_scaled'],self.main_scaled)]
        def bump(t, center):
            x=(t-center)/.3+.5
            if not 0<x<1: return 0.
            v=1/x**2-1/(1-x)**2
            z=math.exp(-abs(v))
            product=z/(1+z)**2
            return product*(2/x**3+2/(1-x)**3)
        def rhs(t, state):
            force=coeff[0]*bump(t,0)+coeff[1]*bump(t,2)
            return [force-slopes[i]*state[i] for i in range(2)]
        flow=solve_ivp(rhs,(start,3),initial,method='DOP853',rtol=2e-13,atol=2e-15,max_step=max_step)
        if not flow.success: raise ArithmeticError('Source integration failed')
        return dict(max_step=max_step,normalized_end_moments=list(map(float,flow.y[:,-1])),
                    maximum_error=float(max(abs(v) for v in flow.y[:,-1])))


def family_increment(order=32, incoming_panels=16):
    """Finite Md=4 amplitude shift with an explicit omitted-bump bound.

    Uses the new lambda/h family and the complete post-pulse swirl. Axial
    bump energy is bounded by 40 lambda^41, giving an amplitude uncertainty
    below 120 lambda^41. This numerical root is a bounded proxy, not the
    exactly corrected amplitude or an admissible Md=4 construction.
    """
    if mp.mp.dps<260:
        raise ValueError('Use at least 260 digits to retain the finite increment')
    prefix = Prefix(order=order,incoming_panels=incoming_panels)
    post = Continuation(prefix,order=order,family=True)
    K, C = pulse_constant(prefix.rule), -mp.expm1(-26)/4
    A0=mp.sqrt(C/K)
    rows=[]
    for eta in (mp.mpf(0),mp.mpf('.5'),mp.mpf(1)):
        mass=post.energy_mass(eta)
        s0=prefix.moments(eta)['s0']
        remainder=prefix.lam*(s0-mp.exp(-26)*mass['total']/2)
        root=mp.sqrt((C-remainder)/K)
        # Rationalization preserves a shift much smaller than the root itself.
        increment=-remainder/(K*(root+A0))
        direct=root-A0
        rows.append(dict(eta=eta,prefix_s0=s0,post_mass=mass,
            amplitude_increment=increment,increment_over_lambda=increment/prefix.lam,
            direct_increment=direct,direct_relative_error=relative(direct,increment),
            omitted_bump_amplitude_error_upper=120*prefix.lam**41,
            binary64_increment=float(root)-float(A0)))
    return dict(T=prefix.T,lambda_value=prefix.lam,h=prefix.lam**2,K=K,
                leading_amplitude=A0,rows=rows)
