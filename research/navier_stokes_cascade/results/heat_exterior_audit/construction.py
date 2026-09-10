"""Computable surrogates for the exact compensated heat reference.

Every Taylor and graded-root truncation has a separate positive error bound.
The continuum argument concerns the exact kernel and its exact small root.
"""
from pathlib import Path
import sys
import mpmath as mp
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.append(str(HERE.parent/'outer_pressure_pilot'))
from schedule import Rule, Schedule, step, step_prime  # noqa: E402


def relative(a,b):
    return abs(a-b)/max(abs(a),abs(b)) if a or b else mp.mpf(0)


def heat_coefficients(h,degree):
    """Coefficient of Z^n in (H-1)/h, with a separate next coefficient."""
    h = mp.mpf(h)
    if not 0<h<mp.mpf('.01') or degree<2:
        raise ValueError('Require 0<h<.01 and Taylor degree>=2')
    b = [mp.mpf(0),-(1+h)]
    for n in range(1,degree+1):
        b.append(-b[-1]*(h+n)*(h+1+n)/(n+1))
    return b


def heat_polynomial(h,z,degree=32):
    """Return H, stable (H-1)/(hZ), and an absolute Taylor remainder."""
    h,z = mp.mpf(h),mp.mpf(z)
    if z<0: raise ValueError('Heat argument must be nonnegative')
    b = heat_coefficients(h,degree)
    normalized = mp.fsum(b[n]*z**(n-1) for n in range(1,degree+1))
    delta = h*z*normalized
    first = h*mp.fsum(n*b[n]*z**(n-1) for n in range(1,degree+1))
    second = h*mp.fsum(n*(n-1)*b[n]*z**(n-2) for n in range(2,degree+1))
    return dict(value=1+delta,delta=delta,normalized_delta=normalized,
        remainder=h*abs(b[degree+1])*z**(degree+1),first_derivative=first,second_derivative=second,
        ode_direct=z*z*second+(1+2*(1+h)*z)*first+h*(1+h)*(1+delta),
        ode_residual=(degree+h)*(degree+1+h)*h*b[degree]*z**degree)


def heat_integral(h,z,order=48):
    """Independent positive Gamma-integral quadrature, normalized first."""
    h,z = mp.mpf(h),mp.mpf(z)
    nodes,weights = mp.gauss_quadrature(order,'glaguerre',alpha=h)
    weights = [w/mp.gamma(1+h) for w in weights]
    def quotient(v):
        if z==0: return -v
        logarithm = mp.log1p(z*v)
        argument = -h*logarithm
        return -logarithm/z*(mp.expm1(argument)/argument if argument else 1)
    return mp.fsum(w*quotient(v) for v,w in zip(nodes,weights))


class TailMoments:
    """Integrate a finite heat Taylor polynomial over the entire tail.

    Three compact pieces use positive quadrature; the infinite power pieces
    are integrated exactly. The angular n=1 term keeps h*(1/h) explicitly.
    """
    def __init__(self,h,log_xtail,order=32,panels=8,degree=32,co='.001'):
        self.h,self.log_xtail,self.co = mp.mpf(h),mp.mpf(log_xtail),mp.mpf(co)
        self.degree,self.order,self.panels = degree,order,panels
        self.b = heat_coefficients(self.h,degree)
        self.rule = Rule(order,panels)
        self.nodes = []
        for left,right in ((mp.mpf('.2'),mp.mpf('.5')),(mp.mpf('.5'),mp.mpf(1)),(mp.mpf(1),mp.mpf(3))):
            for z,w in self.rule.nodes:
                y = left+(right-left)*z
                chi = step((y-mp.mpf('.2'))/mp.mpf('.3'))
                fo = 1-self.co*self.h*step((3-y)/2)
                self.nodes.append((y,(right-left)*w,chi,fo))
        # Cache powers and fixed weights once, instead of nesting quadratures.
        self.integrals = {}
        for name,base,fo_power,chi_power in (
                ('P',1+2*self.h,2,1),('P2',1+2*self.h,2,2),
                ('S',2*self.h,2,1),('S2',2*self.h,2,2),('I',self.h-1,1,1)):
            weighted = [(mp.exp(-y),w*mp.exp(-base*y)*fo**fo_power*chi**chi_power)
                        for y,w,chi,fo in self.nodes]
            powers = [q for q,w in weighted]
            for n in range(1,2*degree+1):
                rate = (n-1)+self.h if name=='I' else n+base
                finite = mp.fsum(v*w for v,(q,w) in zip(powers,weighted))
                if name=='I' and n==1:
                    self.angular_finite = self.h*finite
                    self.angular_leading = mp.exp(-3*self.h)
                    value = self.angular_finite+self.angular_leading
                else:
                    value = finite+mp.exp(-3*rate)/rate
                self.integrals[name,n] = value
                powers = [v*q for v,(q,w) in zip(powers,weighted)]
        self.convolution = {n:mp.fsum(self.b[j]*self.b[n-j]
            for j in range(max(1,n-degree),min(degree,n-1)+1))
            for n in range(2,2*degree+1)}

    def normalized(self,eta):
        """P/(et^2 h/Xtail), S/(et^2 h), I/(sqrt(2) et sqrt(Xtail))."""
        eta = mp.mpf(eta)
        if abs(eta)>1: raise ValueError('eta must lie in [-1,1]')
        d = 1-eta*eta
        terms = {k:[] for k in ('P_linear','P_quadratic','S_linear','S_quadratic','I')}
        derivatives = {k:[] for k in terms}
        for n in range(1,2*self.degree+1):
            power = (2*d)**n
            derivative = -4*eta*n*(2*d)**(n-1)
            scale = mp.exp(-(n-1)*self.log_xtail)
            if n<=self.degree:
                for name,sign in (('P',1),('S',-1)):
                    factor = sign*self.b[n]*scale*self.integrals[name,n]
                    terms[name+'_linear'].append(factor*power)
                    derivatives[name+'_linear'].append(factor*derivative)
                factor = self.b[n]*scale*self.integrals['I',n]*(1 if n==1 else self.h)
                terms['I'].append(factor*power)
                derivatives['I'].append(factor*derivative)
            if n>=2:
                for name,sign in (('P',1),('S',-1)):
                    factor = sign*self.h*self.convolution[n]*scale*self.integrals[name+'2',n]/2
                    terms[name+'_quadratic'].append(factor*power)
                    derivatives[name+'_quadratic'].append(factor*derivative)
        values = {k:mp.fsum(v) for k,v in terms.items()}
        jets = {k:mp.fsum(v) for k,v in derivatives.items()}
        return dict(values=values,derivatives=jets,
            totals=[values['P_linear']+values['P_quadratic'],values['S_linear']+values['S_quadratic'],values['I']],
            eta_totals=[jets['P_linear']+jets['P_quadratic'],jets['S_linear']+jets['S_quadratic'],jets['I']],
            angular_first_finite=-2*d*(1+self.h)*self.angular_finite,
            angular_first_infinite=-2*d*(1+self.h)*self.angular_leading,
            angular_first_infinite_change_from_h0=-2*d*(mp.expm1(-3*self.h)+self.h*mp.exp(-3*self.h)),
            first_order_angular_tail_omitted_at_y1000=-2*d*(1+self.h)*mp.exp(-1000*self.h))

    def target(self,eta,log_xstar,log_hpow_ratio):
        row = self.normalized(eta)
        eta,log_xstar,log_hpow_ratio = map(mp.mpf,(eta,log_xstar,log_hpow_ratio))
        f_inverse,fp_inverse = 1+eta*eta,2*eta
        log_ratio = self.log_xtail-log_xstar
        logs = [mp.log(self.h)+2*log_hpow_ratio-log_xstar-2*log_ratio,
                mp.log(self.h)+2*log_hpow_ratio-self.log_xtail,
                log_hpow_ratio-log_xstar]
        signs,powers = (-1,1,-1),(2,2,1)
        target,derivative = [],[]
        for v,ve,logscale,sign,power in zip(row['totals'],row['eta_totals'],logs,signs,powers):
            scale = sign*mp.exp(logscale)
            target.append(scale*v*f_inverse**power)
            derivative.append(scale*(ve*f_inverse**power+v*power*f_inverse**(power-1)*fp_inverse))
        # A uniform C1 bound on the exact-kernel/Taylor-target difference.
        n = self.degree+1
        common = self.h*abs(self.b[n])*(2*self.degree+3)*mp.exp(n*(mp.log(2)-self.log_xtail))
        rates = [n+1+2*self.h,n+2*self.h,n-1+self.h]
        amplitudes = [36*mp.exp(2*log_hpow_ratio-log_ratio),
                      36*mp.exp(2*log_hpow_ratio),4*mp.exp(log_hpow_ratio+log_ratio)]
        errors = [a*common*mp.exp(-mp.mpf('.2')*k)/k for a,k in zip(amplitudes,rates)]
        return dict(eta=eta,normalized=row,target=target,target_eta=derivative,
            target_log_scales=logs,taylor_target_C1_error=errors)


class Compensation:
    """Exact linear-plus-quadratic moment map for three relative swirl bumps."""
    def __init__(self,lam,rule):
        self.lam,self.rule = mp.mpf(lam),rule
        self.centers = tuple(map(mp.mpf,('.5','2.5','4.5')))
        self.width = mp.mpf('.3')
        rates = (-1-2*self.lam,-2*self.lam,1-self.lam)
        self.B,self.Q = mp.matrix(3,3),mp.matrix(2,3)
        for k,rate in enumerate(rates):
            for j,center in enumerate(self.centers):
                self.B[k,j] = self.width*rule.integrate(lambda z:
                    mp.exp(rate*(center+self.width*(z-mp.mpf('.5'))))*step_prime(z))
                if k<2:
                    self.Q[k,j] = self.width/2*rule.integrate(lambda z:
                        mp.exp(rate*(center+self.width*(z-mp.mpf('.5'))))*step_prime(z)**2)
        self.inverse = self.B**-1

    def quadratic(self,u,v):
        return mp.matrix([mp.fsum(self.Q[k,j]*u[j]*v[j] for j in range(3)) if k<2 else 0 for k in range(3)])

    def newton(self,target,iterations=10):
        target = mp.matrix(target)
        c = self.inverse*target
        for _ in range(iterations):
            residual = self.B*c+self.quadratic(c,c)-target
            jac = self.B.copy()
            for k in range(2):
                for j in range(3): jac[k,j] += 2*self.Q[k,j]*c[j]
            c -= mp.lu_solve(jac,residual)
        return c

    def derivative(self,c,target_eta):
        jac = self.B.copy()
        for k in range(2):
            for j in range(3): jac[k,j] += 2*self.Q[k,j]*c[j]
        return mp.lu_solve(jac,mp.matrix(target_eta))

    def graded(self,target,degree=4):
        """Keep every power of a tiny correction separate; bound the tail."""
        linear = self.inverse*mp.matrix(target)
        scale = max(abs(v) for v in linear)
        if not scale:
            return dict(scale=mp.mpf(0),grades=[],components=[],remainder=mp.mpf(0),row_errors=[])
        grades = [linear/scale]
        for n in range(1,degree):
            rhs = mp.matrix([0,0,0])
            for j in range(n): rhs += self.quadratic(grades[j],grades[n-1-j])
            grades.append(-self.inverse*rhs)
        row_errors = []
        for n,g in enumerate(grades):
            residual = self.B*g-(mp.matrix(target)/scale if n==0 else mp.matrix([0,0,0]))
            if n:
                for j in range(n): residual += self.quadratic(grades[j],grades[n-1-j])
            row_errors.append(max(abs(v) for v in residual)/(1+max(abs(v) for v in self.B*g)))
        ratio = 4*50*100*scale
        if not ratio<mp.mpf('.5'): raise ValueError('Graded root outside the proved small ball')
        return dict(scale=scale,grades=[list(g) for g in grades],
            components=[[scale**(n+1)*v for v in g] for n,g in enumerate(grades)],
            remainder=scale*ratio**degree/(1-ratio),row_errors=row_errors)

    def independent_moments(self,c,cadence=.01):
        """Forward source ODE, normalized before conversion to binary64."""
        scale = max(abs(v) for v in c)
        if not scale: return dict(gaps=[0,0,0],quadratic_omission_upper=0)
        normalized = [v/scale for v in c]
        def rhs(t,state):
            value = mp.fsum(v*step_prime((mp.mpf(t)-center)/self.width+mp.mpf('.5'))
                           for v,center in zip(normalized,self.centers))
            quadratic = scale*value*value/2
            return [float(mp.exp((-1-2*self.lam)*t)*(value+quadratic)),
                    float(mp.exp(-2*self.lam*t)*(value+quadratic)),
                    float(mp.exp((1-self.lam)*t)*value)]
        solution = solve_ivp(rhs,(0,5),[0,0,0],method='DOP853',rtol=2e-12,atol=2e-13,max_step=cadence)
        if not solution.success: raise RuntimeError(solution.message)
        direct = (self.B*mp.matrix(c)+self.quadratic(c,c))/scale
        return dict(gaps=[abs(mp.mpf(a)-b)/(1+abs(b)) for a,b in zip(solution.y[:,-1],direct)],
            quadratic_omission_upper=100*scale,source_final=list(solution.y[:,-1]),formula_final=list(direct))


def selected_fixture(order=32,panels=8):
    """Actual Md4 member of the fixed family; its old axial failure persists."""
    T = mp.exp(4)+10
    lam,h = mp.exp(-4*T),mp.exp(-8*T)
    rule = Rule(order,panels)
    shell = object.__new__(Schedule)
    shell.rule,shell.h = rule,h
    q1 = Schedule.propagate_q(shell,lam,-lam,-1)
    hold = 4*mp.log(1/h)
    before_wait = Schedule.propagate_q(shell,q1+(1-h)*hold,-1,-h)
    co = mp.mpf('.001')
    qp = co*h/(1-co*h)*rule.integrate(lambda z:mp.exp((1-h)*(1+2*z))*step_prime(z))
    wait = mp.log(before_wait/qp)/(1-h)
    tw,tf,tu = 240*T,mp.mpf(1000),120*T
    log_xstar = mp.log(100)+241*T-18
    log_ratio = 20+13/lam+tf+tu+2+hold+wait
    log_rh = -14-mp.log(2)-(lam+h)/2-hold-lam*(20+tf+tu)-h*wait-mp.log1p(-co*h)
    log_estar = -mp.mpf('119.5')*T+mp.mpf('10.3')-lam*(tw-mp.mpf('19.5'))
    return dict(Md=4,T=T,lam=lam,h=h,log_xstar=log_xstar,log_xtail=log_xstar+log_ratio,
        log_hpow_ratio=log_rh,log_estar=log_estar,wait=wait,
        old_axial_failure_retained=True)


def pole_control(h,radius=1000):
    """At eta=1 the edit value vanishes; its derivative repairs angular stress."""
    h,X = mp.mpf(h),mp.mpf(radius)
    D,L = mp.mpf('.5')-h,1-2*h
    a = 2+2*h
    I_eta_over_Hpow = -4*(1+h)
    Q = -D*I_eta_over_Hpow/X
    return dict(heat_value_at_pole=1,heat_value_edit=0,I_eta_over_Hpow=I_eta_over_Hpow,
        I_eta_over_Hpow_correction=-4*h,shear_excess=2*h,
        ps1=X*Q/L,ps1_minus_two=X*Q/L-2,a=a,stress_over_F=X*Q/L-a,
        omitted_derivative_ps1=0,omitted_derivative_stress_over_F=-a)
