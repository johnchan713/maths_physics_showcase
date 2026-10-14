"""Finite, manufactured global extension with nonzero radial memory on I3.

The leading memory is M0=eta*B(X): B rises to one, stays one throughout
the repair patch, and returns to zero only in a later pulse. The piecewise
polynomial B is C4, sufficient for these finite diagnostics, not the actual
analytic leading field or the continuum C-infinity existence argument.
"""
import mpmath as mp
from extension import Repair
from schedule import Rule, step_prime


def integrate_vector(rule, function):
    """Evaluate the vector integrand once per node, with compensated sums."""
    values = [(weight, function(x)) for x, weight in rule.nodes]
    return mp.matrix([mp.fsum(weight*row[i] for weight, row in values)
                      for i in range(len(values[0][1]))])


def smooth_polynomial(x, derivative=0):
    """Degree-nine step with the first four endpoint derivatives zero."""
    coefficients = {5: 126, 6: -420, 7: 540, 8: -315, 9: 70}
    return mp.fsum(mp.mpf(v)*mp.factorial(n)/mp.factorial(n-derivative)*x**(n-derivative)
                   for n, v in coefficients.items() if n >= derivative)


class MemoryFixture:
    def __init__(self, radius='2', h='.005'):
        self.radius, self.h = mp.mpf(radius), mp.mpf(h)
        self.pulse_start = self.radius*mp.exp(12)
        self.support_end = self.radius*mp.exp(13)

    def memory_jets(self, X):
        if X <= 0 or X >= self.support_end:
            return [mp.mpf(0)]*4
        if X < 1:
            return [smooth_polynomial(X, n) for n in range(4)]
        if X <= self.pulse_start:
            return [mp.mpf(1), mp.mpf(0), mp.mpf(0), mp.mpf(0)]
        width = self.support_end-self.pulse_start
        t = (X-self.pulse_start)/width
        return [1-smooth_polynomial(t)] + [-smooth_polynomial(t, n)/width**n for n in range(1, 4)]

    def omega(self, X, eta):
        if X <= 0 or X >= self.support_end:
            return mp.mpf(0)
        B, u, ux, uxx = self.memory_jets(X)
        h, D, L, d = self.h, mp.mpf('.5')-self.h, 1-2*self.h*eta**2, 1-eta**2
        a, ae = 2*eta**2/L, 4*eta/L**2
        V = a*X*u-B
        Vx = (a-1)*u+a*X*ux
        Vxx, Ve = (2*a-1)*ux+a*X*uxx, ae*X*u
        return ((X*Vx+D*eta*Ve)/L + V*(Vx-V/(2*X))
                +eta*u*(-2*eta*X*Vx+d*Ve)/L-2*X*Vxx)

    def source_integrals(self, eta, rule):
        """Integrals of Omega and Omega/X, including the uncut outer pulse."""
        mass = rule.integrate(lambda x: self.omega(x, eta))
        pressure = rule.integrate(lambda x: self.omega(x, eta)/x)
        # On the entire intervening region, U0=0 but V0=-1.
        mass -= mp.log(self.pulse_start)/2
        pressure -= (1-1/self.pulse_start)/2
        width = self.support_end-self.pulse_start
        mass += width*rule.integrate(lambda t: self.omega(self.pulse_start+width*t, eta))
        pressure += width*rule.integrate(lambda t: self.omega(self.pulse_start+width*t, eta)/(self.pulse_start+width*t))
        return mass, pressure

    @staticmethod
    def inner(X, eta):
        if not 0 < X < 1:
            return mp.mpf(0), mp.mpf(0)
        bump = X*(1-X)**6
        return (1-eta/5)*bump, (1+eta/3)*bump

    def incoming(self, eta, rule, omit_omega=False):
        """Actual ordinary integrals of this fixture's tentative extension."""
        f = 1/(1+eta*eta)
        def terms(x):
            F1, U1 = self.inner(x, eta)
            U0 = eta*self.memory_jets(x)[1]
            return [U1, 2*x*F1, 2*f*F1,
                    2*x*(U0*F1+U1*f), 2*U0*U1-2*x*f*F1]
        values = integrate_vector(rule, terms)
        if not omit_omega:
            mass, pressure = self.source_integrals(eta, rule)
            values[2] -= pressure/2
            values[4] += mass/2
        return mp.matrix(values)

    def tail_pressure(self, X, eta, rule):
        """Pi1 outside I3, from zero pressure at infinity, not an artificial cutoff."""
        if not 1 <= X <= self.pulse_start:
            raise ValueError('This diagnostic formula is for the constant-memory region')
        width = self.support_end-self.pulse_start
        outer = width*rule.integrate(lambda t: self.omega(self.pulse_start+width*t, eta)/(self.pulse_start+width*t))
        return (-(1/X-1/self.pulse_start)/2+outer)/2


def direct_patch_changes(coefficients, repair, radius, amplitude, rule):
    """Independent physical integrands in X=rho*exp(y), not a call to B*c."""
    out = mp.matrix(5, 1)
    for index, center in enumerate(repair.centers):
        def values(z):
            y = center+repair.width*(z-mp.mpf('.5'))
            X = radius*mp.exp(y)
            bump = step_prime(z)/repair.width
            U1, E1 = (coefficients[index]*bump, 0) if index < 2 else (0, coefficients[index]*bump)
            E0 = amplitude*mp.exp((-mp.mpf('.5')-repair.lam)*y)
            # Lower-order Omega is unchanged by this FIRST-order variation.
            density = [U1, mp.sqrt(2*X)*E1, E0*E1/X,
                       mp.sqrt(2*X)*E0*U1, -E0*E1]
            return [repair.width*X*v for v in density]
        out += integrate_vector(rule, values)
    return out
