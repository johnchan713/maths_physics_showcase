"""Five moments on the first reserved power-law patch.

Rows are (M, (M-J)/lambda, I, -S, Cp), after the distinct dimensional
normalizations in README. The divided difference stabilizes the inverse but
does NOT remove the lambda^-1 factor in a general incoming discrepancy.
"""
from pathlib import Path
import sys
import mpmath as mp

sys.path.append(str(Path(__file__).resolve().parent.parent/'outer_pressure_pilot'))
from schedule import Rule, step_prime


def divided_weight(lam, y):
    lam, y = mp.mpf(lam), mp.mpf(y)
    if lam < 0:
        raise ValueError('lambda must be nonnegative')
    return y if lam == 0 else -mp.expm1(-lam*y)/lam


class Repair:
    U_CENTERS = ('1', '2')
    E_CENTERS = ('.5', '1.5', '2.5')
    WIDTH = '.1'

    def __init__(self, lam, order=32, panels=8):
        self.lam = mp.mpf(lam)
        if not 0 < self.lam <= mp.mpf('.01'):
            raise ValueError('Exact five-moment repair needs 0<lambda<=.01')
        self.rule = Rule(order, panels)
        self.width = mp.mpf(self.WIDTH)
        self.centers = list(map(mp.mpf,self.U_CENTERS+self.E_CENTERS))
        ordered = sorted(self.centers)
        if ordered[0] <= self.width/2 or ordered[-1]+self.width/2 >= 5 \
                or any(b-a <= self.width for a,b in zip(ordered,ordered[1:])):
            raise ValueError('All five supports must be disjoint and inside the patch')
        self.nodes = [[(c+self.width*(z-mp.mpf('.5')), w*step_prime(z),
                       w*step_prime(z)**2/self.width)
                       for z,w in self.rule.nodes] for c in self.centers]
        self.B = mp.matrix(5,5)
        for j in range(2):
            self.B[0,j] = self.integral(j, lambda y:mp.exp(y))
            self.B[1,j] = self.integral(j, lambda y:mp.exp(y)*divided_weight(self.lam,y))
        for j in range(2,5):
            for row,rate in ((2,mp.mpf('1.5')), (3,mp.mpf('.5')-self.lam),
                             (4,-mp.mpf('.5')-self.lam)):
                self.B[row,j] = self.integral(j, lambda y:mp.exp(rate*y))
        self.inverse = self.B**-1
        self.squares = [self.integral(j,lambda y:mp.exp(y),2) for j in range(5)]
        self.pressure_squares = [self.integral(j,lambda y:1,2) for j in range(2,5)]

    def integral(self, index, weight, power=1):
        if power not in (1,2):
            raise ValueError('Only the exact linear and quadratic moments are used')
        return mp.fsum(row[power]*weight(row[0]) for row in self.nodes[index])

    def quadratic(self, c, d):
        out = mp.matrix(5,1)
        out[3] = -mp.fsum(self.squares[j]*c[j]*d[j] for j in range(2)) \
            +mp.fsum(self.squares[j]*c[j]*d[j]/2 for j in range(2,5))
        out[4] = mp.fsum(self.pressure_squares[j-2]*c[j]*d[j]/2 for j in range(2,5))
        # The U/E supports are disjoint, so J has no quadratic cross term.
        return out

    def value(self, coefficients):
        c = mp.matrix(coefficients)
        return self.B*c+self.quadratic(c,c)

    def solve(self, target):
        """Numerical diagnostic root. Exact existence is the contraction bound."""
        target = mp.matrix(target)
        if max(abs(x) for x in target) > mp.mpf('1e-13'):
            raise ValueError('Target outside this diagnostic small-root domain')
        c = self.inverse*target
        for _ in range(12):
            residual = self.value(c)-target
            jac = self.B.copy()
            for j in range(5):
                unit = mp.matrix(5,1); unit[j] = 1
                column = 2*self.quadratic(c,unit)
                for i in range(5):
                    jac[i,j] += column[i]
            change = mp.lu_solve(jac,residual)
            c -= change
            if max(abs(x) for x in change) < mp.eps*max(abs(x) for x in c):
                break
        return c

    def ordinary_changes(self, coefficients):
        """Return normalized (M,I,J,S,Cp); tiny M-J can round away here.

        Use value()[1] as the stable divided-difference ledger at extreme lambda.
        """
        f = self.value(coefficients)
        return mp.matrix([f[0],f[2],f[0]-self.lam*f[1],-f[3],f[4]])

    def direct_changes(self, coefficients):
        """Independent delta-integrand evaluation in ordinary moment order.

        This control is used only at moderate lambda: unlike the stored map,
        the physical M/J subtraction has no numerical cancellation protection.
        """
        c = mp.matrix(coefficients)
        m,i,j,s,cp = [mp.mpf(0) for _ in range(5)]
        for index in range(5):
            for z,w in self.rule.nodes:
                y = self.centers[index]+self.width*(z-mp.mpf('.5'))
                bump = step_prime(z)/self.width
                du = c[index]*bump if index < 2 else 0
                de = c[index]*bump if index >= 2 else 0
                e0 = mp.exp((-mp.mpf('.5')-self.lam)*y)
                weight = self.width*w
                m += weight*mp.exp(y)*du
                i += weight*mp.exp(mp.mpf('1.5')*y)*de
                j += weight*mp.exp(mp.mpf('1.5')*y)*(e0*du+du*de)
                s += weight*mp.exp(y)*(du*du-e0*de-de*de/2)
                cp += weight*(e0*de+de*de/2)
        return mp.matrix([m,i,j,s,cp])


def precondition_target(ordinary, lam):
    """Transform supplied values; cannot recover an already rounded-away M-J.

    Actual tiny differences require a separate stable/interval moment ledger.
    """
    lam = mp.mpf(lam)
    if lam <= 0:
        raise ValueError('Zero lambda cannot repair independent M and J')
    m,i,j,s,cp = map(mp.mpf, ordinary)
    return mp.matrix([m,(m-j)/lam,i,-s,cp])
