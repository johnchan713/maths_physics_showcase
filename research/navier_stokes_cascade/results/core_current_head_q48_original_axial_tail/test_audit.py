from fractions import Fraction as Q
import json
import math
import unittest
import native_tail as t
import audit


def state():
    # Manufactured lower norm data; construct an actual same-family primary.
    r = t.full.refresh
    lib = r.maps.update.profile.pipeline
    one = lambda j: lib.roots.P(int(j == 0))
    root = lib.roots.RootBudget(H=one, target=one, wave=one, curl=one,
        loghalf=lambda j: lib.roots.P(), K0=1, beta=1, R=1, mass=2, width=Q(1,4))
    hull = dict(l=2, u=3, Xa=1, Xb=8, t1=root.width)
    target = r.maps.update.profile.build(root,
        lambda j: dict(target=1, correction=1, pole=3*j+3),
        Rmin=1, Rmax=8, fixed=lambda j,v: {k:1 for k in lib.fresh.FIXED if k!='rho'},
        compact=lambda j,v:dict(KT=1, KA=1, I=3), slow=lambda j:{k:1 for k in ('H','Kp','Kz','Kg','Kpsi')},
        hull=hull, fast_hull=(1,8), max_output=3, profile=[10000]*17, profile_hull=hull)
    return r.maps.NewPoleMaps(r.maps.update.SignedCovariance(target,0))


def coeffs(v):
    return {j: dict(a.terms) for j,a in v.c.items()}


class Tests(unittest.TestCase):
    def test_original_shell_recovery(self):
        geometry = t.full.refresh.maps.update.profile.pipeline.fresh.geometry.geometry
        for xa, xb in ((Q(4,5),Q(13)), (Q(1,7),Q(25,3))):
            self.assertEqual(t.initial_shell(geometry(2,3,xa,xb,Q(1,5))), (xa,xb))

    def test_shell_divisors(self):
        # Numerical diagnostics for both signs, many nearly resonant directions.
        b = math.sqrt(2)-1
        for K in (1,2,4,8,16,32):
            a = sorted(x-b*y for x in range(-2*K,2*K+1)
                       for y in range(-2*K,2*K+1) if K<=abs(x)+abs(y)<2*K)
            inc = sorted([0]+a)
            self.assertGreater(min(y-x for x,y in zip(inc,inc[1:])), 1/(12*K))
            for p in (1,2,3):
                self.assertLess(sum(abs(x)**(-2*p) for x in a),4*(12*K)**(2*p))

    def test_finite_shifted_mode(self):
        import sympy as s
        u, v = s.symbols('u v', real=True)
        # y Fourier mode one, 2*pi*M*vr.k = pi; u is the output point.
        phi = v**4*(1-v)**4
        e = s.exp(s.I*s.pi*v)
        integral = lambda f:s.integrate(f*e,(v,0,1))*s.exp(-s.I*s.pi*u)
        a = integral(phi)
        self.assertNotEqual(s.simplify(a),0)
        self.assertEqual(s.simplify(s.diff(a,u) - integral(s.diff(phi,v))),0)
        self.assertEqual(s.simplify(a - integral(s.diff(phi,v,2))/(s.I*s.pi)**2),0)

    def test_projected_source_and_geometry(self):
        c = state()
        L = t.ledger(0,1)
        args = dict(d=2, pressure_cutoff=(2,4,[1]*6), axial_cutoff=(3,5,[1,2]),
                    covariance_coverage=6, pole_coverage=90)
        b = c.target.source.haar.initial.budget
        o = t.coefficient(c,0,1,**args)
        D = 7
        W = t.full.convert(b.primary.covariance_R(6)['W'])
        g = 2+6*b.primary.chart.B*28
        fullK = t.full.compact(c,2,2,4,[1]*6,lambda j,v:1,4,90,0)(4,90)['KT']
        pm = fullK*(D*g+D+14)*W
        self.assertEqual(coeffs(pm),coeffs(o['full_pressure']))
        self.assertEqual(coeffs(o['projected_source']),coeffs(2*((D*g+7)*W+D*(W+pm))))
        self.assertEqual(coeffs(o['floor_base']),coeffs(t.EP({0:1,1:3})))
        self.assertEqual(o['floor_power'],1)
        self.assertTrue(o['coefficient'].c)
        # Shrinking only the axial plateau gap worsens the output weight floor.
        a = t.coefficient(c,0,1,**dict(args,axial_cutoff=(Q(9,8),5,[1,2])))
        self.assertNotEqual(coeffs(a['coefficient']),coeffs(o['coefficient']))
        for change in (dict(covariance_coverage=5),dict(pole_coverage=89),
                       dict(pressure_cutoff=(2,4,[1]*5)),dict(axial_cutoff=(3,5,[1]))):
            with self.assertRaises(ValueError):t.coefficient(c,0,1,**dict(args,**change))
        self.assertEqual(L['output_pole'],0)
        self.assertEqual(t.ledger(10)['covariance'],200012)

    def test_evidence(self):
        self.assertEqual(audit.generate(), json.loads((t.HERE/'evidence.json').read_text()))


if __name__ == '__main__':
    unittest.main()
