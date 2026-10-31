import json, unittest
from fractions import Fraction as Q
import flux_tail as q
import audit
t=q.old
def state(m=0):
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
        hull=hull, fast_hull=(1,8), max_output=m+3, profile=[10000]*(m+17), profile_hull=hull)
    return r.maps.NewPoleMaps(r.maps.update.SignedCovariance(target,m))


def coeffs(v):
    return {j: dict(a.terms) for j,a in v.c.items()}


class Tests(unittest.TestCase):
    def test_characteristic_flux_and_surviving_mode(self):
        import sympy as s
        a=s.symbols('a',real=True)
        # One Fourier mode on U in [1,2], d=1, e=1; nonzero connection.
        R=1+a; F=a**3*(1-a)**3; phase=s.exp(s.I*s.pi*a)
        graph=s.diff(F,a)+s.I*s.pi*F+F/R
        self.assertEqual(s.simplify(R*graph*phase-s.diff(R*F*phase,a)),0)
        self.assertEqual(s.simplify(s.integrate(R*graph*phase,(a,0,1))),0)
        I=s.integrate(R*F*phase,(a,0,1))
        self.assertNotEqual(s.simplify(I),0) # surviving slow derivative
        self.assertNotEqual(s.simplify(-s.I*s.pi*I),0) # shift mismatch

    def test_operator(self):
        c=state();p=(2,4,[1]*3);a=(3,6,[1]*2)
        tail=q.InitialTail(c,0,d=2,pressure_cutoff=p,axial_cutoff=a)
        o=tail(0,0); self.assertEqual((o['order'],o['pole']),(0,0))
        W=q.full.convert(tail.b.primary.covariance_R(3)['W'])
        K=q.full.compact(c,2,2,4,[1]*3,lambda n,v:1,1,54,0)(1,54)['KT']
        pm=K*(4*(2+6*64*10)+4+8)*W
        ax,shell=tail.ax
        G=ax.geo;f=G.multipliers(0,1)
        base=(G.Umax-G.Umin)*f['f']*f['k']*ax(0,54)['bounds']['chi_U']
        want=8*base*shell.bounds(54)['edge']*shell.bounds(0)['interior']*(W+pm)
        self.assertEqual(coeffs(o['Ez']),coeffs(want))
        with self.assertRaises(ValueError):tail(1,0)
        with self.assertRaises(ValueError):q.InitialTail(c,0,d=2,pressure_cutoff=(2,4,[1]*2),axial_cutoff=a)

    def test_factory(self):
        c=state(11)
        obj=q.build(c,0,d=1,pressure_cutoff=(2,4,[1]*13),axial_cutoff=(2,4,[1]*12),
                    c0=2,c1=4,shape=[1]*15,time_inverse=lambda n,v:1)
        rows=obj.current.fast.rows
        self.assertIs(rows.initial_A1.b,rows.initial.budget)
        self.assertEqual(coeffs(rows.rows(0)['initial_A1']),coeffs(rows.initial_A1(0,rows.pole)['Ez']))

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((q.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
