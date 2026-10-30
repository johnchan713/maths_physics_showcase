import unittest,json
from fractions import Fraction as Q
import pipeline as p
import audit

def fixture(trace=None):
    P=p.roots.P
    root=p.roots.RootBudget(H=lambda j:P(int(j==0)),target=lambda j:P(int(j==0)),
        wave=lambda j:P(int(j==0)),curl=lambda j:P(),loghalf=lambda j:P(),
        K0=1,beta=1,R=1,mass=2,width=Q(1,4))
    stress=lambda n:dict(target=Q(1),correction=Q(1),pole=3*n+3)
    kw=dict(Rmin=1,Rmax=8,fixed=lambda N,pole:{k:1 for k in p.fresh.FIXED},
        compact=lambda n,pole:dict(KT=1,KA=1,I=6),
        slow=lambda n:{k:1 for k in ('H','Kp','Kz','Kg','Kpsi')},
        hull=dict(l=2,u=3,Xa=1,Xb=8,t1=root.width),fast_hull=(1,8),
        max_output=0,c_rho=lambda n:p.haar.EP(n+1),trace=trace)
    return root,stress,kw

class PipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.trace=set();root,stress,kw=fixture(cls.trace)
        cls.v=p.build(root,stress,**kw);cls.data=cls.v.inputs(0)

    def test_real_recipe_covariance_ceiling(self):
        self.assertEqual(max(n for k,n in self.trace if k=='covariance'),13)
        self.assertEqual(self.v.pole,174)
        self.assertIn(('higher_stress',3),self.trace)

    def test_corrected_background_count_once(self):
        _,stress,_=fixture();n=3
        expected=p.norm.count(n)*p.norm.common.normal.transfer(n,stress(n),1,8)['higher_stress_raw']
        self.assertEqual(self.v.higher(n,self.v.pole).record(),p.haar.EP(expected).record())
        with self.assertRaises(ValueError):self.v.higher(n,1)

    def test_complete_initial_slots_derived(self):
        self.assertEqual(set(self.data),set(p.haar.old.KEYS))
        self.assertEqual(self.data['PiP'].record(),self.v.initial.output(1)['P'].record())
        self.assertEqual(self.data['Ft'].record(),self.v.initial.budget._fast(3)['u'].record())
        self.assertEqual(self.data['U'].record(),self.v.initial.slow_fields(3,True)['U'].record())

    def test_moving_pressure_and_distinct_channels(self):
        kw=dict(pullback=lambda n,e:p.haar.EP(1),bump=lambda n,e:p.haar.EP(1))
        axial=self.v.projected_source(0,1,**kw);angular=self.v.projected_source(0,2,**kw)
        t=self.v.initial.output(1)
        self.assertEqual(axial['moment'].record(),(t['Jz']+2*t['P']).record())
        self.assertEqual(angular['moment'].record(),t['Jtheta'].record())
        self.assertNotEqual(axial['moment'].record(),t['Jz'].record())
        self.assertFalse(axial['signed_primitive_or_root_certified'])

    def test_root_copy_and_no_new_value_data(self):
        root,_,_=fixture();copy=p.copy_root(root)
        for j in range(3):
            for k in ('H','target','wave','curl','loghalf'):
                self.assertEqual(getattr(root,k)(j).c,getattr(copy,k)(j).c)
        self.assertEqual((root.beta,root.width,root.mass),(copy.beta,copy.width,copy.mass))

    def test_guarded_output_and_full_hulls(self):
        with self.assertRaises(ValueError):self.v.inputs(1)
        root,stress,kw=fixture();kw['fast_hull']=(0,8)
        with self.assertRaises(ValueError):p.build(root,stress,**kw)
        kw['fast_hull']=(1,9)
        with self.assertRaises(ValueError):p.build(root,stress,**kw)

    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['complete_residual_or_blowup'])

if __name__=='__main__':unittest.main()
