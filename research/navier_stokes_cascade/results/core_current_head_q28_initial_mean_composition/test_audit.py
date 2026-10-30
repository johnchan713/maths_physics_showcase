import unittest,json
from fractions import Fraction as Q
import sympy as s
import composition as c
import audit

def providers(trace):
    base=c.source.fixture()
    out={k:getattr(base,k) for k in c.KEYS}
    old=out['wave']
    def wave(n):trace.add(n);return old(n)
    out['wave']=wave
    return out

class CompositionTests(unittest.TestCase):
    def test_complete_two_pass_request_depth(self):
        for j in (0,1):
            p,signed=set(),set();first,last=c.compose(providers(p),providers(signed))
            self.assertGreater(last.raw_source(j).evaluate(1),0)
            self.assertEqual(max(p),j+27);self.assertEqual(max(signed),j+15)
            for k in c.source.STATE:
                self.assertEqual(last.initial(j)[k].c,first.final(j)[k].c)
                self.assertGreater(last.initial(j)[k].evaluate(1),0)
                self.assertGreater(last.final(j)[k].evaluate(1),last.initial(j)[k].evaluate(1))

    def test_only_first_pass_starts_zero(self):
        first,last=c.compose(providers(set()),providers(set()))
        self.assertTrue(all(not v.c for v in first.initial(0).values()))
        self.assertTrue(all(v.evaluate(1)>0 for v in last.initial(0).values()))
        for stage in (first,last):
            self.assertGreater(stage.rows(0,True)['pm'].evaluate(1),0)
            self.assertGreater(stage.rows(0,True)['P'].evaluate(1),0)
            self.assertGreater(stage.rows(0,True)['Jz'].evaluate(1),0)

    def test_operator_stages_remain_distinct(self):
        a,b=providers(set()),providers(set());a['KT1']=lambda n:c.P(2);b['KT1']=lambda n:c.P(7)
        first,last=c.compose(a,b)
        self.assertEqual(first.KT1(0).evaluate(1),2);self.assertEqual(last.KT1(0).evaluate(1),7)

    def test_reset_loses_nonzero_terms(self):
        a,b=providers(set()),providers(set());first,last=c.compose(a,b)
        reset=c.source.PrefixBudget(**b,initial=lambda n:{k:c.P() for k in c.source.STATE})
        self.assertGreater(last.source(0)['complete_mean_wave'].evaluate(1),reset.source(0)['complete_mean_wave'].evaluate(1))

    def test_pole_archive_not_only_maximum(self):
        r=c.poles(3,8)
        self.assertEqual(r['pressure_archive'],[6,12,16,32])
        self.assertEqual(r['source'],24);self.assertEqual(r['initial_mean'],6)
        self.assertEqual(c.poles(5,2)['source'],30)
        self.assertEqual(c.poles(0,0)['pressure_archive'],[0])

    def test_projection_fourier_coefficient_and_derivative(self):
        y,z=s.symbols('y z',real=True);E=z*z+3*s.cos(y)+z*s.sin(2*y)
        mean=s.integrate(E,(y,0,2*s.pi))/(2*s.pi)
        for n in range(1,5):self.assertEqual(s.diff(E-mean,y,n),s.diff(E,y,n))
        for k in (1,2):
            mode=s.cos(k*y)
            self.assertEqual(s.integrate(mean*mode,(y,0,2*s.pi)),0)
        # Sup norm projection is not contractive: [1,1,-1] has centered max 4/3.
        values=[Q(1),Q(1),Q(-1)];avg=sum(values)/3
        self.assertGreater(max(abs(v-avg) for v in values),max(map(abs,values)))

    def test_whole_pressure_moment_and_axial_repair(self):
        z=s.symbols('z');rho=1+z*z;P=1+z
        self.assertNotEqual(s.diff(rho*P,z),rho*s.diff(P,z))
        gd=1+z;A=z*z
        self.assertNotEqual(s.expand((gd-A)**2),s.expand(gd**2))

    def test_guards(self):
        p=providers(set());p['initial']=lambda n:None
        with self.assertRaises(ValueError):c.compose(p,providers(set()))
        with self.assertRaises(ValueError):c.poles(-1,2)
        with self.assertRaises(ValueError):c.poles(True,2)

    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['fixed_map_and_background_premises_discharged'])
if __name__=='__main__':unittest.main()
