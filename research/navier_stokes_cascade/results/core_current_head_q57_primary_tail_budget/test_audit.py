import json
import unittest
import sympy as s
import budget
import audit


class Tests(unittest.TestCase):
    def test_both_divisor_directions(self):
        b=s.sqrt(2)-1
        for x,y in [(0,1),(1,0),(29,70),(-70,29),(12,29),(-29,12)]:
            l=x-b*y; t=b*x+y
            self.assertEqual(s.expand(l*t-b*((x+y)**2-2*y*y)),0)
            self.assertGreaterEqual(abs((x+y)**2-2*y*y),1)
            self.assertTrue(bool(abs(l*t)>=b))

    def test_scale_at_crossover(self):
        for n in (0,7):
            for r in (0,1):
                d=budget.plan(n,r)
                self.assertEqual(d['amplitude'],n+199991-r)
                for row in d['shells']:
                    h=row['radial'];q=row['aux']
                    self.assertGreaterEqual(h-r,1)
                    self.assertEqual(-2*h+row['power'],-2*99995)
                    self.assertEqual(h+q+n+1,d['amplitude'])
                self.assertEqual([x['power'] for x in d['shells']],[2,-2])

    def test_real_curl_separation(self):
        A,B,C,D,k=s.symbols('A B C D k',nonzero=True)
        t=s.symbols('t',real=True)
        product=(A*s.cos(t)-B*s.sin(t)/k)*(C*s.cos(t)-D*s.sin(t)/k)
        self.assertEqual(s.simplify(s.integrate(product,(t,0,2*s.pi))/(2*s.pi)-(A*C+B*D/k**2)/2),0)

    def test_frozen(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
