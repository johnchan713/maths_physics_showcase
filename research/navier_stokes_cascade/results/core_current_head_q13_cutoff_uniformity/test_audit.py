import itertools
import json
import math
import unittest
import audit
import cutoff as c
Q=c.Q

class CutoffTests(unittest.TestCase):
    def test_first_derivatives(self):
        for h in (Q(0),Q(1,513),Q(1,10**30),Q(1,100)):
            for e in (Q(-1),Q(0),Q(1,3),Q(1)):
                L=1-2*h*e*e
                for word,expected in [('z',2*e/L),('t',-1/L)]:
                    beta,rows=c.coefficients(word,h)
                    self.assertEqual(rows[0],())
                    self.assertEqual(c.evaluate(rows[1],e)/L**2,expected)
                    self.assertEqual(beta,-(Q(1,2)-h) if word=='z' else -1)
    def test_mixed_partials_commute_as_polynomials(self):
        for h in (Q(0),Q(1,513),Q(1,10**30),Q(1,100)):
            for n in range(2,6):
                for letters in itertools.product('zt',repeat=n):
                    word=''.join(letters)
                    self.assertEqual(c.coefficients(word,h),c.coefficients('z'*word.count('z')+'t'*word.count('t'),h))
    def test_zero_order_and_positive_order_seams(self):
        self.assertEqual(c.coefficients('',Q(1,513)),(0,((Q(1),),)))
        for n in range(1,8):
            for word in ('z'*n,'t'*n):
                _,rows=c.coefficients(word,Q(1,513))
                self.assertEqual(rows[0],())
                self.assertEqual(len(rows),n+1)
    def test_independent_h_zero_monomial_time_derivatives(self):
        for n in range(1,6):
            for d in range(n,8):
                beta,rows=c.coefficients('t'*n,0)
                for e in (Q(-1),Q(1,3),Q(1)):
                    result=sum(c.evaluate(poly,e)*math.factorial(d)//math.factorial(d-k) for k,poly in enumerate(rows))
                    self.assertEqual(result,(-1)**n*math.factorial(d)//math.factorial(d-n))
    def test_independent_h_zero_second_spatial_derivative(self):
        _,rows=c.coefficients('zz',0)
        for d in range(2,8):
            for e in (Q(-1),Q(0),Q(1,3)):
                result=sum(c.evaluate(poly,e)*(math.factorial(d)//math.factorial(d-k)) for k,poly in enumerate(rows))
                self.assertEqual(result,2*d+4*d*(d-1)*e*e)
    def test_coefficient_norm_majorant(self):
        for word in ('z','tt','ztz','zztt'):
            for h in (Q(1,513),Q(1,100)):
                _,rows=c.coefficients(word,h)
                for e in (Q(-1),Q(-1,3),Q(1,2),Q(1)):
                    exact=sum(abs(c.evaluate(poly,e)) for poly in rows)/(1-2*h*e*e)**(2*len(word))
                    self.assertLessEqual(exact,c.coefficient_bound(word,h))
    def test_leibniz_curl_loss(self):
        for m in range(20):
            for k in range(m+2):
                d=m+1-k
                self.assertLessEqual(2*k+d,2*(m+1))
    def test_bad_inputs(self):
        for h in (-1,Q(1,99),0.005,True):
            with self.assertRaises(ValueError):c.coefficients('z',h)
        with self.assertRaises(ValueError):c.coefficients('x',0)
    def test_evidence_and_scope(self):
        r=audit.generate()
        self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['scope']['actual_representative_suprema_proved'])
        self.assertFalse(r['scope']['unconditional_residual_or_blowup'])

if __name__=='__main__':unittest.main()
