import unittest,json,importlib.util
from fractions import Fraction as Q
import provider as p
import audit
s=importlib.util.spec_from_file_location('q42_oracle',p.HERE.parent/'core_current_head_q29_weighted_compact_constants/compact.py')
q29=importlib.util.module_from_spec(s);s.loader.exec_module(q29)
def coeff(v):return {j:dict(c.terms) for j,c in v.c.items()}
def data():return dict(zip(sorted(p.KEYS),range(1,12)))
def raw(n,v):return dict(order=n,pole=v,cutoff_order=n+1,bounds=data())
class Tests(unittest.TestCase):
 def test_exact_Q29_controls(self):
  for n in (0,1,4):
   a=p.constants(n,data());b=q29.constants(n,**data())
   for k in b:self.assertEqual(coeff(a[k]),coeff(p.EP(b[k])))
 def test_polynomial_lower_data(self):
  a=data();a['f']=p.EP({0:2,1:3});b=p.constants(2,a)
  d=data();d['f']=1;c=q29.constants(2,**d)
  for k in c:self.assertEqual(coeff(b[k]),coeff(p.EP({0:2*c[k],1:3*c[k]})))
 def test_exact_pole_archive(self):
  calls=[]
  def source(n,v):
   calls.append((n,v));d=raw(n,v);d['bounds']['kernel']=Q(1,v+1);return d
  c=p.CompactProvider(source,lambda n,v:2,3,8)
  self.assertNotEqual(coeff(c(2,3)['KT']),coeff(c(2,7)['KT']))
  c(2,3);self.assertEqual(calls,[(2,3),(2,7)])
  with self.assertRaises(ValueError):c(4,3)
  with self.assertRaises(ValueError):c(2,9)
 def test_coverage_rejections(self):
  for change in (dict(cutoff_order=2),dict(order=1),dict(pole=5)):
   def bad(n,v):return dict(raw(n,v),**change)
   with self.assertRaises(ValueError):p.CompactProvider(bad,lambda n,v:1,3,8)(2,4)
  with self.assertRaises(ValueError):p.constants(0,{})
 def test_Q40_composition(self):
  from math import comb,factorial
  c=p.CompactProvider(raw,lambda n,v:3,3,8);calls=[]
  def rows(n):calls.append(n);return dict(Etheta=p.EP(2),Ez=p.EP(5))
  n=2;out=p.fast.fast_recipe(n,rows,c,8)
  a=q29.constants(n,**data());b=q29.constants(n+1,**data())
  t=comb(n+5,5)*24*2**(n+4)*3*factorial(n+4)
  r=comb(n+5,5)*24*2**(n+5)*3*factorial(n+5)*b['KT']
  for k,val in dict(u=2*t,z=5*t*(1+a['KA']),x=5*r).items():
   self.assertEqual(coeff(out[k]),coeff(p.EP({1:val})))
  self.assertEqual(calls,[6,7]);self.assertEqual(out['pole'],8)
 def test_evidence_and_scope(self):
  self.assertEqual(audit.generate(),json.loads((p.HERE/'evidence.json').read_text()))
  self.assertFalse(p.record()['initial_A1_tail_full_source_residual_or_blowup_certified'])
if __name__=='__main__':unittest.main()
