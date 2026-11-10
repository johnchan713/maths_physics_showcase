import json,tempfile,unittest
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import q141 as q
import audit

class Tests(unittest.TestCase):
 def fixture(self,T=2):
  seed=q.seed(2,2,1,1,2,2,q.IDENTITY);step=q.switch(Q(1,4),q.IDENTITY)
  traces=q.trace_rows(T,2);geo=q.geometry(2,2,2,2)
  return q.short_taper(seed,step,traces,geo,Q(1,100),Q(1,4),10,10,q.IDENTITY)
 def test_value_based_reference_and_actual_recurrence(self):
  r=self.fixture()
  self.assertFalse(r['reference']['order_dependent_norm_exponentiated'])
  self.assertFalse(r['actual']['order_dependent_norm_exponentiated'])
  self.assertTrue(r['sources']['all_five_incoming_traces_retained'])
  self.assertTrue(r['state']['full_p2_including_Pi_and_Pi_eta'])
  for f in ('a','b','p','E'):
   self.assertGreater(r['state'][f]['P'],0)
 def test_exponential_composition_at_value(self):
  import sympy as s
  x,z=s.symbols('x z');L=(x+z+x*z)/8
  rec=q.exp_comp(q.env(1,2,2),3)
  for i in range(5):
   for k in range(5):
    raw=abs(s.diff(s.exp(L),x,i,z,k).subs({x:0,z:0}))
    cap=q.bound(rec,i,k)['upper']
    self.assertLessEqual(raw,s.Rational(cap.numerator,cap.denominator))
 def test_independent_memory_and_pressure_response(self):
  a=self.fixture(2);b=self.fixture(3)
  self.assertGreater(b['sources']['p1']['P'],a['sources']['p1']['P'])
  self.assertGreater(b['sources']['Xns']['P'],a['sources']['Xns']['P'])
  self.assertGreater(b['state']['p']['P'],a['state']['p']['P'])
  # A pressure trace cannot be replaced by local velocity data.
  import sympy as s
  eta,X,E,h=s.symbols('eta X E h',positive=True);c=eta
  dp=X*(4*(s.Rational(1,2)+h)*eta*c-(1-eta**2)*s.diff(c,eta))/((1-2*h*eta**2)*E)
  self.assertEqual(s.simplify(dp.subs(eta,0)),-X/E)
 def test_old_cap_removed_and_fixed_binary_growth(self):
  r=self.fixture();state=r['state'];common=q.env(state['H'],state['Ay'],state['Be'])
  old=q.previous.q72_log2_lower(Q(1,100),Q(1,100),8)
  new=q.bound(common,8,8,True)['upper']
  self.assertGreater(old,new)
  self.assertIsInstance(q.bound(common,100005,100005,True)['upper'],int)
 def test_Q140_Q139_Q138_binding_and_guards(self):
  r=self.fixture();p=q.previous.previous
  disk=p.phase_disk(p.scalar.parameter_disk(2,8,p.scalar.IDENTITY),p.IDENTITY)
  joined=q.bind(r,1,1,2,8,1,disk,1,8,q.IDENTITY)
  self.assertEqual(joined['joined']['pullback']['N'],8)
  self.assertFalse(joined['actual_order_zero_caps_floors_frequency_and_phase_gates_certified'])
  with self.assertRaises(ValueError):q.short_taper({},q.env(1),{}, {},Q(1,100),Q(1,4),1,1,'different-map')
  with self.assertRaises(ValueError):q.reciprocal(q.env(1),0)
  with self.assertRaises(ValueError):q.bound(q.env(1),100005,0)
 def test_00_record(self):
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'replay.json';frozen=audit.HERE/'evidence.json'
   out.write_text(json.dumps(audit.generate(),indent=2,sort_keys=True)+'\n')
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
