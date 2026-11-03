import unittest,json
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import audit,q119

# Independent normalized series algebra, not the closed chart formula.
def mul(a,b):
 return [sum((a[i]*b[j-i] for i in range(j+1)),Q(0)) for j in range(len(a))]
def expjet(a):
 assert a[0]==0
 e=[Q(1)]
 for j in range(1,len(a)):e.append(sum((i*a[i]*e[j-i] for i in range(1,j+1)),Q(0))/j)
 return e

def full_series(N,a,b):
 Y=[Q(0)]+[sum((Q(32**i*64**(j-i),3) for i in range(j+1)),Q(0)) for j in range(1,N+1)]
 e=expjet(Y);eh=expjet([v/2 for v in Y]);r=[Q(b*b)]+[Q(2*b),Q(1)]+[Q(0)]*max(0,N-2);r=r[:N+1]
 X=[2*v for v in mul(r,e)];X[0]-=2*b*b
 eta=[eh[j]+(4*eh[j-1] if j else 0) for j in range(N+1)];eta[0]-=1
 comp=expjet([x+y+z for x,y,z in zip(Y,X,eta)])
 inv=[Q(a)**(-j-1) for j in range(N+1)]
 return dict(F=[4*x for x in mul(mul(e,inv),comp)],G=[4*x for x in mul(e,comp)],btilde=[2*x for x in mul(eh,comp)])

def budget(cls,trace=None):
 return cls(Q(1,200),leading=lambda m:Q(m+2),profiles=lambda n,m:Q(n+m+1),
            cutoff=lambda m:Q(m+1),Rmin=1,Rmax=8,derivative_basis='euler',trace=trace)

class Tests(unittest.TestCase):
 def test_exact_radius_and_composition_inequalities(self):
  for a,b in ((Q(1,4),3),(1,8),(2,3),(Q(1,2),Q(3,4))):
   r=q119.chart_recipe(0,a,b);t=r['tau']
   Y=(1/((1-32*t)*(1-64*t))-1)/3
   X=144*b*b*t+8*b*t+4*t*t;eta=44*t
   self.assertLessEqual(t,Q(1,1024));self.assertLessEqual(t,Q(a)/2)
   self.assertLessEqual(Y,36*t);self.assertLess(Y,Q(1,16))
   self.assertLess(X,Q(1,4));self.assertLess(eta,Q(1,16))
   self.assertLess(Y+X+eta,1)

 def test_closed_bound_dominates_independent_full_series(self):
  for a,b in ((Q(1,4),Q(3)),(Q(1),Q(8)),(Q(2),Q(3)),(Q(1,2),Q(3,4))):
   table=full_series(9,a,b)
   for m in range(10):
    closed=q119.materialize(q119.chart_recipe(m,a,b))
    for key in ('F','G','btilde'):
     raw=max(factorial(j)*table[key][j] for j in range(m+1))
     self.assertLessEqual(raw,closed[key])

 def test_original_value_comparison_and_low_constants_unchanged(self):
  old=budget(q119.base.BackgroundBudget);new=budget(q119.ClosedBackground)
  for m in range(3):
   self.assertEqual(old.bound(m),new.bound(m))
   self.assertEqual(old.leading_bound(m),new.leading_bound(m))
  args=(Q(1,2),Q(1),Q(1),Q(1),Q(2),Q(1,10))
  self.assertEqual(old.comparison(*args),new.comparison(*args))

 def test_full_tuple_head_tail_and_cutoff_order_retained(self):
  trace=set();new=budget(q119.ClosedBackground,trace);old=budget(q119.base.BackgroundBudget)
  v=new.bound(4)
  self.assertEqual(v['head'],old.bound(4)['head'])
  self.assertIn(('profile',1,4),trace);self.assertIn(('profile',3,4),trace)
  self.assertIn(('cutoff',5),trace);self.assertEqual(v['cutoff_input_order'],5)
  self.assertEqual(v['chart']['kind'],'closed-real-slow-chart-majorant')
  self.assertEqual(v['F'],v['chart']['F']*(Q(6)+v['head']['correction']))
  self.assertFalse(q119.chart_recipe(4,1,8)['actual_profile_inputs_certified'])

 def test_same_current_phase_response_accepts_closed_background(self):
  phase=q119.phase;E=phase.EP
  text=(audit.HERE.parent/'core_current_head_q118_background_phase_binding/test_audit.py').read_text()
  ns=dict(q118=phase,audit=audit,Path=Path,Q=Q,factorial=factorial)
  exec(text[text.index('E=q118.EP'):text.index('\nclass Tests(')],ns)
  old,old_response,_=ns['setup']()
  bg=budget(q119.ClosedBackground)
  native=phase.BackgroundPhase(old.postslow,bg,maximum=3,path=old.recipe.path,Lmin=1,U=Q(1,2),lambda_max=1,gmin=1)
  waves=native.waves(regular_wave=native.state.target.root.wave,inverse_normal_square=1)
  response=phase.BackgroundResponse(phase.native.inputs.PathInputs(phase.native.inputs.source.CurrentSource(waves)))
  out=response.build(0,None,[0,1],None,[1,1],psi_table=[1,1],tangent=1,length=1,entrance_flat_through=0)
  self.assertIs(out['source_state'],old_response.state)
  self.assertEqual(out['background_phase_request']['background_ceiling'],4)
  self.assertEqual(out['complete'][0].record(),(out['principal'][0]+out['curl'][0]).record())
  self.assertFalse(out['phase_provider_inequalities_certified'])
  self.assertEqual(bg.bound(2)['chart'],q119.base.chart_constants(2,1,8))

 def test_deferred_factorials_and_invalid_geometry(self):
  r=q119.chart_recipe(100009,1,8)
  self.assertEqual(r['factorial_argument'],100009)
  with self.assertRaises(ValueError):q119.materialize(r)
  for args in ((-1,1,2),(True,1,2),(0,0,2),(0,3,2)):
   with self.assertRaises(ValueError):q119.chart_recipe(*args)
  with self.assertRaises(ValueError):q119.ClosedBackground(Q(1,200),lambda m:1,lambda n,m:1,lambda m:1,1,8,'frozen')

 def test_record(self):
  # Keep recursive historical imports isolated from the following current graph.
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   output=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(output)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(output.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
