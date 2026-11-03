import unittest,json
from fractions import Fraction as Q
from pathlib import Path
from math import comb
import audit,q117
E=q117.EP
text=(audit.HERE.parent/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text()
ns=dict(q87=type('Module',(),{'full':q117.native.full}),Q=Q)
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']

class Phase:
 def normal(self,j):return E(int(j==0))
 def matrix(self,j):return E(2*int(j==0))
 def moving_normal(self,j):return E(3*int(j==0))

def setup(maximum=4):
 post=fixture(maximum);target=post.current.current.state.target;target.total=lambda j:E(j+1)
 src=q117.native.inputs.source
 w=src.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
 out=q117.CovariantResponse(q117.native.inputs.PathInputs(src.CurrentSource(w)))
 vector=lambda n:{k:E(1) for k in ('x','u','z')}
 post.current.old_full=vector;post.current.fast_unweighted=vector;post.current.slow=lambda n,b:vector(n)
 budget=post.current.initial.budget
 budget._data=lambda n,p:dict(b=E(1),G=E(1),c=E(1),graph=E(2),r=E(3))
 budget.compact=lambda n,p:dict(KT=E(1),KA=E(1),I=E(1))
 return out

def tables(N):return [{k:E(1) for k in ('Eik','bR','Bz','psi')} for _ in range(N+1)]

def constant_value(x):
 assert set(x.c)<= {0}
 if not x.c:return Q(0)
 terms=x.c[0].terms;assert set(terms)<= {Q(0)}
 return terms.get(Q(0),Q(0))

class Tests(unittest.TestCase):
 def test_EP_lift_matches_frozen_rational_recurrence(self):
  for N in (1,2,3):
   G=[Q(j+1) for j in range(N+1)];n=[Q(2)]*(N+1);np=[Q(1)]*(N+1);K=[Q(3)]*(N+1);rate=[Q(1)]*(N+1)
   d=[4*x for x in q117.mixed.convolution(n,n)]
   rational=q117.mixed.coefficients(G,n,np,K,rate,d,q0=1,T=2,L=3,dmin=Q(1,2),entrance_flat_through=N-1)
   for m in (1,-1,2,5):
    out=q117.response(G,n,np,K,rate,1,1,3,damping_floor=Q(1,2),entrance_flat_through=N-1,harmonic=m)
    old=q117.mixed.specialize(rational,m)
    self.assertEqual([constant_value(x) for x in out['velocity_over_residence']],rational['velocity_over_residence'])
    for key,oldkey in (('velocity','velocity'),('pressure','pressure_times_carrier'),('potential','potential_times_carrier')):
     self.assertEqual([constant_value(x) for x in out[key]],old[oldkey])

 def test_source_linearity_and_polynomial_length(self):
  G=[E({0:1,1:1}),E({2:1}),E(3)];H=[E(2),E({1:2}),E(1)]
  n=[1,1,0];np=[1,0,1];K=[2,1,0];rate=[1,0,1]
  def run(G,m=2):return q117.response(G,n,np,K,rate,1,E({1:1}),E({0:1,1:1}),damping_floor=Q(1,4),entrance_flat_through=1,harmonic=m)
  a,b,c=run(G),run(H),run([x+y for x,y in zip(G,H)])
  for key in ('velocity','pressure','potential','velocity_over_residence'):
   for x,y,z in zip(a[key],b[key],c[key]):self.assertEqual((x+y).record(),z.record())
  self.assertEqual(a['new_source_pole_added'],0)
  self.assertEqual(a['residence_upper'].record(),E(Q(4,3)).record())
  self.assertEqual(run(G,1)['residence_upper'].record(),E({0:1,1:1}).record())

 def test_current_source_full_restoration_and_pressure_pole(self):
  out=setup();kw=dict(tangent=E(1),length=E(1),damping_floor=Q(1,4),entrance_flat_through=1)
  v=out.build(1,tables(2),[0,1,2],[1,0,0],[1,1,0],log_poles=[0,2,3],**kw)
  coeff=v['normalized_response'];ev=v['evaluation_factor']
  expected=ev*(coeff['pressure'][1]+coeff['pressure'][0])+coeff['pressure'][0]
  self.assertEqual(v['pressure'][1].record(),expected.record())
  self.assertEqual(v['evaluation_factor'].record(),E({0:1,1:376}).record())
  p0=coeff['potential'][0];p1=ev*(coeff['potential'][1]+p0)+p0
  self.assertEqual(v['curl'][0].record(),(6*p1+3*p0).record())
  self.assertEqual(v['complete'][0].record(),(v['principal'][0]+v['curl'][0]).record())
  data=out.inputs.build(2,tables(2),[0,1,2],[0,2,3])
  bell=q117.native.inputs.quotient([0]*3,[0,1,2],[0,2,3])
  self.assertEqual(v['pole'],data['common_pole']+bell['extra_pole'])
  self.assertEqual((v['source_ceiling'],v['source_factory_reserve']),(2,4))
  self.assertIs(v['source_state'],out.state);self.assertFalse(v['actual_contracts_certified'])
  changed=out.build(1,tables(2),[0,1,2],[1,0,0],[1,0,0],log_poles=[0,2,3],**kw)
  self.assertEqual(v['principal'][0].record(),changed['principal'][0].record())
  self.assertNotEqual(v['curl'][0].record(),changed['curl'][0].record())

 def test_same_state_two_carrier_covariance_uses_new_response(self):
  out=setup();calls=[];original=out.build
  def traced(*a,**kw):calls.append(kw['harmonic']);return original(*a,**kw)
  out.build=traced
  factory=q117.covariance.SourcedCovariance(out)
  v=factory.covariance(0,tables(1),[0,1],[1,0],[1,1],tangent=E(1),length=E(1),damping_floor=1,entrance_flat_through=0)
  self.assertEqual(calls,[1,2]);self.assertIs(v['source_state'],out.state)
  self.assertEqual(v['covariance_pole'],2*v['wave_pole'])
  self.assertEqual(v['change'][0].record(),(v['cross'][0]+v['square'][0]).record())
  self.assertFalse(v['actual_native_contracts_certified']);self.assertFalse(v['means_updated'])
  with self.assertRaises(ValueError):factory.moments(0,v)

 def test_missing_contracts_fail_before_source_evaluation(self):
  out=setup();out.inputs.build=lambda *a,**kw:(_ for _ in ()).throw(AssertionError('Source must not be evaluated'))
  for change in (dict(entrance_flat_through=-1),dict(entrance_flat_through=True),dict(damping_floor=0),dict(damping_floor=1.0),dict(damping_floor=2)):
   kw=dict(tangent=1,length=1,damping_floor=1,entrance_flat_through=0);kw.update(change)
   with self.assertRaises(ValueError):out.build(0,tables(1),[0,1],[1,0],[1,0],**kw)
  out.state=object()
  with self.assertRaises(ValueError):out.build(0,tables(1),[0,1],[1,0],[1,0],tangent=1,length=1,damping_floor=1,entrance_flat_through=0)
  with self.assertRaises(TypeError):q117.CovariantResponse(None)

 def test_same_normal_floor_is_derived_without_new_oracle(self):
  self.assertEqual(q117.damping_lower(Q(3,2)),Q(2,3))
  self.assertEqual(q117.damping_lower(Q(3,2),Q(1,3)),Q(1,3))
  with self.assertRaises(ValueError):q117.damping_lower(Q(3,2),1)
  v=q117.response([1,0],[1,0],[0,0],[0,0],[0,0],2,0,1,entrance_flat_through=0,harmonic=2)
  self.assertEqual(v['normal_reciprocal_derived_floor'],Q(1,2))
  self.assertEqual(v['damping_floor'],Q(1,2))

 def test_full_normal_pressure_and_invalid_inputs(self):
  v=q117.response([1,0],[1,0],[0,0],[0,0],[0,0],1,0,1,damping_floor=1,entrance_flat_through=0,harmonic=2)
  self.assertGreater(constant_value(v['pressure'][0]),0)
  self.assertFalse(v['actual_native_and_entrance_contracts_certified'])
  for change in (dict(harmonic=0),dict(length=0),dict(q0=0),dict(damping_floor=0),dict(entrance_flat_through=-1)):
   kw=dict(q0=1,tangent=0,length=1,damping_floor=1,entrance_flat_through=0);kw.update(change)
   with self.assertRaises(ValueError):q117.response([1,0],[1,0],[0,0],[0,0],[0,0],**kw)

 def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
