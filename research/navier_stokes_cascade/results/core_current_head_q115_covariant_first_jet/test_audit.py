import unittest,json
from fractions import Fraction as Q
import sympy as s
import audit,q115

def leaves(m=1):
 return dict(m=m,T=2,L=5,dmin=Q(1,4),normal=3,normal_prime=2,normal_parameter=4,
             normal_prime_parameter=5,inverse_square=2,K=3,K_parameter=4,
             rate_parameter=5,damping_parameter=2,source=7,source_parameter=11)

class Tests(unittest.TestCase):
 def test_exact_moving_constraint_and_normal_source_identity(self):
  t,a=s.symbols('t a',real=True);I=s.eye(3)
  n=s.Matrix([1,a+t,a*t]);np=n.diff(t);na=n.diff(a);N=(n.T*n)[0]
  Pi=I-n*n.T/N
  K=s.Matrix([[1,a,0],[t,2,1],[0,a*t,3]]);A0=-Pi*K-n*np.T/N
  d=2+a*a+t*t;p=a+t
  u=Pi*s.Matrix([a+t,1,t]);w=Pi*u.diff(a)
  for m in (1,3):
   B=A0-(m*m*d+p)*I
   for normal_force in (0,1):
    g=B*u-u.diff(t)+normal_force*n*(a+t+1)
    F=q115.tangent_forcing(n,np,na,K,K.diff(a),u,g,g.diff(a),d.diff(a),p.diff(a),m)
    residual=w.diff(t)-B*w-F
    J=((n.T*K-np.T)*u+n.T*g)[0]
    Ja=((na.T*K+n.T*K.diff(a)-np.diff(a).T)*u+(n.T*K-np.T)*u.diff(a)+na.T*g+n.T*g.diff(a))[0]
    pressure_residual=s.diff(J/N,a)-Ja/N+2*(n.T*na)[0]*J/N**2
    reconstruction=u.diff(a)-w+n/N*(na.T*u)[0]
    for point in ({t:0,a:0},{t:Q(1,2),a:Q(1,3)},{t:1,a:2}):
     self.assertEqual(residual.subs(point).applyfunc(s.simplify),s.zeros(3,1))
     self.assertEqual(s.simplify(pressure_residual.subs(point)),0)
     self.assertEqual(reconstruction.subs(point).applyfunc(s.simplify),s.zeros(3,1))
    # Purely normal forcing cancels from the velocity identity only when retained
    # in BOTH pressure numerator and derivative of projected forcing.
    if normal_force:
     omitted=F-(Pi*na/N)*(n.T*g)[0]
     self.assertNotEqual((F-omitted).subs({t:1,a:2}),s.zeros(3,1))

 def test_uniform_harmonic_residence_and_full_parameter_bound(self):
  for m in (1,-1,2,-2,3,10,1000):
   out=q115.first_jet(**leaves(m))
   self.assertLessEqual(out['m2_residence'],out['uniform_m2_residence'])
   self.assertEqual(out['full_parameter'],out['tangent_parameter']+out['normal_constraint'])
   self.assertGreater(out['pressure_parameter_times_carrier'],0)
   self.assertFalse(out['actual_source_or_high_jet_contracts_certified'])
  self.assertEqual(q115.residence(1,5,Q(1,4)),5)
  self.assertEqual(q115.residence(2,5,Q(1,4)),Q(4,3))

 def test_normal_direction_blocks_ambient_residence_replacement(self):
  # n=e1,K=0,p=0,d=1: ambient C=diag(0,-m^2,-m^2).
  # With constant normal forcing e1, the zero-entry ambient solution is t*e1.
  t=s.symbols('t',nonnegative=True)
  for m in (2,5,20):
   C=s.diag(0,-m*m,-m*m);y=s.Matrix([t,0,0]);f=s.Matrix([1,0,0])
   self.assertEqual(y.diff(t)-C*y-f,s.zeros(3,1))
   self.assertGreater(Q(5),q115.residence(m,5,1))

 def test_full_normal_pressure_parameter_is_not_velocity(self):
  # Constant n=e1, K=0, g=a*n: u=0 but (pi/P)_a=i/(k*m).
  a=s.symbols('a',real=True);m=3;k=7
  pressure=s.I*a/(k*m)
  self.assertEqual(s.diff(pressure,a),s.I/(k*m))
  v=leaves(m);v.update(normal=1,normal_prime=0,normal_parameter=0,
    normal_prime_parameter=0,inverse_square=1,K=0,K_parameter=0,
    rate_parameter=0,damping_parameter=0,source=0,source_parameter=1)
  out=q115.first_jet(**v)
  self.assertEqual(out['pressure_parameter_times_carrier'],Q(1,m))

 def test_mixed_moving_normal_pressure_term_survives(self):
  a,t=s.symbols('a t',real=True);n=s.Matrix([1,a*t,0]);N=(n.T*n)[0]
  u=(s.eye(3)-n*n.T/N)*s.Matrix([0,1,0])
  numerator=-(n.diff(t).T*u)[0]
  self.assertEqual(s.diff(numerator/N,a).subs({a:0,t:0}),-1)
  self.assertEqual(n.diff(t).subs({a:0,t:0}),s.zeros(3,1))
  self.assertEqual(n.diff(a).subs({a:0,t:0}),s.zeros(3,1))
  self.assertEqual(n.diff(t,a).subs({a:0,t:0}),s.Matrix([0,1,0]))

 def test_invalid_leaves(self):
  for key,val in (('m',0),('m',True),('T',0),('dmin',0),('inverse_square',0),('source',-1),('normal',1.0)):
   with self.assertRaises((TypeError,ValueError)):q115.first_jet(**dict(leaves(),**{key:val}))

 def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
