import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q97
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text();ns={'q87':type('Module',(),{'full':q97.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']
class Phase:
    def normal(self,j):return q97.EP(int(j==0))
    def matrix(self,j):return q97.EP(2*int(j==0))
    def moving_normal(self,j):return q97.EP(3*int(j==0))
class Tests(unittest.TestCase):
    def test_all_six_complex_covariance_entries(self):
        th=s.symbols('theta',real=True)
        A,B,C,D,E,F=map(s.Matrix,([1,2,3],[2,-1,1],[1,3,-2],[-2,1,4],[2,1,-3],[1,4,2]))
        old=A*s.cos(th)+B*s.sin(th);new=C*s.cos(th)+D*s.sin(th)+E*s.cos(2*th)+F*s.sin(2*th)
        exact=((old+new)*(old+new).T-old*old.T).applyfunc(lambda x:s.integrate(s.expand(x),(th,0,2*s.pi))/(2*s.pi))
        expected=(A*C.T+C*A.T+B*D.T+D*B.T+C*C.T+D*D.T+E*E.T+F*F.T)/2
        self.assertEqual(exact,expected)
        self.assertTrue(all(exact[i,j]!=0 for i in range(3) for j in range(i,3)))
        self.assertNotEqual(expected,expected-(B*D.T+D*B.T)/2)
        self.assertNotEqual(expected,expected-(E*E.T+F*F.T)/2)
    def test_global_radial_moment_and_pressure_identities(self):
        R,z=s.symbols('R z');a,b=1,3;h=(R-a)**3*(b-R)**3
        rr=h*(1+z);tt=h*(2+z*z);rz=h*z;zz=h*(3+z);tz=h*(1+z*z)
        gr=-s.diff(rr,R)-rr/R-s.diff(rz,z)+tt/R
        integ=lambda f:s.integrate(f,(R,a,b))
        self.assertEqual(s.simplify(integ(gr)-(integ((tt-rr)/R)-s.diff(integ(rz),z))),0)
        Jz=integ(R*zz)-integ(R*R*gr)/2
        reduced=integ(R*(zz-(rr+tt)/2))+s.diff(integ(R*R*rz),z)/2
        self.assertEqual(s.simplify(Jz-reduced),0)
        # rho=1/(b-a), exact unit mass; pressure has zero traces because gr-rho*P integrates to zero.
        P=integ(gr);pm=s.integrate(gr-P/(b-a),R);pm=pm-pm.subs(R,a)
        crho=integ(R*R/(b-a))/2
        self.assertEqual(s.simplify(pm.subs(R,b)),0)
        self.assertEqual(s.simplify(integ(R*pm)-(-integ(R*R*gr)/2+crho*P)),0)
    def test_raw_covariance_and_whole_moving_product_envelopes(self):
        v=q97.covariance_envelope([1,2,3],[4,5,6])
        self.assertEqual(v['change'][0].record(),q97.EP(24).record())
        self.assertEqual(v['change'][1].record(),q97.EP(66).record())
        out=q97.moment_envelope([1,2,3],1,2,1,[4,5])
        # P0=4,P1=7; Jz1=14; (cP)1=5*4+4*7=48, so whole axial0=62.
        self.assertEqual(out['P'][0].record(),q97.EP(4).record())
        self.assertEqual(out['whole_axial'][0].record(),q97.EP(62).record())
        zero=q97.moment_envelope([1,2,3],1,2,1,[4,0])
        self.assertNotEqual(out['whole_axial'][0].record(),zero['whole_axial'][0].record())
    def test_matches_frozen_zeroth_moment_envelope(self):
        import importlib.util
        spec=importlib.util.spec_from_file_location('q97_old_moments',p/'core_sourced_covariance_moment_join/budgets.py')
        old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
        for D0,D1,a,b,c in ((Q(1),Q(2),Q(1),Q(2),Q(4)),(Q(2,3),Q(3,2),Q(1,2),Q(3),Q(5,4))):
            native=q97.moment_envelope([D0,D1,1],a,b,1,[c,0])
            previous=old.moment_budget(0,D0,D1,a,b,c)
            for key in ('P','Jtheta','Jz','pressure_moment'):
                self.assertEqual(native[key][0].record(),q97.EP(previous[key]).record())
    def test_typed_both_carriers_common_pole_and_state(self):
        post=fixture(3);target=post.current.current.state.target;target.total=lambda j:q97.EP(j+1);EP=q97.EP
        q94=q97.response.inputs.source
        w=q94.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
        response=q97.response.SourcedResponse(q97.response.inputs.PathInputs(q94.CurrentSource(w)))
        vector=lambda n:{k:EP(1) for k in ('x','u','z')}
        post.current.old_full=vector;post.current.fast_unweighted=vector;post.current.slow=lambda n,b:vector(n)
        budget=post.current.initial.budget
        budget._data=lambda n,p:dict(b=EP(1),G=EP(1),c=EP(1),graph=EP(2),r=EP(3))
        budget.compact=lambda n,p:dict(KT=EP(1),KA=EP(1),I=EP(1))
        calls=[];original=response.build
        def traced(*a,**kw):calls.append(kw['harmonic']);return original(*a,**kw)
        response.build=traced
        out=q97.SourcedCovariance(response);table=[{k:EP(1) for k in ('Eik','bR','Bz','psi')} for j in range(2)]
        v=out.covariance(0,table,[0,1],[1,0],[1,1],tangent=EP(1),length=EP(1))
        self.assertEqual(calls,[1,2]);self.assertEqual(v['covariance_pole'],2*v['wave_pole'])
        self.assertEqual(v['change'][0].record(),(v['cross'][0]+v['square'][0]).record())
        self.assertIs(v['source_state'],out.state);self.assertFalse(v['means_updated'])
        self.assertEqual(v['stage'],'post_first_sourced_wave_fixed_means')
        with self.assertRaises(ValueError):out.moments(0,v)
        with self.assertRaises(TypeError):q97.SourcedCovariance(lambda n:1)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
