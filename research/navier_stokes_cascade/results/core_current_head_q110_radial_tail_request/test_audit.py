import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q110

def setup(maximum=1):
    a=q110.axial.angular;p=audit.HERE.parent/'core_current_head_q87_final_angular_adapter/test_audit.py'
    text=p.read_text();ns=dict(q87=type('Module',(),{'full':q110.full}),Q=Q)
    exec(text[text.index('def fixture('):text.index('class Tests')],ns)
    post=ns['fixture'](maximum+12);E=q110.EP;r=q110.radial
    q103=r.post;q102=q103.slow;q101=q102.rows;q100=q101.changes;q99=q100.scalars;q98=q100.pressure;q97=q98.covariance;q94=q97.response.inputs.source
    class Phase:
        def normal(self,j):return E(int(j==0))
        def matrix(self,j):return E(2*int(j==0))
        def moving_normal(self,j):return E(3*int(j==0))
    target=post.current.current.state.target
    wave=q94.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
    cov=q97.SourcedCovariance(q97.response.SourcedResponse(q97.response.inputs.PathInputs(q94.CurrentSource(wave))))
    pole=post.current.fast.rows.pole+2
    pm=q98.ShiftedPressure(cov,c0=2,c1=4,shape=[1]*(maximum+12),time_inverse=lambda n,p:1,maximum=maximum+10,pole=pole)
    total=q101.FinalRows(q100.SourcedRowChanges(q99.ChangedScalars(pm)),maximum=maximum+9,pole=pole,c0=2,c1=4,shape=[1]*(maximum+12),time_inverse=lambda n,p:1)
    fresh=q102.PostfastSlow(total);reset=q103.CurrentPostslow(fresh)
    rad=r.FinalRadial(reset,maximum=maximum,pole=pole,c0=2,c1=4,shape=[1]*(maximum+2),time_inverse=lambda n,p:1)
    projected=a.projected.ProjectedRows(rad)
    obj=q110.RadialAudit(q110.axial.FullAxial(a.FullAngular(projected)))
    bundle=dict(source_state=obj.state,stage='post_first_sourced_wave_fixed_means',means_updated=False,order=maximum+11,covariance_pole=pole,old=[E(1)]*(maximum+12),new=[E(2)]*(maximum+12),change=[E(1)]*(maximum+12))
    return obj,bundle

class Tests(unittest.TestCase):
    def test_exact_radial_split_and_haar_compensation(self):
        R,r,y,Z=s.symbols('R r y Z',real=True);t=R-1;h=t**3*(1-t)**3
        chi=3*t*t-2*t**3+Z*t*t*(1-t)**2;rho=30*t*t*(1-t)**2
        avg=lambda f:s.integrate(s.expand(f),(y,-s.pi,s.pi))/(2*s.pi)
        whole=lambda f:s.integrate(s.expand(f.subs({R:r,y:y+2*(r-R)},simultaneous=True)),(r,1,2))
        A=lambda f:s.diff(chi,R)*whole(f)
        r0=h*(Z+s.cos(y-2*R));P=whole(avg(r0));f=r0-rho*P
        self.assertEqual(s.simplify(avg(A(f))),0)
        self.assertEqual(s.simplify(A(f)-A(r0-avg(r0))),0)
        Er=-rho*P-A(f)
        self.assertEqual(s.simplify(avg(Er)+rho*P),0)
        self.assertNotEqual(s.simplify(A(f).subs(Z,0)),0)
        K=h*(1+s.cos(y-2*R));Dr=lambda g:s.diff(g,R)+2*s.diff(g,y)
        self.assertEqual(s.simplify(A(Dr(K))),0)
        # Removing Dr K under A0 does not remove K from the pressure.
        left=lambda g:s.integrate(s.expand(g.subs({R:r,y:y+2*(r-R)},simultaneous=True)),(r,1,R))
        self.assertEqual(s.simplify(left(Dr(K))-chi*whole(Dr(K))-K),0)
    def test_least_strict_order_and_complete_reserves(self):
        r=q110.request(0)
        self.assertEqual((r['integration_order'],r['raw_source_order']),(50001,100005))
        self.assertEqual((r['mean_fields'],r['rows'],r['covariance']),(100007,100014,100016))
        self.assertEqual((r['mean_factory'],r['target_jets'],r['profile_jets']),(100019,100033,100046))
        self.assertEqual(r['strict_margin'],Q(1,100000));self.assertTrue(r['compensator_sufficient'])
        for k in (Q(1,100000),Q(1,5),Q(249,1000)):
            a=q110.request(2,kappa=k);m=a['integration_order']
            self.assertGreater(a['tail_rate'],a['target_rate'])
            if m>1:self.assertLessEqual(a['source_rate']+(m-1)*k,a['target_rate'])
        k=Q(1,100000);self.assertFalse(q110.request(0,target=3-6*k)['compensator_sufficient'])
        for k in (0,Q(1,4)):
            with self.assertRaises(ValueError):q110.request(0,kappa=k)
    def test_current_source_and_low_archive_gate(self):
        obj,b=setup();E=q110.EP;p=b['covariance_pole'];rad=obj.radial
        calls=[];edge=q110.radial.edge;q110.radial.edge=lambda p,w:E(1)
        rad.fields=lambda n,b:(calls.append(('final_fields',n)) or dict(fields=dict(x=E(1),u=E(2),z=E(3)),pole=p,parts=['all six'],order=n))
        rad.current.initial.budget._data=lambda n,p:{k:E(1) for k in ('b','V','G','c','graph')}
        try:
            out=obj.source(0,b)
            self.assertEqual(out['total_wave_covariance'].record(),E(549).record())
            self.assertEqual(out['coefficient'].record(),E(8489).record())
            self.assertEqual(calls,[('final_fields',2)])
            before=list(calls)
            with self.assertRaises(ValueError):obj.recipe(0,b,geometry=dict(width=Q(1,4),gap=Q(1,8),L=2))
            self.assertEqual(calls,before) # fail before any giant source/norm evaluation
            with self.assertRaises(ValueError):obj.source(0,dict(b,means_updated=True))
        finally:q110.radial.edge=edge
        with self.assertRaises(TypeError):q110.RadialAudit(None)
    def test_moderate_deferred_recipe_keeps_original_raw_profiles(self):
        obj,b=setup(9);geometry=dict(width=Q(1,4),gap=Q(1,8),L=2)
        out=obj.recipe(0,b,geometry=geometry,kappa=Q(1,5))
        self.assertEqual((out['integration_order'],out['raw_source_order']),(3,9))
        self.assertEqual(out['source_factorial'],9)
        self.assertEqual(out['operator'],'A0');self.assertIn('f_0',out['factors'][2])
        self.assertFalse(out['materialized_actual_input_bound']);self.assertFalse(out['radial_cone_closed'])
        self.assertTrue(out['whole_graph_pole_proof_still_required'])
        self.assertTrue(out['original_source_and_output_chart_transfers_required'])
        self.assertFalse(out['chart_factors_materialized'])
        with self.assertRaises(ValueError):obj.recipe(0,dict(b,old=b['old'][:-1]),geometry=geometry,kappa=Q(1,5))
        with self.assertRaises(ValueError):obj.recipe(0,dict(b,order=19),geometry=geometry,kappa=Q(1,5))
        with self.assertRaises(ValueError):obj.recipe(0,dict(b,covariance_pole=obj.radial.pole+1),geometry=geometry,kappa=Q(1,5))
        with self.assertRaises(ValueError):obj.recipe(0,b,geometry=dict(geometry,width=Q(1,8)),kappa=Q(1,5))
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
