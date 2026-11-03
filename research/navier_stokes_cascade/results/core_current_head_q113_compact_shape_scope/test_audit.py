import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q113

text=(audit.HERE.parent/'core_current_head_q110_radial_tail_request/test_audit.py').read_text()
ns=dict(q110=q113.profiles.charts.radial,audit=audit,Path=Path,Q=Q)
exec(text[text.index('def setup('):text.index('\nclass Tests(')],ns);setup=ns['setup']

def scope(maximum=1):
    obj,bundle=setup(maximum)
    return obj,bundle,q113.CompactShapeScope(q113.profiles.KernelProfiles(q113.profiles.charts.ChartBridge(obj)))

class Tests(unittest.TestCase):
    def test_smooth_comparison_has_unbounded_midpoint_derivative(self):
        # An explicit smooth comparison step, not the selected compact cutoff.
        t=s.symbols('t',positive=True)
        g=s.exp(-1/t**2);h=s.exp(-1/(1-t)**2);step=g/(g+h)
        self.assertEqual(step.subs(t,s.Rational(1,2)),s.Rational(1,2))
        self.assertEqual(s.simplify(s.diff(step,t).subs(t,s.Rational(1,2))),8)
        previous=0
        for index in (1,2,7,100):
            v=q113.transition_comparison(Q(2,3),Q(11,3),index)
            w=v['half_width'];m=v['center']
            self.assertGreater(m-w,Q(2,3));self.assertLess(m+w,Q(11,3))
            self.assertEqual(v['derivative_lower_bound'],Q(2*index,3))
            self.assertGreaterEqual(4/w,v['derivative_lower_bound'])
            self.assertGreater(v['derivative_lower_bound'],previous)
            self.assertFalse(v['chosen_actual_shape']);previous=v['derivative_lower_bound']

    def test_independent_mass_and_scaling_diagnostic(self):
        # This polynomial control is only C^2 at the seams, not the C-infinity witness.
        t=s.symbols('t');step=6*t**5-15*t**4+10*t**3
        derivative=s.diff(step,t)
        self.assertEqual(s.factor(derivative),30*t**2*(t-1)**2)
        self.assertEqual(s.integrate(derivative,(t,0,1)),1)
        for index in (1,3,11):
            v=q113.transition_comparison(1,3,index);w=v['half_width']
            self.assertEqual(2*w*v['derivative_lower_bound'],1)
            self.assertGreater(derivative.subs(t,s.Rational(1,2))/(2*w),v['derivative_lower_bound'])

    def test_cutoff_roles_are_distinct_and_substitution_is_rejected(self):
        self.assertEqual(len({r['argument'] for r in q113.ROLES.values()}),3)
        for record in q113.ROLES.values():self.assertTrue((audit.HERE.parent/record['source']).is_file())
        compact=q113.require_compact('compact-X')
        compact['argument']='changed'
        self.assertEqual(q113.require_compact('compact-X')['argument'],'X=R^2/(2s)')
        for role in ('borel','variance','unknown'):
            with self.assertRaises(ValueError):q113.require_compact(role)

    def test_current_finite_order_scope_preserves_large_request(self):
        obj,bundle,proof=scope();out=proof.plan(0)
        self.assertIs(out['source_state'],obj.state)
        self.assertTrue(out['qualitative_fixed_order_band_uniformity'])
        self.assertTrue(out['conditional_on_original_fixed_Cinfinity_choice'])
        for key in ('uniform_over_arbitrary_cutoff_choices','numerical_shape_value_supplied',
                    'numerical_shape_bound_verified','actual_source_or_A0_gain_certified'):
            self.assertFalse(out[key])
        self.assertEqual(out['shape_order'],1)
        self.assertEqual(out['profile_majorants']['source_f0']['factorial_argument'],100005)
        self.assertEqual(out['profile_majorants']['source_factorial_multiplicity_including_f0'],3)
        with self.assertRaises(ValueError):proof.recipe(0,bundle,geometry=dict(width=Q(1,4),gap=Q(1,8),L=2))
        obj.state=object()
        with self.assertRaises(ValueError):proof.plan(0)

    def test_current_deferred_recipe_keeps_all_native_guards(self):
        obj,bundle,proof=scope(9);geometry=dict(width=Q(1,4),gap=Q(1,8),L=2)
        out=proof.recipe(0,bundle,geometry=geometry,kappa=Q(1,5))
        self.assertEqual(out['scope']['profile_majorants']['source_f0']['factorial_argument'],9)
        self.assertIs(out['source_state'],obj.state)
        self.assertTrue(out['actual_high_jet_source_still_required'])
        self.assertFalse(out['full_cone_or_residual_certified'])
        for role in ('borel','variance'):
            with self.assertRaises(ValueError):proof.recipe(0,bundle,geometry=geometry,kappa=Q(1,5),role=role)
        bad=dict(bundle,source_state=object())
        with self.assertRaises(ValueError):proof.recipe(0,bad,geometry=geometry,kappa=Q(1,5))

    def test_invalid_comparison_and_provider(self):
        for args in ((1,1,1),(2,1,1),(1,2,0),(1,2,-1)):
            with self.assertRaises(ValueError):q113.transition_comparison(*args)
        with self.assertRaises(TypeError):q113.CompactShapeScope(None)

    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
