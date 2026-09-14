#!/usr/bin/env python3
"""Compact mixed calculus, source equations, actual constants and scope regressions."""
from copy import deepcopy
from fractions import Fraction as F
import json
import unittest
import mpmath as mp
from jets import FirstJet, SecondJet, pressure_from_moments, shear_from_fields
from bounds import HERE, scalar_bounds, margin_bounds, exponent_ledger, sum_exp, reciprocal_first
from review import manufactured, source_and_jet_review, cutoff_review, shear_sign_review
from audit import (STATUS,FLAGS,RESOLVED,INHERITED_RESOLVED,REMAINING,validate,
                   validate_protocol,inherited_provenance)


class MixedJetTests(unittest.TestCase):
    def test_mixed_product_has_both_cross_terms(self):
        a,b = SecondJet(F(2),F(3),F(5)),SecondJet(F(7),F(11),F(13))
        self.assertEqual((a*b).ye,3*13+5*11)

    def test_raw_second_not_taylor_coefficient(self):
        y = SecondJet(F(1),y=F(1))
        self.assertEqual((y*y).yy,2)

    def test_exact_inverse_identity_all_six_components(self):
        f = SecondJet(*(map(F,(3,2,5,7,11,13))))
        self.assertEqual(f*f.inverse(),SecondJet(F(1)))

    def test_integer_division_does_not_insert_binary_float(self):
        for cls in (FirstJet,SecondJet):
            self.assertEqual(cls(1)/5,cls(F(1,5)))

    def test_angular_derivative_keeps_mixed_and_second(self):
        f = SecondJet(*(map(F,(3,2,5,7,11,13))))
        self.assertEqual(f.partial_first('eta'),FirstJet(F(5),F(11),F(13)))

    def test_radial_derivative_keeps_mixed_and_second(self):
        f = SecondJet(*(map(F,(3,2,5,7,11,13))))
        self.assertEqual(f.partial_first('y'),FirstJet(F(2),F(7),F(11)))

    def test_zero_inverse_rejected(self):
        for cls in (FirstJet,SecondJet):
            with self.assertRaises(ZeroDivisionError):
                cls(F(0)).inverse()

    def test_log_requires_positive_value(self):
        for v in (0,-1):
            with self.assertRaises(ValueError):
                SecondJet(F(v)).log()

    def test_mixed_order_cannot_be_promoted_or_silently_truncated(self):
        with self.assertRaises(TypeError):
            SecondJet(1)+FirstJet(2)
        with self.assertRaises(TypeError):
            FirstJet(1)+SecondJet(2)

    def test_derivative_axis_is_explicit(self):
        with self.assertRaises(ValueError):
            SecondJet(1).partial_first('X')

    def test_shear_uses_logarithmic_radial_derivative(self):
        with mp.workdps(70):
            y,eta = SecondJet(mp.mpf('.2'),y=1),SecondJet(mp.mpf('.3'),e=1)
            E = (-y/5).exp()/(1+eta*eta)
            U = (eta+y*eta*eta)/10
            result = shear_from_fields(E,U)
            self.assertLess(abs(result['a'].v-mp.mpf('1.4')),mp.mpf('1e-60'))
            self.assertLess(abs(result['ts'].v+result['bs'].v/result['a'].v),mp.mpf('1e-60'))


class ReviewTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = source_and_jet_review()
        cls.cutoff = cutoff_review()

    def test_direct_mixed_differentiation(self):
        self.assertTrue(self.source['checks']['mixed_jets_agree_with_direct_differentiation'])

    def test_five_original_radial_primitives(self):
        self.assertTrue(self.source['checks']['all_five_original_radial_identities'])

    def test_original_pressure_source_equations(self):
        self.assertTrue(self.source['checks']['original_source_ODEs_agree'])

    def test_pressure_jets_against_scalar_differentiation(self):
        self.assertTrue(self.source['checks']['pressure_jets_agree_with_scalar_differentiation'])

    def test_first_moment_jets_do_not_suffice(self):
        self.assertTrue(self.source['checks']['missing_curvature_preserves_values'])
        self.assertTrue(self.source['checks']['missing_angular_curvature_detected'])

    def test_axis_datum_cannot_be_dropped(self):
        self.assertTrue(self.source['checks']['changed_axis_datum_detected'])

    def test_pressure_force_cannot_be_dropped(self):
        self.assertTrue(self.source['checks']['dropped_pressure_force_detected'])

    def test_angular_bounds_do_not_control_narrow_radial_jets(self):
        self.assertTrue(self.cutoff['checks']['equal_values_and_angular_derivatives'])
        self.assertTrue(self.cutoff['checks']['second_radial_derivative_has_width_squared_loss'])

    def test_mixed_exponential_product_is_required(self):
        self.assertTrue(self.cutoff['checks']['mixed_and_radial_exponential_rules'])
        self.assertTrue(self.cutoff['checks']['omitted_exponential_product_detected'])

    def test_source_shear_sign(self):
        self.assertTrue(all(shear_sign_review()['checks'].values()))


class BoundTests(unittest.TestCase):
    def test_regional_constants_at_both_precisions(self):
        for digits in (80,110):
            result = scalar_bounds(digits)
            self.assertTrue(all(result['checks'].values()),result['checks'])
            self.assertFalse(result['raw_C_or_activation_width_materialized'])

    def test_original_positive_margins_at_both_precisions(self):
        for digits in (80,110):
            self.assertTrue(all(margin_bounds(digits)['checks'].values()))

    def test_axial_positive_D_term_is_not_lower_bounded_by_its_upper_bound(self):
        # The exact first factor is 1-(2+D)/p + D/2. Dropping its positive
        # term permits a lower bound using Dmax; ADDING Dmax/2 does not.
        p,Dmax = F(40),F(1,10)
        safe = 1-(2+Dmax)/p
        wrong = safe+Dmax/2
        actual_at_D_zero = 1-2/p
        self.assertLess(safe,actual_at_D_zero)
        self.assertGreater(wrong,actual_at_D_zero)

    def test_exponent_ledger_keeps_denominators(self):
        result = exponent_ledger()
        self.assertEqual(result['exponents']['p2_J1'],414)
        self.assertEqual(result['exponents']['vs_J1'],1371)
        self.assertTrue(all(result['checks'].values()))
        self.assertFalse(result['actual_frequency_certified'])

    def test_reciprocal_cost_includes_squared_floor(self):
        self.assertEqual(reciprocal_first(64,4),73)
        self.assertEqual(reciprocal_first(256,22),301)

    def test_bad_exponents_rejected(self):
        for args in ((),(-1,),('64',),(True,)):
            with self.assertRaises(ValueError):
                sum_exp(*args)

    def test_all_frozen_dependencies_unchanged(self):
        self.assertTrue(all(inherited_provenance().values()))


class InputAndScopeTests(unittest.TestCase):
    def data(self):
        y,eta,h = mp.mpf('.2'),mp.mpf('.3'),mp.mpf('.003')
        return y,eta,h,*manufactured(SecondJet(y,y=1),SecondJet(eta,e=1))

    def test_pressure_rejects_missing_moment(self):
        y,e,h,E,U,m,axis = self.data()
        with self.assertRaises(ValueError):
            pressure_from_moments(y,e,h,E,U,m[:4],axis)

    def test_pressure_rejects_first_order_moment_data(self):
        y,e,h,E,U,m,axis = self.data()
        with self.assertRaises(ValueError):
            pressure_from_moments(y,e,h,E,U,[v.first() for v in m],axis)

    def test_pressure_rejects_nonpositive_E(self):
        y,e,h,E,U,m,axis = self.data()
        with self.assertRaises(ValueError):
            pressure_from_moments(y,e,h,SecondJet(0),U,m,axis)

    def test_pressure_rejects_invalid_domain(self):
        y,e,h,E,U,m,axis = self.data()
        for yy,ee,hh in ((mp.inf,e,h),(y,mp.mpf('1.1'),h),(y,e,0)):
            with self.assertRaises(ValueError):
                pressure_from_moments(yy,ee,hh,E,U,m,axis)

    def record(self):
        return deepcopy(dict(status=STATUS,flags=FLAGS,resolved_obligations=RESOLVED,
            inherited_resolved_obligations=INHERITED_RESOLVED,remaining_obligations=REMAINING,gates={'valid':True}))

    def test_only_compact_envelope_newly_closed(self):
        validate(self.record())
        self.assertEqual(RESOLVED,['actual_compact_jet_envelope'])
        self.assertEqual(len(REMAINING),2)

    def test_frequency_and_blowup_promotions_rejected(self):
        for key in ('actual_joined_profile_frequency_selected','full_admissible_stress_realized','blowup_verified'):
            record = self.record()
            record['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(record)

    def test_error_obligations_cannot_disappear(self):
        for key in REMAINING:
            record = self.record()
            del record['remaining_obligations'][key]
            with self.assertRaises(ValueError):
                validate(record)

    def test_failed_bound_rejected(self):
        record = self.record()
        record['gates']['valid'] = False
        with self.assertRaises(ValueError):
            validate(record)

    def test_protocol_derivative_orders_are_locked(self):
        protocol = json.loads((HERE/'protocol.json').read_text())
        validate_protocol(protocol)
        for key,value in (('radial_derivative','d/dX'),('field_and_moment_order',1),
                          ('stress_coordinate_order',2),('envelope','C^512')):
            with self.assertRaises(ValueError):
                validate_protocol(dict(protocol,**{key:value}))


if __name__ == '__main__':
    unittest.main()
