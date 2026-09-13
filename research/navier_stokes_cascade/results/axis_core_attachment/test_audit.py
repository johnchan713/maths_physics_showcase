"""Independent identities, exact controls, and evidence-scope regressions."""
import copy
from fractions import Fraction as Q
import json
from pathlib import Path
import unittest
import mpmath as mp
import audit
from bounds import certificate, elementary_certificate, entry_transfer, lower, upper, point, tail_bound_log, weight
from identities import Polynomial, source_identities, stress_identities, stress_values


class AnalyticCoreTests(unittest.TestCase):
    def test_original_equations_expand_to_the_remainders(self):
        self.assertTrue(all(source_identities().values()))

    def test_polynomial_oracle_does_not_erase_terms(self):
        x, y = Polynomial.variable('x'), Polynomial.variable('y')
        self.assertFalse((x+y)*(x-y)-x*x+y*y)
        self.assertTrue((x+y)*(x-y)-x*x)

    def test_all_index_weight_formulas_at_boundary_indices(self):
        for n in (0,1,7,31):
            for k in (0,1,8,32):
                radial = Q(20*(n+2)**2, (n+1)*(n+k+1))
                angular = Q(20*(n+2)**2*(k+1)**2, (n+1)*(k+2)**2)
                self.assertEqual(weight(n,k)/weight(n+1,k),radial)
                self.assertEqual(weight(n,k+1)/weight(n+1,k),angular)

    def test_exact_weight_convolution_diagnostics(self):
        result=audit.exact_weight_checks(8)
        self.assertTrue(all(result['checks'].values()))
        self.assertFalse(result['exhaustive_all_indices'])

    def test_scalar_comparison_polynomial_and_tail(self):
        t=Q(41,20)
        self.assertEqual(1-t/2+t*t/12-t**3/144,Q(305719,1152000))
        self.assertLess(-Q(1,2)+t/6,0)
        self.assertLess(t/(4*5),1)
        t=Q(99,50)
        self.assertEqual(1-t+t*t/4-t**3/36+t**4/576,-Q(75535511,400000000))
        self.assertLess(-1+Q(2,2)-t*t/12+Q(8,144),0)
        self.assertLess(Q(2,25),1)

    def test_majorant_source_coefficient_budgets(self):
        # Independently bound the fixed analytic coefficients before adding
        # their monomials. invL<=2, |eta|<=2, d<=5, A<=1, D<=1/2.
        self.assertLessEqual(2*(64+Q(1,100)*(1+4*9))+2*64+2*64,386)
        self.assertLessEqual(4+10+Q(8,100)+4+10+10,40)
        self.assertEqual(2*(1+8*9+5*4)+2*64+2*64,442)
        self.assertEqual(8+4+10+10,32)
        self.assertEqual(16+10+8,34)

    def test_interval_core_and_parameter_bounds(self):
        for scale in ('P16','Md4','Md64'):
            result=certificate(scale,80)
            self.assertTrue(all(result['checks'].values()),scale)
        self.assertTrue(all(elementary_certificate().values()))

    def test_old_Lambda_is_outside_this_certificate(self):
        result=certificate('Md64')
        self.assertGreater(lower(result['log_bounds']['map_Lipschitz_numerator']),mp.log(65536))
        # Failure of this sufficient bound is not a nonexistence theorem.

    def test_finite_polynomial_is_not_exact_core(self):
        result=certificate('P16')
        tail=tail_bound_log(result['log_bounds']['ball_norm_bound'],24)
        self.assertGreater(lower(tail),0)

    def test_tail_formula_against_positive_infinite_sum(self):
        mp.iv.dps=80
        with mp.workdps(80):
            for n in (0,4,24):
                actual=mp.nsum(lambda k:mp.mpf('.205')**k/(k+1)**2,[n+1,mp.inf])
                bound=mp.iv.exp(tail_bound_log(point(0),n))
                self.assertLess(actual,lower(bound))

    def test_reject_invalid_parameter_evaluations(self):
        for args in (('unknown',80),('P16',20)):
            with self.assertRaises(ValueError):certificate(*args)
        for n,r in ((-1,'4.1'),(2,'20'),(2,'0')):
            with self.assertRaises(ValueError):tail_bound_log(point(0),n,r)


class ContinuationTests(unittest.TestCase):
    def test_exact_stress_perturbation_identities(self):
        self.assertTrue(all(stress_identities().values()))

    def test_common_kappa_cancels_for_tiny_values(self):
        result=audit.stress_controls()
        self.assertTrue(result['kappa_cancellation'])
        self.assertTrue(result['unactivated_stress_zero'])
        self.assertTrue(result['frozen_E_ratio_changes_Pc'])

    def test_cone_opens_under_shear_reduction(self):
        for q in (0,7,10**30):
            row=stress_values(3,q,Q(1,10**100),1)
            self.assertGreater(row['Pc'],2)
            self.assertLess(row['vs'],1)
            self.assertEqual(row['Jc'],0)
            self.assertGreater(row['Pc']-row['vs'],0)

    def test_zero_stress_cannot_pass_as_strict(self):
        row=stress_values(3,7,1,1)
        self.assertGreater(row['Pc'],2)
        self.assertEqual(row['Pc']-row['vs'],0)

    def test_zero_kappa_is_not_silently_divided(self):
        with self.assertRaises(ValueError):stress_values(3,7,0,1)

    def test_C_independent_shape_and_C_dependent_small_widths(self):
        result=certificate('Md64')
        logs=result['log_bounds']
        self.assertLess(upper(logs['activation_width']),-100)
        self.assertGreater(lower(logs['C']),upper(logs['transition_length_Tsh']))
        self.assertLess(upper(logs['core_ledger_C1_bound']),mp.log(mp.mpf('1e-16')/4))

    def test_enlarging_radius_alone_does_not_bound_pressure_moment(self):
        # Cp=integral E^2/(2X) is invariant under re-labelling X/XR.
        # An actual amplitude decay is required as well as a small xsep.
        e2_over_2x=Q(7,3)
        cp=e2_over_2x*Q(5,2)
        for xr in (10,10**100):
            self.assertEqual(e2_over_2x*xr*(Q(5,2)/xr),cp)

    def test_constant_error_budget_cannot_preserve_j_scale_sign(self):
        floor=Q(1,10**18)
        perturbation=-Q(1,1000)
        self.assertLess(abs(perturbation),Q(1,100))
        self.assertLess(floor+perturbation,0)
        self.assertGreater(floor-floor/100,0)

    def test_actual_entry_bound_transfer(self):
        result=entry_transfer()
        self.assertTrue(all(result['checks'].values()))
        self.assertEqual(len(result['entry_moment_C1_bounds']),5)
        self.assertLess(max(upper(v) for v in result['entry_moment_C1_bounds']),mp.mpf('2.55e-17'))


class RecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=audit.encode(audit.generate())

    def test_all_gates_pass(self):
        self.assertTrue(all(self.record['gates'].values()))

    def test_source_hash_is_exactly_the_frozen_PDF(self):
        self.assertEqual(self.record['protocol']['source_pdf_sha256'],
            '0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f')

    def test_no_full_PDE_or_blowup_promotion(self):
        for name in ('full_admissible_stress_realized','full_PDE_corrections_verified',
                     'smooth_force_verified','blowup_verified'):
            self.assertFalse(self.record['scientific_scope'][name])
        self.assertFalse(self.record['numerical_matched_profile_materialized'])

    def test_refinement_records_cover_selected_large_scale(self):
        self.assertEqual([r['scale'] for r in self.record['refinement']],['P16','Md4','Md64'])
        self.assertEqual(len(self.record['certificates']),6)

    def test_scientific_status_tampering_is_rejected(self):
        changed=copy.deepcopy(self.record)
        changed['scientific_scope']['blowup_verified']=True
        with self.assertRaises(ValueError):audit.compare(changed,self.record)

    def test_missing_gate_is_rejected(self):
        changed=copy.deepcopy(self.record)
        changed['gates'].pop(next(iter(changed['gates'])))
        with self.assertRaises(ValueError):audit.compare(changed,self.record)

    def test_changed_interval_result_is_rejected(self):
        changed=copy.deepcopy(self.record)
        changed['certificates'][0]['log_bounds']['Lambda']['upper']='1'
        with self.assertRaises(ValueError):audit.compare(changed,self.record)

    def test_changed_source_hash_is_rejected(self):
        changed=copy.deepcopy(self.record)
        key=next(iter(changed['source_hashes']))
        changed['source_hashes'][key]='0'*64
        with self.assertRaises(ValueError):audit.compare(changed,self.record)

    def test_frozen_evidence_reproduces_if_present(self):
        path=Path(__file__).with_name('evidence.json')
        if path.exists():audit.compare(json.loads(path.read_text()),self.record)


if __name__=='__main__':
    unittest.main()
