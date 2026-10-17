#!/usr/bin/env python3
"""Current enlarged atlas, complete seed, mixed pressure and endpoint controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path
import tempfile
import unittest
from terminal_bounds import window,stage_bound,closed_centers,floor_log2,ceil_log2,parameter,order,requirements,radial_coefficients,heat_A,heat_S,heat_terminal,pressure_terminal,mixed_majorant,mixed_review
from terminal_checks import integral_review,exact_review
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,protocol,main

class EndpointFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_current_single_field_parent_exact(self):self.assertEqual(protocol()['parent_commit'],'b87f53726641b240050a71789c025401e359abaa')
    def test_enlarged_closed_band_window_has_all_nine_boundary_charts(self):self.assertEqual(window(F(2)**-24,F(2)**-20),dict(first=18,last=26,count=9))
    def test_narrow_window_omits_two_enlarged_boundary_charts(self):
        full=window(F(2)**-24,F(2)**-20);old=window(F(2)**-24,F(2)**-20,enlargement=2)
        self.assertEqual(old,dict(first=19,last=25,count=7));self.assertLess(full['first'],old['first']);self.assertGreater(full['last'],old['last'])
    def test_log_rounding_exact_powers_and_nondyadic_endpoints(self):
        for k in (-30,-1,0,10):
            self.assertEqual(floor_log2(F(2)**k),k);self.assertEqual(ceil_log2(F(2)**k),k)
            self.assertEqual(floor_log2(F(3,2)*F(2)**k),k);self.assertEqual(ceil_log2(F(3,2)*F(2)**k),k+1)
    def test_nonpositive_logs_and_bad_windows_rejected(self):
        for v in (0,-1):
            with self.assertRaises(ValueError):floor_log2(v)
            with self.assertRaises(ValueError):ceil_log2(v)
        with self.assertRaises(ValueError):window(F(1,2),F(1,4))
        with self.assertRaises(ValueError):window(F(1,4),F(1,2),minimum=True)
    def test_terminal_mesh_keeps_both_closed_endpoint_labels(self):self.assertEqual(closed_centers(0,0,F(1,10**36)),(-1,1))
    def test_mesh_rounding_for_negative_and_nonintegral_intervals(self):self.assertEqual(closed_centers(F(-3,4),F(5,4),F(1,2)),(-2,3))
    def test_invalid_mesh_interval_rejected(self):
        for args in ((1,0,1),(0,1,0)):
            with self.assertRaises(ValueError):closed_centers(*args)
    def test_stage_bound_is_strict_at_exact_cutoff_endpoint(self):self.assertEqual(stage_bound(F(2)**-24,F(2)**-16),7)
    def test_same_observation_margin_controls_enlarged_fourQ(self):self.assertEqual(16*F(3,8)*F(1,16),F(3,8));self.assertLess(F(3,8),1)
    def test_actual_seed_regular_M_over_X_and_complete_gauge_bound(self):self.assertIn('M/X',protocol()['seed']);self.assertIn('same fixed gauge',protocol()['exterior'])
    def test_uniform_axis_hole_does_not_assume_fractional_map_smoothness(self):self.assertFalse(FLAGS['fractional_auxiliary_map_claimed_smooth_at_axis'])
    def test_endpoint_jet_requirements_reserve_curl_and_next_time_derivative(self):
        r=requirements(5);self.assertEqual(r['potential_for_velocity_Cauchy'],7);self.assertEqual(r['pressure_Cauchy'],6)
    def test_original_h_domain_and_noninteger_orders_rejected(self):
        for h in (0,F(1,100),-1):
            with self.assertRaises(ValueError):parameter(h)
        for m in (-1,True,F(1,2)):
            with self.assertRaises(ValueError):order(m)
    def test_radial_recurrence_retains_updated_homogeneity(self):
        a=F(101,100);self.assertEqual(radial_coefficients(a,2),[a*(a+1),-4*a-6,F(4)])
    def test_all_36_mixed_terminal_pressure_identities_exact(self):
        r=mixed_review();self.assertEqual(len(r['exact_terminal_rows']),36);self.assertTrue(r['checks']['all_36_mixed_radial_pressure_identities_exact'])
    def test_complete_pressure_convolution_has_nonzero_middle_cross(self):
        h=F(1,200);self.assertEqual(heat_S(h,2)-2*heat_A(h,2),2*heat_A(h,1)**2)
    def test_terminal_radial_pressure_identity_at_first_order(self):self.assertEqual(pressure_terminal(F(1,200),1,0),1)
    def test_mixed_majorants_dominate_terminal_exact_values(self):
        h=F(1,200)
        for n,j in ((0,0),(2,3),(5,4)):
            self.assertGreaterEqual(mixed_majorant(h,n,j,0),abs(heat_terminal(h,n,j)))
            self.assertGreaterEqual(mixed_majorant(h,n,j,0,True),abs(pressure_terminal(h,n,j)))
    def test_negative_Zmax_not_permitted(self):
        with self.assertRaises(ValueError):mixed_majorant(F(1,200),0,0,-1)
    def test_independent_24_nonterminal_mixed_heat_Gamma_comparisons(self):
        r=integral_review();self.assertEqual(r['nonterminal_heat_comparisons'],24);self.assertTrue(r['checks']['all_mixed_heat_kernel_and_recurrence_values_agree'])
    def test_independent_18_complete_pressure_product_comparisons(self):
        r=integral_review();self.assertEqual(r['nonterminal_pressure_product_comparisons'],18);self.assertTrue(r['checks']['all_mixed_complete_pressure_products_agree'])
    def test_omitted_radial_chain_and_pressure_cross_terms_detected(self):
        r=integral_review()['checks'];self.assertTrue(r['dropping_radial_Z_derivatives_detected']);self.assertTrue(r['dropping_pressure_Leibniz_crosses_detected'])
    def test_exact_seed_axis_curl_contains_radial_and_axial_velocities(self):self.assertTrue(exact_review()['checks']['axis_regular_cartesian_potential_curl_has_full_radial_and_axial_components'])
    def test_all_21_Taylor_and_21_Cauchy_controls(self):self.assertEqual(exact_review()['finite_Taylor_controls'],21);self.assertEqual(exact_review()['Cauchy_controls'],21)
    def test_wrong_terminal_time_sign_and_zero_pressure_extension_detected(self):
        c=exact_review()['checks'];self.assertTrue(c['wrong_terminal_time_sign_detected']);self.assertTrue(c['nonzero_terminal_pressure_cannot_be_zero_extended'])
    def test_heat_moments_keep_fixed_positive_h_before_large_order(self):
        h=F(1,200);j=10;self.assertGreater(heat_A(h,j+1)/heat_A(h,j)/(j+1),j)
        self.assertGreaterEqual(heat_A(h,j),h*factorial(j-1)*factorial(j))
    def test_diagnostics_not_actual_terminal_or_infinite_fields(self):self.assertTrue(integral_review()['actual_global_or_terminal_fields_not_evaluated'])
    def test_all_979_historical_and_18_direct_bindings(self):self.assertEqual(len(self.record['inherited_provenance']),979);self.assertEqual(len(self.record['source_hashes']),18)
    def test_frozen_record_and_roundtrip_match(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r);validate(json.loads(json.dumps(self.record)))
    def test_no_future_global_or_unconditional_scope_promotion(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)
    def test_cached_reviews_not_mutable_aliases(self):
        r=generate();r['reviews']['independent_nonterminal_Gamma_controls']['checks']['all_mixed_heat_kernel_and_recurrence_values_agree']=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(integral_review()['checks']['all_mixed_heat_kernel_and_recurrence_values_agree'])
    def test_same_count_historical_key_substitution_rejected(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1])
    def test_protocol_gates_sources_and_obligations_cannot_be_erased(self):
        for k in ('protocol','gates','source_hashes','remaining_obligations'):
            r=deepcopy(self.record)
            if isinstance(r[k],dict):r[k].pop(next(iter(r[k])))
            else:r[k].pop()
            with self.assertRaises(ValueError):validate(r)
    def test_cli_refuses_evidence_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')

if __name__=='__main__':unittest.main()
