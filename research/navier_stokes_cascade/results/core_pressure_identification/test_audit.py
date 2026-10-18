"""Exact scalar inverse, endpoint obstruction and higher-pressure memory checks."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
import audit
from pressure_recovery import (poly,add,scale,mul,power,derivative,value,integral,operator,
    D,compatibility,endpoint_inverse,jet_bounds,stress,inverse_review,obstruction_review,scale_and_memory_review)

class CorePressureIdentificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate()

    def test_live_parent_is_published_canonical_comparison(self):self.assertEqual(audit.protocol()['parent_commit'],'0b16f1b18aa6277dc507e39fd7cecbd9a6e4016a')
    def test_pressure_family_and_leading_identity_are_not_newly_assumed(self):
        self.assertTrue(all(self.record['reviews']['same_reference_core_and_leading_pressure_trace']['checks'].values()))
        self.assertFalse(audit.FLAGS['leading_pressure_identification_still_an_independent_missing_equation'])
    def test_polynomial_derivative_is_exact_fraction(self):self.assertEqual(derivative(poly(1,2,3,4),2),poly(6,24));self.assertIsInstance(derivative(poly(0,0,1))[1],F)
    def test_derivative_orders_rejected(self):
        with self.assertRaises(ValueError):derivative(poly(1),-1)
    def test_polynomial_integral_orientation(self):self.assertEqual(integral(poly(1,2,3),-1,1),4);self.assertEqual(integral(poly(1,2,3),1,-1),-4)
    def test_constant_defect_operator_keeps_eta(self):self.assertEqual(operator(F(51,100),poly(3)),poly(0,F(153,25)))
    def test_operator_rejects_below_core_half_range(self):
        with self.assertRaises(ValueError):operator(F(1,4),poly(1))
    def test_all_20_weighted_endpoint_compatibilities(self):self.assertEqual(inverse_review()['weighted_compatibility_identities'],20)
    def test_all_180_two_sided_endpoint_reconstructions(self):
        r=inverse_review();self.assertEqual(r['endpoint_reconstructions'],180);self.assertTrue(all(a['recovered']==a['expected'] for a in r['rows']))
    def test_all_128_differentiated_operator_identities(self):self.assertEqual(inverse_review()['differentiated_operator_identities'],128)
    def test_both_d_derivatives_are_necessary(self):
        A=F(51,100);delta=poly(1,2,3)
        lhs=operator(A+1,derivative(delta,2))
        omitted_d_second=add(derivative(operator(A,delta),2),scale(derivative(delta),-2*(4*A)))
        self.assertNotEqual(lhs,omitted_d_second)
    def test_plus_endpoint_value(self):
        A=F(1);delta=poly(2,-1,3);g=operator(A,delta)
        self.assertEqual(endpoint_inverse(A,g,1),value(delta,1));self.assertEqual(endpoint_inverse(A,g,1),value(g,1)/(4*A))
    def test_minus_endpoint_has_opposite_sign(self):
        A=F(1);delta=poly(2,-1,3);g=operator(A,delta)
        self.assertEqual(endpoint_inverse(A,g,-1),-value(g,-1)/(4*A))
    def test_constant_source_is_incompatible(self):
        for A in (F(1,2),F(1),F(3,2)):
            self.assertGreater(compatibility(A,poly(1)),0)
            with self.assertRaises(ValueError):endpoint_inverse(A,poly(1),0)
    def test_fixture_inverse_rejects_noninteger_weight(self):
        with self.assertRaises(ValueError):endpoint_inverse(F(51,100),poly(),0)
    def test_endpoint_inverse_rejects_outside_interval(self):
        with self.assertRaises(ValueError):endpoint_inverse(F(1),poly(),2)
    def test_singular_mode_has_zero_operator_but_nonzero_delta(self):
        rows=obstruction_review()['singular_kernel_rows'];self.assertEqual(len(rows),12)
        self.assertTrue(all(F(r['g'])==0 and F(r['singular_delta'])>0 for r in rows))
    def test_singular_mode_grows_towards_endpoint(self):
        rows=[r for r in obstruction_review()['singular_kernel_rows'] if r['A']=='1/2']
        vals=[F(r['singular_delta']) for r in rows];self.assertEqual(vals,sorted(vals));self.assertGreater(vals[-1],50)
    def test_smooth_zero_source_recovers_zero(self):
        for A in (F(1,2),F(1)):
            for eta in (-1,F(-1,2),0,F(1,2),1):self.assertEqual(endpoint_inverse(A,poly(),eta),0)
    def test_jet_bounds_keep_L_first_and_second_derivatives(self):
        G,B=jet_bounds(F(1,100),[1,2,3]);self.assertEqual(G,[F(1),F(51,25),F(16,5)])
        self.assertEqual(B[0],1);self.assertGreater(B[2],G[2])
    def test_jet_bounds_zero_data(self):self.assertEqual(jet_bounds(F(1,100),[0]*8),([F(0)]*8,[F(0)]*8))
    def test_jet_bounds_invalid_inputs_rejected(self):
        for h,n in ((F(-1,100),[1]),(F(1,50),[1]),(F(0),[]),(F(0),[-1])):
            with self.assertRaises(ValueError):jet_bounds(h,n)
    def test_shifted_inverse_cost_is_not_order_uniform(self):
        G,B=jet_bounds(F(1,100),[1]*8);self.assertGreater(B[-1],100*B[1])
    def test_scale_invariance_retains_radial_jacobian(self):
        r=scale_and_memory_review();self.assertEqual(len({a['Cp'] for a in r['scale_rows']}),1)
        self.assertGreater(len({a['wrong_without_inverse_X'] for a in r['scale_rows']}),1)
    def test_pressure_stress_matches_original_ledger(self):
        r=scale_and_memory_review();self.assertEqual(len(r['stress_rows']),36)
        self.assertTrue(all(a['defect']==a['original_ledger'] for a in r['stress_rows']))
    def test_zero_at_eta_zero_does_not_identify_pressure_function(self):
        delta=poly(0,0,1);self.assertEqual(stress(0,F(1,100),1,delta),0)
        self.assertNotEqual(stress(F(1,2),F(1,100),1,delta),0)
    def test_pressure_stress_grows_linearly_in_X(self):
        delta=poly(0,0,1);self.assertEqual(stress(F(1,2),F(1,100),14,delta),2*stress(F(1,2),F(1,100),7,delta))
    def test_pressure_stress_domain_checks(self):
        for eta,h,X in ((2,F(0),1),(0,F(1,50),1),(0,F(0),0)):
            with self.assertRaises(ValueError):stress(eta,h,X,poly(1))
    def test_first_background_pressure_difference_is_nonzero(self):self.assertGreater(F(scale_and_memory_review()['first_pressure_increment']),0)
    def test_second_pressure_constant_is_not_erased(self):self.assertNotEqual(F(scale_and_memory_review()['second_pressure_memory']),0)
    def test_full_corrected_pressure_not_replaced_by_leading_pressure(self):self.assertFalse(audit.FLAGS['full_corrected_NS_pressure_identified_with_leading_centrifugal_pressure'])
    def test_huge_actual_first_band_acceptance_not_claimed(self):
        self.assertFalse(audit.FLAGS['actual_canonical_first_band_value_acceptance_verified'])
        self.assertIn('theta=12CH',audit.protocol()['next_gate'])
    def test_frozen_record_matches_fresh_and_json_roundtrip(self):
        r=json.loads((audit.HERE/'evidence.json').read_text());audit.validate(r);audit.compare_records(self.record,r)
        audit.compare_records(self.record,json.loads(json.dumps(self.record)))
    def test_scope_and_endpoint_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                bad=deepcopy(self.record);bad['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(bad)
    def test_recovered_endpoint_tamper_rejected(self):
        bad=deepcopy(self.record);bad['reviews']['exact_endpoint_inverse_and_fixed_jet_controls']['rows'][0]['recovered']='999'
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_singular_mode_cannot_be_silently_removed(self):
        bad=deepcopy(self.record);bad['reviews']['endpoint_obstructions_and_scope']['singular_kernel_rows'].pop()
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_same_count_provenance_substitution_rejected(self):
        bad=deepcopy(self.record);b=bad['inherited_provenance']['bindings'];b['fake:key']=b.pop(next(iter(b)))
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_source_hash_tamper_rejected(self):
        bad=deepcopy(self.record);bad['source_hashes'][next(iter(bad['source_hashes']))]='0'*64
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_actual_acceptance_gate_cannot_be_erased(self):
        bad=deepcopy(self.record);bad['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_cached_reviews_not_mutable_aliases(self):
        bad=audit.reviews();bad['exact_endpoint_inverse_and_fixed_jet_controls']['rows'][0]['expected']='bad'
        self.assertNotEqual(audit.reviews()['exact_endpoint_inverse_and_fixed_jet_controls']['rows'][0]['expected'],'bad')
    def test_cli_refuses_evidence_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'kept.json';p.write_text('kept')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'kept')
if __name__=='__main__':unittest.main()
