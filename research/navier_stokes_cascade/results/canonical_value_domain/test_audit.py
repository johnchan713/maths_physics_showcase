"""Exact acceptance arithmetic, physical-scale margins and scope regressions."""
from copy import deepcopy
from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path
import tempfile
import unittest
import audit
from acceptance import (bounds,ceil,binary_ceiling,coefficient_exponent,certificate,
    length_margin,accepted_dyadic_band,relative_gate,dyadic_review,domain_review,logarithmic_review)

class ValueDomainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.p=bounds().parent_budget().inputs();cls.cert=certificate(cls.p)

    def test_parent_and_protocol(self):
        self.assertEqual(audit.protocol()['parent_commit'],'edf21536d1268a5d8061d663427128c361ec4970')
        self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_exact_binary_powers(self):
        for n in range(-160,161):
            x=F(2**n) if n>=0 else F(1,2**(-n));self.assertEqual(binary_ceiling(x),n)
    def test_binary_ceiling_independent_rational_grid(self):
        for n in range(1,70):
            for d in range(1,60):
                x=F(n,d);b=binary_ceiling(x);high=F(2**b) if b>=0 else F(1,2**(-b))
                self.assertTrue(high/2<x<=high,(x,b))
    def test_binary_ceiling_one_above_power(self):self.assertEqual(binary_ceiling(F(2**1000+1,2**999)),2)
    def test_binary_ceiling_tiny_nonzero_not_underflow(self):self.assertEqual(binary_ceiling(F(1,2**4000+1)),-4000)
    def test_nonpositive_log_inputs_rejected(self):
        for v in (0,-1):
            with self.assertRaises(ValueError):binary_ceiling(v)
    def test_float_and_bool_inputs_rejected(self):
        for v in (1.0,True,'1'):
            with self.assertRaises(TypeError):binary_ceiling(v)
    def test_rational_ceiling_signs(self):self.assertEqual(ceil(F(-9,2)),-4);self.assertEqual(ceil(F(9,2)),5)
    def test_exact_exponential_coefficient_bound(self):
        self.assertEqual(coefficient_exponent(F(9,8),F(17,3)),7)
    def test_invalid_exponential_bounds_rejected(self):
        for c,e in ((0,1),(1,-1)):
            with self.assertRaises(ValueError):coefficient_exponent(c,e)
    def test_certificate_does_not_mutate_inputs(self):
        p=deepcopy(self.p);certificate(p);self.assertEqual(p,self.p)
    def test_all_dyadic_value_margins(self):
        for r in dyadic_review()['rows']:self.assertGreaterEqual(r['length_margin'],0);self.assertTrue(r['accepted'])
    def test_accepted_half_line_not_single_band(self):
        B=self.cert['accepted_band_log2']
        for j in (0,1,10,100):self.assertTrue(accepted_dyadic_band(self.cert,B+j));self.assertGreaterEqual(length_margin(self.cert,B+j),0)
    def test_preceding_dyadic_band_not_certified(self):self.assertFalse(accepted_dyadic_band(self.cert,self.cert['accepted_band_log2']-1))
    def test_invalid_band_exponents_rejected(self):
        for x in (-1,1.0,True):
            with self.assertRaises(ValueError):accepted_dyadic_band(self.cert,x)
            with self.assertRaises(ValueError):length_margin(self.cert,x)
    def test_radius_sixteen_times_smaller_costs_two_band_log_units(self):
        c=certificate(self.p,radius=F(1,16*6057072));self.assertEqual(c['radius_cost_log2'],self.cert['radius_cost_log2']+4)
        self.assertEqual(c['length_band_log2'],self.cert['length_band_log2']+2)
    def test_radius_not_discarded_as_log_length(self):self.assertGreater(self.cert['length_band_log2'],ceil(F(self.cert['common_length_log2'],2)))
    def test_smaller_h_increases_smallness_threshold(self):self.assertGreater(certificate(self.p,h=F(1,10**5))['smallness_band_log2'],self.cert['smallness_band_log2'])
    def test_original_accepted_floor_retained(self):
        c=certificate(self.p,old_band=2**100000);self.assertGreaterEqual(c['accepted_band_log2'],100000)
    def test_smallness_series_ratio_exact(self):self.assertEqual(F(self.cert['smallness_ratio']),F(945,8192));self.assertLess(F(self.cert['smallness_ratio']),1)
    def test_smallness_constant_arithmetic_independent(self):self.assertLess(3*factorial(8)*4**8,512**4)
    def test_both_cone_lengths_in_common_gate(self):
        for v in self.cert['cone_length_log2'].values():self.assertLessEqual(v,self.cert['common_length_log2'])
    def test_direct_and_comparison_thresholds_are_separate(self):self.assertNotEqual(self.cert['direct_length_log2'],self.cert['comparison_length_log2'])
    def test_comparison_in_common_gate(self):self.assertGreaterEqual(self.cert['common_length_log2'],self.cert['comparison_length_log2'])
    def test_haar_enclosures_change_mass_and_comparison(self):
        c=certificate(self.p,Jmin=F(1,100),Jmax=100);self.assertGreater(c['inverse_mass_log2'],self.cert['inverse_mass_log2']);self.assertGreater(c['comparison_length_log2'],self.cert['comparison_length_log2'])
    def test_smaller_cone_margin_increases_comparison(self):self.assertGreater(certificate(self.p,gamma=F(1,100))['comparison_length_log2'],self.cert['comparison_length_log2'])
    def test_invalid_ranges_rejected(self):
        for kw in ({'h':0},{'h':F(1,100)},{'radius':0},{'radius':F(1,100)},{'gamma':0},{'gamma':2},{'Jmin':2,'Jmax':1},{'old_band':15}):
            with self.assertRaises(ValueError):certificate(self.p,**kw)
    def test_input_floors_revalidated(self):
        p=deepcopy(self.p);p['lam_min']=0
        with self.assertRaises(ValueError):certificate(p)
    def test_no_missing_enclosure_filled_from_manufactured_defaults(self):
        for p in ({},{'Rmin':F(1)},dict(self.p,unreviewed_floor=F(1))):
            with self.assertRaises(ValueError):certificate(p)
    def test_invalid_exact_domain_type_rejected(self):
        with self.assertRaises(TypeError):certificate(self.p,h=.001)
    def test_relative_gate_endpoint_is_half(self):
        r=relative_gate(1,2,5);self.assertTrue(r['accepted']);self.assertEqual(F(r['relative_bound']),F(1,2))
    def test_dropping_plus_one_wrongly_accepts(self):
        r=relative_gate(1,2,4);self.assertLessEqual(F(r['theta']),F(1,4));self.assertFalse(r['accepted']);self.assertGreater(F(r['relative_bound']),F(1,2))
    def test_singular_neumann_case_not_inverted(self):self.assertIsNone(relative_gate(1,2,1)['relative_bound'])
    def test_relative_gate_positive_inputs_required(self):
        with self.assertRaises(ValueError):relative_gate(0,2,5)
    def test_enlarged_chart_both_qstar_and_tail_margins(self):
        rows=domain_review()['chart_rows'];self.assertEqual(len(rows),18);self.assertTrue(all(r['later'] and r['inside'] for r in rows))
    def test_all_independent_logarithmic_majorants(self):self.assertEqual(len(logarithmic_review()['rows']),14);self.assertTrue(all(x['below'] for x in logarithmic_review()['rows']))
    def test_common_acceptance_no_derivative_order_argument(self):
        import inspect
        self.assertNotIn('order',inspect.signature(certificate).parameters)
    def test_actual_first_band_not_claimed(self):self.assertFalse(audit.FLAGS['original_first_band_or_old_localized_force_value_acceptance_verified'])
    def test_no_actual_norm_enclosure_from_inputs_fixture(self):self.assertFalse(audit.FLAGS['actual_outward_norm_enclosures_numerically_verified'])
    def test_historical_evidence_bytes(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_frozen_roundtrip_and_replay(self):
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()));audit.validate(json.loads(json.dumps(self.record)))
    def test_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_budget_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['exact_dyadic_budgets']['rows'][0]['certificate']['accepted_band_log2']-=1
        with self.assertRaises(ValueError):audit.validate(r)
    def test_observation_scope_cannot_be_removed(self):
        r=deepcopy(self.record);r['protocol'].pop('missing')
        with self.assertRaises(ValueError):audit.validate(r)
    def test_cli_requires_fresh_output_and_reproduces(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(out)])
            audit.validate(json.loads(out.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(out)])

if __name__=='__main__':unittest.main()
