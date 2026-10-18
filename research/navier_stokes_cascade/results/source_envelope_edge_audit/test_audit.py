"""Exact mixed normalization, signed half weights, flat edges and cone correction."""
from copy import deepcopy
from fractions import Fraction as F
import inspect,json
from pathlib import Path
import tempfile,unittest
import audit
from source_bounds import (ceil_root,frozen_margin,acceptance_floor,indices,mul,exp_jets,quotient_jets,
    reciprocal_half,edge_logprime,safe_edge_point,edge_absorption,LOGJ,
    acceptance_review,envelope_review,weighted_review)

class SourceEnvelopeEdgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate()
    def test_parent_and_protocol(self):
        self.assertEqual(audit.protocol()['parent_commit'],'ffd0b394d7f0e5df08c9d171937ae3fb1d362d8e')
        self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_integer_roots_independent_rational_grid(self):
        for d in (2,3,6):
            for n in range(1,60):
                for den in range(1,25):
                    x=F(n,den);k=ceil_root(x,d);self.assertLess((k-1)**d,x);self.assertGreaterEqual(k**d,x)
    def test_sixth_root_exact_large_power(self):self.assertEqual(ceil_root(10**120,6),10**20)
    def test_sixth_root_one_above_large_power(self):self.assertEqual(ceil_root(10**120+1,6),10**20+1)
    def test_zero_and_tiny_root(self):self.assertEqual(ceil_root(0,6),0);self.assertEqual(ceil_root(F(1,10**100),6),1)
    def test_root_invalid_inputs_rejected(self):
        for x,d in ((-1,6),(1,0),(1,True)):
            with self.assertRaises(ValueError):ceil_root(x,d)
        with self.assertRaises(TypeError):ceil_root(1.0,6)
    def test_sharp_freezing_half_margin_at_equality(self):
        for g in (F(1,10),F(1,4),F(1,2),F(1)):
            self.assertEqual(frozen_margin(g,g/(4-g)),g/2)
    def test_beyond_sharp_threshold_loses_half_margin(self):
        g=F(1,4);self.assertLess(frozen_margin(g,g/(4-g)+F(1,1000)),g/2)
    def test_padded_eighth_threshold_is_safe(self):
        for g in (F(1,100),F(1,4),F(1)):self.assertGreaterEqual(frozen_margin(g,g/8),g/2)
    def test_freezing_domain_rejected(self):
        for g,t in ((0,0),(2,0),(F(1,4),1),(F(1,4),-1)):
            with self.assertRaises(ValueError):frozen_margin(g,t)
    def test_all_twelve_corrected_certificates(self):
        rows=acceptance_review()['rows'];self.assertEqual(len(rows),12)
        for r in rows:
            self.assertGreaterEqual(F(r['certificate']['freezing_margin']),0)
            self.assertGreaterEqual(F(r['certificate']['frozen_margin_lower']),F(r['gamma0'])/2)
    def test_missing_freezing_control_passes_old_length_and_smallness(self):
        w=acceptance_review()['omitted_freezing_parameter_witness'];self.assertTrue(w['length_passes']);self.assertTrue(w['whole_pulse_smallness_sufficient']);self.assertFalse(w['freezing_passes'])
    def test_freezing_can_dominate_old_displayed_floor(self):
        r=acceptance_floor(25,F(1,200),F(1,6057072),16,10**100,F(1,4),F(1,8),0)
        self.assertGreater(r['corrected_floor'],r['old_displayed_floor'])
    def test_unrelated_larger_length_budget_margin_rejected(self):
        with self.assertRaises(ValueError):acceptance_floor(25,F(1,200),F(1,6057072),16,1,F(1,4),F(1,4),0)
    def test_smaller_length_budget_margin_accepted(self):
        r=acceptance_floor(25,F(1,200),F(1,6057072),16,1,F(1,4),F(1,16),0);self.assertGreaterEqual(F(r['freezing_margin']),0)
    def test_optional_carrier_floor_retained(self):
        r=acceptance_floor(25,F(1,200),F(1,6057072),16,0,F(1,4),F(1,8),10**50)
        self.assertGreaterEqual(F(r['length_margin']),0);self.assertGreater(r['corrected_floor'],10**20)
    def test_invalid_acceptance_inputs_rejected(self):
        args=[25,F(1,200),F(1,6057072),16,1,F(1,4),F(1,8),0]
        for i,x in ((0,0),(1,0),(2,0),(3,15),(4,-1),(5,0),(7,-1)):
            a=args[:];a[i]=x
            with self.assertRaises(ValueError):acceptance_floor(*a)
    def test_new_selector_has_no_jet_order_argument(self):self.assertNotIn('order',inspect.signature(acceptance_floor).parameters)
    def test_mixed_jet_count_is_triangular(self):self.assertEqual(len(indices(6)),28)
    def test_all_83_mixed_quotient_controls(self):self.assertEqual(envelope_review()['mixed_jet_identities'],83)
    def test_exponential_and_inverse_jets_multiply_to_one(self):
        E=exp_jets(LOGJ,5);Q=exp_jets({I:-v for I,v in LOGJ.items()},5);z=mul(E,Q,5)
        self.assertTrue(all(v==(1 if not any(I) else 0) for I,v in z.items()))
    def test_mixed_log_derivative_cost_not_omitted(self):
        E=exp_jets(LOGJ,2);self.assertEqual(E[(1,1)],F(3,5)-6);self.assertNotEqual(E[(1,1)],F(3,5))
    def test_quotient_recovers_independent_source_factor(self):
        E=exp_jets(LOGJ,4);H={I:F(1+I[0]+3*I[1]) for I in indices(4)}
        self.assertEqual(quotient_jets(E,mul(E,H,4),4),H)
    def test_tiny_point_envelope_cancels_without_underflow(self):
        E=exp_jets(LOGJ,3);H={I:F(1+sum(I),3) for I in indices(3)};f=mul(E,H,3);P0=F(1,10**1000)
        self.assertEqual(quotient_jets({I:P0*v for I,v in E.items()},{I:P0*v for I,v in f.items()},3),H)
    def test_quotient_needs_positive_envelope(self):
        with self.assertRaises(ValueError):quotient_jets({(0,0):0},{},2)
    def test_nonnegative_Bell_bounds_dominate_inverse(self):
        B=exp_jets({I:abs(v) for I,v in LOGJ.items()},6);Q=exp_jets({I:-v for I,v in LOGJ.items()},6)
        self.assertTrue(all(abs(Q[I])<=B[I] for I in B))
    def test_empty_log_jets_give_constant_envelope(self):self.assertTrue(all(v==(1 if not any(I) else 0) for I,v in exp_jets({},4).items()))
    def test_reciprocal_half_square_identity(self):
        Y={I:F(1+sum(I),3) for I in indices(4)};Y[(0,0)]=1;R=reciprocal_half(Y,4);z=mul(Y,mul(R,R,4),4)
        self.assertTrue(all(v==(1 if not any(I) else 0) for I,v in z.items()))
    def test_reciprocal_half_rejects_unnormalized_fixture(self):
        with self.assertRaises(ValueError):reciprocal_half({(0,0):2},3)
    def test_all_63_signed_mixed_half_weight_controls(self):self.assertEqual(sum(r['mixed_signed_jets'] for r in weighted_review()['signed_rows']),63)
    def test_no_unweighted_boundary_root_floor_claimed(self):self.assertFalse(audit.FLAGS['uniform_positive_unweighted_reciprocal_root_floor_at_shell_edges'])
    def test_eight_quantitative_edge_derivative_margins(self):
        rows=weighted_review()['edge_rows'];self.assertEqual(len(rows),8);self.assertTrue(all(F(r['log_derivative_margin'])>=0 for r in rows))
    def test_left_primitive_bound_keeps_small_width_cost(self):
        a=F(1,10**40);x=safe_edge_point(a,4,2,11)
        self.assertGreaterEqual(edge_logprime(a,4,2,x,11,20),a/x**3)
    def test_right_edge_is_symmetric_not_left_primitive_without_moment(self):
        b=F(4);d=safe_edge_point(b,1,2,7);x=2-d
        self.assertLessEqual(edge_logprime(1,b,2,x,3,7),-b/d**3)
    def test_nonzero_polynomial_total_detected(self):self.assertEqual(F(weighted_review()['positive_polynomial_total']),F(1,30))
    def test_compact_patching_can_keep_nonzero_total_but_has_remainder(self):
        r=weighted_review()['compact_zero_mode_control'];self.assertEqual(F(r['left_value']),0);self.assertEqual(F(r['right_value']),0);self.assertEqual(F(r['interior_A']),F(1,20))
    def test_compact_primitive_derivative_keeps_nonzero_A_zero_mode(self):self.assertEqual(F(weighted_review()['compact_zero_mode_control']['derivative_mid']),F(1,80))
    def test_flat_absorption_p_zero_and_half_width(self):self.assertEqual(edge_absorption(0,F(1,2)),1);self.assertEqual(edge_absorption(3,F(1,2)),8)
    def test_flat_absorption_small_width_is_not_free(self):self.assertGreater(edge_absorption(10,F(1,1000)),edge_absorption(10,1))
    def test_flat_absorption_invalid_inputs_rejected(self):
        with self.assertRaises(ValueError):edge_absorption(-1,1)
        with self.assertRaises(ValueError):edge_absorption(2,0)
    def test_edge_domain_inputs_rejected(self):
        with self.assertRaises(ValueError):edge_logprime(1,4,2,0,0,0)
        with self.assertRaises(ValueError):safe_edge_point(0,4,2,0)
    def test_recursive_historical_source_bytes(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_frozen_exact_record_and_roundtrip(self):
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()));audit.validate(json.loads(json.dumps(self.record)))
    def test_all_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_missing_freezing_failure_cannot_be_removed(self):
        r=deepcopy(self.record);r['reviews']['corrected_freezing_and_tail_contract']['omitted_freezing_parameter_witness']['freezing_passes']=True
        with self.assertRaises(ValueError):audit.validate(r)
    def test_actual_target_remaining_gate_cannot_be_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_Bell_identity_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['mixed_source_quotient_majorants']['mixed_jet_identities']-=1
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_cli_replay_and_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)])
            audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
