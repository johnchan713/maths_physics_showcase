"""Exact harmonic counterexamples, majorants, full covariance and evidence guards."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE,FLAGS,protocol,generate,validate,compare_records,main,verify_bindings
from mode_bounds import moment,convolution,smoothed,inverse_review,linear_review,pressure_review,cycle_review
from mode_checks import C,P,I,X,curl,divergence,transport,real_field,field_norm,covariance,fixture,add_fields,obstruction,complete_review


class HarmonicTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()
    def test_protocol_current_parent(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()),protocol())
        self.assertEqual(protocol()['parent_commit'],'0c1399dca354f460b80c2b22ef80f63ea285fc95')
    def test_convolution_zero_and_unequal_modes(self):
        self.assertEqual(convolution({1:F(2),-2:F(3)},{-1:F(5),7:F(11)}),{0:10,8:22,-3:15,5:33})
    def test_weight_submultiplicative_arbitrary_signs(self):
        for m in range(-17,18):
            for n in range(-17,18):self.assertLessEqual(1+abs(m+n),(1+abs(m))*(1+abs(n)))
    def test_empty_mode_sums(self):self.assertEqual(moment({},3),0);self.assertEqual(convolution({},{}),{})
    def test_invalid_diagnostic_moment(self):
        for s in (-1,F(1,2),1.5):
            with self.assertRaises(ValueError):moment({1:1},s)
    def test_nonzero_inverse_modes_required(self):
        with self.assertRaises(ValueError):smoothed({0:1},2)
    def test_pressure_velocity_potential_have_different_gains(self):
        r=inverse_review();self.assertTrue(all(r['checks'].values()))
        self.assertEqual({r['smoothing'] for r in r['rows']},{1,2,3})
    def test_high_mode_inverse_weight_gain(self):
        for m in (2,5,17,10000):self.assertLessEqual(F(1,m*m-1),F(4,3*m*m))
    def test_inverse_requires_full_path_envelope_norm(self):
        self.assertIn('full-path envelope-normalized N^W',protocol()['inverse'])
        self.assertTrue(FLAGS['full_path_envelope_normalized_forcing_required_for_inverse'])
    def test_regular_linear_includes_limiting_groups(self):
        r=linear_review();self.assertEqual(min(x['index_gain'] for x in r['rows']),1)
        self.assertTrue(any('pressure' in x['term'] for x in r['rows']))
        self.assertTrue(any('damping' in x['term'] for x in r['rows']))
    def test_normal_pressure_obstructs_three_index_gain(self):self.assertTrue(all(pressure_review()['checks'].values()))
    def test_fraction_complex_arithmetic(self):
        z=C(F(2,3),F(5,7));self.assertEqual(z*z.conj(),C(F(4,9)+F(25,49)))
        self.assertEqual((I*I),C(-1))
    def test_float_diagnostics_rejected(self):
        for x in (0.1,complex(1,2)):
            with self.assertRaises(TypeError):C(x)
    def test_polynomial_product_derivative(self):self.assertEqual((X*X*X).diff(),3*X*X)
    def test_polynomial_negative_powers_rejected(self):
        with self.assertRaises(ValueError):P({-1:1})
    def test_complete_curls_divergence_all_signs(self):
        for m in (-19,-1,1,19):self.assertEqual(divergence(curl([X*X,X*X*X,I*X],m),m),0)
    def test_full_curl_not_principal_transverse(self):
        a=curl([0,I*X,0],1);self.assertNotEqual(a[2],0)
        self.assertEqual(a[0].diff()+a[2]*I,0)
        self.assertNotEqual(a[0].diff(),0)
    def test_direct_vs_repaired_transport(self):
        a,b=fixture();self.assertEqual(transport(a,b),transport(a,b,True))
    def test_mean_input_not_divided_by_zero(self):
        with self.assertRaises(ValueError):transport({0:[X,P(),P()]},{1:[P(),P(1),P()]},True)
    def test_transport_jet_majorants(self):self.assertTrue(complete_review()['checks']['coefficient_transport_majorant'])
    def test_covariance_all_nine_entries_and_symmetric(self):
        a,b=fixture();w=covariance(add_fields(a,b))
        self.assertTrue(all(w[i][j]==w[j][i] for i in range(3) for j in range(3)))
        self.assertTrue(all(w[i][j]!=0 for i in range(3) for j in range(i,3)))
    def test_covariance_without_negative_modes_is_wrong(self):
        a,_=fixture();wrong=covariance({m:x for m,x in a.items() if m>0})
        self.assertNotEqual(wrong,covariance(a))
    def test_drop_high_harmonic_square_detected(self):
        a,b=fixture();full=covariance(add_fields(a,b));bad=covariance(add_fields(a,{m:x for m,x in b.items() if abs(m)!=7}))
        self.assertNotEqual(full,bad)
    def test_signed_covariance_need_not_be_positive(self):
        u=real_field({1:[P(),P(),P(F(1,2))]});zero={}
        self.assertEqual(covariance(zero)[2][2]-covariance(u)[2][2],P(-F(1,2)))
    def test_real_obstruction_sideband_sizes(self):
        for n in (2,11,4096):
            r=obstruction(n);self.assertEqual(r['ab'][n+1][1],P(-F(n,4)));self.assertEqual(r['ab'][n-1][1],P(F(n,4)))
    def test_real_obstruction_inputs_fixed_and_transport_unbounded(self):
        for n in (2,101,10001):
            r=obstruction(n);self.assertEqual(field_norm(r['a'],0,0),2);self.assertEqual(field_norm(r['b'],0,0),1)
            self.assertEqual(field_norm(r['ab'],0,0),n);self.assertEqual(field_norm(r['ba'],0,0),0)
    def test_invalid_obstruction_mode(self):
        for n in (0,1,F(5,2)):
            with self.assertRaises(ValueError):obstruction(n)
    def test_finite_epsilon_cycle_not_promoted_to_mode_summability(self):
        self.assertTrue(all(cycle_review()['checks'].values()))
        self.assertFalse(FLAGS['actual_infinite_nonlinear_weighted_mode_sums_established'])
    def test_pressure_memory_and_original_geometry(self):
        self.assertIn('Pi1,Pi2',protocol()['memory']);self.assertIn('Jz+c_rho*P',protocol()['memory'])
        self.assertIn('C^-20',protocol()['fixed_choices']);self.assertIn('never forward-inverted',protocol()['support'])
    def test_frozen_evidence_and_recursive_sources(self):
        old=json.loads((HERE/'evidence.json').read_text());validate(old);compare_records(self.record,old)
        self.assertGreaterEqual(self.record['inherited_provenance']['records'],75)
        self.assertGreaterEqual(len(self.record['inherited_provenance']['bindings']),1038)
    def test_refuse_scope_promotion(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)
    def test_refuse_missing_or_changed_history(self):
        r=deepcopy(self.record);r['inherited_provenance']['bindings'].pop(next(iter(r['inherited_provenance']['bindings'])))
        with self.assertRaises(ValueError):validate(r)
    def test_refuse_fabricated_transport_success(self):
        r=deepcopy(self.record);r['reviews']['necessary_real_harmonic_loss']['rows'][0]['transport']='0'
        with self.assertRaises(ValueError):validate(r)
    def test_refuse_erased_limits_and_changed_protocol(self):
        for key in ('remaining_obligations','protocol'):
            r=deepcopy(self.record)
            if key=='protocol':r[key]['transport']='no loss'
            else:r[key]=[]
            with self.assertRaises(ValueError):validate(r)
    def test_historical_hash_mismatch_is_error(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'source').write_text('changed');record=p/'evidence.json'
            record.write_text(json.dumps({'source_hashes':{'source':'0'*64}}))
            with self.assertRaises(ValueError):verify_bindings(p,[record])
    def test_binding_outside_project_is_error(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);record=p/'evidence.json';record.write_text(json.dumps({'source_hashes':{'../outside':'0'*64}}))
            with self.assertRaises(ValueError):verify_bindings(p,[record])
    def test_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
    def test_type_corruption_not_equal_boolean(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):validate(r)
    def test_json_roundtrip(self):validate(json.loads(json.dumps(self.record)))
if __name__=='__main__':unittest.main()
