"""Independent sourced-inverse identities, sharp limits and frozen-evidence guards."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from audit import HERE,FLAGS,protocol,generate,validate,compare_records,parent_review,main
from inverse_bounds import mode,residence,jet_majorants,high_mode_scaled_jets,entrance_jets,physical_powers,analytical_review
from inverse_checks import principal_checks,numerical_checks,endpoint_checks


class SourcedInverseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()
    def test_01_current_protocol(self):self.assertEqual(json.loads((HERE/'protocol.json').read_text()),protocol());self.assertEqual(protocol()['parent_commit'],'93103d22e03b2996ad8cf18ef31f77d090c22896')
    def test_02_nonzero_signed_modes(self):
        for m in (1,2,7,1000000):self.assertEqual(mode(m),mode(-m));self.assertEqual(mode(m),m*m-1)
    def test_03_invalid_modes(self):
        for m in (0,True,1.5,'2'):
            with self.assertRaises(ValueError):mode(m)
    def test_04_resonant_residence(self):
        self.assertEqual(residence(1,20,F(1,2)),20);self.assertEqual(residence(-1,20,F(1,2)),20)
    def test_05_high_residence(self):self.assertEqual(residence(2,20,F(1,2)),F(2,3));self.assertEqual(residence(100,20,F(1,2)),F(2,9999))
    def test_06_short_pulse_residence(self):self.assertEqual(residence(2,F(1,100),F(1,2)),F(1,100))
    def test_07_invalid_floors(self):
        for length,d in ((0,1),(1,0),(-1,1)):
            with self.assertRaises(ValueError):residence(2,length,d)
    def test_08_exact_all_fixed_parameter_recursion(self):self.assertTrue(analytical_review()['checks']['parameter_jets_have_high_mode_inverse_square_gain'])
    def test_09_differentiated_damping_not_omitted(self):
        a=jet_majorants(10,30,F(1,2),2,(1,0),(0,0),(0,1))
        b=jet_majorants(10,30,F(1,2),2,(1,0),(0,0),(0,0))
        self.assertGreater(a[1],0);self.assertEqual(b[1],0)
    def test_10_binomial_jet_multiplicities(self):
        z=high_mode_scaled_jets(1,1,(1,0,0),(0,0,0),(0,1,0))
        self.assertEqual(z,[1,1,2])
    def test_11_negative_jet_bounds_rejected(self):
        with self.assertRaises(ValueError):jet_majorants(2,1,1,1,(1,-1),(0,0),(0,0))
        with self.assertRaises(ValueError):high_mode_scaled_jets(1,1,(1,),(0,0),(0,))
    def test_12_generic_entrance_derivative_has_no_inverse_square_gain(self):
        for m in (2,7,1000000):self.assertEqual(entrance_jets(m,1,(1,0,0))[1],1)
    def test_13_full_flat_entrance_supplies_all_fixed_zero_data(self):
        for m in (1,2,10):self.assertEqual(entrance_jets(m,1,(0,)*12),[F(0)]*13)
    def test_14_finite_compatibility_is_not_infinite_flatness(self):
        for m in (2,10):self.assertEqual(entrance_jets(m,1,(0,0,0,1)),[0,0,0,0,1])
    def test_15_all_rational_moving_principal_components(self):self.assertEqual(principal_checks()['exact_ODE_real_components'],144);self.assertTrue(principal_checks()['checks']['all_manufactured_principal_equations_and_pressure_exact'])
    def test_16_both_normal_terms_detected(self):
        r=principal_checks()['checks'];self.assertTrue(r['moving_normal_omission_detected']);self.assertTrue(r['source_normal_omission_detected'])
    def test_17_complex_negative_harmonics(self):self.assertTrue(principal_checks()['checks']['negative_harmonic_reality'])
    def test_18_pressure_has_sharp_direct_source_loss(self):self.assertTrue(principal_checks()['checks']['pressure_cannot_inherit_velocity_inverse_square_gain'])
    def test_19_noncommuting_propagator_factorization(self):
        r=numerical_checks();self.assertTrue(r['checks']['noncommuting_matrix_fixture']);self.assertLess(r['max_factor_absolute_error'],2e-8)
    def test_20_energy_and_residence_controls(self):
        r=numerical_checks()['checks'];self.assertTrue(r['energy_propagator_majorant']);self.assertTrue(r['forced_value_residence_bound'])
    def test_21_scalar_parameter_jets(self):self.assertEqual(len(numerical_checks()['parameter_rows']),21);self.assertTrue(numerical_checks()['checks']['parameter_jet_inverse_square_gain'])
    def test_22_uncut_exit_and_full_path_support(self):self.assertTrue(all(endpoint_checks()['checks'].values()));self.assertIn('uncut exit may be nonzero',protocol()['supports'])
    def test_23_physical_transfer_keeps_pressure_half_gain(self):
        h=F(1,200);a=F(3,2);p=physical_powers(h,a,2,3)
        self.assertEqual(p['pressure']-p['potential'],-1-h)
        self.assertEqual(p['velocity']-p['potential'],-F(1,2)-h/2)
    def test_24_invalid_physical_orders(self):
        for h,a,b in ((0,1,1),(F(1,100),1,1),(F(1,200),True,1),(F(1,200),1,-1)):
            with self.assertRaises(ValueError):physical_powers(h,1,a,b)
    def test_25_no_unconditional_or_infinite_harmonic_promotion(self):
        self.assertTrue(FLAGS['explicit_finite_sourced_ODE_existence_uniqueness_from_reviewed_frame']);self.assertFalse(FLAGS['unconditional_all_actual_pulse_family_verified']);self.assertFalse(FLAGS['infinite_harmonic_physical_summability_verified'])
    def test_26_frozen_record_and_source_bindings(self):
        frozen=json.loads((HERE/'evidence.json').read_text());validate(frozen);compare_records(self.record,frozen)
        self.assertEqual(len(self.record['source_hashes']),13);self.assertGreater(len(self.record['inherited_provenance']),3000)
    def test_27_parent_and_import_path_preserved(self):
        before=sys.path[:];self.assertTrue(all(parent_review()[1].values()));self.assertEqual(sys.path,before)
    def test_28_reject_scope_promotions(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)
    def test_29_reject_erased_history_or_obligation(self):
        for key in ('source_hashes','inherited_provenance','remaining_obligations','gates'):
            r=deepcopy(self.record)
            if isinstance(r[key],dict):r[key].pop(next(iter(r[key])))
            else:r[key].pop()
            with self.assertRaises(ValueError):validate(r)
    def test_30_reject_pressure_and_entrance_condition_tampering(self):
        for key in ('pressure','pulse_jet_condition'):
            r=deepcopy(self.record);r['protocol'][key]='zero'
            with self.assertRaises(ValueError):validate(r)
    def test_31_reject_failed_exact_control(self):
        r=deepcopy(self.record);r['reviews']['exact_moving_principal_pressure_controls']['checks']['source_normal_omission_detected']=False
        with self.assertRaises(ValueError):validate(r)
    def test_32_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
    def test_33_fresh_process_selects_local_audit(self):
        code='import sys;sys.path.insert(0,'+repr(str(HERE))+');before=sys.path[:];import audit;assert str(audit.HERE)=='+repr(str(HERE))+';assert before==sys.path'
        subprocess.run([sys.executable,'-c',code],check=True,capture_output=True,text=True)
    def test_34_all_gates_and_roundtrip(self):self.assertTrue(all(self.record['gates'].values()));validate(json.loads(json.dumps(self.record)))
if __name__=='__main__':unittest.main()
