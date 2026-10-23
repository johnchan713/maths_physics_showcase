"""Exact mixed-jet controls and frozen analytical-scope regressions."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review,manufactured,rational_box,original_fixture,x
from bounds import P,inverse,CanonicalEntranceBudget,from_original,root_wave_join

class CanonicalEntranceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        vals=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(vals,key);self.assertTrue(all(vals),key)
    def test_exact_constraint_and_full_homogeneous_equation(self):self.control('moving_normal_constraint_exact');self.control('manufactured_homogeneous_equation_exact')
    def test_ambient_normal_response_and_tangent_agreement(self):self.control('ambient_normal_response');self.control('ambient_extension_agrees')
    def test_homogeneous_pressure_is_nonzero(self):self.control('homogeneous_pressure_is_nonzero')
    def test_moving_normal_pressure_cannot_be_discarded(self):self.control('dropping_moving_normal_changes_pressure')
    def test_potential_formula(self):self.control('exact_potential')
    def test_zero_entrance_inverse_cannot_replace_primary(self):self.control('zero_source_zero_entry_inverse')
    def test_canonical_reference_entrance_is_preserved(self):self.control('reference_entrance_mixing');self.control('canonical_entrance_is_er')
    def test_actual_wa_family_is_distinct(self):self.control('inherited_actual_wa_entrance_not_identified')
    def test_manufactured_tangent_majorant(self):self.control('tangent_growth_exponent');self.control('conservative_tangent_constant')
    def test_entrance_Gaussian_derivatives_not_dropped(self):self.control('entrance_Gaussian_first_derivative');self.control('entrance_Gaussian_second_derivative');self.control('first_initial_jet_contains')
    def test_Bell_and_initial_Leibniz_recurrences(self):self.control('entrance_Bell');self.control('initial_Leibniz')
    def test_velocity_all_mixed_jets_on_whole_box(self):self.control('velocity_complete_mixed_jets_bound')
    def test_pressure_all_mixed_jets_on_whole_box(self):self.control('pressure_complete_mixed_jets_bound')
    def test_potential_all_mixed_jets_on_whole_box(self):self.control('potential_complete_mixed_jets_bound')
    def test_path_derivative_part_of_induction(self):self.control('pulse_derivative_uses_Ahat')
    def test_varying_tangential_inverse_sqrt(self):self.control('inverse_sqrt_differential_identity');self.control('varying_tangential_inverse_sqrt_bound')
    def test_ambient_floor_cannot_replace_tangential_floor(self):self.control('full_normal_floor_alone')
    def test_production_factory_retains_original_tail(self):self.control('original_normal_and_tangent_gates');self.control('original_same_band_and_length')
    def test_production_tangential_and_full_floors_have_same_proof(self):self.control('original_refined_floor_proved')
    def test_primary_source_zero_does_not_change_forced_pressure(self):self.control('primary_source_exactly_zero');self.assertTrue(audit.FLAGS['complete_forced_pressure_source_and_global_memory_obligations_retained'])
    def test_root_callbacks_and_value_domain_unchanged(self):self.control('root_H_target_curl_loghalf');self.control('root_value_gaps_mass_radius_width')
    def test_root_cutoff_jets_are_retained(self):self.control('root_wave_slot_is_full_cutoff_jet');self.control('cutoff_derivative_contributes')
    def test_derived_wave_drives_root_and_signed_recipes(self):self.control('joined_root_amplitude_and_signed_difference')
    def test_invalid_denominator_control_is_rejected(self):
        with self.assertRaises(ValueError):rational_box(1/(1-x))
    def test_invalid_primary_and_tangential_floor_are_rejected(self):
        base=manufactured()[-2]
        for args in ((0,1,0),(1,0,0),(1,1,-1)):
            with self.assertRaises(ValueError):CanonicalEntranceBudget(base,tangential_inverse_square=args[0],cmax=args[1],u=args[2])
        with self.assertRaises(ValueError):CanonicalEntranceBudget(inverse.fixture(harmonic=2),tangential_inverse_square=1,cmax=1,u=0)
        with self.assertRaises(TypeError):CanonicalEntranceBudget(base,tangential_inverse_square=0.5,cmax=1,u=0)
    def test_invalid_orders_dimensions_and_joins_are_rejected(self):
        h=manufactured()[-1]
        for j in (-1,0.5):
            with self.assertRaises((ValueError,TypeError)):h.velocity(j)
        with self.assertRaises(ValueError):h.wave(0,lambda j:P(1),0)
        with self.assertRaises(TypeError):root_wave_join(object(),h,lambda j:P(1))
    def test_unaccepted_band_is_not_reselected(self):
        p,sp,ell,kw=original_fixture()
        with self.assertRaises(ValueError):from_original(p,sp,16,**kw)
    def test_all_controls_and_historical_sources_pass(self):self.assertTrue(all(self.record['gates'].values()));self.assertTrue(all(self.r['checks'].values()))
    def test_eight_historical_records_and_frozen_local_sources(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),8)
        self.assertEqual(self.record['source_hashes'],audit.sources())
    def test_scope_does_not_promote_actual_suprema_or_blowup(self):
        for k in ('canonical_and_inherited_physical_primaries_identified','actual_background_suprema_Haar_column_jets_or_all_generated_sources_certified','unconditional_pulse_complete_residual_infinite_field_or_blowup_proved'):
            self.assertIs(self.record['flags'][k],False)
        self.assertEqual(len(self.record['remaining_obligations']),4)
    def test_tampered_recipe_scope_or_provenance_is_rejected(self):
        for field in ('flags','source_hashes','review'):
            r=deepcopy(self.record);r[field]['tampered']=True
            with self.assertRaises(ValueError):audit.validate(r)
    def test_frozen_record_matches_current_sources(self):
        audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()))
    def test_cli_replays_and_refuses_evidence_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            dest=Path(d)/'replay.json';audit.main(['--verify-record',str(Path(__file__).parent/'evidence.json'),'--output',str(dest)])
            self.assertEqual(json.loads(dest.read_text()),self.record)
            with self.assertRaises(SystemExit):audit.main(['--output',str(dest)])
if __name__=='__main__':unittest.main()
