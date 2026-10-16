#!/usr/bin/env python3
"""Physical derivative losses, graph terms, domain margins, and scope controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from audit import HERE,FLAGS,compare_records,generate,main,validate
from transfer import losses,primary_powers,common_loss,first_band,est


class PhysicalPulseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=generate()

    def test_exact_graph_and_physical_power_controls(self):
        r=self.record['reviews']['exact_graph_and_physical_power_ledger']
        self.assertEqual(r['exact_power_checks'],210)
        self.assertEqual(r['exact_cover_checks'],9)
        self.assertTrue(all(r['checks'].values()))

    def test_spatial_and_time_harmonic_losses(self):
        l=losses(Q(1,200))
        self.assertEqual(l['spatial'],Q(201,400))
        self.assertEqual(l['time'],Q(403,400))

    def test_harmonic_losses_dominate_all_graph_and_frame_losses(self):
        for h in (Q(1,10**12),Q(1,200),Q(99,10000)):
            l=losses(h)
            self.assertGreater(l['spatial'],l['coefficient_radial'])
            self.assertGreater(l['spatial'],l['coefficient_axial'])
            self.assertGreater(l['spatial'],l['frame'])
            self.assertGreater(l['time'],l['coefficient_time'])

    def test_primary_velocity_pressure_potential_zero_order_powers(self):
        self.assertEqual(primary_powers(Q(1,200),0,0),dict(velocity=Q(-201,400),pressure=Q(-201,200),potential=Q(0)))

    def test_primary_mixed_physical_derivative_powers(self):
        self.assertEqual(primary_powers(Q(1,200),2,1),dict(velocity=Q(-503,200),pressure=Q(-1207,400),potential=Q(-161,80)))

    def test_common_potential_loss_has_no_stage_argument(self):
        self.assertEqual(common_loss(Q(1,200),3),Q(101,100)+4*Q(403,400))

    def test_invalid_h_and_derivative_orders_rejected(self):
        for h in (0,-1,Q(1,100)):
            with self.assertRaises(ValueError):
                losses(h)
        for a,b in ((-1,0),(0,-1),(True,0),(0,1.5)):
            with self.assertRaises(ValueError):
                primary_powers(Q(1,200),a,b)

    def test_first_band_requires_a_real_background_domain(self):
        for qstar in (0,-1,2,float('inf')):
            with self.assertRaises(ValueError):
                first_band(Q(1,200),qstar,25)

    def test_short_pulse_invalid_geometry_and_radius_rejected(self):
        for length,geo,r in ((24,16,est.RADIUS),(25,15,est.RADIUS),(25,16,Q(1,25))):
            with self.assertRaises(ValueError):
                first_band(Q(1,200),'.001',length,geo,r)

    def test_all_later_domain_length_and_phase_checks(self):
        r=self.record['reviews']['common_band_domain_and_enlargement']
        self.assertEqual(r['diagnostic_acceptance_checks'],18)
        self.assertTrue(all(r['checks'].values()))

    def test_enlarged_domain_retains_four_Q_margin(self):
        self.assertEqual(self.record['protocol']['enlarged_background_margin'],'4Q<qstar')
        with mp.workdps(60):
            qstar=mp.power(2,-1000000)
            row=first_band(Q(1,200),qstar,25,10**6)
            self.assertLess(4*mp.power(2,-row['first_band']),qstar)

    def test_merely_Q_below_qstar_can_leave_the_domain(self):
        qstar=Q(1,2**20);Qband=3*qstar/4
        self.assertLess(Qband,qstar)
        self.assertGreater(2*Qband,qstar)

    def test_working_domain_does_not_invoke_omitted_early_bands(self):
        for row in self.record['reviews']['common_band_domain_and_enlargement']['rows']:
            self.assertEqual(row['working_log2_q'],-row['first_band']-2)

    def test_operator_control_is_independent_direct_pullback(self):
        r=self.record['reviews']['independent_physical_graph_and_phase_derivatives']
        self.assertLess(r['maximum_graph_and_harmonic_derivative_error'],1e-60)
        self.assertTrue(all(r['checks'].values()))

    def test_radial_graph_cancels_only_in_phase(self):
        r=self.record['reviews']['independent_physical_graph_and_phase_derivatives']
        self.assertLess(r['maximum_radial_phase_cancellation_error'],1e-60)
        self.assertTrue(all(row['omitted_radial_graph_error']>1 for row in r['rows']))

    def test_fast_auxiliary_time_not_discarded(self):
        r=self.record['reviews']['independent_physical_graph_and_phase_derivatives']
        self.assertTrue(all(row['omitted_fast_time_error']>1 for row in r['rows']))

    def test_operator_fixtures_are_not_accepted_actual_pulses(self):
        self.assertTrue(self.record['reviews']['independent_physical_graph_and_phase_derivatives']['checks']['short_operator_fixture_bands_not_claimed_as_accepted_pulses'])
        self.assertTrue(self.record['protocol']['manufactured_derivative_lift_translation_not_actual_map_change'])

    def test_discrete_choices_and_Q_not_differentiated(self):
        self.assertTrue(self.record['protocol']['fixed_labels_carriers_covers_Q_epsilon_S_during_differentiation'])

    def test_physical_Gaussian_tails_absorb_all_fixed_derivative_losses(self):
        r=self.record['reviews']['physical_flat_tails_and_closed_supports']
        self.assertEqual([row['derivative_order'] for row in r['rows']],[0,2,6])
        self.assertTrue(all(row['maximum_log_tail_over_q_seventh']<0 for row in r['rows']))

    def test_finite_endpoint_envelope_does_not_supply_zero_extension(self):
        self.assertTrue(self.record['reviews']['physical_flat_tails_and_closed_supports']['checks']['endpoint_envelope_is_positive_not_a_zero_extension'])
        self.assertTrue(FLAGS['sourced_inverse_requires_input_class_and_smooth_endpoint_extension'])

    def test_fixed_potential_and_pressure_cut_before_curl(self):
        self.assertTrue(self.record['reviews']['physical_flat_tails_and_closed_supports']['checks']['potentials_and_pressures_cut_before_curl'])

    def test_closed_supports_and_one_sided_boundary_preserved(self):
        r=self.record['reviews']['physical_flat_tails_and_closed_supports']
        self.assertEqual(r['supports']['transverse'],'abs(xi)<=r0/2')
        self.assertEqual(r['supports']['pulse'],'abs(v/L-1/2)<=2/5')
        self.assertTrue(r['checks']['one_sided_tau_zero_jets_not_zero_extended'])

    def test_primary_covariance_not_absolutely_decaying(self):
        self.assertLess(Q(self.record['reviews']['exact_graph_and_physical_power_ledger']['physical_covariance_remainder_power']),0)

    def test_full_nonlinear_source_estimates_remain_open(self):
        self.assertFalse(FLAGS['all_nonlinear_generated_source_class_estimates_verified'])
        self.assertFalse(self.record['reviews']['physical_flat_tails_and_closed_supports']['full_nonlinear_generated_source_estimates_verified'])

    def test_both_mean_pressure_tails_retained(self):
        self.assertTrue(FLAGS['first_two_mean_pressure_tails_and_complete_covariance_obligation_retained'])

    def test_parent_background_and_580_inherited_hashes_replayed(self):
        self.assertTrue(self.record['gates']['audited_background_inputs_replayed'])
        self.assertEqual(len(self.record['inherited_provenance']),580)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_frozen_record_reproduces(self):
        frozen=json.loads((HERE/'evidence.json').read_text())
        validate(frozen);compare_records(self.record,frozen)

    def test_unverified_source_mean_residual_and_blowup_promotions_rejected(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):
                    validate(r)

    def test_analytic_dependency_cannot_be_removed(self):
        r=deepcopy(self.record);r['flags']['relies_on_replayed_background_geometry_and_pulse_lemmas']=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_changed_physical_normalization_rejected(self):
        r=deepcopy(self.record);r['protocol']['physical_rescalings'][1]='Q^-A pressure'
        with self.assertRaises(ValueError):
            validate(r)

    def test_missing_gate_or_obligation_rejected(self):
        for section in ('gates','remaining_obligations'):
            r=deepcopy(self.record)
            if section=='gates':
                r[section].pop(next(iter(r[section])))
            else:
                r[section].pop()
            with self.assertRaises(ValueError):
                validate(r)

    def test_changed_historical_and_current_sources_rejected(self):
        for section in ('source_hashes','inherited_provenance'):
            r=deepcopy(self.record);r[section][next(iter(r[section]))]=False
            with self.assertRaises(ValueError):
                validate(r)

    def test_nonfinite_diagnostic_rejected(self):
        r=deepcopy(self.record);r['reviews']['independent_physical_graph_and_phase_derivatives']['maximum_graph_and_harmonic_derivative_error']=float('nan')
        with self.assertRaises(ValueError):
            compare_records(self.record,r)

    def test_CLI_does_not_overwrite_frozen_record(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':
    unittest.main()
