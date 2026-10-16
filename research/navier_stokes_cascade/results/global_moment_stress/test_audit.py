#!/usr/bin/env python3
"""Equation, tail, signed stress, frozen provenance and scope regressions."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE, FLAGS, compare_records, generate, main, protocol, validate
from ledger import exterior_defects, heat_factor
import mpmath as mp


class GlobalMomentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=generate()

    def test_exact_five_basis_and_combined_subtractions(self):
        r=self.record['reviews']['exact_exterior_defects']
        self.assertEqual(r['rational_equation_checks'],90)
        self.assertTrue(all(r['checks'].values()))

    def test_zero_functional_defects_restore_exterior(self):
        self.assertEqual(exterior_defects(Q(1),Q(1,100),Q(100),Q(2),[(Q(0),Q(0))]*5),(0,0,0))

    def test_M_memory_changes_radial_velocity_even_when_U_zero(self):
        ds=[(Q(1),Q(0))]+[(Q(0),Q(0))]*4
        self.assertNotEqual(exterior_defects(Q(1,2),Q(1,100),Q(100),Q(2),ds)[0],0)

    def test_M_memory_term_in_angular_stress_must_not_be_deleted(self):
        ds=[(Q(1),Q(0))]+[(Q(0),Q(0))]*4
        self.assertNotEqual(exterior_defects(Q(1,2),Q(1,100),Q(100),Q(2),ds)[1],0)

    def test_renormalized_I_memory_keeps_sign_at_pole(self):
        ds=[(Q(0),Q(0))]*5;ds[1]=(Q(1),Q(0))
        self.assertGreater(exterior_defects(Q(1),Q(1,100),Q(100),Q(2),ds)[1],0)

    def test_J_derivative_memory_survives_at_eta_zero(self):
        ds=[(Q(0),Q(0))]*5;ds[2]=(Q(0),Q(1))
        self.assertEqual(exterior_defects(Q(0),Q(1,100),Q(100),Q(2),ds)[1],-1)

    def test_S_derivative_memory_survives_at_eta_zero(self):
        ds=[(Q(0),Q(0))]*5;ds[3]=(Q(0),Q(1))
        self.assertEqual(exterior_defects(Q(0),Q(1,100),Q(100),Q(2),ds)[2],-1)

    def test_pressure_value_match_at_one_angle_is_insufficient(self):
        ds=[(Q(0),Q(0))]*5;ds[4]=(Q(0),Q(1,10**30))
        self.assertEqual(exterior_defects(Q(0),Q(1,100),Q(100),Q(2),ds)[2],Q(-1,10**28))

    def test_pressure_memory_is_linear_in_weighted_axial_radius_X(self):
        ds=[(Q(0),Q(0))]*5;ds[4]=(Q(1),Q(0))
        a=exterior_defects(Q(1,2),Q(1,100),Q(100),Q(2),ds)[2]
        b=exterior_defects(Q(1,2),Q(1,100),Q(200),Q(2),ds)[2]
        self.assertEqual(b,2*a)

    def test_finite_h_is_not_dropped_in_energy_memory(self):
        ds=[(Q(0),Q(0))]*5;ds[3]=(Q(1),Q(0))
        self.assertGreater(exterior_defects(Q(1,2),Q(1,100),Q(100),Q(2),ds)[2],0)
        self.assertEqual(exterior_defects(Q(1,2),Q(0),Q(100),Q(2),ds)[2],0)

    def test_signed_heat_compensation_restores_three_totals(self):
        r=self.record['reviews']['signed_edit_ledger']
        self.assertEqual(r['heat_plus_compensation'],['0']*5)
        self.assertTrue(all(r['checks'].values()))

    def test_raw_I_is_never_used_as_finite_global_total(self):
        self.assertIn('H-Hpow',protocol()['renormalized_I'])
        self.assertEqual(protocol()['global_rows'][1],'renormalized_I')

    def test_power_subtraction_and_nonuniform_heat_integral_are_explicit(self):
        r=self.record['reviews']['integrable_renormalized_tail']
        self.assertEqual(r['h_to_zero_integrated_heat_defect'],'-2*k*(1-eta^2), not zero')
        self.assertTrue(all(r['checks'].values()))

    def test_second_angular_tail_derivative_keeps_both_terms(self):
        r=self.record['reviews']['integrable_renormalized_tail']
        self.assertGreater(Q(r['normalized_I_tail_bounds'][2]),Q(r['normalized_I_tail_bounds'][1]))

    def test_physical_scaling_and_minus_forward_plus_backward_sign(self):
        r=self.record['reviews']['physical_signed_scaling']
        self.assertEqual(r['rational_scaling_checks'],18)
        self.assertTrue(all(r['checks'].values()))

    def test_heat_at_zero_uses_finite_derivatives_not_division_by_d(self):
        with mp.workdps(60):
            h=mp.mpf('0.125')
            self.assertEqual(heat_factor(mp.mpf(0),h),1)
            self.assertEqual(heat_factor(mp.mpf(0),h,1),-h*(1+h))
            self.assertEqual(heat_factor(mp.mpf(0),h,2),h*(1+h)**2*(2+h))

    def test_independent_gamma_heat_factor_and_zero_stress_diagnostic(self):
        r=self.record['reviews']['labelled_heat_and_Fubini_diagnostics']
        self.assertTrue(all(r['checks'].values()))
        self.assertLess(max(x['zero_heat_stress_error'] for x in r['high_precision']['rows']),1e-30)

    def test_heat_pole_value_zero_does_not_erase_angular_derivative(self):
        row=self.record['reviews']['labelled_heat_and_Fubini_diagnostics']['high_precision']['rows'][-1]
        self.assertEqual(mp.mpf(row['tail_I']),0)
        self.assertEqual(row['pole_I_eta_over_Hpow'],-4.5)

    def test_Taylor_and_quadrature_errors_remain_labelled(self):
        r=self.record['reviews']['labelled_heat_and_Fubini_diagnostics']
        self.assertLess(r['high_precision']['integrated_Taylor_error_upper'],1e-55)
        self.assertLess(r['refinement_error'],1e-25)
        self.assertEqual(protocol()['diagnostic_quadrature_orders'],[24,36])

    def test_pressure_Fubini_sign_is_independently_checked(self):
        r=self.record['reviews']['labelled_heat_and_Fubini_diagnostics']
        self.assertLess(r['high_precision']['manufactured_Fubini_error'],1e-60)

    def test_exact_root_dependencies_and_521_inherited_hashes(self):
        self.assertTrue(self.record['gates']['exact_root_dependency_acceptance'])
        self.assertEqual(len(self.record['inherited_provenance']),521)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_frozen_record_replays(self):
        frozen=json.loads((HERE/'evidence.json').read_text())
        validate(frozen);compare_records(self.record,frozen)

    def test_full_profile_support_pulse_residual_and_blowup_promotion_rejected(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):
                    validate(r)

    def test_compact_input_hypotheses_cannot_be_removed(self):
        r=deepcopy(self.record);r['flags']['compact_input_hypotheses_remain_dependencies']=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_omitted_global_row_rejected(self):
        r=deepcopy(self.record);r['protocol']['global_rows'].pop()
        with self.assertRaises(ValueError):
            validate(r)

    def test_missing_signed_stress_gate_rejected(self):
        r=deepcopy(self.record);r['gates'].pop('physical_signed_scaling')
        with self.assertRaises(ValueError):
            validate(r)

    def test_failing_pressure_memory_control_rejected(self):
        r=deepcopy(self.record)
        r['reviews']['exact_exterior_defects']['checks']['one_angle_pressure_value_does_not_restore_axial_stress']=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_changed_frozen_provenance_rejected(self):
        r=deepcopy(self.record);r['inherited_provenance'][next(iter(r['inherited_provenance']))]=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_nonfinite_diagnostic_rejected(self):
        r=deepcopy(self.record);r['reviews']['labelled_heat_and_Fubini_diagnostics']['refinement_error']=float('nan')
        with self.assertRaises(ValueError):
            compare_records(self.record,r)

    def test_missing_remaining_support_obligation_rejected(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(1)
        with self.assertRaises(ValueError):
            validate(r)

    def test_CLI_cannot_overwrite_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':
    unittest.main()
