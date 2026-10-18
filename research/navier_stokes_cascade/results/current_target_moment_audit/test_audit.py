"""Independent exact finite controls and explicit scope/provenance guards."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import target_algebra as a

class CurrentTargetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']['exact_current_state']
    def identity(self,key):self.assertTrue(self.r['identities'][key],key)
    def test_current_divergence(self):self.identity('actual_divergence')
    def test_angular_velocity_moment(self):self.identity('angular_velocity_moment')
    def test_axial_velocity_moment(self):self.identity('axial_velocity_moment')
    def test_moving_pressure_bump_unit_integral(self):self.identity('pressure_bump_unit')
    def test_pressure_obstruction_closed_formula(self):self.identity('current_P_closed_formula')
    def test_radial_second_moment_closed_formula(self):self.identity('current_Qg_closed_formula')
    def test_current_pressure_derivative(self):self.identity('pressure_derivative')
    def test_pressure_compact_at_both_edges(self):
        self.identity('pressure_left_compact');self.identity('pressure_right_compact')
    def test_complete_pressure_moment(self):self.identity('full_pressure_moment')
    def test_angular_viscous_flux_decomposition(self):self.identity('angular_flux_decomposition')
    def test_axial_viscous_flux_decomposition(self):self.identity('axial_flux_decomposition')
    def test_current_angular_integrated_residual(self):self.identity('angular_integrated_residual')
    def test_current_axial_integrated_residual(self):self.identity('axial_integrated_residual')
    def test_compensated_transports_have_zero_total(self):
        self.identity('angular_compensated_transport_zero_total');self.identity('axial_compensated_transport_zero_total')
    def test_signed_angular_flux_plus_primitive(self):self.identity('angular_target_flux_plus_primitive')
    def test_signed_axial_flux_plus_primitive(self):self.identity('axial_target_flux_plus_primitive')
    def test_signed_targets_compact_at_both_edges(self):
        for k in ('angular_target_left_compact','angular_target_right_compact','axial_target_left_compact','axial_target_right_compact'):self.identity(k)
    def test_second_current_source_linear_and_quadratic_expansion(self):self.identity('current_source_full_linear_quadratic_difference')
    def test_second_current_pressure_obstruction_difference(self):self.identity('current_pressure_obstruction_difference')
    def test_second_current_angular_flux_difference(self):self.identity('current_angular_flux_difference')
    def test_second_current_axial_flux_difference(self):self.identity('current_axial_flux_difference')
    def test_moving_bump_derivative_not_zero(self):self.assertTrue(self.r['controls']['nonzero_moving_pressure_bump_term'])
    def test_product_of_haar_means_is_not_mean_product(self):self.assertTrue(self.r['controls']['old_mean_correlation_not_product_of_Haar_means'])
    def test_both_old_new_correlations_survive(self):
        self.assertTrue(self.r['controls']['old_new_angular_correlation_nonzero']);self.assertTrue(self.r['controls']['old_new_axial_correlation_nonzero'])
    def test_all_six_covariances_and_both_targets_nonzero(self):
        for k in ('full_covariance_six_entries_nonzero','angular_target_nonzero','axial_target_nonzero'):self.assertTrue(self.r['controls'][k])
    def test_each_covariance_entry_has_nonzero_complete_equation_response(self):
        for k in ('nonzero_radial_covariance_pressure','nonzero_centrifugal_covariance_pressure',
                  'nonzero_axial_covariance_transport','nonzero_angular_axial_covariance_flux',
                  'nonzero_axial_square_covariance_flux','nonzero_radial_angular_covariance_target'):
            self.assertTrue(self.r['controls'][k],k)
    def test_independent_time_viscous_swirl_control_at_two_scales(self):
        for row in self.record['reviews']['swirl_time_viscous_controls']:
            self.assertTrue(all(v for k,v in row.items() if k!='epsilon'))
    def test_control_does_not_exclude_order_at_equality(self):self.assertTrue(a.ceiling(F(2),F(1))['order_not_excluded_by_this_control'])
    def test_higher_order_requires_extra_cancellation(self):self.assertFalse(a.ceiling(F(5,2),F(1))['order_not_excluded_by_this_control'])
    def test_later_canonical_rows_exceed_raw_ceiling(self):
        rows=self.record['reviews']['raw_state_ceiling'];self.assertTrue(rows[0]['order_not_excluded_by_this_control'])
        self.assertTrue(all(F(r['exponent_margin'])<0 for r in rows[1:]))
    def test_invalid_preserved_mean_order_rejected(self):
        for mu in (0,-1):
            with self.assertRaises(ValueError):a.ceiling(2,mu)
    def test_smooth_control_primitive_identity_finite_analogue(self):
        psi=a.R**3*a.F;f=a.s.diff(psi,a.R)/a.R**2
        self.assertTrue(a.zero(a.primitive(f,2)-psi/a.R**2));self.assertTrue(a.zero(a.integ(a.R**2*f)))
    def test_protocol_and_correct_parent(self):
        self.assertEqual(audit.protocol()['parent_commit'],'08caced6f10105d98faf59af4e251af26c3a2791')
        self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_all_historical_source_bytes_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_frozen_exact_replay(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_every_unproved_scope_promotion_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_current_identity_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['exact_current_state']['identities']['full_pressure_moment']=False
        with self.assertRaises(ValueError):audit.validate(r)
    def test_cancellation_obligation_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_integer_substituted_for_boolean_rejected(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_cli_replay_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)])
            audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
