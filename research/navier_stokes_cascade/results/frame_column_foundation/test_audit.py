"""Exact projected-frame, phase, column and scientific-scope regression checks."""
from copy import deepcopy
from fractions import Fraction as F
import json
from math import isfinite
from pathlib import Path
import tempfile
import unittest
from audit import HERE,FLAGS,protocol,generate,validate,compare_records,main
from frame_budget import inputs,primitive,budgets,log_cone_budget,budget_review,all_band_review
from frame_identity import matrices,identity_review,eigen_review,column_review,gauge_review,mul,sub
from phase_checks import nearest_nonzero,phase_data,phase_review


class FrameColumnTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()
    def test_current_parent_and_protocol(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()),protocol())
        self.assertEqual(protocol()['parent_commit'],'255ae2e44292419235dd9b7e1a340f20e690dc39')
    def test_positive_floor_inputs(self):
        for kw in ({'gmin':0},{'bmin':0},{'cmin':0},{'u':0},{'Rmax':F(1,2)}):
            with self.assertRaises(ValueError):inputs(**kw)
    def test_unknown_input_and_inexact_budget_rejected(self):
        with self.assertRaises(ValueError):inputs(extra=1)
        for value in (0.1,True):
            with self.assertRaises(TypeError):inputs(Db=value)
    def test_primitive_rounding_term_is_retained(self):
        p=inputs(Db=0,Dbox=0,Mz=0);b=primitive(p)
        self.assertGreaterEqual(b['Dt'],1/p['Rmin']);self.assertEqual(b['DF'],0)
    def test_radial_reference_rounding_gradient_cost(self):
        p=inputs(Db=0);b=primitive(p);self.assertEqual(b['De'],p['gmax']/p['Rmin'])
    def test_operator_bound_includes_connections(self):
        p=inputs();b=budgets(p)
        self.assertEqual(b['CE'],b['kM']*b['DH']+b['S0']*b['S1']+p['lam_max']*b['Ds'])
        self.assertGreater(b['CE'],b['DH'])
    def test_damping_bound_keeps_calibration_factor_four(self):
        p=inputs();b=budgets(p);dn=b['Dr']+b['Dt']
        self.assertEqual(b['Cd'],4*dn*(2*p['bmax']*(1+b['Sref'])+dn))
    def test_value_floors_suffice_without_angular_inverse(self):
        p=inputs();b=budgets(p);self.assertEqual(b['beta'],p['bmin']/2)
        self.assertEqual(b['normal_acceptance'],max(1,2*b['Dt']/p['bmin'],b['Ds']))
    def test_Riccati_length_is_explicit(self):
        b=budgets(inputs());self.assertEqual(b['Riccati_acceptance'],8*b['CE']/b['mu'])
    def test_complete_column_ratio_contribution_not_zero(self):
        p=inputs();b=budgets(p);self.assertGreater(b['D0'],b['Ds']+b['Sref']*b['Da'])
        self.assertEqual(b['D1'],p['u']*(1+p['cmax']))
    def test_logarithmic_cone_length_is_finite(self):
        r=log_cone_budget(inputs());self.assertTrue(all(isfinite(float(x)) for x in r.values()))
        self.assertGreater(r['log_direction_constant'],1000)
    def test_invalid_cone_margin_rejected(self):
        for gamma in (0,2):
            with self.assertRaises(ValueError):log_cone_budget(inputs(),gamma)
    def test_exact_projected_identity_all_cases(self):
        r=identity_review();self.assertEqual(r['exact_frame_cases'],24);self.assertTrue(all(r['checks'].values()))
    def test_exact_original_reference_eigenvalues(self):self.assertTrue(all(eigen_review()['checks'].values()))
    def test_reference_mixing_operator_budget(self):
        p=inputs();b=budgets(p)
        self.assertEqual(b['CEref'],b['kMref']*(b['DH']+2*p['Fmax']*b['Ds'])+b['Sref']*p['u'])
    def test_reference_mixing_exact_generator_conjugacy(self):
        self.assertTrue(gauge_review()['checks']['exact_time_dependent_generator_conjugacy'])
    def test_entrance_must_be_transported(self):
        r=gauge_review();self.assertTrue(r['checks']['transported_entrance_preserves_physical_pulse'])
        self.assertTrue(r['checks']['resetting_growing_coordinate_changes_initial_physical_pulse'])
        self.assertNotEqual(r['transported_repository_entrance'][1],['0'])
    def test_original_tilt_obstructs_equal_scales_before_rounding(self):
        self.assertGreater(F(gauge_review()['pre_rounding_entrance_squared_scale_gap']),0)
    def test_reference_and_actual_frame_forcing_controls(self):self.assertTrue(all(gauge_review()['checks'].values()))
    def test_protocol_disallows_silent_frame_replacement(self):
        self.assertIn('not reset eplus',protocol()['frame_convention'])
        self.assertFalse(FLAGS['canonical_growing_entrances_in_both_frames_identified_as_same_pulse'])
    def test_moving_normal_omission_detected(self):self.assertTrue(identity_review()['checks']['omitting_moving_normal_pressure_detected'])
    def test_frame_prime_omission_detected(self):self.assertTrue(identity_review()['checks']['omitting_Bprime_detected'])
    def test_invalid_exact_frame_rejected(self):
        with self.assertRaises(ValueError):matrices((1,1),0,1,-1,1,(1,1),0,0)
    def test_left_inverse_is_exact(self):
        r=matrices((F(3,5),F(-4,5)),F(3,4),F(5,4),F(-1,2),1,(F(-8,5),F(-6,5)),F(1,50),F(1,100))
        self.assertEqual(mul(r['left'],r['B']),[[1,0],[0,1]])
        self.assertEqual(r['direct'],r['reduced'])
    def test_complete_principal_column_identity(self):self.assertTrue(all(column_review()['checks'].values()))
    def test_ratio_bound_is_sharp_at_negative_one_eighth(self):
        q=F(-1,8);self.assertEqual(abs((1-q)/(1+q)-1),F(16,7)*abs(q))
    def test_original_nonzero_rounding_error(self):
        self.assertNotEqual(nearest_nonzero(0,10),0);self.assertEqual(abs(nearest_nonzero(0,10)),F(1,10))
    def test_invalid_carrier_rejected(self):
        for k in (0,-1,F(3,2)):
            with self.assertRaises(ValueError):nearest_nonzero(1,k)
    def test_full_phase_calibration_and_accepted_domain(self):
        d=phase_data(10**14,1,F(1,2));self.assertEqual(d['epsilon']*d['k']**2,1)
        self.assertLessEqual(d['S']**2*(d['epsilon']+d['epsilon']**2+F(1,d['k'])),1)
        self.assertGreaterEqual(d['L'],d['budget']['Riccati_acceptance'])
    def test_full_axial_normal_and_derivative_kept(self):
        d=phase_data(10**14,1,1);self.assertNotEqual(d['n'][2],d['pz']);self.assertNotEqual(d['nprime'][2],0)
    def test_phase_endpoints_and_signs(self):
        r=phase_review();self.assertEqual(len(r['rows']),12);self.assertTrue(all(r['checks'].values()))
        self.assertEqual({x['pulse_fraction'] for x in r['rows']},{'0','1/2','1'})
    def test_no_manufactured_schedule_promoted(self):self.assertTrue(all_band_review()['checks']['actual_schedule_uses_its_existing_acceptance_not_manufactured_epsilon'])
    def test_full_pressure_memory_and_fixed_geometry(self):
        self.assertIn('Pi1/Pi2',protocol()['memory']);self.assertIn('never forward-inverted',protocol()['memory'])
        self.assertIn('C_setup^-20',protocol()['fixed_choices']);self.assertIn('Bprime',protocol()['connection'])
    def test_frozen_evidence_and_all_sources(self):
        old=json.loads((HERE/'evidence.json').read_text());validate(old);compare_records(self.record,old)
        self.assertEqual(len(old['source_hashes']),19);self.assertEqual(old['inherited_provenance']['records'],77)
    def test_scope_promotion_rejected(self):
        for key,value in FLAGS.items():
            if not value:
                r=deepcopy(self.record);r['flags'][key]=True
                with self.assertRaises(ValueError):validate(r)
    def test_erased_history_and_limits_rejected(self):
        for key in ('inherited_provenance','remaining_obligations'):
            r=deepcopy(self.record)
            if key=='inherited_provenance':r[key]['bindings'].pop(next(iter(r[key]['bindings'])))
            else:r[key]=[]
            with self.assertRaises(ValueError):validate(r)
    def test_numerical_identity_corruption_not_hidden_by_float_tolerance(self):
        r=deepcopy(self.record);r['reviews']['original_full_phase_high_precision_controls']['maximum_direct_reduced_identity_error']=1e-20
        with self.assertRaises(ValueError):validate(r)
    def test_algebra_and_protocol_corruption_rejected(self):
        r=deepcopy(self.record);r['protocol']['column']='q=0'
        with self.assertRaises(ValueError):validate(r)
    def test_failed_phase_gate_rejected(self):
        r=deepcopy(self.record);r['reviews']['original_full_phase_high_precision_controls']['checks']['explicit_damping_budget_in_phase_controls']=False
        with self.assertRaises(ValueError):validate(r)
    def test_frozen_output_not_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
    def test_json_roundtrip(self):validate(json.loads(json.dumps(self.record)))
if __name__=='__main__':unittest.main()
