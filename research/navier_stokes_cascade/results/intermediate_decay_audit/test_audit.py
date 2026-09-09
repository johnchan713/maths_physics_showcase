#!/usr/bin/env python3
"""Stress bounds, independent source checks and controls for lost small terms."""
import copy
from fractions import Fraction
import json
import math
import unittest
import mpmath as mp
from bounds import decay_bound,finite_cone_bound,elementary_checks,axial,HERE
from reference import Reference
from audit import validate,compare_record,relative,pole_control


class BoundsTests(unittest.TestCase):
    def setUp(self):
        mp.mp.dps=110

    def test_elementary_inequalities(self):
        self.assertTrue(all(elementary_checks().values()))

    def test_uniform_decay_certificate(self):
        b=decay_bound(64)
        self.assertLess(b['universal_lambda_w_squared_upper'],mp.mpf('1e-42'))
        self.assertGreater(b['second_A24_margin_lower'],mp.mpf('1.99'))
        self.assertLess(b['second_A24_margin_lower'],2)

    def test_parameter_order_and_reserved_patches(self):
        b=decay_bound(64)
        self.assertLess(b['log_h_interval'][1],b['log_lambda_interval'][0])
        self.assertLess(b['log_h_interval'][1],-2*b['Td_interval'][1])
        self.assertGreater(b['Tw_interval'][0],25)

    def test_explicit_radius_floor(self):
        b=finite_cone_bound(axial.axial_bound(64,256,50))
        self.assertTrue(b['all_pass'])
        self.assertGreater(b['axial_transformed_gap_lower'],mp.mpf('.38'))
        self.assertGreater(b['intermediate_quadratic_gap_lower'],mp.mpf('.3'))

    def test_transformed_axial_cone_identity(self):
        for b,w in [(Fraction(-2,3),Fraction(7,5)),(Fraction(1,100),Fraction(1000)),(Fraction(0),Fraction(99))]:
            self.assertEqual(b*b/2*(w+b/2)**2,(b*w+b*b/2)**2/2)

    def test_intermediate_quadratic_cone_test(self):
        for lam in [mp.mpf('.01'),mp.mpf('1e-60')]:
            ps1=mp.mpf(5);a=2+2*lam;w=mp.sqrt(mp.mpf('1e-42')/lam)
            self.assertGreater(ps1,a)
            self.assertLess((a-2)*(ps1*w)**2,2*(ps1-a)**2)

    def test_future_tail_power_bound(self):
        for lam in map(mp.mpf,['.001','.1','.5','.99']):
            self.assertLess(mp.exp(-13/lam),lam**8)

    def test_wide_exponent_still_loses_added_shear(self):
        log_lam=-4*(mp.exp(64)+10)
        lam=mp.exp(log_lam)
        self.assertGreater(lam,0)
        self.assertEqual(math.exp(float(log_lam)),0)
        self.assertEqual(2+2*lam-2,0)
        self.assertTrue(mp.isfinite(mp.log(2)+log_lam))

    def test_invalid_certificate_parameters(self):
        for args in [(3,60),(4.5,60),(64,40)]:
            with self.assertRaises(ValueError):decay_bound(*args)


class ReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=320
        cls.r=Reference()

    def setUp(self):
        mp.mp.dps=320

    def test_energy_scale_from_composed_profile(self):
        r=self.r;e=mp.mpf('.5')
        for y in (0,3*r.T,r.Tw):
            row=r.state('power',y,e)
            # Independently compose P*, the first ramp, axial decay, and U=0 stages.
            E=mp.exp(r.T+1)*mp.exp(mp.mpf('-.2'))*mp.exp(-r.T/2)
            E*=mp.exp(-(1+y)/2-r.lam*(y+mp.mpf('.5')))/(1+e*e)
            self.assertLess(relative(row['E_squared'],E*E),mp.mpf('1e-200'))

    def test_omitted_axial_energy_decay_is_detectable(self):
        row=self.r.state('power',0,mp.sqrt(self.r.delta))
        self.assertLess(row['lambda_w_squared'],mp.mpf('1e-42'))
        self.assertGreater(row['lambda_w_squared']/self.r.delta,mp.mpf('1e-42'))

    def test_stage_boundary_continuity(self):
        for eta in (0,mp.sqrt(self.r.delta),mp.mpf('.5'),1):
            a=self.r.state('ramp',1,eta);b=self.r.state('power',0,eta)
            for k in ('Q','v','E_squared'):
                self.assertLess(relative(a[k],b[k]),mp.mpf('1e-190'))

    def test_source_equations_at_late_and_small_angle_points(self):
        for y in (3*self.r.T,self.r.Tw):
            for e in (0,mp.sqrt(self.r.lam),1):
                a=self.r.state('power',y,e);b=self.r.source_power(y,e)
                for k in ('Q','v'):
                    self.assertLess(relative(a[k],b[k]),mp.mpf('1e-100'))

    def test_independent_pressure_and_flow_IVPs(self):
        e=mp.sqrt(self.r.delta)
        for row in self.r.independent_ramp(e,[0,.5,1]):
            actual=self.r.state('ramp',row['x'],e)
            self.assertLess(relative(actual['Q']/row['scale'],row['Q_over_scale']),mp.mpf('1e-9'))
            self.assertLess(relative(actual['v'],row['v']),mp.mpf('1e-9'))
            self.assertLess(relative(self.r.ramp_factors(row['x'])['Z'],row['Z']),mp.mpf('1e-9'))

    def test_nonzero_pole_increment(self):
        row=pole_control(self.r)
        self.assertLess(row['increment'],0)
        self.assertLess(abs(row['increment_over_h_y']+2),mp.mpf('1e-30'))
        self.assertLess(row['relative_gap'],mp.mpf('1e-40'))

    def test_omitted_h_and_binary64_erase_increment(self):
        row=pole_control(self.r)
        self.assertEqual(row['omitted_h_increment'],0)
        self.assertEqual(row['binary64_increment'],0)

    def test_axis_memory_is_retained(self):
        for stage,y in [('ramp',0),('power',3*self.r.T),('power',self.r.Tw)]:
            row=self.r.state(stage,y,0)
            self.assertGreater(row['Q'],0)
            self.assertEqual(row['lambda_w_squared'],0)

    def test_eta_symmetry(self):
        for e in (mp.sqrt(self.r.lam),mp.mpf('.5')):
            a=self.r.state('power',3*self.r.T,e);b=self.r.state('power',3*self.r.T,-e)
            for k in ('Q','v','lambda_w_squared'):self.assertEqual(a[k],b[k])

    def test_positive_future_error_budget(self):
        row=self.r.state('power',self.r.Tw,1)
        self.assertGreater(row['v_tail_error_upper'],0)
        self.assertLess(row['v_tail_error_upper'],mp.mpf('1e-800'))

    def test_reject_unsafe_precision_and_unsupported_scale(self):
        with mp.workdps(100):
            with self.assertRaises(ValueError):Reference()
        with self.assertRaises(ValueError):Reference(md=64)

    def test_coordinate_domains(self):
        for args in [('other',0,0),('ramp',-1,0),('power',self.r.Tw+1,0),('ramp',0,2)]:
            with self.assertRaises(ValueError):self.r.state(*args)


class RecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=json.loads((HERE/'evidence.json').read_text())

    def test_pinned_record(self):validate(self.record)

    def test_scope_promotion_is_rejected(self):
        row=copy.deepcopy(self.record);row['status']='global-solution-proved'
        with self.assertRaises(ValueError):validate(row)

    def test_relabelled_gate_is_rejected(self):
        row=copy.deepcopy(self.record);row['gates']['fake_gate']=row['gates'].pop('integrated_energy_scale')
        with self.assertRaises(ValueError):validate(row)

    def test_pole_increment_erasure_is_rejected(self):
        row=copy.deepcopy(self.record);row['pilot']['pole']['increment']='0'
        with self.assertRaises(ValueError):validate(row)

    def test_old_failure_erasure_is_rejected(self):
        row=copy.deepcopy(self.record);row['old_lambda_failure']['lambda_w_squared']='.1'
        with self.assertRaises(ValueError):validate(row)

    def test_provenance_change_is_rejected(self):
        row=copy.deepcopy(self.record);row['provenance'][next(iter(row['provenance']))]='0'*64
        with self.assertRaises(ValueError):validate(row)

    def test_tiny_increment_comparison_has_no_unit_floor(self):
        with self.assertRaises(ValueError):compare_record('1e-220','0','/pilot/pole/increment')


if __name__=='__main__':unittest.main()
