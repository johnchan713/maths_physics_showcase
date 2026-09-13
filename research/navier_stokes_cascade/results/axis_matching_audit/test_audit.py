#!/usr/bin/env python3
"""Tests of identities, nonlinear corrections, angular data, and scope limits."""
import copy
import json
from pathlib import Path
import unittest
import mpmath as mp
from moments import Jet,MomentMap,fixture,zeta,relative
from annulus import (Annulus,ideal,accept_entry_bounds,physical_to_rows,
                     rows_to_physical,radius_requirements)
from bounds import certificate

HERE=Path(__file__).resolve().parent


class BoundsTests(unittest.TestCase):
    def test_outward_constants(self):
        c=certificate(60)
        self.assertTrue(all(c['checks'].values()))
        self.assertGreater(c['Pc_lower'],2)
        self.assertLess(c['old_XR100_Pc_upper'],2)

    def test_actual_axis_is_not_inferred(self):
        c=certificate(60)
        self.assertFalse(c['actual_axis_entry_verified'])
        self.assertFalse(c['full_admissible_cone_realized'])
        self.assertFalse(c['blowup_verified'])

    def test_uniform_entry_bounds_required(self):
        self.assertTrue(accept_entry_bounds(['1e-16']*5,'1e-16'))
        self.assertFalse(accept_entry_bounds(['1e-15']*5,'1e-16'))
        self.assertFalse(accept_entry_bounds(['nan']*5,0))
        with self.assertRaises(ValueError):accept_entry_bounds([0]*4,0)

    def test_radius_alone_does_not_bound_the_core(self):
        r=radius_requirements(120,mp.log(16),1000)
        self.assertTrue(r['radius_large_enough'])
        self.assertTrue(r['transition_before_restoration'])
        self.assertFalse(r['core_moments_verified'])
        self.assertFalse(accept_entry_bounds([0]*5,'1e-10'))

    def test_transition_length_is_fixed_before_amplitude(self):
        a=radius_requirements(50,mp.log(16),1000)
        b=radius_requirements(120,mp.log(16),1000)
        self.assertFalse(a['transition_before_restoration'])
        self.assertTrue(b['transition_before_restoration'])

    def test_shrinking_j_alone_loses_axis_alternative(self):
        j=mp.mpf('1e-18');old_sigma=mp.mpf('.002')
        self.assertLess(j*j/(j*j+old_sigma**2),mp.mpf('.99'))
        sigma=j/100
        self.assertGreater(j*j/(j*j+sigma**2),mp.mpf('.99'))
        # eta=0 belongs to |Z|<=j P^2/100 for the P=100 control.
        self.assertLess(mp.mpf('4.505')*j,j*100**2/100)


class ConstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=100
        cls.matrix=MomentMap(32,8)
        cls.eta=mp.mpf('.5');cls.P=mp.mpf(16)
        cls.data=fixture(cls.eta)
        cls.solution=cls.matrix.solve(cls.data,cls.eta,cls.P)
        cls.annulus=Annulus(cls.matrix,'.005',cls.P)

    def test_value_derivative_product_rule(self):
        e=Jet('.5',1);v=(1+e**2)**-2
        self.assertAlmostEqual(float(v.d),float(-4*e.v/(1+e.v**2)**3),places=14)

    def test_physical_moment_adapter_keeps_cross_derivatives(self):
        rows=self.data['entry']
        for radius in (mp.mpf(10000),mp.mpf('1e30')):
            physical=rows_to_physical(rows,self.eta,self.P,radius)
            recovered=physical_to_rows(physical,self.eta,self.P,radius)
            self.assertLess(max(relative(a.v,b.v) for a,b in zip(rows,recovered)),mp.mpf('1e-90'))
            self.assertLess(max(relative(a.d,b.d) for a,b in zip(rows,recovered)),mp.mpf('1e-90'))

    def test_adapter_rejects_missing_derivatives(self):
        with self.assertRaises(ValueError):physical_to_rows([0]*5,.5,16,10000)

    def test_disjoint_supports_make_block_structure(self):
        for t in map(mp.mpf,['.15','.3','.5','.7','.85']):
            self.assertEqual(sum(self.matrix.bump(j,t)!=0 for j in range(5)),1)
        B=self.matrix.B
        self.assertTrue(all(B[i,j]==0 for i in range(2) for j in range(2,5)))
        self.assertTrue(all(B[i,j]==0 for i in range(2,5) for j in range(2)))

    def test_exact_bump_mass_controls_weighted_integrals(self):
        self.assertLess(abs(self.matrix.integral(0,0)-mp.mpf('.08')),mp.mpf('1e-30'))
        for j in range(5):
            r=mp.mpf('1.6');c=self.matrix.centers[j];w=self.matrix.width
            value=self.matrix.integral(j,r)
            self.assertGreater(value,w*mp.exp(r*(c-w/2)))
            self.assertLess(value,w*mp.exp(r*(c+w/2)))

    def test_inverse_matches_outward_bound(self):
        norm=max(mp.fsum(abs(v) for v in self.matrix.inverse[i,:]) for i in range(5))
        self.assertLess(norm,1000)
        self.assertLess(mp.norm(self.matrix.inverse*self.matrix.B-mp.eye(5)),mp.mpf('1e-90'))

    def test_duplicate_bumps_rejected(self):
        class Degenerate(MomentMap):U_CENTERS=('.3','.3')
        with self.assertRaises(ValueError):Degenerate(16,1)

    def test_nonlinear_root_closes_each_row(self):
        r=self.matrix.ledger(self.data,self.solution,self.eta,self.P,1)
        scale=max(abs(v) for v in self.solution['target'])
        self.assertLess(max(abs(v.v) for v in r)/scale,mp.mpf('1e-85'))
        self.assertLess(max(abs(v.d) for v in r)/scale,mp.mpf('1e-85'))

    def test_linear_root_has_nonzero_pressure_error(self):
        c=self.matrix.inverse*mp.matrix(self.solution['target'])
        q=self.matrix.quadratic(c,c,self.eta,self.P)
        self.assertGreater(q[4],0)
        residual=self.matrix.B*c+q-mp.matrix(self.solution['target'])
        self.assertLess(relative(residual[4],q[4]),mp.mpf('1e-60'))

    def test_quadratic_grade_and_remainder_are_retained(self):
        g=self.matrix.graded(self.solution['target'],self.eta,self.P,4)
        self.assertGreater(max(abs(v) for v in g['components'][1]),0)
        self.assertGreater(g['remainder'],0)
        summed=[mp.fsum(row[j] for row in g['components']) for j in range(5)]
        self.assertLess(max(abs(a-b) for a,b in zip(summed,self.solution['root'])),g['remainder'])

    def test_implicit_eta_derivative(self):
        dx=mp.mpf('1e-4')
        samples=[self.matrix.solve(fixture(self.eta+k*dx),self.eta+k*dx,self.P)['root'] for k in (-2,-1,1,2)]
        expected=[(samples[0][i]-8*samples[1][i]+8*samples[2][i]-samples[3][i])/(12*dx) for i in range(5)]
        self.assertLess(max(relative(a,b) for a,b in zip(expected,self.solution['root_eta'])),mp.mpf('1e-12'))

    def test_zero_value_does_not_remove_eta_derivative(self):
        d=fixture(0);s=self.matrix.solve(d,0,self.P)
        self.assertTrue(all(v==0 for v in s['root']))
        self.assertTrue(any(v!=0 for v in s['root_eta']))
        correct=self.annulus.state('-5',0,d,s)
        frozen=self.annulus.state('-5',0,d,s,True)
        self.assertLess(abs(correct['Q_defect'])/d['epsilon'],mp.mpf('1e-80'))
        self.assertGreater(abs(frozen['Q_defect'])/d['epsilon'],mp.mpf('.001'))

    def test_endpoint_fields_restore_exactly(self):
        f=self.annulus.fields('-5',self.eta,self.data,self.solution)
        self.assertEqual(f['U'].v,4*self.eta)
        self.assertEqual(f['epsilon'].v,0)
        self.assertEqual(f['epsilon_t'],0)
        self.assertEqual(f['Ut'],0)

    def test_closed_ideal_formula_matches_five_primitives(self):
        data=dict(entry=[Jet(0)]*5,g=Jet(0),epsilon=mp.mpf('1e-16'))
        s=dict(root=[mp.mpf(0)]*5,root_eta=[mp.mpf(0)]*5)
        for eta in map(mp.mpf,[-1,0,1]):
            row=self.annulus.state('-7',eta,data,s)
            ref=ideal(eta,self.annulus.h,self.P,mp.exp(-7))
            self.assertLess(abs(row['Q']-ref['Q']),mp.mpf('1e-90'))
            self.assertLess(abs(row['N_over_P2']-ref['N_over_P2']),mp.mpf('1e-90'))

    def test_missing_finite_h_term_is_detected(self):
        eta=mp.mpf('.5');h=mp.mpf('.005')
        correct=ideal(eta,h,self.P,mp.exp(-7))['Q']
        f=1/(1+eta**2);D=mp.mpf('.5')-h;d=1-eta**2
        wrong=mp.mpf(9)/8-5*h/8-h*eta**2/2+mp.mpf(5)/4*eta**2*(D+4*d)*f
        self.assertLess(abs((correct-wrong)/h-mp.mpf('.625')),mp.mpf('1e-90'))

    def test_full_patch_relaxed_inequalities(self):
        for y in map(mp.mpf,['-8','-7.5','-5.72','-5.48','-5.13','-5']):
            row=self.annulus.state(y,self.eta,self.data,self.solution)
            self.assertTrue(row['strict_relaxed_cone'])
            self.assertFalse(row['full_admissible_cone_claimed'])

    def test_old_radius_fails_joining_pressure_criterion(self):
        d=dict(entry=[Jet(0)]*5,g=Jet(0),epsilon=mp.mpf('1e-16'))
        s=dict(root=[mp.mpf(0)]*5,root_eta=[mp.mpf(0)]*5)
        small=Annulus(self.matrix,'.005',16,100).state('-8',0,d,s)
        self.assertLess(small['Pc'],2)
        self.assertFalse(small['strict_relaxed_cone'])

    def test_original_radial_sources_detect_the_small_edit(self):
        for y in ['-7.5','-5.72','-5.48']:
            r=self.annulus.source_check(y,self.eta,self.data,self.solution)
            self.assertLess(max(r['normalized_source_gaps']),mp.mpf('1e-30'))

    def test_independent_physical_density_ODE(self):
        r=self.matrix.independent_integrals(self.solution['root'],self.eta,self.P,.005)
        self.assertLess(max(r['gaps']),mp.mpf('1e-9'))

    def test_independent_quadratic_integration(self):
        r=self.matrix.tanh_sinh_quadratic(self.solution['root'],self.eta,self.P)
        self.assertLess(max(r['gaps']),mp.mpf('1e-30'))
        self.assertTrue(all(r[k]>0 for k in ('U_square','E_energy_square','E_pressure_square')))


class RecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from audit import validate
        cls.validate=staticmethod(validate)
        cls.record=json.loads((HERE/'evidence.json').read_text())

    def test_record_valid(self):self.validate(self.record)

    def test_missing_gate_rejected(self):
        r=copy.deepcopy(self.record);r['gates'].pop(next(iter(r['gates'])))
        with self.assertRaises(ValueError):self.validate(r)

    def test_faked_axis_entry_rejected(self):
        r=copy.deepcopy(self.record);r['certificate']['actual_axis_entry_verified']=True
        with self.assertRaises(ValueError):self.validate(r)

    def test_full_cone_promotion_rejected(self):
        r=copy.deepcopy(self.record);r['certificate']['full_admissible_cone_realized']=True
        with self.assertRaises(ValueError):self.validate(r)

    def test_negative_control_cannot_be_erased(self):
        r=copy.deepcopy(self.record);r['controls']['old_radius_Pc']='3'
        with self.assertRaises(ValueError):self.validate(r)

    def test_protocol_or_source_changes_rejected(self):
        for key in ('protocol','provenance'):
            r=copy.deepcopy(self.record);r[key]={}
            with self.assertRaises(ValueError):self.validate(r)

    def test_inherited_failure_cannot_be_erased(self):
        r=copy.deepcopy(self.record);r['inherited']['old_Md_failure']={}
        with self.assertRaises(ValueError):self.validate(r)


if __name__=='__main__':unittest.main()
