#!/usr/bin/env python3
"""Local regression, nonlinear-equation, and non-promotion tests."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from inner import (Parameters,add,mul,inv,derivative,axis_series,construct,
                   diagnostics,direct_heat_join)
from physical import Jet,physical_point,residual,finite_difference,cartesian_finite_difference
from run import HERE,encode,validate,comparison,moment_quadrature,main,compare_record

PROTOCOL=json.loads((HERE/'protocol.json').read_text())


class AlgebraTests(unittest.TestCase):
    def setUp(self):
        mp.mp.dps=60

    def test_convolution(self):
        self.assertEqual(mul([1,2],[3,4],2),[3,10,8])

    def test_inverse(self):
        result=mul([mp.mpf(1),mp.mpf(2)],inv([mp.mpf(1),mp.mpf(2)],10),10)
        self.assertEqual(result,[1]+[0]*10)

    def test_zero_inverse_rejected(self):
        with self.assertRaises(ValueError):
            inv([0],5)

    def test_eta_derivative(self):
        self.assertEqual(derivative([1,2,3,4]),[2,6,12])

    def test_jet_product_hessian(self):
        x,y=Jet.variable(2,0),Jet.variable(3,1)
        value=x*x*y
        self.assertEqual(value.v,12)
        self.assertEqual(value.g[:2],[12,4])
        self.assertEqual(value.H[0][:2],[6,4])

    def test_jet_reciprocal(self):
        x=Jet.variable(2,0)
        value=1/x
        self.assertEqual(value.v,mp.mpf('.5'))
        self.assertEqual(value.g[0],mp.mpf('-.25'))
        self.assertEqual(value.H[0][0],mp.mpf('.25'))

    def test_invalid_parameters(self):
        for key,value in [('h','0'),('j','0'),('sigma','nan'),('Lambda','0'),('pressure_scale','-1')]:
            changed=dict(PROTOCOL['parameters']);changed[key]=value
            with self.assertRaises(ValueError):
                Parameters.from_dict(changed)


class ProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=90
        cls.p=Parameters.from_dict(PROTOCOL['parameters'])
        cls.profile=construct(cls.p,mp.mpf('.5'),24)
        cls.center=construct(cls.p,cls.p.hzero(),24)

    def setUp(self):
        mp.mp.dps=90

    def test_axis_data(self):
        row=diagnostics(self.profile,mp.mpf(0))
        self.assertEqual(row['Phi'],1)
        self.assertEqual(row['axis_U_error'],0)
        self.assertTrue(all(v==0 for v in row['moments']))

    def test_first_radial_coefficients(self):
        s=axis_series(self.p,self.profile.eta0,4)
        phi1=(-s['chi'][0]+(s['W'][0]+self.p.h*(1-2*self.profile.eta0*s['U'][0]))/(self.p.lam*s['L'][0]))/4
        u1=-s['Z'][0]/(2*s['L'][0])
        self.assertLess(abs(self.profile.phi[1][0]-phi1),mp.mpf('1e-80'))
        self.assertLess(abs(self.profile.u[1][0]-u1),mp.mpf('1e-80'))

    def test_seed_pressure_monotonicity(self):
        for e in map(mp.mpf,['-.9','-.1','.1','.9']):
            self.assertGreater(e*(4*self.p.pressure*e/(1+e*e)**3),0)

    def test_nonzero_tiny_swirl(self):
        g=mp.exp(self.profile.log_g0)
        self.assertGreater(g,0)
        self.assertEqual(float(g),0.)
        self.assertLess(self.profile.log_g0,-self.p.lam)

    def test_local_equations(self):
        for prof in [self.profile,self.center]:
            row=diagnostics(prof,mp.mpf('4.1'))
            for key in ('angular_residual','axial_residual','pressure_residual'):
                self.assertLess(row[key],mp.mpf('1e-12'))
            self.assertGreater(row['Phi'],mp.mpf('.25'))

    def test_not_just_comparison_profile(self):
        Y=mp.mpf(4);eta=self.profile.eta0
        chi=self.p.H(eta)**2/(self.p.H(eta)**2+self.p.sigma**2)
        comparison_value=2*mp.besselj(1,mp.sqrt(2*Y*chi))/mp.sqrt(2*Y*chi)
        self.assertGreater(abs(self.profile.value(self.profile.phi,Y)-comparison_value),mp.mpf('1e-8'))

    def test_degree_convergence(self):
        lower=construct(self.p,self.profile.eta0,18)
        self.assertLess(comparison(lower,self.profile,[mp.mpf(4)]),mp.mpf('1e-12'))

    def test_missing_transport_detected(self):
        self.assertGreater(diagnostics(self.profile,mp.mpf(4))['omitted_eta_transport_gap'],mp.mpf('1e-6'))

    def test_missing_pressure_detected(self):
        self.assertGreater(diagnostics(self.profile,mp.mpf(4))['omitted_axial_pressure_gap'],mp.mpf('1e-6'))

    def test_five_moments(self):
        self.assertLess(moment_quadrature(self.profile,mp.mpf(4))['maximum_error'],mp.mpf('1e-40'))

    def test_direct_join_rejected(self):
        result=direct_heat_join(self.profile)
        self.assertFalse(result['matched'])
        self.assertGreater(result['U_gap'],mp.mpf('1e-10'))
        self.assertGreater(result['M_gap'],mp.mpf('1e-10'))

    def test_cartesian_full_residual_is_nonzero(self):
        row=residual(self.profile,physical_point(self.profile,mp.mpf(2)))
        self.assertLess(row['divergence'],mp.mpf('1e-40'))
        self.assertLess(max(row['normalized_leading']),mp.mpf('1e-12'))
        self.assertGreater(max(row['normalized_full']),mp.mpf('.1'))

    def test_component_separated_physical_stencil(self):
        point=physical_point(self.profile,mp.mpf(2));exact=residual(self.profile,point)
        approximation=finite_difference(self.profile,point,mp.mpf('.005'))
        self.assertLess(max(abs(a-b)/s for a,b,s in zip(approximation,exact['total'],exact['scales'])),mp.mpf('1e-6'))

    def test_mixed_cartesian_absorption_failure_retained(self):
        point=physical_point(self.profile,mp.mpf(2));exact=residual(self.profile,point)
        bad=cartesian_finite_difference(self.profile,point,mp.mpf('.005'))
        self.assertGreater(abs(bad[1]-exact['total'][1])/exact['scales'][1],1)

    def test_invalid_construction(self):
        for eta,n in [(mp.mpf(2),12),(mp.mpf(0),1),(mp.mpf(0),41)]:
            with self.assertRaises(ValueError):
                construct(self.p,eta,n)

    def test_invalid_physical_point(self):
        with self.assertRaises(ValueError):
            physical_point(self.profile,mp.mpf(0))


class RecordTests(unittest.TestCase):
    def record(self):
        mp.mp.dps=90
        return json.loads((HERE/'evidence.json').read_text())

    def test_record_passes_only_local_gates(self):
        validate(self.record(),PROTOCOL)

    def test_missing_gate_rejected(self):
        record=self.record();record['checks'].pop('axis-data')
        with self.assertRaises(ValueError):
            validate(record,PROTOCOL)

    def test_changed_threshold_rejected(self):
        record=self.record();record['checks']['axis-data']['limit']='1'
        with self.assertRaisesRegex(ValueError,'threshold'):
            validate(record,PROTOCOL)

    def test_nonfinite_raw_evidence_rejected(self):
        record=self.record();record['samples'][0]['refinements'][0]['samples'][0]['Phi']='nan'
        with self.assertRaisesRegex(ValueError,'Nonfinite'):
            validate(record,PROTOCOL)

    def test_unverified_obligations_preserved(self):
        record=self.record();record['unverified']=[]
        with self.assertRaises(ValueError):
            validate(record,PROTOCOL)

    def test_proof_promotion_rejected(self):
        for name in ('full_proof_verified','new_candidate','matched_paper_profile'):
            record=self.record();record[name]=True
            with self.assertRaises(ValueError):
                validate(record,PROTOCOL)

    def test_relabelled_pressure_seed_rejected(self):
        record=self.record();record['pressure_seed_not_paper_schedule']=False
        with self.assertRaises(ValueError):
            validate(record,PROTOCOL)

    def test_changed_measurement_rejected(self):
        record=self.record();changed=copy.deepcopy(record)
        changed['samples'][0]['refinements'][0]['samples'][1]['Phi']='2'
        with self.assertRaisesRegex(ValueError,'measurements'):
            compare_record(record,changed)

    def test_join_promotion_rejected(self):
        record=self.record();record['direct_join_gates'][0]['passed']=True
        with self.assertRaises(ValueError):
            validate(record,PROTOCOL)

    def test_missing_physical_samples_rejected(self):
        record=self.record();record['physical'].pop()
        with self.assertRaises(ValueError):
            validate(record,PROTOCOL)

    def test_nonfinite_encoding_rejected(self):
        for value in [mp.inf,mp.nan,-mp.inf]:
            with self.assertRaises(ValueError):
                encode(value)

    def test_existing_output_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'existing.json'
            output.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output',str(output)])
            self.assertEqual(output.read_text(),'preserve')

    def test_wrong_paper_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            paper=Path(directory)/'wrong.pdf';paper.write_text('not the paper')
            with self.assertRaisesRegex(ValueError,'Paper hash'):
                main(['--output',str(Path(directory)/'new.json'),'--paper',str(paper)])


if __name__=='__main__':
    unittest.main()
