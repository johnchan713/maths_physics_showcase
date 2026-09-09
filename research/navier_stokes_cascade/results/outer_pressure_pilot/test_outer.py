#!/usr/bin/env python3
"""Source-equation, pressure-tail, moment, precision and provenance regressions."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from schedule import Schedule,step,step_prime,power_jet
from coupling import CoupledParameters,SeedDatum,Parameters,construct,compare,diagnostics,axis_series,PREVIOUS
from audit import HERE,validate,compare_record,main,encode

PROTOCOL=json.loads((HERE/'protocol.json').read_text())


class ConstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=160
        cls.s=Schedule(PROTOCOL['parameters'],32)
        cls.p=CoupledParameters.from_schedule(cls.s,PROTOCOL['parameters'])
        cls.profile=construct(cls.p,mp.mpf('.5'),24)

    def setUp(self):
        mp.mp.dps=160

    def test_step_symmetry_and_derivative(self):
        for x in map(mp.mpf,['.1','.37','.5','.81']):
            self.assertLess(abs(step(x)+step(1-x)-1),mp.mpf('1e-150'))
            self.assertLess(abs(mp.diff(step,x)-step_prime(x)),mp.mpf('1e-140'))
        self.assertEqual(step(-1),0)
        self.assertEqual(step(2),1)
        self.assertEqual(step_prime(0),0)

    def test_pressure_jet_against_differentiation(self):
        eta=mp.mpf('.37');exponent=mp.mpf('1.3')
        jet=power_jet(eta,exponent,10)
        for k in (1,3,8,10):
            self.assertLess(abs(jet[k]*mp.factorial(k)-mp.diff(lambda e:(1+e*e)**-exponent,eta,k)),mp.mpf('1e-130'))

    def test_explicit_parameter_order(self):
        self.assertGreater(self.s.log_p,self.s.td)
        self.assertLess(self.s.h,mp.exp(-self.s.td))
        self.assertLess(self.s.h,self.s.lam)
        self.assertGreater(self.s.tw,25)
        self.assertGreater(self.s.wait,0)

    def test_invalid_schedule_rejected(self):
        data=dict(PROTOCOL['parameters']);data['lambda']='0'
        with self.assertRaises(ValueError):
            Schedule(data,16)

    def test_complete_schedule_and_pressure_mass(self):
        self.assertEqual(len(self.s.stages),13)
        self.assertEqual(len(self.s.components),14)
        self.assertEqual(self.s.components[0]['nodes'][0],(mp.mpf(5),mp.mpf(1)))
        for c in self.s.components:
            self.assertTrue(all(w>0 and 0<=theta<=1 for w,theta in c['nodes']))

    def test_boundary_values_and_slopes(self):
        for left,right in zip(self.s.stages[:-1],self.s.stages[1:]):
            for e in map(mp.mpf,['0','.5','1']):
                actual=left.log_amplitude+self.s.log_shape(left,left.length,e)
                expected=right.log_amplitude+self.s.log_shape(right,mp.mpf(0),e)
                self.assertLess(abs(actual-expected),mp.mpf('1e-130'))
                self.assertLess(abs(self.s.slope(left,left.length,e)-self.s.slope(right,mp.mpf(0),e)),mp.mpf('1e-130'))

    def test_terminal_wait_hits_prescribed_Q(self):
        end=self.s.q_before_wait*mp.exp(-(1-self.s.h)*self.s.wait)
        self.assertLess(abs(end/self.s.qp-1),mp.mpf('1e-130'))

    def test_terminal_slope_stays_in_paper_interval(self):
        stage=self.s.stages[-2]
        for y in [mp.mpf(k)/20 for k in range(61)]:
            relative=(self.s.slope(stage,y)+self.s.h)/self.s.h
            self.assertGreaterEqual(relative,0)
            self.assertLess(relative,mp.mpf('.25'))

    def test_pressure_tail_retained_with_arbitrary_exponent(self):
        part=self.s.pressure_parts(mp.mpf('.5'),0)[-1]
        mass=-mp.exp(part['log_scale'])*part['jet'][0]
        self.assertGreater(mass,0)
        self.assertEqual(float(mass),0.)

    def test_pressure_sign_bound_and_evenness(self):
        for e in map(mp.mpf,['0','.1','.5','1']):
            p=self.s.pressure(e,True)
            self.assertLess(p,-mp.mpf('2.5')/(1+e*e)**2)
            self.assertEqual(p,self.s.pressure(-e,True))
            if e:
                self.assertGreater(self.s.pressure_jet(e,1,True)[1],0)

    def test_pressure_scale_is_derived(self):
        jet=self.s.pressure_jet(mp.mpf('.5'),2)
        normalized=self.s.pressure_jet(mp.mpf('.5'),2,True)
        for a,b in zip(jet,normalized):
            self.assertLess(abs(a/(b*mp.exp(2*self.s.log_p))-1),mp.mpf('1e-140'))

    def test_actual_angular_correction_and_nonzero_endpoint_memory(self):
        for e in (mp.mpf('.5'),mp.mpf(1)):
            row=self.s.angular_correction(e)
            self.assertTrue(all(v!=0 for v in row['coefficients']))
            self.assertLess(max(row['relative_errors']),mp.mpf('1e-100'))
            self.assertLess(row['smallness_bound'],1)
            self.assertLess(row['maximum_relative_edit_bound'],mp.mpf('1e-8'))
            self.assertGreater(row['negative_control_pressure_gap'],mp.mpf('1e-6'))

    def test_first_radial_coefficients_use_scheduled_pressure(self):
        p=self.p;e=self.profile.eta0;b=axis_series(p,e,4)
        self.assertLess(abs(self.profile.u[1][0]/(-b['Z'][0]/(2*b['L'][0]))-1),mp.mpf('1e-130'))
        expected=(-b['chi'][0]+(b['W'][0]+p.h*(1-2*e*b['U'][0]))/(p.lam*b['L'][0]))/4
        self.assertLess(abs(self.profile.phi[1][0]-expected),mp.mpf('1e-130'))

    def test_nonlinear_leading_balance(self):
        row=diagnostics(self.profile,mp.mpf('4.1'))
        self.assertLess(max(row[k] for k in ('angular_residual','axial_residual','pressure_residual')),mp.mpf('1e-10'))
        self.assertGreater(row['Phi'],mp.mpf('.25'))
        self.assertGreater(row['normalized_axial_pressure_contribution'],mp.mpf('.9'))

    def test_second_assembly_matches_previous_seed_solver(self):
        from inner import construct as old_construct
        data=json.loads((PREVIOUS/'protocol.json').read_text())['parameters']
        p=Parameters.from_dict(data)
        cp=CoupledParameters(p.h,p.j,p.sigma,p.lam,mp.mpf(0),SeedDatum(p.pressure))
        e=p.hzero()
        self.assertLess(compare(old_construct(p,e,18),construct(cp,e,18),[mp.mpf(4)]),mp.mpf('1e-100'))


class RecordTests(unittest.TestCase):
    def setUp(self):
        mp.mp.dps=160

    def record(self):
        return json.loads((HERE/'evidence.json').read_text())

    def test_record_has_limited_scope(self):
        validate(self.record(),PROTOCOL)

    def test_no_proof_or_matching_promotion(self):
        for key in ('full_proof_verified','new_candidate','matched_paper_profile','global_thresholds_verified'):
            r=self.record();r[key]=True
            with self.assertRaisesRegex(ValueError,'scientific'):
                validate(r,PROTOCOL)

    def test_changed_limits_rejected(self):
        r=self.record();r['checks']['original-leading-equations']['limit']='1'
        with self.assertRaisesRegex(ValueError,'limit'):
            validate(r,PROTOCOL)

    def test_changed_measurements_rejected_even_with_old_pass_flag(self):
        r=self.record();r['inner_samples'][0]['refinements'][-1]['samples'][-1]['Phi']='0'
        with self.assertRaisesRegex(ValueError,'measurements'):
            validate(r,PROTOCOL)

    def test_nonfinite_evidence_rejected(self):
        r=self.record();r['pressure_samples'][0]['normalized_pressure']='nan'
        with self.assertRaisesRegex(ValueError,'Nonfinite'):
            validate(r,PROTOCOL)

    def test_zero_tail_rejected(self):
        r=self.record();r['pressure_decomposition']['positive_tail_mass_normalized']='0'
        with self.assertRaisesRegex(ValueError,'tail'):
            validate(r,PROTOCOL)

    def test_reproduction_checks_tiny_tail_relatively(self):
        r=self.record();changed=copy.deepcopy(r)
        changed['pressure_decomposition']['positive_tail_mass_normalized']='0'
        with self.assertRaisesRegex(ValueError,'differs'):
            compare_record(r,changed)

    def test_missing_obligations_rejected(self):
        r=self.record();r['unverified']=[]
        with self.assertRaisesRegex(ValueError,'limitations'):
            validate(r,PROTOCOL)

    def test_changed_refinement_rejected(self):
        r=self.record();r['inner_samples'][0]['refinements'].pop()
        with self.assertRaisesRegex(ValueError,'refinements'):
            validate(r,PROTOCOL)

    def test_existing_record_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'existing.json';output.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output',str(output)])
            self.assertEqual(output.read_text(),'preserve')

    def test_wrong_pdf_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            paper=Path(directory)/'wrong.pdf';paper.write_text('wrong')
            with self.assertRaisesRegex(ValueError,'Paper hash'):
                main(['--paper',str(paper),'--output',str(Path(directory)/'unused.json')])

    def test_nonfinite_encoding_rejected(self):
        with self.assertRaises(ValueError):
            encode(mp.inf)


if __name__=='__main__':
    unittest.main()
