import json
import unittest
import audit
import ceiling as c
Q=c.Q
cert=c.load('q14_certified','core_certified_finite_stage_domain_join/certified_recurrence.py')

class CeilingTests(unittest.TestCase):
    def test_interface_not_primitive(self):
        self.assertEqual(c.q12.state_ceiling(0),500038)
        self.assertEqual(c.expand(c.q12.state_ceiling(0)),900047)
        self.assertEqual(c.expand(c.q12.state_ceiling(1)),1400076)
    def test_compressed_states_match_certified_provider(self):
        for j in (0,1,4,9):
            for n in (0,3,100):
                old=cert.state_join(j,n)
                self.assertEqual(c.q12.state_ceiling(j,n),old['first_real_interface_ceiling'])
                self.assertEqual(c.expand(c.q12.state_ceiling(j,n)),old['fully_expanded_primitive_ceiling'])
    def test_enlarged_physical_request(self):
        for h in (Q(1,200),Q(1,10**30)):
            for j,m,N in ((0,0,0),(1,2,3),(5,7,11)):
                row=c.request(h,j,m,N);old=cert.physical_join(j,m,N,h)
                self.assertEqual(row['enlarged_physical_primitive'],old['fully_expanded_primitive_ceiling']+2)
                self.assertEqual(row['combined_primitive'],row['combined_interface']+400009)
                self.assertEqual(row['base_strict_margins'],old['strict_margins'])
    def test_large_stage_without_recursion(self):
        h=Q(1,10**100)
        j=c.q12.rational_plan(h)['stage']
        row=c.request(h,j)
        self.assertEqual(row['state_primitive'],900047+99994*j*j+400035*j)
        self.assertEqual(row['combined_primitive'],max(row['state_primitive'],row['enlarged_physical_primitive']))
    def test_pole_increment_with_dominating_branches(self):
        for I in (0,1,500038,10**40):
            for original in (0,1,10**50):
                for retained in (0,30,10**60):
                    delta=c.pole(c.expand(I),original,retained)-c.pole(I,original,retained)
                    self.assertTrue(0<=delta<=2400054)
        self.assertEqual(c.pole(c.expand(500038))-c.pole(500038),2400054)
    def test_archive_retains_distinct_domains(self):
        row=c.archive([900047,1400076,900047])
        self.assertEqual(len(row['distinct']),2)
        self.assertEqual(row['by_request'][1],row['by_request'][2])
        for q,ops in zip(row['distinct'],row['operators']):
            self.assertEqual(ops['pressure'],(6*q,12*q))
            self.assertEqual(ops['moment'],(6*q,12*q))
            self.assertEqual(ops['oriented'],6*q)
    def test_invalid_indices(self):
        for n in (-1,True,Q(1,2)):
            with self.assertRaises(ValueError):c.expand(n)
        with self.assertRaises(ValueError):c.pole(0,-1)
    def test_provenance_and_scope(self):
        r=audit.generate()
        self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['scope']['representative_suprema_proved'])
        self.assertFalse(r['scope']['physical_field_changed'])

if __name__=='__main__':unittest.main()
