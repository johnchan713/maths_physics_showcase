import json
import unittest
from fractions import Fraction as F
import sympy as sp

import audit
import oscillation as o


class Tests(unittest.TestCase):
    def test_cover_preserves_area_with_all_deck_copies(self):
        # Lift two spanning vectors back by each actual transition matrix.
        J = sp.Matrix([[3,1],[1,5]])
        V = sp.Matrix([[1,sp.sqrt(2)-1],[-(sp.sqrt(2)-1),1]])
        for d in range(4):
            B = J**d
            lifted = B.inv()*V
            self.assertEqual(sp.simplify(B.det()*lifted.det()-V.det()), 0)
        r = F(1,6057072)
        self.assertLess(o.union_area(3*2*3**3,r), 1)
        self.assertEqual(o.union_area(2,F(1,8)),1)  # no gap follows

    def test_physical_normalization_on_overlapping_bands(self):
        q,s,h = sp.symbols('q s h', positive=True)
        ref = q/s
        physical = q**(-1-h)
        ratio = ref**(1+2*h)*physical/(ref**h*s**(-1-h))
        self.assertEqual(sp.powsimp(sp.expand_power_base(ratio, force=True),force=True).simplify(),1)

    def test_signed_losses_and_nontrivial_tail_threshold(self):
        b,c = F(7,13),F(97)
        self.assertEqual(o.normalized_floor(b,b/2,c,b/(4*c)), b/(4*c))
        self.assertEqual(o.normalized_floor(b,b/2,c,b/c),0)
        self.assertEqual(o.normalized_floor(b,2*b,c,0),0)

    def test_frozen(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))


if __name__ == '__main__':
    unittest.main()
