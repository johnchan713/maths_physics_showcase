#!/usr/bin/env python3
"""Regression checks for unit-scale roundoff in a tiny removed-energy fraction."""
import unittest

import numpy as np

from audit import close


class ArithmeticComparisonTests(unittest.TestCase):
    def test_ci_subtraction_roundoff_is_accepted_only_with_absolute_budget(self):
        observed, recorded = 2.486736954576685e-6, 2.4867369541325957e-6
        with self.assertRaises(ValueError):
            close(observed, recorded)
        close(observed, recorded, absolute_tolerance=64 * np.finfo(float).eps)

    def test_material_fraction_difference_is_rejected(self):
        with self.assertRaises(ValueError):
            close(2.49e-6, 2.48e-6, absolute_tolerance=64 * np.finfo(float).eps)

    def test_zero_uses_only_explicit_absolute_budget(self):
        close(0., 0.)
        with self.assertRaises(ValueError):
            close(0., 1e-15)
        close(0., 1e-15, absolute_tolerance=64 * np.finfo(float).eps)

    def test_nonfinite_values_and_tolerances_are_rejected(self):
        for value in (float('nan'), float('inf'), -float('inf')):
            for arguments in ((value, 1.), (1., value), (1., 1., value), (1., 1., 1e-10, value)):
                with self.assertRaises(ValueError):
                    close(*arguments)

    def test_negative_tolerances_are_rejected(self):
        for arguments in ((1., 1., -1.), (1., 1., 1e-10, -1.)):
            with self.assertRaises(ValueError):
                close(*arguments)


if __name__ == '__main__':
    unittest.main()
