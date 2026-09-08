#!/usr/bin/env python3
"""Analytic Fourier and adversarial bookkeeping checks; no trajectory runs."""
import unittest

import numpy as np

from analyze import (compare_fields, direct_curl, embed, low_pass, peak_sum,
                     projection_fractions, require_field, sample_curl, spectral_norms)


def shear(n=8, amplitude=1.):
    field = np.zeros((n, n, n, 3), dtype=complex)
    field[0, 1, 0, 0] = field[0, -1, 0, 0] = amplitude / 2
    return field


class FourierTests(unittest.TestCase):
    def test_analytic_shear_curl_and_parseval(self):
        field = shear()
        omega = sample_curl(field, 16)
        expected = np.sin(2 * np.pi * np.arange(16) / 16)
        np.testing.assert_allclose(omega[..., 2], np.broadcast_to(expected[None, :, None], (16, 16, 16)), atol=1e-14)
        self.assertEqual(float(np.max(np.abs(omega[..., :2]))), 0.)
        self.assertAlmostEqual(peak_sum(omega)["maximum"], 1.)
        norms = spectral_norms(field)
        for value in norms.values():
            self.assertAlmostEqual(value, .5)
        self.assertAlmostEqual(float(np.mean(np.sum(omega**2, axis=-1))), norms["vorticity_l2_squared"])

    def test_embedding_preserves_physical_field(self):
        field = shear()
        np.testing.assert_allclose(sample_curl(field, 32), sample_curl(embed(field, 16), 32), atol=1e-14)
        self.assertEqual(spectral_norms(field), spectral_norms(embed(field, 16)))

    def test_direct_sum_matches_analytic_values(self):
        points = [[0, 4, 0], [3, 12, 7], [0, 0, 0]]
        values = direct_curl(shear(), points, 16)
        np.testing.assert_allclose(values, [[0, 0, 1], [0, 0, -1], [0, 0, 0]], atol=1e-14)

    def test_filter_separates_known_modes(self):
        field = shear()
        field[2, 0, 0, 1] = field[-2, 0, 0, 1] = .25
        np.testing.assert_array_equal(low_pass(field, 1), shear())
        np.testing.assert_array_equal(low_pass(field, 2), field)
        np.testing.assert_array_equal(low_pass(field, 0), np.zeros_like(field))

    def test_nonfinite_field_is_rejected(self):
        field = shear()
        field[0, 0, 0, 0] = np.nan
        with self.assertRaises(ValueError):
            require_field(field)

    def test_noncubic_field_is_rejected(self):
        with self.assertRaises(ValueError):
            require_field(np.zeros((8, 8, 16, 3)))

    def test_downsampling_is_rejected(self):
        with self.assertRaises(ValueError):
            sample_curl(shear(16), 8)
        with self.assertRaises(ValueError):
            embed(shear(16), 8)

    def test_nonnested_sampling_is_rejected(self):
        with self.assertRaises(ValueError):
            sample_curl(shear(), 24)

    def test_invalid_filter_is_rejected(self):
        for cutoff in (-1, 3):
            with self.assertRaises(ValueError):
                low_pass(shear(), cutoff)


class DecompositionTests(unittest.TestCase):
    def test_known_low_error_and_high_tail(self):
        coarse, fine = shear(8, .8), shear(16)
        fine[2, 0, 0, 1] = fine[-2, 0, 0, 1] = .25
        result = compare_fields(coarse, fine, 1, 32, [1, 2, 5])
        for key, value in (("coarse", .8), ("fine", 2.), ("filtered_fine", 1.),
                           ("extra", 1.), ("shared_error", .2), ("total_error", 1.2)):
            self.assertAlmostEqual(result["peaks"][key]["maximum"], value)
        self.assertAlmostEqual(result["high_mode_fractions_of_fine_squared_norms"]["vorticity_l2_squared"], .5)
        self.assertAlmostEqual(result["shared_mode_fractions_of_squared_error"]["vorticity_l2_squared"], .02 / .52)
        self.assertLess(max(result["spectral_split_relative_errors"].values()), 1e-12)
        self.assertLess(result["sampled_curl_identity_relative_error"], 1e-12)
        self.assertLess(max(result["direct_fourier_relative_errors"].values()), 1e-12)
        fractions = result["pointwise"][1]["projected_error_fractions"]
        self.assertAlmostEqual(fractions["shared"], 1 / 6)
        self.assertAlmostEqual(fractions["extra"], 5 / 6)

    def test_sum_of_maxima_is_not_used(self):
        a = np.zeros((8, 8, 8, 3))
        a[..., 0] = 1.
        self.assertEqual(peak_sum(a)["maximum"], 1.)
        self.assertEqual(peak_sum(-a)["maximum"], 1.)
        self.assertEqual(peak_sum(a, -a)["maximum"], 0.)

    def test_signed_projections_allow_cancellation(self):
        fractions = projection_fractions([-1., 0., 0.], [2., 0., 0.])
        self.assertEqual(fractions, {"shared": -1., "extra": 2.})

    def test_zero_vector_error_has_no_projection_fraction(self):
        self.assertIsNone(projection_fractions([1., 0., 0.], [-1., 0., 0.]))


if __name__ == "__main__":
    unittest.main()
