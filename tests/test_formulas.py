"""Numerical smoke tests for formulas represented in the figures."""

from __future__ import annotations

import unittest

import numpy as np

from scripts.brownian_contour import heat_kernel, upper_hitting_probability
from scripts.srw_convergence import biased_increment_cf, brownian_cf
from scripts.transform_tightness import cantor_cdf, uniform_shift_cf


class FormulaTests(unittest.TestCase):
    def test_heat_kernel_normalizes(self) -> None:
        x = np.linspace(-12.0, 12.0, 100_001)
        for t in (0.2, 1.0, 3.0):
            self.assertTrue(
                np.isclose(np.trapezoid(heat_kernel(x, t), x), 1.0, atol=1e-10)
            )

    def test_upper_hitting_probability_is_increasing_in_drift(self) -> None:
        mu = np.linspace(-3.0, 3.0, 301)
        probability = upper_hitting_probability(mu)
        self.assertTrue(np.all(np.diff(probability) > 0.0))
        self.assertTrue(np.isclose(probability[len(probability) // 2], 0.5))

    def test_arw_characteristic_function_converges_on_compact_set(self) -> None:
        u = np.linspace(-4.0, 4.0, 801)
        mu, sigma, horizon = 0.6, 1.1, 1.0
        errors = []
        target = brownian_cf(u, horizon, mu, sigma)
        for dt in (0.1, 0.025, 0.00625):
            discrete = biased_increment_cf(u, dt, mu, sigma) ** round(horizon / dt)
            errors.append(np.max(np.abs(discrete - target)))
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])

    def test_cantor_cdf_and_escaping_uniform_cf_endpoints(self) -> None:
        x = np.array([-1.0, 0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, 2.0])
        values = cantor_cdf(x)
        self.assertTrue(np.all((0.0 <= values) & (values <= 1.0)))
        self.assertTrue(np.allclose(values[[0, 1, 4, 5]], [0.0, 0.0, 1.0, 1.0]))
        self.assertTrue(
            np.isclose(uniform_shift_cf(np.array([0.0]), n=20)[0], 1.0 + 0.0j)
        )


if __name__ == "__main__":
    unittest.main()
