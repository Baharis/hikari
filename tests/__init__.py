"""
Set environment variable "HIKARI_TESTS_RAPID=True" to skip the slowest tests.
"""

from unittest import TestCase

import numpy as np


class NumpyTestCase(TestCase):
    def assertAllClose(self, a, b, **kwargs):
        if kwargs.get('atol', None) is None:
            kwargs['atol'] = 1e-15
        self.assertIsNone(np.testing.assert_allclose(a, b, **kwargs))
