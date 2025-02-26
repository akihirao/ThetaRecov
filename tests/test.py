import unittest

import numpy as np
import pandas as pd

from cyvcf2 import VCF

import ThetaRecov
from ThetaRecov.core import calc_tajimaD_overall
from ThetaRecov.core import calc_tajimaD_windows


class TestTajimadOverall(unittest.TestCase):
	"""test class of calc_tajimaD_overall
	"""

	def test_tajimaD_overall(self):
		"""test method for calc_tajimaD_overall
		"""
		actual = ThetaRecov.core.calc_tajimaD_overall("test.vcf.gz", output_csv="test_tajimaD_overall.csv")
		expected = pd.read_csv("tajimaD_overall.csv")

        assert_frame_equal(actual, expected)

if __name__ == "__main__":
    unittest.main()