import unittest
import filecmp

import numpy as np
import pandas as pd

from pandas.testing import assert_frame_equal

from cyvcf2 import VCF

import ThetaRecov
from ThetaRecov.core import calc_tajimaD_overall
from ThetaRecov.core import calc_tajimaD_windows


class Test_calc_tajimaD(unittest.TestCase):
	"""
	test class of calc_tajimaD_overall/calc_tajimaD_windows
	"""

	def test_tajimaD_overall(self):
		"""
		test method for calc_tajimaD_overall
		"""
		actual = ThetaRecov.core.calc_tajimaD_overall("test.vcf.gz", output_csv="test_tajimaD_overall.csv")
		expected = pd.read_csv("tajimaD_overall.csv")
		assert_frame_equal(actual, expected)

	def test_tajimaD_windows(self):
		"""
		test method for calc_tajimaD_windows
		"""
		actual = ThetaRecov.core.calc_tajimaD_windows("test.vcf.gz", windows_size=1000,output_csv="test_tajimaD_windows.csv")
		expected = pd.read_csv("tajimaD_windows.csv")
		actual_TajimaD = pd.DataFrame({'Tajima_D': actual['Tajima_D'].astype(float)})
		expected_TajimaD = pd.DataFrame({'Tajima_D': expected['Tajima_D'].astype(float)})
		assert_frame_equal(actual_TajimaD,expected_TajimaD)


if __name__ == "__main__":
    unittest.main()