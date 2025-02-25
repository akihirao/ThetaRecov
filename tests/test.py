import sys

import numpy as np
import pandas as pd

from cyvcf2 import VCF

import ThetaRecov
from ThetaRecov.core import calc_tajimaD_overall
from ThetaRecov.core import calc_tajimaD_windows

# Computation of Tajima's D across overall sites to write a csv file.
#tajimaD_overall input_vcf output_csv
ThetaRecov.cor.tajimaD_overall test.vcf.gz test_TajimaD_overall.csv

# Computation of Tajima's D in sliding windows to write a csv file.
#tajimaD_windows input_vcf window size output_csv
ThetaRecov.cor.tajimaD_windows test.vcf.gz 1000 test_TajimaD_windows.csv
