import sys

import numpy as np
import pandas as pd

from cyvcf2 import VCF

from ThetaRecov import ThetaRecov
from ThetaRecov import calc

# Computation of Tajima's D across overall sites to write a csv file.
#tajimaD_overall(vcf_path, output_path)
calc.tajimaD_overall("test.vcf.gz","TajimaD_overall.csv")

# Computation of Tajima's D in sliding windows to write a csv file.
#tajimaD_windows(vcf_path, window size, output_path)
calc.tajimaD_windows("test.vcf.gz",1000,"TajimaD_windows.csv")
