import sys

import numpy as np
import pandas as pd

from cyvcf2 import VCF

import ThetaRecov



def exec_calc_tajimaD_windows():
    """
    Computation and output of parameters, theta and Tajima's D, in sliding windows
    Parameters:
        argv[0]: path of a vcf file
        argv[1]: size of sliding windows
        argv[2]: name of an output csv file
    Returns:
        pandas data frame of estimated parameters 
    """
	vcf_path = sys.argv[0]
	windows_size = sys.argv[1]
	output_csv = sys.argv[2]
	
	if not len(sys.argv)==3:
		
	ThetaRecov.calc_tajimaD_windows(vcf_path, windows_size, output_csv)
	
	
def exec_calc_tajimaD_overall():
    """
    Computation and output of parameters, theta and Tajima's D, across overall
    Parameters:
        argv[0]: path of a vcf file
        argv[1]: name of an output csv file
    Returns:
        pandas data frame of estimated parameters 
    """
	vcf_path = sys.argv[0]
	output_csv = sys.argv[1]
	
	ThetaRecov.calc_tajimaD_overall(vcf_path, output_csv)

