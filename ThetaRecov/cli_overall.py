import argparse

import numpy as np
import pandas as pd

from cyvcf2 import VCF

import ThetaRecov
from ThetaRecov.core import calc_tajimaD_overall


def main():
    """
    Computation and output of parameters, theta and Tajima's D, across overall genome
    Parameters:
        vcf_file: path of a vcf file
        output_file: name of an output csv file
    Returns:
        pandas data frame of estimated parameters 
    """
    parser = argparse.ArgumentParser(description="Compute theta and Tajima's D across overall genome from a vcf file.")
    parser.add_argument("vcf_gz_file", type=str, help="Input VCF file")
    parser.add_argument("output_file", type=str, help="Output results file")
    args = parser.parse_args()
    
	ThetaRecov.corecalc_tajimaD_overall(args.vcf_gz_file, 
        output_csv = args.output_file)
	

if __name__ == "__main__":
    main()
