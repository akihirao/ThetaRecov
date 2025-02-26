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
        input_vcf: an input VCF/VCF.gz file
        output_csv: an output file as csv
    Returns:
        pandas data frame of estimated parameters
    """
    parser = argparse.ArgumentParser(description="Compute theta and Tajima's D across overall genome.\nBases: bases sequenced\nS: number of seqregating sites\nTheta_Watterson: Watterson's theta\nTheta pi\nTajima_D: Tajima's D")
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument('--input_vcf', '-in', type=str, nargs = '?', help = 'Input a VCF/VCF.gz file', required = True)
    optional.add_argument('--output_csv', '-out', type=str, nargs = '?', help = 'Output a csv file (default: tajimaD_overall.csv)', default = "tajimaD_overall.csv")
    
    args = parser.parse_args()
    
    ThetaRecov.core.calc_tajimaD_overall(args.input_vcf,
        output_csv = args.output_csv)


if __name__ == "__main__":
    main()
