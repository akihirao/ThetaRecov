import argparse

import numpy as np
import pandas as pd

from cyvcf2 import VCF

import ThetaRecov
from ThetaRecov.core import calc_inbreed

def main():
    """
    Computation and output of parameters, theta and Tajima's D, in sliding windows
    Parameters:
        input_vcf: path of a input VCF/VCF.gz file
        output_csv: name of an output csv file
    Returns:
        pandas data frame of estimated parameters 
    """
    parser = argparse.ArgumentParser(description="Compute theta and Tajima's D in sliding windows.\nChromosome: chromosome or scaffold\nWindwos_Start: Window start position\nBases: bases sequenced\nS: number of seqregating sites\nTheta_Watterson: Watterson's theta\nTheta pi\nTajima_D: Tajima's D")
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument('--input_vcf', '-in', type=str, nargs = '?', help = 'Input a VCF/VCF.gz file', required = True)
    optional.add_argument('--output_csv', '-out', type=str, nargs = '?', help = 'Output a csv file (default: inbreed.csv)', default = "inbreed.csv")
    
    args = parser.parse_args()
        
    ThetaRecov.core.calc_inbreed(args.input_vcf,
        output_csv = args.output_csv)
    

if __name__ == "__main__":
    main()
