import argparse

import numpy as np
import pandas as pd

from cyvcf2 import VCF

import ThetaRecov
from ThetaRecov.core import calc_tajimaD_windows

def main():
    """
    Computation and output of parameters, theta and Tajima's D, in sliding windows
    Parameters:
        input_vcf: path of a input VCF/VCF.gz file
        windows_size: size of sliding windows (bp)
        output_csv: name of an output csv file
    Returns:
        pandas data frame of estimated parameters 
    """
    parser = argparse.ArgumentParser(description="Compute theta and Tajima's D in sliding windows.\nChromosome: chromosome or scaffold\nWindwos_Start: Window start position\nBases: bases sequenced\nS: number of seqregating sites\nTheta_Watterson: Watterson's theta\nTheta pi\nTajima_D: Tajima's D")
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument('--input_vcf', '-in', type=str, nargs = '?', help = 'Input a VCF/VCF.gz file', required = True)
    optional.add_argument('--output_csv', '-out', type=str, nargs = '?', help = 'Output a csv file (default: tajimaD_windows.csv)', default = "tajimaD_windows.csv")
    optional.add_argument('--windows_size', '-size', type=int, nargs = '?', help = 'Size of sliding windows (bp) (default: 1000)', default = 1000)
    
    args = parser.parse_args()
        
    ThetaRecov.core.calc_tajimaD_windows(args.input_vcf,
        windows_size = args.windows_size,
        output_csv = args.output_csv)
    

if __name__ == "__main__":
    main()
