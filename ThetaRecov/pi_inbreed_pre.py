import argparse

import numpy as np
import pandas as pd

import numexpr as ne

from cyvcf2 import VCF

import ThetaRecov
from ThetaRecov.core import calc_inbreed_pre

def main():
    """
    Computation and output of pi_within, pi_among, and 1 - pi_within/pi_among as the deviation from Hardy-Weinberg equilibrium
    Parameters:
        input_vcf: path of a input VCF/VCF.gz file
        output_csv: name of an output csv file
    Returns:
        pandas data frame of estimated parameters 
    """

    start_time = time.time() #stamp start time

    parser = argparse.ArgumentParser(description="Compute the deviance frow H-W equibrium")
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument('--input_vcf', '-in', type=str, nargs = '?', help = 'Input a VCF/VCF.gz file', required = True)
    optional.add_argument('--output_csv', '-out', type=str, nargs = '?', help = 'Output a csv file (default: inbreed.csv)', default = "inbreed.csv")
    
    args = parser.parse_args()
        
    ThetaRecov.core.calc_inbreed_pre(args.input_vcf,
        output_csv = args.output_csv)
    
    end_time = time.time() #stamp end time
    elasped_time = end_time - start_time # (sec)
    print(f"Function execution time: {elasped_time:.1f} seconds")


if __name__ == "__main__":
    main()
