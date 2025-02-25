import argparse

import numpy as np
import pandas as pd

from cyvcf2 import VCF

import ThetaRecov
form ThetaRecov.core import calc_tajimaD_windows

def main():
    """
    Computation and output of parameters, theta and Tajima's D, in sliding windows
    Parameters:
        vcf_file: path of a vcf file
        windows_size: size of sliding windows (bp)
        output_file: name of an output csv file
    Returns:
        pandas data frame of estimated parameters 
    """
    parser = argparse.ArgumentParser(description="Compute theta and Tajima's D in sliding windows from a vcf file.")
    parser.add_argument("vcf_gz_file", type=str, help="Input VCF file")
    parser.add_argument("windows_size", type=int, help="size of windows")
    parser.add_argument("output_file", type=str, help="Output results file")

    args = parser.parse_args()
        
    ThetaRecov.core.calc_tajimaD_windows(args.vcf_gz_file,
        windows_size = args.windows_size,
        output_csv = args.output_file)
    

if __name__ == "__main__":
    main()
