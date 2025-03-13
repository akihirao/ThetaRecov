import argparse

import numpy as np
import pandas as pd

from itertools import combinations
from multiprocessing import Pool
from functools import partial

from cyvcf2 import VCF

import ThetaRecov

from ThetaRecov.core import vcf2gt_matrix
from ThetaRecov.core import diff_count_within
from ThetaRecov.core import calc_pi_within_elements_indiv_i
from ThetaRecov.core import calc_pi_among_elements_indiv_ij

def main():
    """
    Computation and output of pi_within, pi_among, and 1 - pi_within/pi_among as the deviation from Hardy-Weinberg equilibrium
    Parameters:
        input_vcf: path of a input VCF/VCF.gz file
        output_csv: name of an output csv file
    Returns:
        pandas data frame of estimated parameters 
    """
    parser = argparse.ArgumentParser(description="Compute the deviance frow H-W equibrium")
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument('--input_vcf', '-in', type=str, nargs = '?', help = 'Input a VCF/VCF.gz file', required = True)
    optional.add_argument('--output_csv', '-out', type=str, nargs = '?', help = 'Output a csv file (default: inbreed.csv)', default = "inbreed.csv")
    optional.add_argument('--threads', '-t', type=int, nargs = '?', help = 'Number of threads (default: 2)', default = 2)
   
    args = parser.parse_args()

    

    IN_VCF = args.input_vcf
    OUT_CSV = args.output_csv
    num_threads = args.threads

    results = []

    gt_matrix = vcf2gt_matrix(IN_VCF)
    

    #vcf_reader = VCF(IN_VCF)

    
    L = gt_matrix.shape[1] #length of sequences
    #samples = vcf_reader.samples #list of samples
    num_indiv = int(gt_matrix.shape[0] / 2)

    i_series = list(range(num_indiv))
    
    pairs = list(combinations(range(num_indiv), 2))

    
    result_within = []

    # Count difference in nucleotide between two sequences within individuals
    print("processing count within individuals")
    diff_within, count_within = diff_count_within(gt_matrix)

    print(f"Number fo sample pairs: {len(pairs)}") # check how much pairs
    print("processing count among individuals")
    with Pool(num_threads) as pool:
        #result_among =  pool.map(partial(ThetaRecov.core.calc_pi_among_elements_indiv_ij, IN_VCF), pairs[:20])
        result_among =  pool.map(partial(ThetaRecov.core.calc_pi_among_elements_indiv_ij, IN_VCF), pairs)
    diff_count_among = np.array(result_among).sum(axis=0)

    #print(f"diff_count_within: {diff_count_within}") #for debug
    print(f"diff_count_among: {diff_count_among}") #for debug

    pi_overall = (diff_within + diff_count_among[0])/(count_within + diff_count_among[1])
    pi_within = diff_within/count_within
    pi_among = diff_count_among[0]/diff_count_among[1]
    homo_deviance = 1 - pi_within/pi_among

    results.append([pi_overall,pi_within,pi_among,homo_deviance])
    df = pd.DataFrame(results, columns=["pi_overall","pi_within","pi_among","homo_deviance"])
    df.to_csv(args.output_csv, sep=",", index=False)
    
    #results.append([pi_within])
    #df = pd.DataFrame(results, columns=["pi_within"])
    #df.to_csv(args.output_csv, sep=",", index=False)
    

if __name__ == "__main__":
    main()
