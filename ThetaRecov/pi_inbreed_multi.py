import argparse

import numpy as np
import pandas as pd

from itertools import combinations
from multiprocessing import Pool
from multiprocessing import shared_memory

from functools import partial

from cyvcf2 import VCF

import ThetaRecov

from ThetaRecov.core import vcf2gt_matrix
from ThetaRecov.core import calc_pi_within_elements_indiv_i
from ThetaRecov.core import calc_pi_among_elements_indiv_ij


def worker(name, shape, dtype):
    """共有メモリを開いて処理"""
    shm = shared_memory.SharedMemory(name=name)
    array = np.ndarray(shape, dtype=dtype, buffer=shm.buf)

    shm.array




def pi_within_elements_indiv_i(index, name, shape, dtype):
    """共有メモリを開いて, gt_matrixからpi_withinを算出"""

    shm = shared_memory.SharedMemory(name=name)
    matrix = np.ndarray(shape, dtype=dtype, buffer=shm.buf)

    matrix_n_2_m = matrix.reshape(-1,2,matrix.shape[1])
    
    #row_1 = matrix[2 * index + 1, :]
    #row_2 = matrix[2 * index + 2, :]

    #target_indiv_matrix = np.vstack((row_1, row_2))
    #nan_mask = np.isnan(target_indiv_matrix).any(axis=0)
    #filtered_target_indiv_matrix = target_indiv_matrix[~nan_mask]


    target_indiv_matrix = matrix_n_2_m[index]
    mask = ~np.isnan(target_indiv_matrix).any(axis=0)
    filtered_target_indiv_matrix = target_indiv_matrix[:, mask]
    
    diff_within = np.sum(np.abs(np.diff(filtered_target_indiv_matrix, axis=0)))
    count_within = filtered_target_indiv_matrix.shape[1]

    shm.close() #メモリを閉じる

    #return count_within
    return diff_within, count_within



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

    results = []

    IN_VCF = args.input_vcf
    OUT_CSV = args.output_csv
    num_threads = args.threads


    #vcf_reader = VCF(IN_VCF)

    gt_matrix = ThetaRecov.core.vcf2gt_matrix(IN_VCF) #matrix of 2n samples by m loci 
    
    #共有メモリを作成
    shm = shared_memory.SharedMemory(create=True, size = gt_matrix.nbytes)

    #共有メモリをndarry配列にリンク
    shared_array = np.ndarray(gt_matrix.shape, dtype=gt_matrix.dtype, buffer=shm.buf)
    shared_array = gt_matrix[:] #データをコピー

    #L = vcf_reader.seqlens[0] #length of sequences
    #samples = vcf_reader.samples #list of samples
    num_samples = int(gt_matrix.shape[0] / 2)

    i_series = list(range(num_samples))
    
    print(i_series)
    pairs = list(combinations(range(num_samples), 2))

    
    with Pool(num_threads) as pool:
        #result_within = []
        #for res_within in pool.imap_unordered(partial(ThetaRecov.core.calc_pi_within_elements_indiv_i, gt_matrix), i_series):
        #    result_within.append(res_within)
        result_within = pool.starmap(pi_within_elements_indiv_i, [(i, shm.name, gt_matrix.shape, gt_matrix.dtype) for i in i_series])

        #result_among = []
        #for res_among in pool.imap_unordered(partial(ThetaRecov.core.calc_pi_among_elements_indiv_ij, gt_matrix), pairs):
        #    result_among.append(res_among)
    
    diff_count_within = np.array(result_within).sum(axis=0)
    #diff_count_among = np.array(result_among).sum(axis=0)
    print(f"result_within: {result_within}") #for debug
    print(f"diff_count_within: {diff_count_within}") #for debug

    #print(f"diff_count_among: {diff_count_among}") #for debug

    #pi_overall = (diff_count_within[0] + diff_count_among[0])/(diff_count_within[1] + diff_count_among[1])
    #pi_within = diff_count_within[0]/diff_count_within[1]
    #pi_among = diff_count_among[0]/diff_count_among[1]
    #homo_deviance = 1 - pi_within/pi_among

    #results.append([pi_within])
    #df = pd.DataFrame(results, columns=["pi_within"])
    #df.to_csv(OUT_CSV, sep=",", index=False)
    
    #results.append([pi_overall,pi_within,pi_among,homo_deviance])
    #df = pd.DataFrame(results, columns=["pi_overall","pi_within","pi_among","homo_deviance"])
    #df.to_csv(OUT_CSV, sep=",", index=False)

    shm.close()
    shm.unlink()
    
if __name__ == "__main__":
    main()
