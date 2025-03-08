import numpy as np
import pandas as pd

from itertools import combinations

from cyvcf2 import VCF



#=======================================
def vcf2gt_matrix(input_vcf_file):
    """
    Generate gt_matrix as ndarray object from a vcf file
    Parameters:
        input_vcf_file:
    Returns:
        gt_matrix: a ndarray object that stores genotypic information
    """
    vcf_reader = VCF(input_vcf_file)

    L = vcf_reader.seqlens[0] #length of sequences
    samples = vcf_reader.samples #list of samples
    num_samples = len(samples)
    num_seq = num_samples * 2

    gt_matrix = [] # Set matrix of genotypes

    for record in vcf_reader:
        gts = record.genotypes
        gts_list = []
        for gt in gts:
            gts_list.append(gt[0])
            gts_list.append(gt[1])
        gt_matrix.append(gts_list)

    gt_matrix = np.array(gt_matrix).T
       
    #if the source vcf file includes missing genotypes, substitute "nan".
    if np.isnan(gt_matrix).any():
        gt_matrix = np.where(gt_matrix == -1, np.nan, gt_matrix)
        warnings.warn("The input vcf includes missing genotypes.")
    
    return gt_matrix

    

#=========================================================
def parse_genotypes(variants):

    all_genotypes = []

    for v in variants:
        raw_genotypes = v.genotypes
        processed_genotypes = []

        for gt in raw_genotypes:
            alleles = gt[:2]
            alleles = [np.nan if a == -1 else a for a in alleles]
            processed_genotypes.extend(alleles)

        all_genotypes.append(processed_genotypes)

    return np.array(all_genotypes).T


#=======================================
# define function for caluclating Tajima's D
def calc_tajimaD(S, theta_w_region, theta_pi_region,n):
    if S == 0:
        tajima_d = np.nan

    if n > 1:
        a1 = sum(1 / i for i in range(1, n))
        a2 = sum (1 / (i**2) for i in range(1,n))
        b1 = (n + 1) / (3 * (n -1))
        b2 = (2 * (n ** 2 + n + 3)) / (9 * n * (n - 1))
        c1 = b1 - 1 / a1
        c2 = b2 - (n + 2) / (a1 * n ) + a2 / (a1 ** 2)
        e1 = c1 / a1
        e2 = c2 / (a1 ** 2 + a2)

        variance = (e1 * S) + (e2 * S * (S - 1))
    
        if not variance == 0:
            tajima_d = (theta_pi_region - theta_w_region) / np.sqrt(variance)
        else:
            tajima_d = np.nan #void division by zero
    else:
        tajima_d = np.nan

    return tajima_d



#=======================================
def calc_gt_matrix2tajimaD(gt_matrix):
    """
    compulation theta from gt_matrix as an ndarray object
    Parameters:
        gt_matrix: an ndarray object that stores genotypic information
    Returns:
        summary_statistics
    
    """
    num_allel_exp, L_obs = gt_matrix.shape
    num_samples = int(num_allel_exp /2)

    num_ref = np.sum(gt_matrix == 0, axis=0)
    num_alt = np.sum(gt_matrix == 1, axis=0)
    
    num_allel = np.sum(np.isin(gt_matrix, [0,1]), axis=0).astype(int)
    
    ave_num_allel = np.nanmean(num_allel)
    ave_num_allel = np.rint(ave_num_allel).astype(int)
    num_comparison_ave_num_allel = ave_num_allel * (ave_num_allel - 1) /2
    
    num_diff = num_ref * num_alt
    num_comparisons = num_allel * (num_allel - 1) /2
    num_comparisons_complete = num_allel_exp * (num_allel_exp - 1) /2
    
    freq_ref = num_ref / num_allel
    freq_alt = num_alt / num_allel
    
    def make_het(freq_ref, freq_alt, num_allel):
        freq_ref_hom = freq_ref**2
        freq_alt_hom = freq_alt**2
        return  (1 - freq_ref_hom - freq_alt_hom ) * (num_allel) /(num_allel -1)
    
    het = np.where(num_allel > 1, make_het(freq_ref, freq_alt, num_allel), np.nan)
    
    a1_site = np.array([np.sum(1 / np.arange(1,i)) if i > 1 else 0 for i in num_allel])
    S_site = (het > 0).astype(int) / a1_site
    
    S = np.sum(het > 0).astype(int)

    theta_pi = np.sum(num_diff/num_comparisons_complete)/L_obs
    theta_pi = float(theta_pi)
    theta_pi_unweight = np.sum(num_diff/num_comparisons)/L_obs
    theta_pi_unweight = float(theta_pi_unweight)
    theta_pixy = np.sum(num_diff) / (np.sum(num_comparisons))
    theta_pixy = float(theta_pixy)

    ##theta_pi_unweigted is derived from the equation (2) in Samuk (2022)
    ##theta_pi_region is derived from the equation (11) in Tajima (1989)
    # theta_pi_region / L_obs = theta_pi_he = theta_pi_unweigted
    
    theta_pi_region = float(np.nansum(het))
    num_no_count_sites = np.sum(np.isnan(het))

    theta_pi_he = float(theta_pi_region / (L_obs - num_no_count_sites))

    a1 = sum(1 / i for i in range(1, num_allel_exp))
    a1_effective = sum(1 / i for i in range(1, ave_num_allel))
    
    theta_w_region = float(S/ a1)
    theta_w_region_corr = np.nansum(S_site)
    theta_w_region_corr = float(theta_w_region_corr)
    
    theta_w = float(S / a1 / (L_obs - num_no_count_sites))
    theta_w_corr = float(S / a1_effective / (L_obs - num_no_count_sites))
    
    tajima_D = calc_tajimaD(S, theta_w_region_corr, theta_pi_region, ave_num_allel)
    tajima_D = float(tajima_D)
    
    summary_statistics = {"L": L_obs,\
              "S": S,\
              "num_seq": num_allel_exp,\
              "effective_num_seq": ave_num_allel,\
              #"theta_w": theta_w,\
              #"theta_w_corr": theta_w_corr,\
              #"theta_w_region": theta_w_region,\
              "theta_w_region_corr": theta_w_region_corr,\
              "theta_pi_region": theta_pi_region,\
              #"theta_pi": theta_pi, \
              #"theta_pi_unweight": theta_pi_unweight,\
              "theta_pixy": theta_pixy,\
              "theta_pi_he": theta_pi_he,\
              "tajima_D": tajima_D,\
             }
    
    return summary_statistics



#==========================================================================
def calc_tajimaD_windows(vcf_path, windows_size = 1_000_000, output_csv = "tajimaD_windows.csv"):
    """
    calculation of Tajima's D in sliding windows from a vcf file
    Parameters:
        vcf_path: path of a vcf file
        windows_size: size of sliding windows (default: 1000000)
        output_csv: an ouput file name (defualt: "tajimaD_windows.csv")
    Returns:
        a csv file contains positons and estimates as follows: Chromosome, Windows position, Bases, S, Pi, Pixy, Tajima's D
    
    """

    vcf = VCF(vcf_path)
    tajimaD_windows_results = []
    windows = {}
    
    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS
        window_start = ((pos - 1) // windows_size) * windows_size #Window start position
        windows.setdefault((chrom, window_start), []).append(variant)

    
    for (chrom, start), variants in windows.items():
        gt_matrix_unfiltered_nan = parse_genotypes(variants)
        nan_cols = np.all(np.isnan(gt_matrix_unfiltered_nan), axis=0) #Determine if all are NaN for each column
        gt_matrix = gt_matrix_unfiltered_nan[:, ~nan_cols] # Delete NaN-only columns

        base_sequenced = gt_matrix.shape[1]
        start_lab = start + 1
        summary_statistics = calc_gt_matrix2tajimaD(gt_matrix)
        S = summary_statistics.get("S")
        num_seq = summary_statistics.get("num_seq")
        effective_num_seq = summary_statistics.get("effective_num_seq")
        theta_w_region = summary_statistics.get("theta_w_region_corr")
        theta_pi_region = summary_statistics.get("theta_pi_region")
        pi_h = summary_statistics.get("theta_pi_he")
        pixy = summary_statistics.get("theta_pixy")
        tajima_D = summary_statistics.get("tajima_D")
        tajimaD_windows_results.append([chrom,start_lab,base_sequenced,S,num_seq,effective_num_seq,theta_w_region,theta_pi_region,tajima_D])

    df = pd.DataFrame(tajimaD_windows_results, columns=["Chromosome","Windows_Start","Bases","S","N_seq","Ne_seq","Theta_Watterson","Theta_pi","Tajima_D"])
    df.to_csv(output_csv, sep=",", index=False)
    
    return df


#==========================================================================
def calc_tajimaD_overall(vcf_path, output_csv = "tajimaD_overall.csv"):
    """
    calculation of Tajima's D based on overall sequenced sites in a vcf file
    Parameters:
        vcf_path: path of a vcf file
        output_csv: an ouput file name (defualt: "tajimaD_overall.csv")
    Returns:
        a csv file contains positons and estimates as follows: Bases, S, Pi, Pixy, Tajima's D
    
    """

    tajimaD_overall_results = []
    
    gt_matrix_unfiltered_nan = vcf2gt_matrix(vcf_path)

    nan_cols = np.all(np.isnan(gt_matrix_unfiltered_nan), axis=0) #Determine if all are NaN for each column
    gt_matrix = gt_matrix_unfiltered_nan[:, ~nan_cols] # Delete NaN-only columns

    base_sequenced = gt_matrix.shape[1]
    summary_statistics = calc_gt_matrix2tajimaD(gt_matrix)
    S = summary_statistics.get("S")
    num_seq = summary_statistics.get("num_seq")
    effective_num_seq = summary_statistics.get("effective_num_seq")
    theta_w_region = summary_statistics.get("theta_w_region_corr")
    theta_pi_region = summary_statistics.get("theta_pi_region")
    pi_h = summary_statistics.get("theta_pi_he")
    pixy = summary_statistics.get("theta_pixy")
    tajima_D = summary_statistics.get("tajima_D")
    tajimaD_overall_results.append([base_sequenced,S,num_seq,effective_num_seq,theta_w_region,theta_pi_region,tajima_D])

    df = pd.DataFrame(tajimaD_overall_results, columns=["Bases","S","N_seq","Ne_seq","Theta_Watterson","Theta_Pi","Tajima_D"])
    df.to_csv(output_csv, sep=",", index=False)
    
    return df




#==========================================================================
def calc_inbreed(vcf_path, output_csv = "inbreed.csv"):
    """
    calculation of Tajima's D based on overall sequenced sites in a vcf file
    Parameters:
        vcf_path: path of a vcf file
        output_csv: an ouput file name (defualt: "inbreed.csv")
    Returns:
        a csv file contains positons and estimates as follows: pi_within, pi_among, homo_deviance
    
    """

    inbreed_results = []
    
    gt_matrix = vcf2gt_matrix(vcf_path)

    if np.isnan(gt_matrix).any():
        warnings.warn("Missing genotypes must be imputed!")

    base_sequenced = gt_matrix.shape[1]
    gt_matrix_n_2_m = gt_matrix.reshape(-1,2,gt_matrix.shape[1])

    num_indiv = gt_matrix_n_2_m.shape[0]
    
    pi_within = np.sum(np.abs(gt_matrix_n_2_m[:, 0, :] - gt_matrix_n_2_m[:, 1, :]))/num_indiv/base_sequenced


    diff_among = 0
    num_diff_comb = num_indiv * (num_indiv - 1) // 2

    # (i, j) の組み合わせを全て列挙 (i < j のみ)
    for i, j in combinations(range(gt_matrix_n_2_m.shape[0]), 2):
        # 4つのペアの絶対差を求める
        diff_11 = np.abs(gt_matrix_n_2_m[i, 0] - gt_matrix_n_2_m[j, 0])  # (iの1行目, jの1行目)
        diff_12 = np.abs(gt_matrix_n_2_m[i, 0] - gt_matrix_n_2_m[j, 1])  # (iの1行目, jの2行目)
        diff_21 = np.abs(gt_matrix_n_2_m[i, 1] - gt_matrix_n_2_m[j, 0])  # (iの2行目, jの1行目)
        diff_22 = np.abs(gt_matrix_n_2_m[i, 1] - gt_matrix_n_2_m[j, 1])  # (iの2行目, jの2行目)

        # 4ペア分の総和を計算
        diff_among += np.sum(diff_11) + np.sum(diff_12) + np.sum(diff_21) + np.sum(diff_22)


    pi_among = diff_among/4/num_diff_comb/base_sequenced

    homo_deviance = 1 - pi_within/pi_among

    inbreed_results.append([pi_within,pi_among,homo_deviance])

    df = pd.DataFrame(inbreed_results, columns=["pi_within","pi_among","homo_deviance"])
    df.to_csv(output_csv, sep=",", index=False)
    
    return df