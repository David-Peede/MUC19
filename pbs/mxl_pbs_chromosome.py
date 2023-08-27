# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###


# Define a function to load genotyope and positions arrays.
def load_gt_pos(prefix, chrom):
    # Intialize the file path.
    path = '../vcf_data/zarr_data/{0}_chr{1}.zarr'.format(prefix, chrom)
    # Load the zarr array.
    zarr_array = zarr.open_group(path, mode='r')
    # Extract the genotype callset.
    callset = zarr_array['{0}/calldata/GT'.format(chrom)]
    # Convert the genotype callset to an array.
    gt = allel.GenotypeArray(callset)
    # Load the positions.
    pos = allel.SortedIndex(zarr_array['{0}/variants/POS'.format(chrom)])
    return gt, pos

# Define a function to calculate PBS.
def calc_pbs(gt, pop_a, pop_b, pop_c):
    # Determine allele counts.
    a_ac = gt.take(pop_a, axis=1).count_alleles()
    b_ac = gt.take(pop_b, axis=1).count_alleles()
    c_ac = gt.take(pop_c, axis=1).count_alleles()
    # Calculate the numerator and denominator for Hudson's Fst estimator.
    a_b_num, a_b_den = allel.hudson_fst(a_ac, b_ac)
    a_c_num, a_c_den = allel.hudson_fst(a_ac, c_ac)
    c_b_num, c_b_den = allel.hudson_fst(c_ac, b_ac)
    # Calculate Fst.
    a_b_raw_fst = a_b_num / a_b_den
    a_c_raw_fst = a_c_num / a_c_den
    c_b_raw_fst = c_b_num / c_b_den
    # Correct for Fst values of 1 that will lead to inf.
    a_b_fst_raw = np.where(a_b_raw_fst == 1, 0.99999, a_b_raw_fst)
    a_c_fst_raw = np.where(a_c_raw_fst == 1, 0.99999, a_c_raw_fst)
    c_b_fst_raw = np.where(c_b_raw_fst == 1, 0.99999, c_b_raw_fst)
    # Correct for negative Fst values.
    a_b_fst = np.where(a_b_fst_raw < 0, 0, a_b_fst_raw)
    a_c_fst = np.where(a_c_fst_raw < 0, 0, a_c_fst_raw)
    c_b_fst = np.where(c_b_fst_raw < 0, 0, c_b_fst_raw)
    # Calculate PBS.
    pbs = (
        ((np.log(1.0 - a_b_fst) * -1.0) +\
         (np.log(1.0 - a_c_fst) * -1.0) -\
         (np.log(1.0 - c_b_fst) * -1.0)) / float(2)
    )
    return pbs, np.nanmean(pbs)

# Define a function to calculate pbs per snp per chromosome.
def pbs_chromosome(chrom):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {
        'MXL_NAT': np.loadtxt('../meta_data/mxl_nat_idx.csv', delimiter=',', dtype=int),
        'MXL_NOT': np.loadtxt('../meta_data/mxl_not_idx.csv', delimiter=',', dtype=int),
    }
    # For every population...
    for pop in ['MXL', 'CEU', 'CHB']:
        # Append the dictionary with sample indicies.
        samp_idx_dicc[pop] = tgp_meta_df[tgp_meta_df['POP'] == pop].index.values
    # Extract the genotype callset and positions.
    gt, _ = load_gt_pos('tgp_mod_arc_anc', chrom)
    # Compute PBS for the entire population.
    pop_pbs, _  = calc_pbs(
        gt=gt, pop_a=samp_idx_dicc['MXL'],
        pop_b=samp_idx_dicc['CHB'], pop_c=samp_idx_dicc['CEU'],
    )
    # Export the results.
    np.savetxt(
        f'./chroms/mxl_pop_chb_ceu_pbs_chr{chrom}.csv.gz',
        [pop_pbs], fmt='%1.15f', delimiter=',', newline='\n',
    )
    # For each ancestry partition...
    for anc in ['NAT', 'NOT']:
        # Compute PBS.
        anc_pbs, _  = calc_pbs(
            gt=gt, pop_a=samp_idx_dicc['MXL'+'_'+anc],
            pop_b=samp_idx_dicc['CHB'],
            pop_c=samp_idx_dicc['CEU'],
        )
        # Export the results.
        np.savetxt(
            f'./chroms/mxl_{anc.lower()}_chb_ceu_pbs_chr{chrom}.csv.gz',
            [anc_pbs], fmt='%1.15f', delimiter=',', newline='\n',
        )
    return


# Calculate pbs.
pbs_chromosome(chrom=int(sys.argv[1]))
