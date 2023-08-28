# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###
### sys.argv[2] = window size ###


# Define a function to load genotyope and positions arrays.
def load_callset_pos(prefix, chrom):
    # Intialize the file path.
    path = '../vcf_data/zarr_data/{0}_chr{1}.zarr'.format(prefix, chrom)
    # Load the zarr array.
    zarr_array = zarr.open_group(path, mode='r')
    # Extract the genotype callset.
    callset = zarr_array['{0}/calldata/GT'.format(chrom)]
    # Load the positions.
    pos = allel.SortedIndex(zarr_array['{0}/variants/POS'.format(chrom)])
    return callset, pos

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

# Define a function to calculate a windowed average pbs.
def pbs_windows(chromosome, window_size):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {}
    # For every population...
    for pop in ['MXL', 'CEU', 'CHB']:
        # Append the dictionary with sample indicies.
        samp_idx_dicc[pop] = tgp_meta_df[tgp_meta_df['POP'] == pop].index.values
    # Load the windows data frame.
    qc_windows_df = pd.read_csv('../windowing/tgp/{0}kb_nonoverlapping_variant_windows.csv.gz'.format(window_size))
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chromosome]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of windows.
    n_windows = chr_qc_windows_df.shape[0]
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('tgp_mod_arc_anc', chromosome)
    # Intialize a list to store the results.
    pbs_list = []
    # For every window...
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # Compute the windowed average PBS for the entire population.
        _, pop_pbs  = calc_pbs(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=samp_idx_dicc['MXL'],
            pop_b=samp_idx_dicc['CHB'], pop_c=samp_idx_dicc['CEU'],
        )
        # Append the results.
        pbs_list.append(pop_pbs)
    # Export the results.
    np.savetxt(
        f'./windows/mxl_pop_chb_ceu_pbs_chr{chromosome}_{window_size}kb.csv.gz',
        np.array(pbs_list), fmt='%1.15f', delimiter=',', newline='\n',
    )
    return


# Calculate pbs per non-overlapping window.
pbs_windows(chromosome=int(sys.argv[1]), window_size=int(sys.argv[2]))
