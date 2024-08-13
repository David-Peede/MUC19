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
    path = f'../zarr_data/{prefix}_chr{chrom}.zarr'
    # Load the zarr array.
    zarr_array = zarr.open_group(path, mode='r')
    # Extract the genotype callset.
    callset = zarr_array[f'{chrom}/calldata/GT']
    # Load the positions.
    pos = allel.SortedIndex(zarr_array[f'{chrom}/variants/POS'])
    return callset, pos

# Define a function to compute Hudson's estimator of Fst per site
def calc_fst_per_site(ac1, ac2):
    # Compute.
    num, den = allel.hudson_fst(ac1, ac2)
    raw_fst = num / den
    # Correct for when the denominator is zero.
    raw_fst = np.where(den == 0, 0, raw_fst)
    # Correct for negative Fst values.
    raw_fst = np.where(raw_fst < 0, 0, raw_fst)
    # Correct for Fst values of 1 that will lead to inf.
    raw_fst = np.where(raw_fst == 1, 0.99999, raw_fst)
    return raw_fst

# Define a function to calculate PBS.
def calc_pbs_per_site(gt, pop_a, pop_b, pop_c):
    # Determine the sample list.
    samp_list = np.concatenate([pop_a, pop_b, pop_c])
    # Determine if any of the sites are segregating.
    seg_mask = gt.take(samp_list, axis=1).count_alleles().is_segregating()
    # Determine allele counts.
    a_ac = gt.take(pop_a, axis=1).count_alleles()
    b_ac = gt.take(pop_b, axis=1).count_alleles()
    c_ac = gt.take(pop_c, axis=1).count_alleles()
    # Calculate Fst with corrections.
    a_b_fst = calc_fst_per_site(a_ac, b_ac)
    a_c_fst = calc_fst_per_site(a_ac, c_ac)
    c_b_fst = calc_fst_per_site(c_ac, b_ac)
    # Calculate PBS.
    pbs = (
        ((np.log(1.0 - a_b_fst) * -1.0) +\
         (np.log(1.0 - a_c_fst) * -1.0) -\
         (np.log(1.0 - c_b_fst) * -1.0)) / 2.0
    )
    # Set the sites that are not segregating to undefined.
    pbs[~seg_mask] = np.nan
    return pbs

# Define a function to load the genome-wide PBS values and compute the critical threshold.
def load_gw_crit_thresh_amr_chb_ceu_pbs(critical_threshold):
    # Intialize a list and dictionary.
    all_pbs = []
    pbs_99th_percentile = {}
    # For all chromosomes...
    for chrom in range(1, 23):
        # Load and append the pbs data.
        all_pbs.append(
            np.loadtxt(
                f'../muc19_results/tgp_mod_no_aa/amr_chb_ceu_pbs_chr{chrom}.txt.gz',
        ))
    # Concatenate the chromosome results.
    pbs_chroms = np.concatenate(all_pbs)
    # For every population.
    for i, pop in enumerate(['PEL', 'CLM', 'PUR']):
        # Set negative values to zero.
        cleaned_pbs = np.where(pbs_chroms[:, i] < 0, 0, pbs_chroms[:, i])
        # Update the dictionary.
        pbs_99th_percentile[pop] = np.nanpercentile(cleaned_pbs, critical_threshold)
    return pbs_99th_percentile


# Define a function to determine outlier PBS_{AMR:CHB:CEU} clusters in non-overlapping windows.
def pbs_outlier_props_per_window(chrom, window_size):
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {}
    # For every population...
    for pop in ['PEL', 'CLM', 'PUR', 'CEU', 'CHB']:
        # Append the dictionary with sample indicies.
        samp_idx_dicc[pop] = meta_df[meta_df['POP'] == pop].index.values
    # Compute the critical thrersholds.
    crit_dicc = load_gw_crit_thresh_amr_chb_ceu_pbs(99.95)
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('tgp_mod_no_aa', chrom)
    # Load the windows data frame.
    qc_windows_df = pd.read_csv(f'../windowing/tgp_mod_no_aa/{window_size}kb_nonoverlapping_variant_windows.csv.gz')
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chrom]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of winds.
    n_winds = chr_qc_windows_df.shape[0]
    # Intialize a dictionary to store the results.
    results_dicc = {}
    # For every ancestry partition.
    for amr in ['PEL', 'CLM', 'PUR']:
        # Intialize a matrix to store the results.
        results_dicc[amr] = np.empty((n_winds, 3))
    # For ever window.
    for wind in range(n_winds):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # For every ancestry partition.
        for i, amr in enumerate(['PEL', 'CLM', 'PUR']):
            # Compute PBS.
            pbs_per_site = calc_pbs_per_site(
                gt=allel.GenotypeArray(callset[wind_loc]), pop_a=samp_idx_dicc[amr],
                pop_b=samp_idx_dicc['CHB'], pop_c=samp_idx_dicc['CEU'],
            )
            # Mask undefined PBS values.
            cleaned_pbs_per_site = pbs_per_site[~np.isnan(pbs_per_site)]
            # Determine the number of defined PBS values.
            n_pbs_per_site = cleaned_pbs_per_site.size
            # If there is at least one defined pbs value.
            if n_pbs_per_site > 0:
                # Determine the number of sites above the critical threshold.
                n_crit_pbs = (cleaned_pbs_per_site > crit_dicc[amr]).sum()
                # Determine the proportion of PBS values above the critical threshold.
                n_crit_prop = n_crit_pbs / n_pbs_per_site
                # Update the results.
                results_dicc[amr][wind, :] = np.array([n_pbs_per_site, n_crit_pbs, n_crit_prop])
            # Else, there are no defined PBS vlaues.
            else:
                # Update the results.
                results_dicc[amr][wind, :] = np.array([0, 0, np.nan])
    # For every ancestry partition.
    for amr in ['PEL', 'CLM', 'PUR']:
        # Export the results.
        np.savetxt(
            f'../muc19_results/tgp_mod_no_aa/{amr.lower()}_chb_ceu_pbs_outlier_props_chr{chrom}_{window_size}kb.txt.gz',
            results_dicc[amr], fmt='%1.15f',
        )
    return


# Idenitify the PBS values critical thresholds per window.
pbs_outlier_props_per_window(chrom=int(sys.argv[1]), window_size=int(sys.argv[2]))