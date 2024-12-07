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

# Define a function to calculate alternative allele frequencies for a single individual.
def calc_ind_alt_freqs(gt):
    # Compute the frequency for each site.
    raw_freqs = np.nansum(gt, axis=2).flatten() / 2
    # Set missing data to Nan
    alt_freqs = np.where(raw_freqs == -1, np.nan, raw_freqs)
    return alt_freqs

# Define a function to calculate pap scores.
def psuedo_ancestry_painting(gt, target, source_1, source_2):
    # Determine the alternative allele frequencies for the target individual.
    t_aaf = calc_ind_alt_freqs(gt.take(target, axis=1))
    # Create a mask for the heterozygous sites.
    het_mask = t_aaf == 0.5
    # If there are no sites heterozygous sites/
    if het_mask.sum() == 0:
        pap = np.array([np.nan, 0])
    # Determine the alternative allele frequencies for the sites of interest.
    t_aaf = t_aaf[het_mask]
    s1_aaf = calc_ind_alt_freqs(gt.take(source_1, axis=1).compress(het_mask, axis=0))
    s2_aaf = calc_ind_alt_freqs(gt.take(source_2, axis=1).compress(het_mask, axis=0))
    # Determine the sites that can be painted.
    s1_ref_s2_alt = (s1_aaf == 0) & (s2_aaf == 1) & (t_aaf == 0.5)
    s1_alt_s2_ref = (s1_aaf == 1) & (s2_aaf == 0) & (t_aaf == 0.5)
    # Compute the pap sites and total possible sites.
    pap = np.array([(s1_ref_s2_alt | s1_alt_s2_ref).sum(), het_mask.sum()])
    return pap

# Define a function to compute the number of pap sites in non-overlapping windows.
def tgp_pap_windows(chromosome, window_size):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary of focal individuals.
    pap_idx_dicc = {
        'ALT': np.array([2347]), 'CHA': np.array([2348]),
        'VIN': np.array([2349]), 'DEN': np.array([2350]),
        'MXL': tgp_meta_df[tgp_meta_df['IND'] == 'NA19664'].index.values,
        'YRI': tgp_meta_df[tgp_meta_df['IND'] == 'NA19190'].index.values,
    }
    # Intialize the pap configurations.
    pap_configs = [
        ('CHA', 'MXL', 'YRI'), ('VIN', 'MXL', 'YRI'),
        ('DEN', 'MXL', 'YRI'), ('ALT', 'MXL', 'YRI'),
    ]
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('tgp_arcs_masked_no_aa', chromosome)
    # Load the windows data frame.
    qc_windows_df = pd.read_csv(f'../windowing/tgp_arcs_masked_no_aa/{window_size}kb_nonoverlapping_variant_windows.csv.gz')
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chromosome]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of windows.
    n_windows = chr_qc_windows_df.shape[0]
    # Intialize a results dictionary.
    pap_dicc = {}
    # For every pap configuration.
    for config in pap_configs:
        # Intialize the results matrix.
        pap_dicc[config] = np.empty((n_windows, 2))
    # For every window...
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # For every pap configuration.
        for config in pap_configs:
            # Unpack.
            t, s1, s2 = config
            # Update the results.
            pap_dicc[config][wind, :] = psuedo_ancestry_painting(
                allel.GenotypeArray(callset[wind_loc]), 
                pap_idx_dicc[t], pap_idx_dicc[s1], pap_idx_dicc[s2],
            )
    # For every pap configuration.
    for config in pap_configs:
        # Unpack.
        t, s1, s2 = config
        # Export the the results matrix.
        np.savetxt(
            f'../muc19_results/tgp_arcs_masked_no_aa/{t.lower()}_{s1.lower()}_{s2.lower()}_pap_counts_chr{chromosome}_{window_size}kb.txt.gz',
            pap_dicc[config], fmt='%d',
        )
    return

# Calculate pap sites in windows.
tgp_pap_windows(chromosome=int(sys.argv[1]), window_size=int(sys.argv[2]))