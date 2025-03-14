# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###
### sys.argv[2] = window size ###
### sys.argv[3] = archaic ###


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

# Define a function to compute pairwise differences.
def pwd_per_site(px, py):
    return ((px * (1 - py)) + (py * (1 - px)))

# Define a function to calculate sequence divergence between archaic diplotypes and modern human haplotypes.
def calc_hum_hap_v_arc_dip_diffs(gt):
    # Load the meta data.
    meta_data = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store the results.
    pwd_dicc = {'hap_1': np.array([]), 'hap_2': np.array([])}
    # Intialize the archaic's allele frequency.
    arc_freq = calc_ind_alt_freqs(gt.take([2347], axis=1))
    # For every human sample.
    for samp in range(meta_data.shape[0]):
        # Extract the sample's haplotypes.
        hap_1 = gt[:, samp, 0]
        hap_2 = gt[:, samp, 1]
        # Compute the distance from the archaic diplotype.
        pwd_1 = pwd_per_site(hap_1, arc_freq)
        pwd_2 = pwd_per_site(hap_2, arc_freq)
        # Update dictionaries.
        pwd_dicc['hap_1'] = np.append(pwd_dicc['hap_1'], np.nansum(pwd_1))
        pwd_dicc['hap_2'] = np.append(pwd_dicc['hap_2'], np.nansum(pwd_2))
    return pwd_dicc

# Define a function to calculate haplotype divergence in non-overlapping windows.
def tgp_hap_v_arc_dip_diffs_windows(chromosome, window_size, archaic):
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos(f'tgp_{archaic.lower()}_masked_no_aa', chromosome)
    # Load the windows data frame.
    qc_windows_df = pd.read_csv(f'../windowing/tgp_{archaic.lower()}_masked_no_aa/{window_size}kb_nonoverlapping_variant_windows.csv.gz')
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chromosome]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of windows.
    n_windows = chr_qc_windows_df.shape[0]
    # Intialize a dictionary to store the results.
    pw_diffs_dicc = {
        'hap_1': np.empty((n_windows, meta_df.shape[0])),
        'hap_2': np.empty((n_windows, meta_df.shape[0])),
    }
    # For every window.
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # Compute divergence.
        pw_diffs = calc_hum_hap_v_arc_dip_diffs(allel.GenotypeArray(callset[wind_loc]))
        # For every haplotype.
        for hap in ['hap_1', 'hap_2']:
            # Append the haplotype results.
            pw_diffs_dicc[hap][wind, :] = pw_diffs[hap]
    # For every haplotype.
    for hap in ['hap_1', 'hap_2']:
        # Export the haplotype distances.
        np.savetxt(
            f'../muc19_results/tgp_{archaic.lower()}_masked_no_aa/{archaic.lower()}_{hap}_pw_diffs_chr{chromosome}_{window_size}kb.txt.gz',
            pw_diffs_dicc[hap], fmt='%1.15f',
        )
    return

# Calculate haplotype divergence in non-overlapping windows.
tgp_hap_v_arc_dip_diffs_windows(
    chromosome=int(sys.argv[1]),
    window_size=int(sys.argv[2]),
    archaic=str(sys.argv[3]),
)