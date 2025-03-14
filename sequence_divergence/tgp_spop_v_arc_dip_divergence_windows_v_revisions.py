# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###
### sys.argv[2] = window size ###
### sys.argv[3] = super population ###
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

# Define a function to calculate alternative allele frequencies.
def calc_alt_freqs(gt):
    # If there are no altenative alleles...
    if (gt.count_alleles().shape[1] == 1):
        # Calculate alternative allele frequencies.
        alt_freqs = gt.count_alleles().to_frequencies()[:, 0] - 1
    # Else...
    else:
        # Calculate alternative allele frequencies.
        alt_freqs = gt.count_alleles().to_frequencies()[:, 1]
    return alt_freqs

# Define a function to compute pairwise differences.
def pwd_per_site(px, py):
    return ((px * (1 - py)) + (py * (1 - px)))

# Define a function to calculate haplotype divergence in non-overlapping windows.
def tgp_spop_v_arc_dip_diffs_windows(chromosome, window_size, spop, archaic):
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the super population indicies.
    spop_idx = meta_df[meta_df['SUPERPOP'] == spop].index.values
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
    # Intialize an array to store the results.
    avg_pwd = np.empty(n_windows)
    # For every window.
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # Subset the genotype matrix.
        wind_gt = allel.GenotypeArray(callset[wind_loc])
        # Compute alternative allele frequencies.
        spop_aaf = calc_alt_freqs(wind_gt.take(spop_idx, axis=1))
        arc_aaf = calc_ind_alt_freqs(wind_gt.take([2347], axis=1))
        # Compute divergence and update the results.
        avg_pwd[wind] = np.nansum(pwd_per_site(spop_aaf, arc_aaf))
    # Export the haplotype distances.
    np.savetxt(
        f'../muc19_results/tgp_{archaic.lower()}_masked_no_aa/{spop.lower()}_{archaic.lower()}_avg_pw_diffs_chr{chromosome}_{window_size}kb.txt.gz',
        [avg_pwd], fmt='%1.15f',
    )
    return

# Calculate sequence divergence in non-overlapping windows.
tgp_spop_v_arc_dip_diffs_windows(
    chromosome=int(sys.argv[1]),
    window_size=int(sys.argv[2]),
    spop=str(sys.argv[3]),
    archaic=str(sys.argv[4]),
)