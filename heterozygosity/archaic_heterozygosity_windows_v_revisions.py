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

# Define a function to compute per site heterozygosity for a single individual.
def ind_per_site_heterozygosity(gt):
    # Compute the alternative allele frequency.
    p = calc_ind_alt_freqs(gt)
    # Compute the per site heterozygosity.
    per_site_h = (1 - ((p ** 2) + ((1 - p) ** 2)))
    return per_site_h

# Define a function to calculate the number of hertozygous sites per individual and the heterozygosity.
def archaic_het_windows(chromosome, window_size, archaic):
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos(f'{archaic.lower()}_masked_no_aa', chromosome)
    # Load the windows data frame.
    qc_windows_df = pd.read_csv(f'../windowing/{archaic.lower()}_masked_no_aa/{window_size}kb_nonoverlapping_variant_windows.csv.gz')
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chromosome]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of windows.
    n_windows = chr_qc_windows_df.shape[0]
    # Intialize a results matrices to store the results.
    heterozygosity_mat = np.empty((n_windows, 2))
    # For every window.
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # Extract the genotype matrix.
        wind_gt = allel.GenotypeArray(callset[wind_loc])
        # Determine the number of heterozygous sites per individual.
        heterozygosity_mat[wind, 0] = wind_gt.count_het(axis=0)
        # Determine the heterozygosity for the individual.
        heterozygosity_mat[wind, 1] = np.nansum(ind_per_site_heterozygosity(wind_gt))
    # Export the the heterozygosity matrix.
    np.savetxt(
        f'../muc19_results/{archaic.lower()}_masked_no_aa/archaic_het_sites_heterozygosity_chr{chromosome}_{window_size}kb.txt.gz',
        heterozygosity_mat, fmt='%1.15f',
    )
    return

# Calculate heterozygosity in windows.
archaic_het_windows(chromosome=int(sys.argv[1]), window_size=int(sys.argv[2]), archaic=str(sys.argv[3]))