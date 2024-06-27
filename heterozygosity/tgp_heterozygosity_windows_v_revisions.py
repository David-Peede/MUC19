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

# Define a function to compute per site heterozygosity.
def per_site_heterozygosity(gt):
    # Determine the number of sequences.
    n = gt.shape[1] * 2
    # Compute the alternative allele frequency.
    p = calc_alt_freqs(gt)
    # Compute the per site heterozygosity.
    per_site_h = (n / (n - 1)) * (1 - ((p ** 2) + ((1 - p) ** 2)))
    return per_site_h

# Define a function to calculate the number of hertozygous sites per individual and the heterozygosity.
def tgp_het_windows(chromosome, window_size):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary of sample indicies.
    samp_idx_dicc = {
        'AFR': tgp_meta_df[tgp_meta_df['SUPERPOP'] == 'AFR'].index.values,
        'HET': np.loadtxt('../meta_data/72kb_all_het_int_idx.csv', delimiter=',', dtype=int),
        'HOM': np.loadtxt('../meta_data/72kb_all_hom_int_idx.csv', delimiter=',', dtype=int),
    }
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('tgp_mod_no_aa', chromosome)
    # Load the windows data frame.
    qc_windows_df = pd.read_csv(f'../windowing/tgp_mod_no_aa/{window_size}kb_nonoverlapping_variant_windows.csv.gz')
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chromosome]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of windows.
    n_windows = chr_qc_windows_df.shape[0]
    # Intialize a results matrix and dictionary to store the results.
    heterozygosity_mat = np.empty((n_windows, 3))
    # Intialize a dictionary to store the results.
    het_sites = {}
    # For every focal set.
    for focal_set in samp_idx_dicc:
        # Intialize a results matrix.
        het_sites[focal_set] = np.empty((n_windows, samp_idx_dicc[focal_set].size))
    # For every window.
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # For every focal set.
        for i, focal_set in enumerate(het_sites):
            # Extract the genotype matrix.
            wind_gt = allel.GenotypeArray(callset[wind_loc]).take(samp_idx_dicc[focal_set], axis=1)
            # Determine the number of heterozygous sites per individual.
            het_sites[focal_set][wind, :] = wind_gt.count_het(axis=0)
            # Determine the heterozygosity for the focal set.
            heterozygosity_mat[wind, i] = np.nansum(per_site_heterozygosity(wind_gt))
    # Export the the heterozygosity matrix.
    np.savetxt(
        f'../muc19_results/tgp_mod_no_aa/afr_het_hom_heterozygosity_chr{chromosome}_{window_size}kb.txt.gz',
        heterozygosity_mat, fmt='%1.15f',
    )
    # For every focal set.
    for focal_set in het_sites:
        # Export the number of heterozygous sites per individual.
        np.savetxt(
            f'../muc19_results/tgp_mod_no_aa/{focal_set.lower()}_het_sites_chr{chromosome}_{window_size}kb.txt.gz',
            het_sites[focal_set], fmt='%d',
        )
    return

# Calculate heterozygosity in windows.
tgp_het_windows(chromosome=int(sys.argv[1]), window_size=int(sys.argv[2]))