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

# Define a function to calculate per indvidual heterozygosity in windows.
def tgp_het_windows(chromosome, window_size):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('tgp_mod_arc_anc', chromosome)
    # Load the windows data frame.
    qc_windows_df = pd.read_csv('../windowing/tgp/{0}kb_nonoverlapping_variant_windows.csv.gz'.format(window_size))
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chromosome]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of windows.
    n_windows = chr_qc_windows_df.shape[0]
    # Intialize a results matrix to store the results.
    results_mat = np.empty((n_windows, tgp_meta_df.index.values.size+4))
    # For every window...
    for wind in range(n_windows):
            # Locate the window.
            wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
            # Determine the site pattern counts.
            het_counts = allel.GenotypeArray(callset[wind_loc]).count_het(axis=0)
            # Append the result matrix.
            results_mat[wind, :] = het_counts[:-1] # Don't care about the ancestral genotype as it is always homozygous.
    # Compile the results file name.
    results_file = 'tgp_arc_het_counts_chr{0}_{1}kb.csv.gz'.format(chromosome, window_size)
    # Export the the results matrix.
    np.savetxt(
        './tgp/windows/'+results_file,
        results_mat, fmt='%d', delimiter=',', newline='\n',
    )
    return


# Calculate heterozygosity in windows.
tgp_het_windows(chromosome=int(sys.argv[1]), window_size=int(sys.argv[2]))
