# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###
### sys.argv[2] = window size ###
### sys.argv[3] = den/den-like ind idx ###


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

def afr_alt_den_div(gt, tgp_df, den_like_idx):
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = { 
        'AFR': tgp_df[tgp_df['SUPERPOP'] == 'AFR'].index.values,
        'ALT': np.array([2347]), 'DEN': np.array([den_like_idx]),
    }
    # Concatenate the sample list.
    samp_list = np.concatenate([
        samp_idx_dicc['AFR'], samp_idx_dicc['ALT'], samp_idx_dicc['DEN'],
    ])
    # Determine the indicies where all samples are called.
    called_mask = (gt.take(samp_list, axis=1).is_called() == True).all(axis=1)
    # If there are no sites called between all three samples...
    if (called_mask.sum() == 0):
        # Set the results to 0 since we are iterating over QC'ed regions.
        result = np.zeros(2)
    # Else...
    else:
        # Determine the indicies where we have varibale sites.
        var_mask = gt.take(samp_list, axis=1).compress(called_mask, axis=0).count_alleles().is_variant()
        # If there are no variable sites...
        if (var_mask.sum() == 0):
            # Set the results to 0 since we are iterating over QC'ed regions.
            result = np.zeros(2)
        # Else...
        else:
            # Calculate alternative allele frequencies.
            afr_alt_freq = calc_alt_freqs(gt.take(samp_idx_dicc['AFR'], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            alt_alt_freq = calc_alt_freqs(gt.take(samp_idx_dicc['ALT'], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            den_alt_freq = calc_alt_freqs(gt.take(samp_idx_dicc['DEN'], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            anc_freq = calc_alt_freqs(gt.take([-1], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            # Polarize the samples.
            afr_der_freq = np.where(anc_freq == 1, np.abs(afr_alt_freq - 1), afr_alt_freq)
            alt_der_freq = np.where(anc_freq == 1, np.abs(alt_alt_freq - 1), alt_alt_freq)
            den_der_freq = np.where(anc_freq == 1, np.abs(den_alt_freq - 1), den_alt_freq)
            # Determine the number of sites where the derived allele is segregating at 95%
            # frequncy or higher in AFR, fixed in ALT, and absent in DEN.
            der_div = np.where((afr_der_freq >= 0.95) & (alt_der_freq == 1.0) & (den_der_freq == 0.0))[0].size
            # Determine the number of sites where the derived allele is segregating at 5%
            # frequncy or less in AFR, absent in ALT, and fixed in DEN.
            anc_div = np.where((afr_der_freq <= 0.05) & (alt_der_freq == 0.0) & (den_der_freq == 1.0))[0].size
            # Compile the results.
            result = np.array([der_div, anc_div])
    return result

# Define a function to calculate (AFR, NEA)-DEN divergence in windows.
def afr_alt_den_divergence_windows(chromosome, window_size, den_like_ind):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod_arc.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract all samples and populations.
    tgp_ind = tgp_meta_df['IND'].values
    tgp_pop = tgp_meta_df['POP'].values
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
    # Intialize a list to store results.
    results = np.empty((n_windows, 2))
    # For every window...
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # Append the results matricies.
        results[wind, :] = afr_alt_den_div(
            gt=allel.GenotypeArray(callset[wind_loc]),
            tgp_df=tgp_meta_df,
            den_like_idx=den_like_ind,
        )
    # Compile the results file name.
    results_file = 'afr_alt_{0}_{1}_chr{2}_{3}kb.csv.gz'.format(
        tgp_ind[den_like_ind].lower(),
        tgp_pop[den_like_ind].lower(),
        chromosome,
        window_size,
    )
    # Export the the results matrix.
    np.savetxt(
        './afr_alt_den/windows/'+results_file,
        results, fmt='%1.15f', delimiter=',', newline='\n',
    )
    return

# Calculate divergence in windows.
afr_alt_den_divergence_windows(
    chromosome=int(sys.argv[1]),
    window_size=int(sys.argv[2]),
    den_like_ind=int(sys.argv[3]),
)
