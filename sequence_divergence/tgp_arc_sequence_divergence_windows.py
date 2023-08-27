# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###
### sys.argv[2] = window size ###
### sys.argv[3] = tgp population ###
### sys.argv[4] = archaic ###


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

# Define a function to calculate sequence divergence.
def tgp_arc_divergence(gt, tgp_idx, arc_idx):
    # Determine the sample list.
    samp_list = np.concatenate([tgp_idx, arc_idx])
    # Determine the indicies where all samples are called.
    called_mask = (gt.take(samp_list, axis=1).is_called() == True).all(axis=1)
    # If there are no sites called between all three samples...
    if (called_mask.sum() == 0):
        # Set the results to 0 since we are iterating over QC'ed regions.
        abs_div = np.zeros(tgp_idx.size)
    # Else...
    else:
        # Determine the indicies where we have varibale sites.
        var_mask = gt.take(samp_list, axis=1).compress(called_mask, axis=0).count_alleles().is_variant()
        # If there are no variable sites...
        if (var_mask.sum() == 0):
            # Set the results to 0 since we are iterating over QC'ed regions.
            abs_div = np.zeros(tgp_idx.size)
        # Else...
        else:
            # Intialize lists to store results.
            abs_div = []
            # Determine the alternate allele frequency for the archaic.
            arc_freq = calc_alt_freqs(gt.take(arc_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            # For every individual in the target population...
            for idx in tgp_idx:
                # Determine the alternative allele frequency.
                ind_freq = calc_alt_freqs(gt.take([idx], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
                # Compute the number of pairwise differences.
                pw_diffs = ((arc_freq * (1 - ind_freq)) + (ind_freq * (1 - arc_freq)))
                # Compute the total distance.
                seq_dist = np.nansum(pw_diffs)
                # Append the results.
                abs_div.append(seq_dist)
            # Convert the results to numpy arrays.
            abs_div = np.array(abs_div)
    return abs_div

# Define a function to calculate sequence divergence in windows.
def tgp_arc_divergence_windows(chromosome, window_size, tgp_pop, arc_pop):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize an ordered population list.
    tgp_pop_list = [
        'LWK', 'GWD', 'MSL', 'ESN', 'YRI', # AFR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
        'PEL', 'MXL', 'CLM', 'PUR', # AMR.
    ]
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {
        'ALT': np.array([2347]), 'CHA': np.array([2348]),
        'VIN': np.array([2349]), 'DEN': np.array([2350]),
    }
    # For every population...
    for pop in tgp_pop_list:
        # Append the dictionary with sample indicies.
        samp_idx_dicc[pop] = tgp_meta_df[tgp_meta_df['POP'] == pop].index.values
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
    abs_dist_mat = np.empty((n_windows, samp_idx_dicc[tgp_pop].size))
    # For every window...
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # Append the results matricies.
        abs_dist_mat[wind, :] = tgp_arc_divergence(
            gt=allel.GenotypeArray(callset[wind_loc]),
            tgp_idx=samp_idx_dicc[tgp_pop], arc_idx=samp_idx_dicc[arc_pop],
        )
    # Compile the results file name.
    results_file = '{0}_{1}_chr{2}_{3}kb.csv.gz'.format(
        tgp_pop.lower(), arc_pop.lower(), chromosome, window_size,
    )
    # Export the the results matrix.
    np.savetxt(
        './tgp/windows/pw_diffs_'+results_file,
        abs_dist_mat, fmt='%1.15f', delimiter=',', newline='\n',
    )
    return

# Calculate sequence divergence in non-overlapping windows.
tgp_arc_divergence_windows(
    chromosome=int(sys.argv[1]),
    window_size=int(sys.argv[2]),
    tgp_pop=str(sys.argv[3]),
    arc_pop=str(sys.argv[4]),
)