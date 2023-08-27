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
def arc_distance(gt):
    # Intialize a list of archaics and there indicies.
    arc_idx_dicc = {
        'DEN': 3, 'ALT': 0,
        'CHA': 1, 'VIN': 2,
    }
    # Intialize a list of all combinations.
    arc_combos = [('DEN', 'ALT'), ('DEN', 'CHA'), ('DEN', 'VIN')]
    # Intialize a dictionary for distance results.
    dist_dicc = {
        'DEN': {'ALT': {}, 'CHA': {}, 'VIN': {}},
    }
    # For every combo...
    for idx in range(len(arc_combos)):
        # Grab the two archaics.
        arc_1 = arc_combos[idx][0]
        arc_2 = arc_combos[idx][1]
        # Determine the sample list.
        samp_list = [arc_idx_dicc[arc_1], arc_idx_dicc[arc_2]]
        # Determine the indicies where all samples have variable called sites.
        called_mask = (gt.take(samp_list, axis=1).is_called() == True).all(axis=1)
        # If there are no variable sites called between the samples...
        if (called_mask.sum() == 0):
            dist_dicc[arc_1][arc_2]['ABS'] = 0
        # Else...
        else:
            # Determine the alternative allele frequencies per archaic.
            arc_1_freq = calc_alt_freqs(gt.take([arc_idx_dicc[arc_1]], axis=1).compress(called_mask, axis=0))
            arc_2_freq = calc_alt_freqs(gt.take([arc_idx_dicc[arc_2]], axis=1).compress(called_mask, axis=0))
            # Compute the number of pairwise differences.
            pw_diffs = ((arc_1_freq * (1 - arc_2_freq)) + (arc_2_freq * (1 - arc_1_freq)))
            # Compute the total distance.
            seq_dist = np.nansum(pw_diffs)
            # Fill the dictionary.
            dist_dicc[arc_1][arc_2]['ABS'] = seq_dist
    return dist_dicc

# Define a function to calculate sequence divergence in windows.
def den_nea_div_windows(chromosome, window_size):
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('arc_anc', chromosome)
    # Load the windows data frame.
    qc_windows_df = pd.read_csv('../windowing/arc/{0}kb_nonoverlapping_variant_windows.csv.gz'.format(window_size))
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chromosome]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of windows.
    n_windows = chr_qc_windows_df.shape[0]
    # Intialize a results matrix to store the results.
    results_mat = np.empty((n_windows, 3))
    # For every window...
    for wind in range(n_windows):
            # Locate the window.
            wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
            # Determine the site pattern counts.
            all_arcs_dist = arc_distance(allel.GenotypeArray(callset[wind_loc]))
            # Compile the divergence results.
            arc_div = []
            for key_1 in all_arcs_dist.keys():
                for key_2 in all_arcs_dist[key_1].keys():
                    arc_div.append(all_arcs_dist[key_1][key_2]['ABS'])
            # Append the result matrix.
            results_mat[wind, :] = np.array(arc_div)
    # Compile the results file name.
    results_file = 'den_nea_div_chr{0}_{1}kb.csv.gz'.format(chromosome, window_size)
    # Export the the results matrix.
    np.savetxt(
        './arc/windows/'+results_file,
        results_mat, fmt='%1.15f', delimiter=',', newline='\n',
    )
    return

# Calculate sequence divergence in windows.
den_nea_div_windows(chromosome=int(sys.argv[1]), window_size=int(sys.argv[2]))