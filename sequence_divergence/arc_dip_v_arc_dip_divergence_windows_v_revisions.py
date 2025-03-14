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

# Define a function to compute pairwise differences.
def pwd_per_site(px, py):
    return ((px * (1 - py)) + (py * (1 - px)))

# Define a function to calculate sequence divergence between archaic diplotypes.
def calc_arc_dip_v_arc_dip_diffs(gt):
    # Intialize a list of archaics and there indicies.
    arc_idx_dicc = {
        'DEN': 3, 'ALT': 0,
        'CHA': 1, 'VIN': 2,
    }
    # Intialize a list of all combinations.
    arc_combos = [
        ('DEN', 'ALT'), ('DEN', 'CHA'), ('DEN', 'VIN'),
        ('ALT', 'CHA'), ('ALT', 'VIN'), ('CHA', 'VIN'),
    ]
    # Intialize data structures to store the results.
    pwd_dicc = {}
    # For every archaic combination.
    for combo in arc_combos:
        # Unpack the two archaics.
        arc_1, arc_2 = combo
        # Determine the sample list.
        samp_list = [arc_idx_dicc[arc_1], arc_idx_dicc[arc_2]]
        # Determine the indicies where all samples have variable called sites.
        called_mask = (gt.take(samp_list, axis=1).is_called() == True).all(axis=1)
        # If there are no variable sites called between the samples...
        if (called_mask.sum() == 0):
            # Fill the dictionary.
            pwd_dicc[combo] = np.nan
        # Else...
        else:
            # Determine the alternative allele frequencies per archaic.
            arc_1_freq = calc_ind_alt_freqs(gt.take([arc_idx_dicc[arc_1]], axis=1).compress(called_mask, axis=0))
            arc_2_freq = calc_ind_alt_freqs(gt.take([arc_idx_dicc[arc_2]], axis=1).compress(called_mask, axis=0))
            # Compute the number of pairwise differences.
            pw_diffs = pwd_per_site(arc_1_freq, arc_2_freq)
            # Fill the dictionary.
            pwd_dicc[combo] = np.nansum(pw_diffs)
    return pwd_dicc

# Define a function to calculate sequence divergence in windows.
def arc_dip_v_arc_dip_diffs_windows(chromosome, window_size):
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('arcs_masked_no_aa', chromosome)
    # Load the windows data frame.
    qc_windows_df = pd.read_csv(f'../windowing/arcs_masked_no_aa/{window_size}kb_nonoverlapping_variant_windows.csv.gz')
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chromosome]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of windows.
    n_windows = chr_qc_windows_df.shape[0]
    # Intialize a list of all combinations.
    arc_pw_combos = [
        ('DEN', 'ALT'), ('DEN', 'CHA'), ('DEN', 'VIN'),
        ('ALT', 'CHA'), ('ALT', 'VIN'), ('CHA', 'VIN'),
    ]
    # Intialize a results matricies to store the results.
    pwd_mat_dicc = {combo: np.empty(n_windows) for combo in arc_pw_combos}
    # For every window...
    for wind in range(n_windows):
            # Locate the window.
            wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
            # Determine the site pattern counts.
            pw_diffs = calc_arc_dip_v_arc_dip_diffs(allel.GenotypeArray(callset[wind_loc]))
            # For every archaic combination.
            for combo in arc_pw_combos:
                # Update the results.
                pwd_mat_dicc[combo][wind] = pw_diffs[combo]
    # For every archaic combination.
    for combo in arc_pw_combos:
        # Unpack.
        arc1, arc2 = combo
        # Export.
        np.savetxt(
            f'../muc19_results/arcs_masked_no_aa/{arc1.lower()}_v_{arc2.lower()}_pw_diffs_chr{chromosome}_{window_size}kb.txt.gz',
                [pwd_mat_dicc[combo]], fmt='%1.15f',
            )
    return

# Calculate sequence divergence in windows.
arc_dip_v_arc_dip_diffs_windows(
    chromosome=int(sys.argv[1]),
    window_size=int(sys.argv[2]),
)