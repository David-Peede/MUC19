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

# Define a site pattern function.
def site_patterns(p1, p2, p3):
    # Calculate site pattern counts.
    abba = np.sum((1 - p1) * (p2) * (p3))
    baba = np.sum((p1) * (1 - p2) * (p3))
    bbaa = np.sum((p1) * (p2) * (1 - p3))
    baaa = np.sum((p1) * (1 - p2) * (1 - p3))
    abaa = np.sum((1 - p1) * (p2) * (1 - p3))
    aaba = np.sum((1 - p1) * (1 - p2) * (p3))
    return abba, baba, bbaa, baaa, abaa, aaba

# Define a function to calculate alternative allele frequencies for a single individual.
def calc_ind_alt_freqs(gt):
    # Compute the frequency for each site.
    raw_freqs = np.nansum(gt, axis=2).flatten() / 2
    # Set missing data to Nan
    alt_freqs = np.where(raw_freqs == -1, np.nan, raw_freqs)
    return alt_freqs

# Define a function to filter and calculate site pattern counts for the archaics.
def arc_site_patterns(gt, p1, p2, p3):
    # Intialize a dictionary of indicies.
    idx_dicc = {
        'ALT': 0, 'CHA': 1,
        'VIN': 2, 'DEN': 3,
    }
    # Determine the sample list.
    samp_list = [idx_dicc[p1], idx_dicc[p2], idx_dicc[p3]]
    # Determine the indicies where all samples are called.
    called_mask = (gt.take(samp_list, axis=1).is_called() == True).all(axis=1)
    # If there are no sites called between all three samples...
    if (called_mask.sum() == 0):
        # Set the results to 0 since we are iterating over QC'ed regions.
        results = np.zeros(6)
    # Else...
    else:
        # Determine the indicies where we have varibale sites.
        var_mask = gt.take(samp_list, axis=1).compress(called_mask, axis=0).count_alleles().is_segregating()
        # If there are no variable sites...
        if (var_mask.sum() == 0):
            # Set the results to 0 since we are iterating over QC'ed regions.
            results = np.zeros(6)
        # Else...
        else:
            # Calculate the alternative allele frequencies.
            p1_alt_freqs = calc_ind_alt_freqs(gt.take([idx_dicc[p1]], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p2_alt_freqs = calc_ind_alt_freqs(gt.take([idx_dicc[p2]], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p3_alt_freqs = calc_ind_alt_freqs(gt.take([idx_dicc[p3]], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            anc_freqs = calc_ind_alt_freqs(gt.take([-1], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            # Polarize the samples.
            p1_der_freqs = np.where(anc_freqs == 1, np.abs(p1_alt_freqs - 1), p1_alt_freqs)
            p2_der_freqs = np.where(anc_freqs == 1, np.abs(p2_alt_freqs - 1), p2_alt_freqs)
            p3_der_freqs = np.where(anc_freqs == 1, np.abs(p3_alt_freqs - 1), p3_alt_freqs)
            # Calculate the site pattern counts.
            abba, baba, bbaa, baaa, abaa, aaba = site_patterns(
                p1_der_freqs, p2_der_freqs, p3_der_freqs,
            )
            # Intialize a results array.
            results = np.array([abba, baba, bbaa, baaa, abaa, aaba])
    return results

# Define a function to calculate archaic site patterns in windows.
def archaic_site_patterns_windows(chromosome, window_size):
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('arcs_masked_aa', chromosome)
    # Load the windows data frame.
    qc_windows_df = pd.read_csv(f'../windowing/arcs_masked_aa/{window_size}kb_nonoverlapping_variant_windows.csv.gz')
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chromosome]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of windows.
    n_windows = chr_qc_windows_df.shape[0]
    # Intialize a list of site pattern configurations.
    config_list = [
        ('ALT', 'CHA', 'DEN'),
        ('ALT', 'VIN', 'DEN'),
    ]
    # Intialize a dictionary to store the results.
    sp_dicc = {}
    # For every site pattern configuration.
    for config in config_list:
        # Intialize a results matrix to store the results.
        sp_dicc[config] = np.empty((n_windows, 6))
    # For every window.
    for wind in range(n_windows):
            # Locate the window.
            wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
            # For every site pattern configuration.
            for config in config_list:
                # Unpack.
                p1_arc, p2_arc, p3_arc = config
                # Append the result matrix.
                sp_dicc[config][wind, :] = arc_site_patterns(
                    gt=allel.GenotypeArray(callset[wind_loc]),
                    p1=p1_arc, p2=p2_arc, p3=p3_arc,
                )
    # For every site pattern configuration.
    for config in config_list:
        # Unpack.
        p1_arc, p2_arc, p3_arc = config
        # Export the the results matrix.
        np.savetxt(
            f'../muc19_results/arcs_masked_aa/{p1_arc.lower()}_{p2_arc.lower()}_{p3_arc.lower()}_chr{chromosome}_{window_size}kb.txt.gz',
            sp_dicc[config], fmt='%1.15f',
        )
    return

# Calculate site pattern counts in windows.
archaic_site_patterns_windows(
    chromosome=int(sys.argv[1]),
    window_size=int(sys.argv[2]),
)