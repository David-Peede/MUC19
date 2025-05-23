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

# Define a function to calculate alternative allele frequencies for a single individual.
def calc_ind_alt_freqs(gt):
    # Compute the frequency for each site.
    raw_freqs = np.nansum(gt, axis=2).flatten() / 2
    # Set missing data to Nan
    alt_freqs = np.where(raw_freqs == -1, np.nan, raw_freqs)
    return alt_freqs

# Define a function to calculate site patterns between the TGP and archaics.
def tgp_arc_site_patterns(
    gt,
    p1_idx, p2_idx, p3_idx,
):
    # Determine the sample list.
    samp_list = np.concatenate([p1_idx, p2_idx, p3_idx])
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
            p1_alt_freqs = calc_alt_freqs(gt.take(p1_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p2_alt_freqs = calc_alt_freqs(gt.take(p2_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p3_alt_freqs = calc_ind_alt_freqs(gt.take(p3_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
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

# Define a function to calculate site patterns in windows.
def mxl_archaic_site_patterns_windows(chromosome, window_size, archaic):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a list of focal individuals.
    p2_list = ['NA19664']
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {
        f'{archaic}': np.array([2347]),
        'YRI': tgp_meta_df[tgp_meta_df['POP'] == 'YRI'].index.values,
    }
    # For every p2 population.
    for ind in p2_list:
        # Append the dictionary with sample indicies.
        samp_idx_dicc[ind] = tgp_meta_df[tgp_meta_df['IND'] == ind].index.values
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos(f'tgp_{archaic.lower()}_masked_aa', chromosome)
    # Load the windows data frame.
    qc_windows_df = pd.read_csv(f'../windowing/tgp_{archaic.lower()}_masked_aa/{window_size}kb_nonoverlapping_variant_windows.csv.gz')
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chromosome]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of windows.
    n_windows = chr_qc_windows_df.shape[0]
    # Intialize a dictionary to store the results.
    sp_dicc = {}
    # For ever P2.
    for p2 in p2_list:
        # Intialize a a results matrix.
        sp_dicc[p2] = np.empty((n_windows, 6))
    # For every window.
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # For every P2 population.
        for p2 in p2_list:
            # Update the results matrix.
            sp_dicc[p2][wind, :] = tgp_arc_site_patterns(
                gt=allel.GenotypeArray(callset[wind_loc]),
                p1_idx=samp_idx_dicc['YRI'], p2_idx=samp_idx_dicc[p2], p3_idx=samp_idx_dicc[f'{archaic}'],
            )
    # For every P2 population.
    for p2 in p2_list:
        # Export the the results matrix.
        np.savetxt(
            f'../muc19_results/tgp_{archaic.lower()}_masked_aa/yri_{p2.lower()}_{archaic.lower()}_chr{chromosome}_{window_size}kb.txt.gz',
            sp_dicc[p2], fmt='%1.15f',
        )
    return

# Calculate site pattern counts.
mxl_archaic_site_patterns_windows(
    chromosome=int(sys.argv[1]),
    window_size=int(sys.argv[2]),
    archaic=str(sys.argv[3]),
)