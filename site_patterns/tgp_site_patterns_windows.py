# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###
### sys.argv[2] = window size ###
### sys.argv[3] = p1 ###
### sys.argv[4] = p2 ###
### sys.argv[5] = p3 ###


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
        var_mask = gt.take(samp_list, axis=1).compress(called_mask, axis=0).count_alleles().is_variant()
        # If there are no variable sites...
        if (var_mask.sum() == 0):
            # Set the results to 0 since we are iterating over QC'ed regions.
            results = np.zeros(6)
        # Else...
        else:
            # Calculate the alternative allele frequencies.
            p1_alt_freqs = calc_alt_freqs(gt.take(p1_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p2_alt_freqs = calc_alt_freqs(gt.take(p2_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p3_alt_freqs = calc_alt_freqs(gt.take(p3_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            anc_freqs = calc_alt_freqs(gt.take([-1], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
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
def tgp_arc_site_patterns_windows(chromosome, window_size, p1_pop, p2_pop, p3_pop):
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
    results_mat = np.empty((n_windows, 6))
    # For every window...
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # Append the result matrix.
        results_mat[wind, :] = tgp_arc_site_patterns(
            gt=allel.GenotypeArray(callset[wind_loc]),
            p1_idx=samp_idx_dicc[p1_pop], p2_idx=samp_idx_dicc[p2_pop], p3_idx=samp_idx_dicc[p3_pop],
        )
    # Compile the results file name.
    results_file = '{0}_{1}_{2}_chr{3}_{4}kb.csv.gz'.format(
        p1_pop.lower(), p2_pop.lower(), p3_pop.lower(), chromosome, window_size,
    )
    # Export the the results matrix.
    np.savetxt(
        './tgp/windows/'+results_file,
        results_mat, fmt='%1.15f', delimiter=',', newline='\n',
    )
    return

# Calculate site pattern counts in windows.
tgp_arc_site_patterns_windows(
    chromosome=int(sys.argv[1]),
    window_size=int(sys.argv[2]),
    p1_pop=str(sys.argv[3]),
    p2_pop=str(sys.argv[4]),
    p3_pop=str(sys.argv[5]),
)
