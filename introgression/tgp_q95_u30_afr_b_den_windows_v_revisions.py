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

# Define a funct to compute Q95_ABC(w, y).
def Q95_ABC(gt, A, B, C, w, y):
    # Calculate alternative allele frequencies.
    A_alt_freq = calc_alt_freqs(gt.take(A, axis=1))
    B_alt_freq = calc_alt_freqs(gt.take(B, axis=1))
    C_alt_freq = calc_ind_alt_freqs(gt.take(C, axis=1))
    # Intialize conditions.
    A_ref = (1 - A_alt_freq) < w
    A_alt = A_alt_freq < w
    C_ref = (1 - C_alt_freq) == y
    C_alt = C_alt_freq == y
    Q_ref = A_ref & C_ref
    Q_alt = A_alt & C_alt
    # If there are no archaic sites.
    if (Q_ref | Q_alt).sum() == 0:
        return 0
    # Else there are archaic sites.
    else:
        # Subset the archaic alleles.
        Q_freqs = np.concatenate([B_alt_freq[Q_alt], (1 - B_alt_freq)[Q_ref]])
        # Find the 95% percentile.
        Q95_thresh = np.percentile(Q_freqs, 95)
        return Q95_thresh

# Define a funct to compute U_ABC(w, x, y).
def U_ABC(gt, A, B, C, w, x, y):
    # Calculate alternative allele frequencies.
    A_alt_freq = calc_alt_freqs(gt.take(A, axis=1))
    B_alt_freq = calc_alt_freqs(gt.take(B, axis=1))
    C_alt_freq = calc_ind_alt_freqs(gt.take(C, axis=1))
    # Intialize conditions.
    A_ref = (1 - A_alt_freq) < w
    A_alt = A_alt_freq < w
    B_ref = (1 - B_alt_freq) > x
    B_alt = B_alt_freq > x
    C_ref = (1 - C_alt_freq) == y
    C_alt = C_alt_freq == y
    # Determine the U sites.
    U_ref = A_ref & B_ref & C_ref
    U_alt = A_alt & B_alt & C_alt
    U_sites = U_ref | U_alt
    return U_sites.sum()


# Define a function to compute Q95_{AFR,B,ARC}(1%, 100%) and U_{AFR,B,ARC}(1%, 20%/30%, 100%) statistics in non-overlapping windows.
def tgp_Q95_U20_U30_windows(chromosome, window_size, archaic):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {f'{archaic}': np.array([2347])}
    # Intialize a list of focal populations.
    ooa_list = [
        'MXL', 'PEL', 'CLM', 'PUR', # AMR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
    ]
    # Update the dictionary with the AFR indicies.
    samp_idx_dicc['AFR'] = tgp_meta_df[tgp_meta_df['SUPERPOP'] == 'AFR'].index.values
    # For every OOA population.
    for pop in ooa_list:
        # Update the dictionary.
        samp_idx_dicc[pop] = tgp_meta_df[tgp_meta_df['POP'] == pop].index.values
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos(f'tgp_{archaic.lower()}_masked_no_aa', chromosome)
    # Load the windows data frame.
    qc_windows_df = pd.read_csv(f'../windowing/tgp_{archaic.lower()}_masked_no_aa/{window_size}kb_nonoverlapping_variant_windows.csv.gz')
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chromosome]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of windows.
    n_windows = chr_qc_windows_df.shape[0]
    # Intialize a results matrix to store the results.
    Q95_U_dicc = {pop: np.empty((n_windows, 3)) for pop in ooa_list}
    # For every window.
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # For every OOA population.
        for pop in ooa_list:
            # Compute the Q95 stats.
            Q95 = Q95_ABC(
                gt=allel.GenotypeArray(callset[wind_loc]),
                A=samp_idx_dicc['AFR'], B=samp_idx_dicc[pop], C=samp_idx_dicc[f'{archaic}'],
                w=0.01, y=1,
            )
            # Compute the U_ABC stats.
            U20 = U_ABC(
                gt=allel.GenotypeArray(callset[wind_loc]),
                A=samp_idx_dicc['AFR'], B=samp_idx_dicc[pop], C=samp_idx_dicc[f'{archaic}'],
                w=0.01, x=0.2, y=1,
            )
            U30 = U_ABC(
                gt=allel.GenotypeArray(callset[wind_loc]),
                A=samp_idx_dicc['AFR'], B=samp_idx_dicc[pop], C=samp_idx_dicc[f'{archaic}'],
                w=0.01, x=0.3, y=1,
            )
            # Update the results.
            Q95_U_dicc[pop][wind, :] = np.array([Q95, U20, U30])
    # For every OOA population.
    for pop in ooa_list:
        # Export the the results matrix.
        np.savetxt(
            f'../muc19_results/tgp_{archaic.lower()}_masked_no_aa/q95_u20_u30_afr_{pop.lower()}_{archaic.lower()}_chr{chromosome}_{window_size}kb.txt.gz',
            Q95_U_dicc[pop], fmt='%1.15f',
        )
    return

# Compute the Q95 and U-statistics.
tgp_Q95_U20_U30_windows(
    chromosome=int(sys.argv[1]),
    window_size=int(sys.argv[2]),
    archaic=str(sys.argv[3]),
)