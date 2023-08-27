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

# Define a function to calculate alternative allele frequencies for the archaics
def calc_arc_alt_freqs(arc_gt):
    # Compute the frequency for each site.
    raw_freqs = np.nansum(arc_gt, axis=2).flatten() / 2
    # Set missing data to Nan
    alt_freqs = np.where(raw_freqs == -1, np.nan, raw_freqs)
    return alt_freqs

# Define a function to compute sequence divergence.
def sequence_divergence(px, py):
    # Calculate sequence divergence.
    seq_div = np.nansum(((px * (1 - py)) + (py * (1 - px))))
    return seq_div

# Define a function to polarize a genotype matrix.
def polarize_gt(gt):
    # Extract and determine the ancestral sequence.
    anc_seq = np.mean(gt[:, -1, :], axis=1)
    # Intialize an empty genotype matrix.
    p_gt = np.empty_like(gt[:, :-1, :])
    # For every sample...
    for samp in range(gt[:, :-1, :].shape[1]):
        # Extract the haplotypes.
        hap_1 = gt[:, samp, 0]
        hap_2 = gt[:, samp, 1]
        # Polarize.
        p_hap_1 = np.where(anc_seq == 1, np.abs(hap_1 - 1), hap_1)
        p_hap_2 = np.where(anc_seq == 1, np.abs(hap_2 - 1), hap_2)
        # Build the polarized genotype matrix.
        p_gt[:, samp, 0] = np.where(p_hap_1 == 2, -1, p_hap_1)
        p_gt[:, samp, 1] = np.where(p_hap_2 == 2, -1, p_hap_2)
    # Intialize the polarized genotype matrix into sci-kit allele.
    polarized_gt = allel.GenotypeArray(p_gt)
    return polarized_gt

# Define a function to calculate sequence divergence between archaic diplotypes and tgp haplotypes.
def tgp_haplotype_distances(p_gt):
    # Load the meta data.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize archaic information dictionaries.
    freq_dicc = {
        'DEN': calc_arc_alt_freqs(p_gt.take([2350], axis=1)),
        'ALT': calc_arc_alt_freqs(p_gt.take([2347], axis=1)),
        'CHA': calc_arc_alt_freqs(p_gt.take([2348], axis=1)),
        'VIN': calc_arc_alt_freqs(p_gt.take([2349], axis=1)),
    }
    het_array = np.array([])
    pwd_dicc = {
        'DEN': {'hap_1': np.array([]), 'hap_2': np.array([])},
        'ALT': {'hap_1': np.array([]), 'hap_2': np.array([])},
        'CHA': {'hap_1': np.array([]), 'hap_2': np.array([])},
        'VIN': {'hap_1': np.array([]), 'hap_2': np.array([])},
    }
    fxd_dicc = {
        'DEN': {'hap_1': np.array([]), 'hap_2': np.array([])},
        'ALT': {'hap_1': np.array([]), 'hap_2': np.array([])},
        'CHA': {'hap_1': np.array([]), 'hap_2': np.array([])},
        'VIN': {'hap_1': np.array([]), 'hap_2': np.array([])},
    }
    # Extract haplotype arrays.
    tgp_hap_1 = p_gt[:, :-4, 0]
    tgp_hap_2 = p_gt[:, :-4, 1]
    # For every archaic individual...
    for arc in ['DEN', 'ALT', 'CHA', 'VIN']:
        # Create a mask for heterozygous sites.
        het_mask = (freq_dicc[arc] == 0.5)
        # Update the number of heterozygous sites in this window.
        het_array = np.append(het_array, het_mask.sum())
        # For every tgp sample...
        for samp in range(p_gt[:, :-4, :].shape[1]):
            # Extract the sample's haplotypes.
            hap_1 = tgp_hap_1[:, samp]
            hap_2 = tgp_hap_2[:, samp]
            # Compute the distance from the archaic diplotype.
            pwd_1 = sequence_divergence(hap_1, freq_dicc[arc])
            pwd_2 = sequence_divergence(hap_2, freq_dicc[arc])
            fxd_1 = sequence_divergence(hap_1[~het_mask], freq_dicc[arc][~het_mask])
            fxd_2 = sequence_divergence(hap_2[~het_mask], freq_dicc[arc][~het_mask])
            # Update dictionaries.
            pwd_dicc[arc]['hap_1'] = np.append(pwd_dicc[arc]['hap_1'], pwd_1)
            pwd_dicc[arc]['hap_2'] = np.append(pwd_dicc[arc]['hap_2'], pwd_2)
            fxd_dicc[arc]['hap_1'] = np.append(fxd_dicc[arc]['hap_1'], fxd_1)
            fxd_dicc[arc]['hap_2'] = np.append(fxd_dicc[arc]['hap_2'], fxd_2)
    return het_array, pwd_dicc, fxd_dicc

# Define a function to calculate haplotype divergence in windows.
def hap_divergence_windows(chromosome, window_size):
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
    # Intialize a list of archaics.
    arc_list = ['DEN', 'ALT', 'CHA', 'VIN']
    # Intialize a results matrix to store the number of heterozygous sites.
    het_mat = np.empty((n_windows, len(arc_list)))
    # Intialize dictionaries to store the results.
    pw_diffs_dicc = {}
    fx_diffs_dicc = {}
    # For every archaic...
    for arc in arc_list:
        # Intialize the subdictionaries.
        pw_diffs_dicc[arc] = {}
        fx_diffs_dicc[arc] = {}
        # For every haplotype.
        for hap in ['hap_1', 'hap_2']:
            # Intialize the results matrix.
            pw_diffs_dicc[arc][hap] = np.empty((n_windows, tgp_meta_df.shape[0]))
            fx_diffs_dicc[arc][hap] = np.empty((n_windows, tgp_meta_df.shape[0]))
    # For every window...
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # Polarize the genotype matrix.
        p_gt = polarize_gt(gt=allel.GenotypeArray(callset[wind_loc]))
        # Compute haplotype distances.
        hets, pw_diffs, fx_diffs = tgp_haplotype_distances(p_gt=p_gt)
        # Append the het sites results.
        het_mat[wind, :] = hets # [DEN, ALT, CHA, VIN]
        # For every archaic...
        for arc in arc_list:
            # For every haplotype.
            for hap in ['hap_1', 'hap_2']:
                # Append the haplotype results.
                pw_diffs_dicc[arc][hap][wind, :] = pw_diffs[arc][hap]
                fx_diffs_dicc[arc][hap][wind, :] = fx_diffs[arc][hap]
    # Export the the het sites results matrix.
    np.savetxt(
        f'./tgp/windows/arc_het_sites_chr{chromosome}_{window_size}kb.csv.gz',
        het_mat, fmt='%d', delimiter=',', newline='\n',
    )
    # For every archaic...
    for arc in arc_list:
        # For every haplotype.
        for hap in ['hap_1', 'hap_2']:
            # Export the haplotype distances.
            np.savetxt(
                f'./tgp/windows/{arc.lower()}_{hap}_pw_diffs_chr{chromosome}_{window_size}kb.csv.gz',
                pw_diffs_dicc[arc][hap], fmt='%1.15f', delimiter=',', newline='\n',
            )
            np.savetxt(
                f'./tgp/windows/{arc.lower()}_{hap}_fx_diffs_chr{chromosome}_{window_size}kb.csv.gz',
                fx_diffs_dicc[arc][hap], fmt='%1.15f', delimiter=',', newline='\n',
            )
    return

# Calculate haplotype divergence in non-overlapping windows.
hap_divergence_windows(
    chromosome=int(sys.argv[1]),
    window_size=int(sys.argv[2]),
)