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

# Define a function to calculate alternative allele frequencies for the archaics
def calc_arc_alt_freqs(arc_gt):
    # Compute the frequency for each site.
    raw_freqs = np.nansum(arc_gt, axis=2).flatten() / 2
    # Set missing data to Nan
    alt_freqs = np.where(raw_freqs == -1, np.nan, raw_freqs)
    return alt_freqs

# Define a function to compute sequence divergence between haplotype groups.
def hap_mat_diffs(hap_mat_x, hap_mat_y):
    # Extract the derived allele frequency arrays.
    px = hap_mat_x.sum(axis=1) / hap_mat_x.shape[1]
    # If the y group is an archaic...
    if hap_mat_y.ndim == 1:
        py = hap_mat_y
    # Else...
    else:
        py = hap_mat_y.sum(axis=1) / hap_mat_y.shape[1]
    # Compute the number of pairwise differences. 
    pw_diffs = np.nansum(((px * (1 - py)) + (py * (1 - px))))
    return pw_diffs


# Define a function to compute haplo-group differences.
def haplo_group_diffs(p_gt, haplo_groups):
    # Extract the denisovan-like haplotypes.
    arc_hap_mat = np.concatenate((
        p_gt[:, haplo_groups['ARC']['hap_1'], 0],
        p_gt[:, haplo_groups['ARC']['hap_2'], 1],
    ), axis=1)
    # Extract the human-like haplotypes.
    hum_hap_mat = np.concatenate((
        p_gt[:, haplo_groups['HUM']['hap_1'], 0],
        p_gt[:, haplo_groups['HUM']['hap_2'], 1],
    ), axis=1)
    # Extract the recombinant haplotypes.
    rec_hap_mat = np.concatenate((
        p_gt[:, haplo_groups['REC']['hap_1'], 0],
        p_gt[:, haplo_groups['REC']['hap_2'], 1],
    ), axis=1)
    # Create a dictionary of haplotype matricies.
    hap_mat_dicc = {
        'ARC': arc_hap_mat,
        'HUM': hum_hap_mat,
        'REC': rec_hap_mat,
    }
    # Intialize archaic diplotype dictionary.
    arc_dip_dicc = {
        'DEN': calc_arc_alt_freqs(p_gt.take([2350], axis=1)),
        'ALT': calc_arc_alt_freqs(p_gt.take([2347], axis=1)),
        'CHA': calc_arc_alt_freqs(p_gt.take([2348], axis=1)),
        'VIN': calc_arc_alt_freqs(p_gt.take([2349], axis=1)),
    }
    # Intialize a list to store the results.
    diff_list = []
    # For every haplotype group...
    for hap in hap_mat_dicc.keys():
        # For every archaic...
        for arc in arc_dip_dicc.keys():
            # Compute the number of pairwise differences.
            diff_list.append(hap_mat_diffs(hap_mat_dicc[hap], arc_dip_dicc[arc]))
    # Intialize an ordered list of pairwise comparisons.
    combos = [
        ('ARC', 'HUM'),
        ('ARC', 'REC'),
        ('HUM', 'REC'),
        ('ARC', 'ARC'),
        ('HUM', 'HUM'),
        ('REC', 'REC'),
    ]
    # For every combination...
    for combo in combos:
        # Extract the haplotype groups.
        key_1, key_2 = combo
        # Compute the number of pairwise differences.
        diff_list.append(hap_mat_diffs(hap_mat_dicc[key_1], hap_mat_dicc[key_2]))
    return np.array(diff_list)

# Define a function to calculate haplo-group divergence in windows.
def haplo_group_pw_diffs_windows(chromosome, window_size):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary of haplo-groups.
    hap_dicc = {
        'ARC': {
            'hap_1': np.loadtxt(f'../meta_data/{window_size}kb_arc_like_hap_1_idx.csv', delimiter=',', dtype=int),
            'hap_2': np.loadtxt(f'../meta_data/{window_size}kb_arc_like_hap_2_idx.csv', delimiter=',', dtype=int),
        },
        'HUM': {
            'hap_1': np.loadtxt(f'../meta_data/{window_size}kb_hum_like_hap_1_idx.csv', delimiter=',', dtype=int),
            'hap_2': np.loadtxt(f'../meta_data/{window_size}kb_hum_like_hap_2_idx.csv', delimiter=',', dtype=int),
        },
        'REC': {
            'hap_1': np.loadtxt('../meta_data/72kb_rec_like_hap_1_idx.csv', delimiter=',', dtype=int),
            'hap_2': [np.loadtxt('../meta_data/72kb_rec_like_hap_2_idx.csv', delimiter=',', dtype=int)],
        },
    }
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
    # Intialize a results matrix.
    dist_mat = np.empty((n_windows, 18))
    # For every window...
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # Polarize the genotype matrix.
        p_gt = polarize_gt(gt=allel.GenotypeArray(callset[wind_loc]))
        # Compute haplotype distances.
        dist_mat[wind, :] = haplo_group_diffs(p_gt=p_gt, haplo_groups=hap_dicc)
    # Export the the results matrix.
    np.savetxt(
        f'./tgp/windows/haplo_group_pw_diffs_chr{chromosome}_{window_size}kb.csv.gz',
        dist_mat, fmt='%1.15f', delimiter=',', newline='\n',
    )
    return

# Calculate haplotype divergence in non-overlapping windows.
haplo_group_pw_diffs_windows(
    chromosome=int(sys.argv[1]),
    window_size=int(sys.argv[2]),
)