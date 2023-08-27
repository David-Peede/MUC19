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

# Define a function to calculate alternative allele frequencies for the archaics
def calc_arc_alt_freqs(arc_gt):
    # Compute the frequency for each site.
    raw_freqs = np.nansum(arc_gt, axis=2).flatten() / 2
    # Set missing data to Nan
    alt_freqs = np.where(raw_freqs == -1, np.nan, raw_freqs)
    return alt_freqs

# Define a function to count the number of archaic specific derived snps.
def count_arc_der_sites(gt, pop_dicc):
    # Calculate alternative allele frequencies.
    afr_alt_freq = calc_alt_freqs(gt.take(pop_dicc['AFR'], axis=1))
    ooa_alt_freq = calc_alt_freqs(gt.take(pop_dicc['OOA'], axis=1))
    alt_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['ALT'], axis=1))
    cha_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['CHA'], axis=1))
    vin_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['VIN'], axis=1))
    den_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['DEN'], axis=1))
    anc_freq = calc_alt_freqs(gt.take([-1], axis=1))
    # Polarize the samples.
    afr_der_freq = np.where(anc_freq == 1, np.abs(afr_alt_freq - 1), afr_alt_freq)
    ooa_der_freq = np.where(anc_freq == 1, np.abs(ooa_alt_freq - 1), ooa_alt_freq)
    alt_der_freq = np.where(anc_freq == 1, np.abs(alt_alt_freq - 1), alt_alt_freq)
    cha_der_freq = np.where(anc_freq == 1, np.abs(cha_alt_freq - 1), cha_alt_freq)
    vin_der_freq = np.where(anc_freq == 1, np.abs(vin_alt_freq - 1), vin_alt_freq)
    den_der_freq = np.where(anc_freq == 1, np.abs(den_alt_freq - 1), den_alt_freq)
    # Intialize conditions.
    c_hum = (afr_der_freq < 0.01) & (ooa_der_freq > 0)
    c_den = (den_der_freq == 1)
    c_not_den = (den_der_freq != 1)
    c_alt = (alt_der_freq == 1)
    c_not_alt = (alt_der_freq != 1)
    c_cha = (cha_der_freq == 1)
    c_not_cha = (cha_der_freq != 1)
    c_vin = (vin_der_freq == 1)
    c_not_vin = (vin_der_freq != 1)
    c_nea = c_alt | c_cha | c_vin
    c_not_nea = c_not_alt & c_not_cha & c_not_vin
    # Determine the archaic specific sites.
    den_only_idx = np.where(c_hum & c_not_nea & c_den)[0]
    nea_only_idx = np.where(c_hum & c_not_den & c_nea)[0]
    return [den_only_idx.size, nea_only_idx.size]

# Define a function to count the number of archaic specific ancestral snps.
def count_arc_anc_sites(gt, pop_dicc):
    # Calculate alternative allele frequencies.
    afr_alt_freq = calc_alt_freqs(gt.take(pop_dicc['AFR'], axis=1))
    ooa_alt_freq = calc_alt_freqs(gt.take(pop_dicc['OOA'], axis=1))
    alt_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['ALT'], axis=1))
    cha_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['CHA'], axis=1))
    vin_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['VIN'], axis=1))
    den_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['DEN'], axis=1))
    anc_freq = calc_alt_freqs(gt.take([-1], axis=1))
    # Polarize the samples.
    afr_der_freq = np.where(anc_freq == 1, np.abs(afr_alt_freq - 1), afr_alt_freq)
    ooa_der_freq = np.where(anc_freq == 1, np.abs(ooa_alt_freq - 1), ooa_alt_freq)
    alt_der_freq = np.where(anc_freq == 1, np.abs(alt_alt_freq - 1), alt_alt_freq)
    cha_der_freq = np.where(anc_freq == 1, np.abs(cha_alt_freq - 1), cha_alt_freq)
    vin_der_freq = np.where(anc_freq == 1, np.abs(vin_alt_freq - 1), vin_alt_freq)
    den_der_freq = np.where(anc_freq == 1, np.abs(den_alt_freq - 1), den_alt_freq)
    # Intialize conditions.
    c_hum = ((1 - afr_der_freq) < 0.01) & ((1 - ooa_der_freq) > 0)
    c_den = (den_der_freq == 0)
    c_not_den = (den_der_freq != 0)
    c_alt = (alt_der_freq == 0)
    c_not_alt = (alt_der_freq != 0)
    c_cha = (cha_der_freq == 0)
    c_not_cha = (cha_der_freq != 0)
    c_vin = (vin_der_freq == 0)
    c_not_vin = (vin_der_freq != 0)
    c_nea = c_alt | c_cha | c_vin
    c_not_nea = c_not_alt & c_not_cha & c_not_vin
    # Determine the archaic specific sites.
    den_only_idx = np.where(c_hum & c_not_nea & c_den)[0]
    nea_only_idx = np.where(c_hum & c_not_den & c_nea)[0]
    return [den_only_idx.size, nea_only_idx.size]

# Define a function to compute archaic snp denisty.
def tgp_arc_snp_denisty(chromosome, window_size):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {
        'ALT': np.array([2347]), 'CHA': np.array([2348]),
        'VIN': np.array([2349]), 'DEN': np.array([2350]),
    }
    # Append the dictionary with sample AFR and OOA indicies.
    samp_idx_dicc['AFR'] = tgp_meta_df[tgp_meta_df['SUPERPOP'] == 'AFR'].index.values
    samp_idx_dicc['OOA'] = tgp_meta_df[tgp_meta_df['SUPERPOP'] != 'AFR'].index.values
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
    results_mat = np.empty((n_windows, 4))
    # For every window...
    for wind in range(n_windows):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # Compute the number of derived and ancestral archaic snps.
        c_der = count_arc_der_sites(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_dicc=samp_idx_dicc,
        )
        c_anc = count_arc_anc_sites(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_dicc=samp_idx_dicc,
        )
        # Append the results matrix
        results_mat[wind, :] = np.array(c_der+c_anc)
    # Compile the results file name.
    results_file = f'chr{chromosome}_{window_size}kb.csv.gz'
    # Export the the results matrix.
    np.savetxt(
        './tgp/windows/'+results_file,
        results_mat, fmt='%d', delimiter=',', newline='\n',
    )
    return

# Compute the archaic allele density.
tgp_arc_snp_denisty(
    chromosome=int(sys.argv[1]),
    window_size=int(sys.argv[2]),
)
