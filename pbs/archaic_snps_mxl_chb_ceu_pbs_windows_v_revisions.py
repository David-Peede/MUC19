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

# Define a  function to compute Hudson's Fst estimator for a region.
def calc_fst_per_region(ac1, ac2):
    # Compute.
    num, den = allel.hudson_fst(ac1, ac2)
    # Account for the denominator being zero.
    fst = np.nansum(num) / np.nansum(den) if np.nansum(den) != 0 else 0
    # Correct for negative Fst values.
    return max(fst, 0)

# Define a function to calculate PBS per region.
def calc_pbs_per_region(gt, pop_a, pop_b, pop_c):
    # Determine the sample list.
    samp_list = np.concatenate([pop_a, pop_b, pop_c])
    # Determine if any of the sites are segregating.
    seg_mask = gt.take(samp_list, axis=1).count_alleles().is_segregating()
    # If no sites are segregating.
    if seg_mask.sum() == 0:
        # Return undefinded.
        return np.nan
    # Else there are segreagting sites.
    else:
        # Determine allele counts.
        a_ac = gt.take(pop_a, axis=1).count_alleles()
        b_ac = gt.take(pop_b, axis=1).count_alleles()
        c_ac = gt.take(pop_c, axis=1).count_alleles()
        # Calculate Fst.
        a_b_fst = calc_fst_per_region(a_ac, b_ac)
        a_c_fst = calc_fst_per_region(a_ac, c_ac)
        c_b_fst = calc_fst_per_region(c_ac, b_ac)
        # Correct for Fst values of 1 that will lead to inf.
        a_b_fst = min(a_b_fst, 0.99999)
        a_c_fst = min(a_c_fst, 0.99999)
        c_b_fst = min(c_b_fst, 0.99999)
        # Calculate PBS.
        pbs = (
            ((np.log(1.0 - a_b_fst) * -1.0) +\
             (np.log(1.0 - a_c_fst) * -1.0) -\
             (np.log(1.0 - c_b_fst) * -1.0)) / 2.0
        )
    return pbs

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

# Define a function to find the number of archaic snps.
def find_arc_alleles(gt, pop_dicc, ooa_pop, ooa_freq):
    # Calculate alternative allele frequencies.
    afr_alt_freq = calc_alt_freqs(gt.take(pop_dicc['AFR'], axis=1))
    ooa_alt_freq = calc_alt_freqs(gt.take(pop_dicc[ooa_pop], axis=1))
    alt_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['ALT'], axis=1))
    cha_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['CHA'], axis=1))
    vin_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['VIN'], axis=1))
    den_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['DEN'], axis=1))
    # Intialize the conditions.
    c_ref_hum = ((1 - afr_alt_freq) < 0.01) & ((1 - ooa_alt_freq) > ooa_freq)
    c_ref_den = (den_alt_freq == 0)
    c_ref_not_den = (den_alt_freq != 0)
    c_ref_alt = (alt_alt_freq == 0)
    c_ref_not_alt = (alt_alt_freq != 0)
    c_ref_cha = (cha_alt_freq == 0)
    c_ref_not_cha = (cha_alt_freq != 0)
    c_ref_vin = (vin_alt_freq == 0)
    c_ref_not_vin = (vin_alt_freq != 0)
    c_ref_nea = c_ref_alt | c_ref_cha | c_ref_vin
    c_ref_not_nea = c_ref_not_alt & c_ref_not_cha & c_ref_not_vin
    c_ref_shared = c_ref_nea & c_ref_den
    c_alt_hum = (afr_alt_freq < 0.01) & (ooa_alt_freq > ooa_freq)
    c_alt_den = (den_alt_freq == 1)
    c_alt_not_den = (den_alt_freq != 1)
    c_alt_alt = (alt_alt_freq == 1)
    c_alt_not_alt = (alt_alt_freq != 1)
    c_alt_cha = (cha_alt_freq == 1)
    c_alt_not_cha = (cha_alt_freq != 1)
    c_alt_vin = (vin_alt_freq == 1)
    c_alt_not_vin = (vin_alt_freq != 1)
    c_alt_nea = c_alt_alt | c_alt_cha | c_alt_vin
    c_alt_not_nea = c_alt_not_alt & c_alt_not_cha & c_alt_not_vin
    c_alt_shared = c_alt_nea & c_alt_den
    # Determine the archaic specific masks.
    den_only_ref_mask = c_ref_hum & c_ref_not_nea & c_ref_den
    nea_only_ref_mask = c_ref_hum & c_ref_not_den & c_ref_nea
    shared_ref_mask = c_ref_hum & c_ref_shared
    arc_ref_mask = den_only_ref_mask | nea_only_ref_mask | shared_ref_mask
    den_only_alt_mask = c_alt_hum & c_alt_not_nea & c_alt_den
    nea_only_alt_mask = c_alt_hum & c_alt_not_den & c_alt_nea
    shared_alt_mask = c_alt_hum & c_alt_shared
    arc_alt_mask = den_only_alt_mask | nea_only_alt_mask | shared_alt_mask
    # Construct a dictionary.
    arc_dicc = {
        'DEN': (den_only_ref_mask | den_only_alt_mask),
        'NEA': (nea_only_ref_mask | nea_only_alt_mask),
        'SHR': (shared_ref_mask | shared_alt_mask),
        'ARC': (arc_ref_mask | arc_alt_mask),
    }
    return arc_dicc

# Define a function to calculate PBS_{MXL:CHB:CEU} in non-overlapping windows per archaic snp partition.
def archaic_snps_pbs_windows(chrom, window_size):
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize the snp-type dictionary.
    snp_type_dicc = {
        'DEN': 'denisovan_specific_snps',
        'NEA': 'neanderthal_specific_snps',
        'SHR': 'shared_archaic_snps',
        'ARC': 'archaic_specific_snps',
    }
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {
        'MXL_NAT': np.loadtxt('../amr_lai/anc_props/mxl_nat_idx.txt.gz', dtype=int),
        'MXL_NOT': np.loadtxt('../amr_lai/anc_props/mxl_not_idx.txt.gz', dtype=int),
        'AFR': meta_df[meta_df['SUPERPOP'] == 'AFR'].index.values,
        'OOA': meta_df[meta_df['SUPERPOP'] != 'AFR'].index.values,
        'ALT': np.array([2347]), 'CHA': np.array([2348]),
        'VIN': np.array([2349]), 'DEN': np.array([2350]),
    }
    # For every population...
    for pop in ['MXL', 'CEU', 'CHB']:
        # Append the dictionary with sample indicies.
        samp_idx_dicc[pop] = meta_df[meta_df['POP'] == pop].index.values
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('tgp_arcs_masked_no_aa', chrom)
    # Load the windows data frame.
    qc_windows_df = pd.read_csv(f'../windowing/tgp_arcs_masked_no_aa/{window_size}kb_nonoverlapping_variant_windows.csv.gz')
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chrom]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of winds.
    n_winds = chr_qc_windows_df.shape[0]
    # Intialize a dictionary.
    pbs_dicc = {}
    # For every archaic snp type.
    for snp_type in snp_type_dicc:
        # Intialize a matrix to store the results.
        pbs_dicc[snp_type] = np.empty((n_winds, 3))
    # For ever window.
    for wind in range(n_winds):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # Generate the archaic specific snp masks.
        arc_masks = find_arc_alleles(
            gt=allel.GenotypeArray(callset[wind_loc]),
            pop_dicc=samp_idx_dicc, ooa_pop='MXL', ooa_freq=0.01,
        )
        # For every archaic snp type.
        for snp_type in arc_masks:
            # If there are no archaic snps.
            if (arc_masks[snp_type].sum() == 0):
                # Update the results.
                pbs_dicc[snp_type][wind, :] = np.array([np.nan, np.nan, np.nan])
            # Else, there are archaic snps to perform computations on.
            else:
                # For every ancestry partition.
                for i, mxl in enumerate(['MXL', 'MXL_NAT', 'MXL_NOT']):
                    # Compute PBS.
                    pbs_dicc[snp_type][wind, i] = calc_pbs_per_region(
                        gt=allel.GenotypeArray(callset[wind_loc]).compress(arc_masks[snp_type], axis=0),
                        pop_a=samp_idx_dicc[mxl], pop_b=samp_idx_dicc['CHB'], pop_c=samp_idx_dicc['CEU'],
                    )
    # For every archaic snp type.
    for snp_type in pbs_dicc:
        # Export the results.
        np.savetxt(
            f'../muc19_results/tgp_arcs_masked_no_aa/{snp_type_dicc[snp_type]}_mxl_chb_ceu_pbs_chr{chrom}_{window_size}kb.txt.gz',
            pbs_dicc[snp_type], fmt='%1.15f',
        )
    return


# Calculate pbs.
archaic_snps_pbs_windows(chrom=int(sys.argv[1]), window_size=int(sys.argv[2]))