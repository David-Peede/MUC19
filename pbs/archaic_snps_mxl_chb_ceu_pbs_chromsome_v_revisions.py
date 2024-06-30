# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###


# Define a function to load genotyope and positions arrays.
def load_gt_pos(prefix, chrom):
    # Intialize the file path.
    path = f'../zarr_data/{prefix}_chr{chrom}.zarr'
    # Load the zarr array.
    zarr_array = zarr.open_group(path, mode='r')
    # Extract the genotype callset.
    callset = zarr_array[f'{chrom}/calldata/GT']
    # Load the positions.
    pos = allel.SortedIndex(zarr_array[f'{chrom}/variants/POS'])
    # Convert the genotype callset to an array.
    gt = allel.GenotypeArray(callset)
    # Load the positions.
    pos = allel.SortedIndex(zarr_array['{0}/variants/POS'.format(chrom)])
    return gt, pos

# Define a function to compute Hudson's estimator of Fst per site
def calc_fst_per_site(ac1, ac2):
    # Compute.
    num, den = allel.hudson_fst(ac1, ac2)
    raw_fst = num / den
    # Correct for when the denominator is zero.
    raw_fst = np.where(den == 0, 0, raw_fst)
    # Correct for negative Fst values.
    raw_fst = np.where(raw_fst < 0, 0, raw_fst)
    # Correct for Fst values of 1 that will lead to inf.
    raw_fst = np.where(raw_fst == 1, 0.99999, raw_fst)
    return raw_fst

# Define a function to calculate PBS.
def calc_pbs_per_site(gt, pop_a, pop_b, pop_c):
    # Determine the sample list.
    samp_list = np.concatenate([pop_a, pop_b, pop_c])
    # Determine if any of the sites are segregating.
    seg_mask = gt.take(samp_list, axis=1).count_alleles().is_segregating()
    # Determine allele counts.
    a_ac = gt.take(pop_a, axis=1).count_alleles()
    b_ac = gt.take(pop_b, axis=1).count_alleles()
    c_ac = gt.take(pop_c, axis=1).count_alleles()
    # Calculate Fst with corrections.
    a_b_fst = calc_fst_per_site(a_ac, b_ac)
    a_c_fst = calc_fst_per_site(a_ac, c_ac)
    c_b_fst = calc_fst_per_site(c_ac, b_ac)
    # Calculate PBS.
    pbs = (
        ((np.log(1.0 - a_b_fst) * -1.0) +\
         (np.log(1.0 - a_c_fst) * -1.0) -\
         (np.log(1.0 - c_b_fst) * -1.0)) / 2.0
    )
    # Set the sites that are not segregating to undefined.
    pbs[~seg_mask] = np.nan
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

# Define a function to calculate pbs per archaic snp partition per chromosome.
def archaic_snps_pbs_chromosome(chrom):
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
    gt, pos = load_gt_pos('tgp_arcs_masked_no_aa', chrom)
    # Intialize a dictionary.
    pbs_dicc = {}
    # For every archaic snp type.
    for snp_type in snp_type_dicc:
        # Intialize a matrix to store the results.
        pbs_dicc[snp_type] = np.empty((pos.size, 3))
    # Generate the archaic specific masks.
    arc_masks = find_arc_alleles(gt=gt, pop_dicc=samp_idx_dicc, ooa_pop='MXL', ooa_freq=0.01)
    # For every archaic snp type.
    for snp_type in snp_type_dicc:
        # Intialize the results matrix.
        results_mat = np.empty((pos[arc_masks[snp_type]].size, 3))
        # For every ancestry partition.
        for i, mxl in enumerate(['MXL', 'MXL_NAT', 'MXL_NOT']):
            # Compute PBS.
            results_mat[:, i] = calc_pbs_per_site(
                gt=gt.compress(arc_masks[snp_type], axis=0),
                pop_a=samp_idx_dicc[mxl], pop_b=samp_idx_dicc['CHB'], pop_c=samp_idx_dicc['CEU'],
            )
        # Export the results.
        np.savetxt(
            f'../muc19_results/tgp_arcs_masked_no_aa/{snp_type_dicc[snp_type]}_mxl_chb_ceu_pbs_chr{chrom}.txt.gz',
            results_mat, fmt='%1.15f',
        )
    return


# Calculate pbs.
archaic_snps_pbs_chromosome(chrom=int(sys.argv[1]))