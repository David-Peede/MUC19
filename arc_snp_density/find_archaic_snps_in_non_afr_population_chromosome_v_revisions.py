# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###
### sys.argv[2] = non-afr population ###


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

# Define a function to identify the archaic specific snps for a non-AFR population.
def identify_archaic_snps_chromosome(chrom, pop):
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
        f'{pop}': meta_df[meta_df['POP'] == pop].index.values,
        'AFR': meta_df[meta_df['SUPERPOP'] == 'AFR'].index.values,
        'OOA': meta_df[meta_df['SUPERPOP'] != 'AFR'].index.values,
        'ALT': np.array([2347]), 'CHA': np.array([2348]),
        'VIN': np.array([2349]), 'DEN': np.array([2350]),
    }
    # Extract the genotype callset and positions.
    gt, pos = load_gt_pos('tgp_arcs_masked_no_aa', chrom)
    # Generate the archaic specific masks.
    arc_masks = find_arc_alleles(gt=gt, pop_dicc=samp_idx_dicc, ooa_pop=pop, ooa_freq=0.01)
    # For every archaic snp type.
    for snp_type in snp_type_dicc:
        # Export the results.
        np.savetxt(
            f'../muc19_results/tgp_arcs_masked_no_aa/{pop.lower()}_{snp_type_dicc[snp_type]}_chr{chrom}.txt.gz',
            [pos[arc_masks[snp_type]]], fmt='%d',
        )
    return


# Determine what poistions are archaic snps.
identify_archaic_snps_chromosome(chrom=int(sys.argv[1]), pop=str(sys.argv[2]))