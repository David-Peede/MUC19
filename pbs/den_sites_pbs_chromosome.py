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
    path = '../vcf_data/zarr_data/{0}_chr{1}.zarr'.format(prefix, chrom)
    # Load the zarr array.
    zarr_array = zarr.open_group(path, mode='r')
    # Extract the genotype callset.
    callset = zarr_array['{0}/calldata/GT'.format(chrom)]
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

# Define a function to calculate alternative allele frequencies for the archaics
def calc_arc_alt_freqs(arc_gt):
    # Compute the frequency for each site.
    raw_freqs = np.nansum(arc_gt, axis=2).flatten() / 2
    # Set missing data to Nan
    alt_freqs = np.where(raw_freqs == -1, np.nan, raw_freqs)
    return alt_freqs

# Define a function to exctract a archaic specific snps.
def extract_arc_der_sites(gt, pos, pop_dicc):
    # Calculate alternative allele frequencies.
    afr_alt_freq = calc_alt_freqs(gt.take(pop_dicc['AFR'], axis=1))
    ooa_alt_freq = calc_alt_freqs(gt.take(pop_dicc['OOA'], axis=1))
    alt_alt_freq = calc_alt_freqs(gt.take(pop_dicc['ALT'], axis=1))
    cha_alt_freq = calc_alt_freqs(gt.take(pop_dicc['CHA'], axis=1))
    vin_alt_freq = calc_alt_freqs(gt.take(pop_dicc['VIN'], axis=1))
    den_alt_freq = calc_alt_freqs(gt.take(pop_dicc['DEN'], axis=1))
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
    # Construct the masks.
    den_mask = np.in1d(np.arange(pos.size), den_only_idx)
    nea_mask = np.in1d(np.arange(pos.size), nea_only_idx)
    # Compile a dictionary with the results.
    arc_sites_dicc = {
        'DEN': {'IDX': den_only_idx, 'MASK': den_mask},
        'NEA': {'IDX': nea_only_idx, 'MASK': nea_mask},
    }
    return arc_sites_dicc

# Define a function to exctract a archaic specific snps.
def extract_arc_anc_sites(gt, pos, pop_dicc):
    # Calculate alternative allele frequencies.
    afr_alt_freq = calc_alt_freqs(gt.take(pop_dicc['AFR'], axis=1))
    ooa_alt_freq = calc_alt_freqs(gt.take(pop_dicc['OOA'], axis=1))
    alt_alt_freq = calc_alt_freqs(gt.take(pop_dicc['ALT'], axis=1))
    cha_alt_freq = calc_alt_freqs(gt.take(pop_dicc['CHA'], axis=1))
    vin_alt_freq = calc_alt_freqs(gt.take(pop_dicc['VIN'], axis=1))
    den_alt_freq = calc_alt_freqs(gt.take(pop_dicc['DEN'], axis=1))
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
    # Construct the masks.
    den_mask = np.in1d(np.arange(pos.size), den_only_idx)
    nea_mask = np.in1d(np.arange(pos.size), nea_only_idx)
    # Compile a dictionary with the results.
    arc_sites_dicc = {
        'DEN': {'IDX': den_only_idx, 'MASK': den_mask},
        'NEA': {'IDX': nea_only_idx, 'MASK': nea_mask},
    }
    return arc_sites_dicc

# Define a function to calculate PBS.
def calc_pbs(gt, pop_a, pop_b, pop_c):
    # Determine allele counts.
    a_ac = gt.take(pop_a, axis=1).count_alleles()
    b_ac = gt.take(pop_b, axis=1).count_alleles()
    c_ac = gt.take(pop_c, axis=1).count_alleles()
    # Calculate the numerator and denominator for Hudson's Fst estimator.
    a_b_num, a_b_den = allel.hudson_fst(a_ac, b_ac)
    a_c_num, a_c_den = allel.hudson_fst(a_ac, c_ac)
    c_b_num, c_b_den = allel.hudson_fst(c_ac, b_ac)
    # Calculate Fst.
    a_b_raw_fst = a_b_num / a_b_den
    a_c_raw_fst = a_c_num / a_c_den
    c_b_raw_fst = c_b_num / c_b_den
    # Correct for Fst values of 1 that will lead to inf.
    a_b_fst_raw = np.where(a_b_raw_fst == 1, 0.99999, a_b_raw_fst)
    a_c_fst_raw = np.where(a_c_raw_fst == 1, 0.99999, a_c_raw_fst)
    c_b_fst_raw = np.where(c_b_raw_fst == 1, 0.99999, c_b_raw_fst)
    # Correct for negative Fst values.
    a_b_fst = np.where(a_b_fst_raw < 0, 0, a_b_fst_raw)
    a_c_fst = np.where(a_c_fst_raw < 0, 0, a_c_fst_raw)
    c_b_fst = np.where(c_b_fst_raw < 0, 0, c_b_fst_raw)
    # Calculate PBS.
    pbs = (
        ((np.log(1.0 - a_b_fst) * -1.0) +\
         (np.log(1.0 - a_c_fst) * -1.0) -\
         (np.log(1.0 - c_b_fst) * -1.0)) / float(2)
    )
    return pbs, np.nanmean(pbs)

# Define a function to calculate pbs per snp per chromosome.
def pbs_chromosome(chrom):
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
    # For every population...
    for pop in ['MXL', 'CEU', 'CHB']:
        # Append the dictionary with sample indicies.
        samp_idx_dicc[pop] = tgp_meta_df[tgp_meta_df['POP'] == pop].index.values
    # Extract the genotype callset and positions.
    gt, pos = load_gt_pos('tgp_mod_arc_anc', chrom)
    # Determine which sites are denisovan specific.
    arc_der_dicc = extract_arc_der_sites(gt=gt, pos=pos, pop_dicc=samp_idx_dicc)
    arc_anc_dicc = extract_arc_anc_sites(gt=gt, pos=pos, pop_dicc=samp_idx_dicc)
    # Create the Denisovan specific mask.
    den_mask = arc_der_dicc['DEN']['MASK'] | arc_anc_dicc['DEN']['MASK']
    # Calculate pbs.
    mxl_pbs_all, _  = calc_pbs(
        gt=gt.compress(den_mask, axis=0), pop_a=samp_idx_dicc['MXL'],
        pop_b=samp_idx_dicc['CHB'], pop_c=samp_idx_dicc['CEU'],
    )
    # Export the the results matrix.
    np.savetxt(
        './chroms/mxl_chb_ceu_den_sites_pbs_all_chr{0}.csv.gz'.format(chrom),
        [mxl_pbs_all], fmt='%1.15f', delimiter=',', newline='\n',
    )
    return


# Calculate pbs.
pbs_chromosome(chrom=int(sys.argv[1]))