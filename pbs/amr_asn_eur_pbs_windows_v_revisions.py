# Import packages.
import allel
from itertools import product
import numpy as np
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

# Define a function to calculate PBS_{AMR:SAS/EAS:EUR} in non-overlapping windows.
def pbs_windows(chrom, window_size):
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the asian, european, and admixed american populations.
    asn_pops = np.unique(meta_df[
        (meta_df['SUPERPOP'] == 'SAS') | (meta_df['SUPERPOP'] == 'EAS')
    ]['POP'].values)
    eur_pops = np.unique(meta_df[meta_df['SUPERPOP'] == 'EUR']['POP'].values)
    amr_pops = np.unique(meta_df[meta_df['SUPERPOP'] == 'AMR']['POP'].values)
    # Find all unique combinations.
    amr_asn_eur_combos = list(product(amr_pops, asn_pops, eur_pops))
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {}
    # For every population...
    for pop in np.concatenate((amr_pops, asn_pops, eur_pops)):
        # Append the dictionary with sample indicies.
        samp_idx_dicc[pop] = meta_df[meta_df['POP'] == pop].index.values
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('tgp_mod_no_aa', chrom)
    # Load the windows data frame.
    qc_windows_df = pd.read_csv(f'../windowing/tgp_mod_no_aa/{window_size}kb_nonoverlapping_variant_windows.csv.gz')
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chrom]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of winds.
    n_winds = chr_qc_windows_df.shape[0]
    # Intialize a dictionary to store the results.
    results_dicc = {combo: np.empty((n_winds)) for combo in amr_asn_eur_combos}
    # For ever window.
    for wind in range(n_winds):
        # Locate the window.
        wind_loc = all_pos.locate_range(chr_starts[wind], chr_stops[wind])
        # For every combination of population b and c.
        for a_pop, b_pop, c_pop in amr_asn_eur_combos:
            # Compute PBS.
            results_dicc[(a_pop, b_pop, c_pop)][wind] = calc_pbs_per_region(
                gt=allel.GenotypeArray(callset[wind_loc]), pop_a=samp_idx_dicc[a_pop],
                pop_b=samp_idx_dicc[b_pop], pop_c=samp_idx_dicc[c_pop],
            )
    # For every combination of population b and c.
    for a_pop, b_pop, c_pop in amr_asn_eur_combos:
        # Export the results.
        np.savetxt(
            f'../muc19_results/tgp_mod_no_aa/{a_pop.lower()}_{b_pop.lower()}_{c_pop.lower()}_pbs_chr{chrom}_{window_size}kb.txt.gz',
            [results_dicc[(a_pop, b_pop, c_pop)]], fmt='%1.15f',
        )
    return


# Calculate pbs.
pbs_windows(chrom=int(sys.argv[1]), window_size=int(sys.argv[2]))