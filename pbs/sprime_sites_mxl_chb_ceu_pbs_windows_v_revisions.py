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


# Define a function to calculate pbs for all sprime sites in a window.
def sprime_sites_pbs_windows(chrom, window_size):
    # Intialize a list of sprime site types.
    sprime_sites = ['all_arc']
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('tgp_mod_no_aa', chrom)
    # Intilize a dictionary.
    sprime_dicc = {}
    # For every sprime site.
    for site in sprime_sites:
        # Load the information.
        sprime_pos = np.loadtxt(f'../meta_data/mxl_sprime_{site}_sites_chr{chrom}.txt.gz', dtype=int)
        # Determine the sites in the dataset.
        dataset_mask = np.isin(sprime_pos, all_pos)
        # Fill the sites of interest.
        sprime_dicc[site] = sprime_pos[dataset_mask]
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {}
    # For every population...
    for pop in ['MXL', 'CHB', 'CEU']:
        # Append the dictionary with sample indicies.
        samp_idx_dicc[pop] = meta_df[meta_df['POP'] == pop].index.values
    # Load the windows data frame.
    qc_windows_df = pd.read_csv(f'../windowing/tgp_mod_no_aa/{window_size}kb_nonoverlapping_variant_windows.csv.gz')
    # Subset the the windows for the chromosome.
    chr_qc_windows_df = qc_windows_df[qc_windows_df['CHR'] == chrom]
    # Extract the start and stop positions.
    chr_starts = chr_qc_windows_df['START'].values
    chr_stops = chr_qc_windows_df['STOP'].values
    # Determine the number of winds.
    n_winds = chr_qc_windows_df.shape[0]
    # Intialize a dictionary.
    results_dicc = {arc_type: np.empty((n_winds, 2)) for arc_type in sprime_sites}
    # For ever gene.
    for wind in range(n_winds):
        # For every archaic snp type.
        for arc_type in sprime_sites:
            # Create a mask of the archaic sites.
            arc_mask = ((chr_starts[wind] <= sprime_dicc[arc_type]) & (sprime_dicc[arc_type] <= chr_stops[wind]))
            # If there are no archaic alleles.
            if (arc_mask.sum() == 0):
                # Update the results.
                results_dicc[arc_type][wind, :] = np.array([np.nan, 0])
            # Else there are archaic alleles to do computations on.
            else:
                # Generate a mask.
                arc_pos_mask = np.isin(all_pos, sprime_dicc[arc_type][arc_mask])
                # Compute PBS.
                wind_pbs = calc_pbs_per_region(
                    gt=allel.GenotypeArray(callset).compress(arc_pos_mask, axis=0),
                    pop_a=samp_idx_dicc['MXL'], pop_b=samp_idx_dicc['CHB'], pop_c=samp_idx_dicc['CEU'],
                )
                # Update the results.
                results_dicc[arc_type][wind, :] = np.array([wind_pbs, arc_pos_mask.sum()])
    # For every archaic site type.
    for arc_type in sprime_sites:
        # Export the results.
        np.savetxt(
            f'../muc19_results/tgp_mod_no_aa/mxl_chb_ceu_pbs_sprime_{arc_type.lower()}_chr{chrom}_{window_size}kb.txt.gz',
            results_dicc[arc_type], fmt='%1.15f',
        )
    return


# Calculate pbs.
sprime_sites_pbs_windows(chrom=int(sys.argv[1]), window_size=int(sys.argv[2]))