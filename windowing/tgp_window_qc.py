# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = window size ###
### sys.argv[2] = chromosome ###


# Define a function to compute adjusted chromosome lengths.
def chr_seq_len(window_size):
    # Intialize the original assembly size for hg19.
    chr_dicc = {
        1: 249250621, 2: 243199373, 3: 198022430,
        4: 191154276, 5: 180915260, 6: 171115067,
        7: 159138663, 8: 146364022, 9: 141213431,
        10: 135534747, 11: 135006516, 12: 133851895,
        13: 115169878, 14: 107349540, 15: 102531392,
        16: 90354753, 17: 81195210, 18: 78077248,
        19: 59128983, 20: 63025520, 21: 48129895,
        22: 51304566,
    }
    # Initialize new empty dictionary to store the output.
    new_chr_dicc = {}
    # Iterate through every chromosome.
    for key in chr_dicc :
        # Floor divide every chromosome length by the window size and multiply
        # to find the new chromosome length.
        chr_len = chr_dicc[key]
        new_chr_len = (chr_len//window_size)*window_size
        # Refill dictionary with the new chromosome lengths.
        new_chr_dicc[key] = new_chr_len
    return new_chr_dicc

# Define a function to break up a chromosome into windows.
def window_info(positions, window_size, sequence_length):
    # Intialize a dicctionary with the start and stop position for each window.
    windows = {}
    start_stop = {}
    index = 0
    for window_start in range(1, int(sequence_length), int(window_size)):
        windows[index] = np.where(((window_start <= positions) & (positions < (window_start+window_size))))[0]
        start_stop[index] = [window_start, (window_start+window_size)]
        index += 1
    return windows, start_stop

# Define a function to calculate the number of variable sites.
def variable_site_check(gt, idx_dicc):
    # Intialize a list to store the number of varibale sites.
    var_sites = []
    # For every archaic...
    for key in idx_dicc.keys():
        # Determine the indicies where all samples are called.
        called_mask = (gt.take(idx_dicc[key], axis=1).is_called() == True).all(axis=1)
        # If there are no sites called between all samples...
        if (called_mask.sum() == 0):
            # Append the results.
            var_sites.append(0)
        # Else:
        else:
            # Determine the number of variable sites.
            var_bool = gt.take(idx_dicc[key], axis=1).compress(called_mask, axis=0).count_alleles().is_variant()
            # Append the results.
            var_sites.append(var_bool.sum())
    return np.array(var_sites)

# Define a function to calculate the number of segregating sites.
def count_seg_sites(gt):
    # Count the number of segregsting sites.
    seg_sites = gt.take(np.arange(gt.shape[1]-1), axis=1).count_alleles().count_segregating()
    return seg_sites

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

# Define a function to qc non-overlapping windows.
def windows_qc(window_size, chromosome):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Identify all sample indicies.
    tgp_idx = tgp_meta_df.index.values
    # Intialize sample indicies dictionaries.
    var_check_dicc = {
        'ALT': np.append(tgp_idx, np.array([2347])),
        'CHA': np.append(tgp_idx, np.array([2348])),
        'VIN': np.append(tgp_idx, np.array([2349])),
        'DEN': np.append(tgp_idx, np.array([2350])),
    }
    # Initialize the window length.
    window_length = window_size * 1_000
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('tgp_mod_arc_anc', chromosome)
    # Calculate the adjusted chromosome lengths.
    chr_len_dicc = chr_seq_len(window_length)
    # Construct the window dictionaries.
    wind_dicc, start_stop = window_info(
        positions=all_pos,
        window_size=window_length,
        sequence_length=chr_len_dicc[chromosome],
    )
    # Load the all site information.
    all_sites_mat = np.loadtxt(
        '../vcf_data/vcf_bookkeeping/tgp_all_archaics_merged_all_sites_qc/chr{0}_all_sites_dedupped_report.csv'.format(chromosome),
        delimiter=',', dtype=int,
    )
    # Determine the number of windows.
    n_windows = len(wind_dicc.keys())
    # Intialize a results matrix to store the results.
    results_mat = np.empty((n_windows, 17))
    # For every window...
    for wind in wind_dicc.keys():
        # Extract the start and stop positions.
        start = start_stop[wind][0]
        stop = start_stop[wind][1] - 1 # Sci-kit allele uses [inclusive, inclsuive] indexing.
        # Extract the variable sites information.
        var_idx = wind_dicc[wind]
        # Extract the effective sequence length information.
        esl_idx = np.where(((start <= all_sites_mat[:, 1]) & (all_sites_mat[:, 1] <= stop)))[0]
        # Determine the effective sequence length for the tgp.
        tgp_esl = esl_idx.size
        # If the window has variant sites to do calculations on...
        if (var_idx.size > 0):
            # Locate the window.
            wind_loc = all_pos.locate_range(start, stop)
            # Determine the number of variable sites.
            wind_var_sites = variable_site_check(
                gt=allel.GenotypeArray(callset[wind_loc]),
                idx_dicc=var_check_dicc,
            )
            # Determine the number of segregating sites.
            s = count_seg_sites(gt=allel.GenotypeArray(callset[wind_loc]))
            # Determine the effective sequence length for each sample.
            eff_seq_len = all_sites_mat[esl_idx, 2:].sum(axis=0)
            # Determine the effective sequence length per pairwise comparison.
            den_alt = np.count_nonzero((all_sites_mat[esl_idx, 2] == 1) & (all_sites_mat[esl_idx, 5] == 1))
            den_cha = np.count_nonzero((all_sites_mat[esl_idx, 3] == 1) & (all_sites_mat[esl_idx, 5] == 1))
            den_vin = np.count_nonzero((all_sites_mat[esl_idx, 4] == 1) & (all_sites_mat[esl_idx, 5] == 1))
            alt_cha = np.count_nonzero((all_sites_mat[esl_idx, 3] == 1) & (all_sites_mat[esl_idx, 2] == 1))
            alt_vin = np.count_nonzero((all_sites_mat[esl_idx, 4] == 1) & (all_sites_mat[esl_idx, 2] == 1))
            cha_vin = np.count_nonzero((all_sites_mat[esl_idx, 4] == 1) & (all_sites_mat[esl_idx, 3] == 1))
            pw_list = [
                den_alt, den_cha, den_vin,
                alt_cha, alt_vin, cha_vin,
            ]
            # If all archaics have variant site information.
            if (wind_var_sites > 0).all():
                # Append the results.
                results_mat[wind, :] = np.array(
                    [wind, start, stop]+eff_seq_len.tolist()+pw_list+[tgp_esl, s]+[1], dtype=object,
                )
            # Else-if all archaics passed QC at at least one site.
            elif (np.array(eff_seq_len.tolist()+pw_list) > 0).all():
                # Append the results.
                results_mat[wind, :] = np.array(
                    [wind, start, stop]+eff_seq_len.tolist()+pw_list+[tgp_esl, s]+[1], dtype=object,
                )
            # Else.
            else:
                # Append the results.
                results_mat[wind, :] = np.array(
                    [wind, start, stop]+eff_seq_len.tolist()+pw_list+[tgp_esl, s]+[0], dtype=object,
                )
        # Else-if there sites that passed qc for this window.
        elif (esl_idx.size > 0):
            # Determine the effective sequence length for each sample.
            eff_seq_len = all_sites_mat[esl_idx, 2:].sum(axis=0)
            # Determine the effective sequence length per pairwise comparison.
            den_alt = np.count_nonzero((all_sites_mat[esl_idx, 2] == 1) & (all_sites_mat[esl_idx, 5] == 1))
            den_cha = np.count_nonzero((all_sites_mat[esl_idx, 3] == 1) & (all_sites_mat[esl_idx, 5] == 1))
            den_vin = np.count_nonzero((all_sites_mat[esl_idx, 4] == 1) & (all_sites_mat[esl_idx, 5] == 1))
            alt_cha = np.count_nonzero((all_sites_mat[esl_idx, 3] == 1) & (all_sites_mat[esl_idx, 2] == 1))
            alt_vin = np.count_nonzero((all_sites_mat[esl_idx, 4] == 1) & (all_sites_mat[esl_idx, 2] == 1))
            cha_vin = np.count_nonzero((all_sites_mat[esl_idx, 4] == 1) & (all_sites_mat[esl_idx, 3] == 1))
            pw_list = [
                den_alt, den_cha, den_vin,
                alt_cha, alt_vin, cha_vin,
            ]
            # If all archaics passed QC at at least one site.
            if (np.array(eff_seq_len.tolist()+pw_list) > 0).all():
                # Append the results.
                results_mat[wind, :] = np.array(
                    [wind, start, stop]+eff_seq_len.tolist()+pw_list+[tgp_esl, 0]+[1], dtype=object,
                )
            # Else.
            else:
                # Append the results.
                results_mat[wind, :] = np.array(
                    [wind, start, stop]+eff_seq_len.tolist()+pw_list+[tgp_esl, 0]+[0], dtype=object,
                )
        # Else.
        else:
            # Append the results.
            results_mat[wind, :] = np.array(
                [wind, start, stop]+np.zeros(12).tolist()+[0], dtype=object,
            )
    # Export the the results matrix as a pandas dataframe.
    wind_df = pd.DataFrame(
        results_mat,
        columns=[
            'IDX',
            'START', 'STOP',
            'ALT', 'CHA',
            'VIN', 'DEN',
            'DEN-ALT', 'DEN-CHA',
            'DEN-VIN', 'ALT-CHA',
            'ALT-VIN', 'CHA-VIN',
            'TGP', 'S', 'QC',
        ],
    )
    wind_df = wind_df.astype({
        'IDX': 'int',
        'START': 'int', 'STOP': 'int',
        'ALT': 'int', 'CHA': 'int',
        'VIN': 'int', 'DEN': 'int',
        'DEN-ALT': 'int', 'DEN-CHA': 'int',
        'DEN-VIN': 'int', 'ALT-CHA': 'int',
        'ALT-VIN': 'int', 'CHA-VIN': 'int',
        'TGP': 'int', 'S': 'int',
        'QC': 'int',
    })
    wind_df.to_csv('./tgp/{0}kb_window_summary_chr{1}.csv'.format(window_size, chromosome), index=False)
    return

# Conduct the QC.
windows_qc(window_size=int(sys.argv[1]), chromosome=int(sys.argv[2]))
