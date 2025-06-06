# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = window size ###
### sys.argv[2] = chromosome ###
### sys.argv[3] = dataset prefix ###


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

# Define a function to break up a chromosome into non-overlapping windows.
def window_info(positions, window_size, sequence_length):
    # Intialize a dicctionary with the start and stop position for each window.
    windows = {}
    start_stop = {}
    index = 0
    # For every non-overlapping window.
    for window_start in range(1, int(sequence_length), int(window_size)):
        # Update the position indicies of the window.
        windows[index] = np.where(((window_start <= positions) & (positions < (window_start+window_size))))[0]
        # Update the start and stop positions.
        start_stop[index] = [window_start, (window_start+window_size)]
        index += 1
    return windows, start_stop

# Define a function to calculate the number of segregating sites.
def count_seg_sites(gt):
    # Count the number of segregsting sites.
    seg_sites = gt.count_alleles().count_segregating()
    return seg_sites

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

# Define a function to qc non-overlapping windows.
def windows_qc(window_size, chromosome, prefix):
    # Intialize a dictionary to store the all sites path information.
    path_dicc = {
        'den_masked_no_aa': '../vcf_data/bookkeeping/den_masked_no_aa_calls_all_sites',
        'alt_masked_no_aa': '../vcf_data/bookkeeping/alt_masked_no_aa_calls_all_sites',
        'cha_masked_no_aa': '../vcf_data/bookkeeping/cha_masked_no_aa_calls_all_sites',
        'vin_masked_no_aa': '../vcf_data/bookkeeping/vin_masked_no_aa_calls_all_sites',
    }
    # Intitialize the all sites path, window path, and column variable.
    arc_all_sites_path = f'{path_dicc[prefix]}_chr{chromosome}.txt.gz'
    arc_window_path = f'./{prefix}/{window_size}kb_window_summary_chr{chromosome}.csv.gz'
    # Extract the archaic.
    archaic = prefix.split('_')[0]
    # Initialize the window length.
    window_length = window_size * 1_000
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos(prefix, chromosome)
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
        arc_all_sites_path,
        usecols=(0, 1),
        dtype=int,
    )
    # Determine the number of windows.
    n_windows = len(wind_dicc.keys())
    # Intialize a results matrix to store the results.
    results_mat = np.empty((n_windows, 6))
    # For every window...
    for wind in wind_dicc.keys():
        # Extract the start and stop positions.
        start = start_stop[wind][0]
        stop = start_stop[wind][1] - 1 # Sci-kit allele uses [inclusive, inclsuive] indexing.
        # Extract the variable sites information.
        var_idx = wind_dicc[wind]
        # Extract the effective sequence length information.
        esl_idx = np.where(((start <= all_sites_mat[:, 1]) & (all_sites_mat[:, 1] <= stop)))[0]
        # Determine the effective sequence length.
        arc_esl = esl_idx.size
        # If the window has variant sites to do calculations on...
        if (var_idx.size > 0):
            # Locate the window.
            wind_loc = all_pos.locate_range(start, stop)
            # Determine the number of segregating sites.
            s = count_seg_sites(gt=allel.GenotypeArray(callset[wind_loc]))
            # Append the results.
            results_mat[wind, :] = np.array(
                [wind, start, stop, arc_esl, s, 1], dtype=object,
            )
        # Else-if there are sites that passed qc for this window.
        elif (esl_idx.size > 0):
            # Append the results.
            results_mat[wind, :] = np.array(
                [wind, start, stop, arc_esl, 0, 1], dtype=object,
            )
        # Else.
        else:
            # Append the results.
            results_mat[wind, :] = np.array(
                [wind, start, stop, 0, 0, 0], dtype=object,
            )
    # Export the the results matrix as a pandas dataframe.
    wind_df = pd.DataFrame(
        results_mat,
        columns=[
            'IDX',
            'START', 'STOP',
            f'{archaic.upper()}', 'S', 'QC',
        ],
    )
    wind_df = wind_df.astype({
        'IDX': 'int',
        'START': 'int', 'STOP': 'int',
        f'{archaic.upper()}': 'int', 'S': 'int',
        'QC': 'int',
    })
    wind_df.to_csv(arc_window_path, index=False)
    return

# Conduct the QC.
windows_qc(window_size=int(sys.argv[1]), chromosome=int(sys.argv[2]), prefix=str(sys.argv[3]))