# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###
### sys.argv[2] = papuan prefix ###


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

# Define a function to calculate alternative allele frequencies for a single individual.
def calc_ind_alt_freqs(gt):
    # Compute the frequency for each site.
    raw_freqs = np.nansum(gt, axis=2).flatten() / 2
    # Set missing data to Nan
    alt_freqs = np.where(raw_freqs == -1, np.nan, raw_freqs)
    return alt_freqs

# Define a function to compute pairwise differences.
def pwd_per_site(px, py):
    return ((px * (1 - py)) + (py * (1 - px)))

# Define a function to calculate haplotype divergence for papuan tracts.
def pap_v_den_seq_div_intro_tracts(chromosome, pap):
    # Intialize paths and dataset variables.
    meta_path = '../meta_data/sgdp.txt'
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        meta_path, sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Grab the papuan index.
    pap_idx = np.where(meta_df.IND.values == pap)[0][0]
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('sgdp_den_masked_no_aa', chromosome)
    # For each haplotype.
    for i, hap in enumerate(['hap1', 'hap2']):
        # Intialize the path to the qced tracts.
        tracts_path = f'../windowing/sgdp_den_masked_no_aa/{pap}_{hap}_den_intro_variant_tracts.csv.gz'
        # Load the tracts data frame.
        qc_tracts_df = pd.read_csv(tracts_path)
        # Subset the the tracts for the chromosome.
        chr_qc_tracts_df = qc_tracts_df[qc_tracts_df['CHR'] == chromosome]
        # Extract the start and stop positions.
        chr_starts = chr_qc_tracts_df['START'].values
        chr_stops = chr_qc_tracts_df['STOP'].values
        # Extract the effective sequence lengths.
        chr_esl = chr_qc_tracts_df['DEN'].values
        # Determine the number of tracts.
        n_tracts = chr_qc_tracts_df.shape[0]
        # Intialize a matrix to store the results.
        div_mat = np.empty((n_tracts, 2))
        # For every tract...
        for tract in range(n_tracts):
            # Locate the window.
            tract_loc = all_pos.locate_range(chr_starts[tract], chr_stops[tract])
            # Extract the genotype matrix.
            tract_gt = allel.GenotypeArray(callset[tract_loc])
            # Intialize the denisovan allele frequencies.
            den_freqs = calc_ind_alt_freqs(tract_gt.take([278], axis=1))
            # Extract the Papuan haplotype.
            pap_hap = tract_gt[:, pap_idx, i]
            # Compute the number of pairwise differences.
            pw_diffs = np.nansum(pwd_per_site(pap_hap, den_freqs))
            # Compute the sequence divergence.
            seq_div = pw_diffs / chr_esl[tract]
            # Update the results.
            div_mat[tract, :] = np.array([pw_diffs, seq_div])
        # Export the results.
        np.savetxt(
            f'../muc19_results/sgdp_den_masked_no_aa/{pap}_{hap}_den_pw_diffs_seq_div_at_den_intro_tracts_chr{chromosome}.txt.gz',
            div_mat, fmt='%1.15f',
        )
    return

# Calculate divergence for papuan tracts.
pap_v_den_seq_div_intro_tracts(
    chromosome=int(sys.argv[1]),
    pap=str(sys.argv[2]),
)