# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###


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

# Define a funct to compute U_ABC(w, x, y).
def U_ABC(gt, A, B, C, w, x, y):
    # Calculate alternative allele frequencies.
    A_alt_freq = calc_alt_freqs(gt.take(A, axis=1))
    B_alt_freq = calc_alt_freqs(gt.take(B, axis=1))
    C_alt_freq = calc_ind_alt_freqs(gt.take(C, axis=1))
    # Intialize conditions.
    A_ref = (1 - A_alt_freq) < w
    A_alt = A_alt_freq < w
    B_ref = (1 - B_alt_freq) > x
    B_alt = B_alt_freq > x
    C_ref = (1 - C_alt_freq) == y
    C_alt = C_alt_freq == y
    # Determine the U sites.
    U_ref = A_ref & B_ref & C_ref
    U_alt = A_alt & B_alt & C_alt
    U_sites = U_ref | U_alt
    return U_sites.sum()


# Define a function to compute U_{AFR,B,DEN}(1%, 30%, 100%) statistics for every ncbi refseq gene.
def tgp_U30_genes(chromosome):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {'DEN': np.array([2347])}
    # Intialize a list of focal populations.
    ooa_list = [
        'MXL', 'PEL', 'CLM', 'PUR', # AMR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
    ]
    # Update the dictionary with the AFR indicies.
    samp_idx_dicc['AFR'] = tgp_meta_df[tgp_meta_df['SUPERPOP'] == 'AFR'].index.values
    # For every OOA population.
    for pop in ooa_list:
        # Update the dictionary.
        samp_idx_dicc[pop] = tgp_meta_df[tgp_meta_df['POP'] == pop].index.values
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos('tgp_den_masked_no_aa', chromosome)
    # Load the genes data frame.
    qc_genes_df = pd.read_csv(f'../annotations/tgp_den_masked_no_aa/ncbi_refseq_variant_genes.csv.gz')
    # Subset the the genes for the chromosome.
    chr_qc_genes_df = qc_genes_df[qc_genes_df['CHR'] == chromosome]
    # Extract the start and stop positions.
    chr_starts = chr_qc_genes_df['START'].values
    chr_stops = chr_qc_genes_df['STOP'].values
    # Determine the number of genes.
    n_genes = chr_qc_genes_df.shape[0]
    # Intialize a results matrix to store the results.
    U_AB_den_mat = np.empty((n_genes, 19))
    # For every gene.
    for gene in range(n_genes):
        # Locate the gene.
        gene_loc = all_pos.locate_range(chr_starts[gene], chr_stops[gene])
        # Intialize a list to store the reseults.
        U_AB_den_list = []
        # For every OOA population.
        for pop in ooa_list:
            # Compute the U_ABC stats.
            U_AB_den = U_ABC(
                gt=allel.GenotypeArray(callset[gene_loc]),
                A=samp_idx_dicc['AFR'], B=samp_idx_dicc[pop], C=samp_idx_dicc['DEN'],
                w=0.01, x=0.3, y=1,
            )
            # Update the list.
            U_AB_den_list.append(U_AB_den)
        # Append the results matrix
        U_AB_den_mat[gene, :] = np.array(U_AB_den_list)
    # Export the the results matrix.
    np.savetxt(
        f'../muc19_results/tgp_den_masked_no_aa/u30_afr_b_den_per_gene_chr{chromosome}.txt.gz',
        U_AB_den_mat, fmt='%d',
    )
    return

# Compute the U-statistics.
tgp_U30_genes(
    chromosome=int(sys.argv[1]),
)