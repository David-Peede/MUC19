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
    # Convert the genotype callset to an array.
    gt = allel.GenotypeArray(callset)
    # Load the positions.
    pos = allel.SortedIndex(zarr_array[f'{chrom}/variants/POS'])
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

# Define a function to compute per site heterozygosity.
def per_site_heterozygosity(gt):
    # Determine the number of sequences.
    n = gt.shape[1] * 2
    # Compute the alternative allele frequency.
    p = calc_alt_freqs(gt)
    # Compute the per site heterozygosity.
    per_site_h = (n / (n - 1)) * (1 - ((p ** 2) + ((1 - p) ** 2)))
    return per_site_h

# Define a function to compute the number of heterozygous sites and heterozygoisty for a chromosome.
def tgp_het_chromosome(chromosome):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
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
    # Load the qc site information.
    qc_sites_mat = np.loadtxt(
        f'../vcf_data/bookkeeping/tgp_no_aa_calls_failed_qc_chr{chromosome}.txt.gz',
        usecols=(1), dtype=int,
    )
    qced_sites = np.unique(qc_sites_mat).size
    # Extract the genotypes.
    gt, pos = load_gt_pos('tgp_mod_no_aa', chromosome)
    # Count the number of heterozygous sites per individual.
    het_counts = gt.count_het(axis=0)
    # Determine the heterozygosity in afr.
    afr_heterozygosity = np.nansum(per_site_heterozygosity(gt.take(tgp_meta_df[tgp_meta_df['SUPERPOP'] == 'AFR'].index.values, axis=1)))
    # Export the number of heterozygous sites per individual.
    np.savetxt(
        f'../muc19_results/tgp_mod_no_aa/tgp_het_sites_chr{chromosome}.txt.gz',
        [het_counts], fmt='%d',
    )
    # Export the afr heterozygosity results.
    np.savetxt(
        f'../muc19_results/tgp_mod_no_aa/afr_heterozygosity_eff_seq_len_chr{chromosome}.txt.gz',
        [np.array([afr_heterozygosity, chr_dicc[chromosome] - qced_sites])], fmt='%1.15f',
    )
    return

# Compute heterozygosity.
tgp_het_chromosome(chromosome=int(sys.argv[1]))