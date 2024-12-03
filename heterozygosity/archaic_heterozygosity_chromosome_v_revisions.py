# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###
### sys.argv[2] = archaic ###

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

# Define a function to calculate alternative allele frequencies for a single individual.
def calc_ind_alt_freqs(gt):
    # Compute the frequency for each site.
    raw_freqs = np.nansum(gt, axis=2).flatten() / 2
    # Set missing data to Nan
    alt_freqs = np.where(raw_freqs == -1, np.nan, raw_freqs)
    return alt_freqs

# Define a function to compute per site heterozygosity for a single individual.
def ind_per_site_heterozygosity(gt):
    # Compute the alternative allele frequency.
    p = calc_ind_alt_freqs(gt)
    # Compute the per site heterozygosity.
    per_site_h = (1 - ((p ** 2) + ((1 - p) ** 2)))
    return per_site_h

# Define a function to compute the number of heterozygous sites and heterozygoisty for a chromosome.
def archaic_het_chromosome(chromosome, archaic):
    # Extract the genotypes.
    gt, _ = load_gt_pos(f'{archaic.lower()}_masked_no_aa', chromosome)
    # Load the all sites matrix.
    arc_all_sites_info = np.loadtxt(
        f'../vcf_data/bookkeeping/{archaic.lower()}_masked_no_aa_calls_all_sites_chr{chromosome}.txt.gz',
        usecols=(1), dtype=int,
    )
    # Determine the effective sequence length for this chromosome.
    arc_chr_esl = arc_all_sites_info.size
    # Count the number of heterozygous sites for the individual.
    het_counts = gt.count_het(axis=0)[0]
    # Compute heterozygosity for the individual.
    heterozygosity = np.nansum(ind_per_site_heterozygosity(gt=gt))
    # Export the results.
    np.savetxt(
        f'../muc19_results/{archaic.lower()}_masked_no_aa/archaic_het_sites_heterozygosity_eff_seq_len_chr{chromosome}.txt.gz',
        [np.array([het_counts, heterozygosity, arc_chr_esl])], fmt='%1.15f',
    )
    return

# Compute heterozygosity.
archaic_het_chromosome(chromosome=int(sys.argv[1]), archaic=str(sys.argv[2]))