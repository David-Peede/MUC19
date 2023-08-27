# Import packages.
import allel
import gzip
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

# Define a function to load genotyope and positions arrays.
def load_callset_pos(prefix, chrom):
    # Intialize the file path.
    path = './zarr_data/{0}_chr{1}.zarr'.format(prefix, chrom)
    # Load the zarr array.
    zarr_array = zarr.open_group(path, mode='r')
    # Extract the genotype callset.
    callset = zarr_array['{0}/calldata/GT'.format(chrom)]
    # Load the positions.
    pos = allel.SortedIndex(zarr_array['{0}/variants/POS'.format(chrom)])
    return callset, pos

# Load the duplicates dataframe.
dup_df = pd.read_csv('../meta_data/tgp_dup_sites.csv')

# Intialize the chromosome.
chrom = int(sys.argv[1])
# Load the all site information.
all_sites_mat = np.loadtxt(
    f'./vcf_bookkeeping/tgp_all_archaics_merged_all_sites_qc/chr{chrom}_all_sites_report.txt',
    delimiter='\t', dtype=int,
)
# Extract the positions for the joint vcf.
_, vcf_pos = load_callset_pos('tgp_mod_arc_anc', chrom)
# Subset the dataframe.
sub_df = dup_df[dup_df['CHR'] == chrom]
# Extract the duplicate poistions.
dup_pos = sub_df['POS'].values
# Find the duplicate records in the vcf.
vcf_dups = np.intersect1d(vcf_pos, dup_pos)
# Create a mask for the duplicated sites.
dup_mask = np.in1d(all_sites_mat[:, 1], dup_pos)
# Mask the matrix for exporting.
export_mat = all_sites_mat[~dup_mask]
# Export the matrix.
np.savetxt(
    f'./vcf_bookkeeping/tgp_all_archaics_merged_all_sites_qc/chr{chrom}_all_sites_dedupped_report.csv',
    export_mat, fmt='%d', delimiter=',', newline='\n',
)
