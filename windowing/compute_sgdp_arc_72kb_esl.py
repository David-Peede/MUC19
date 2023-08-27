# Import packages.
import numpy as np

# Load the all sites information for chromosome 12.
sgdp_all_sites_chr12 = np.loadtxt(
    './vcf_data/vcf_bookkeeping/sgdp_all_archaics_merged_all_sites_qc/chr12_all_sites_report.txt',
    delimiter='\t', dtype=int,
)

# Intialize the 72kb coordinates.
coord_s = 40758000
coord_e = 40830000

# Determine how the indicies of the sites that fall in this window.
sgdp_sites_idx = np.where((coord_s <= sgdp_all_sites_chr12[:, 1]) & (sgdp_all_sites_chr12[:, 1] <= coord_e))[0]

# Find the effective sequence lengths.
## [ALT, CHA, VIN, DEN] ##
sgdp_coord_eff_seq_len = sgdp_all_sites_chr12[sgdp_sites_idx, 2:].sum(axis=0)
# Determine the effective sequence length per pairwise comparison.
sgdp_den_alt = np.count_nonzero((sgdp_all_sites_chr12[sgdp_sites_idx, 2] == 1) & (sgdp_all_sites_chr12[sgdp_sites_idx, 5] == 1))
sgdp_den_cha = np.count_nonzero((sgdp_all_sites_chr12[sgdp_sites_idx, 3] == 1) & (sgdp_all_sites_chr12[sgdp_sites_idx, 5] == 1))
sgdp_den_vin = np.count_nonzero((sgdp_all_sites_chr12[sgdp_sites_idx, 4] == 1) & (sgdp_all_sites_chr12[sgdp_sites_idx, 5] == 1))
sgdp_alt_cha = np.count_nonzero((sgdp_all_sites_chr12[sgdp_sites_idx, 3] == 1) & (sgdp_all_sites_chr12[sgdp_sites_idx, 2] == 1))
sgdp_alt_vin = np.count_nonzero((sgdp_all_sites_chr12[sgdp_sites_idx, 4] == 1) & (sgdp_all_sites_chr12[sgdp_sites_idx, 2] == 1))
sgdp_cha_vin = np.count_nonzero((sgdp_all_sites_chr12[sgdp_sites_idx, 4] == 1) & (sgdp_all_sites_chr12[sgdp_sites_idx, 3] == 1))
sgdp_pw_list = [
    sgdp_den_alt, sgdp_den_cha, sgdp_den_vin,
    sgdp_alt_cha, sgdp_alt_vin, sgdp_cha_vin,
]

# Intialize an array to store the results.
esl_array = np.array(sgdp_coord_eff_seq_len.tolist()+sgdp_pw_list+[sgdp_sites_idx.size], dtype=int)

# Export the results
np.savetxt(
    f'../meta_data/sgdp_72kb_eff_seq_len.csv',
    [esl_array], fmt='%d', delimiter=',', newline='\n',
)