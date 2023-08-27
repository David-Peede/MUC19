# Import packages.
import numpy as np
import sys

### sys.argv[1] = chromosome ###
### sys.argv[2] = region start ###
### sys.argv[3] = region end ###
### sys.argv[4] = label for output ###


# Define a function to find the effective sequence length for a set of coordinates.
def region_esl(chrom, start, end, label):
    # Load the all sites information.
    arc_all_sites_info = np.loadtxt(
        f'../vcf_data/vcf_bookkeeping/all_archaics_merged_all_sites_qc/chr{chrom}_all_sites_report.txt',
        delimiter='\t', dtype=int,
    )
    tgp_all_sites_info = np.loadtxt(
        f'../vcf_data/vcf_bookkeeping/tgp_all_archaics_merged_all_sites_qc/chr{chrom}_all_sites_dedupped_report.csv',
        delimiter=',', dtype=int,
    )
    # Intialize the coordinates.
    coord_s = start
    coord_e = end
    # Determine how the indicies of the sites that fall in this window.
    arc_sites_idx = np.where((coord_s <= arc_all_sites_info[:, 1]) & (arc_all_sites_info[:, 1] <= coord_e))[0]
    tgp_sites_idx = np.where((coord_s <= tgp_all_sites_info[:, 1]) & (tgp_all_sites_info[:, 1] <= coord_e))[0]
    # Find the effective sequence lengths.
    ## [ALT, CHA, VIN, DEN] ##
    arc_coord_eff_seq_len = arc_all_sites_info[arc_sites_idx, 2:].sum(axis=0)
    tgp_coord_eff_seq_len = tgp_all_sites_info[tgp_sites_idx, 2:].sum(axis=0)
    # Determine the effective sequence length per pairwise comparison.
    arc_den_alt = np.count_nonzero((arc_all_sites_info[arc_sites_idx, 2] == 1) & (arc_all_sites_info[arc_sites_idx, 5] == 1))
    arc_den_cha = np.count_nonzero((arc_all_sites_info[arc_sites_idx, 3] == 1) & (arc_all_sites_info[arc_sites_idx, 5] == 1))
    arc_den_vin = np.count_nonzero((arc_all_sites_info[arc_sites_idx, 4] == 1) & (arc_all_sites_info[arc_sites_idx, 5] == 1))
    arc_alt_cha = np.count_nonzero((arc_all_sites_info[arc_sites_idx, 3] == 1) & (arc_all_sites_info[arc_sites_idx, 2] == 1))
    arc_alt_vin = np.count_nonzero((arc_all_sites_info[arc_sites_idx, 4] == 1) & (arc_all_sites_info[arc_sites_idx, 2] == 1))
    arc_cha_vin = np.count_nonzero((arc_all_sites_info[arc_sites_idx, 4] == 1) & (arc_all_sites_info[arc_sites_idx, 3] == 1))
    arc_pw_list = [
        arc_den_alt, arc_den_cha, arc_den_vin,
        arc_alt_cha, arc_alt_vin, arc_cha_vin,
    ]
    tgp_den_alt = np.count_nonzero((tgp_all_sites_info[tgp_sites_idx, 2] == 1) & (tgp_all_sites_info[tgp_sites_idx, 5] == 1))
    tgp_den_cha = np.count_nonzero((tgp_all_sites_info[tgp_sites_idx, 3] == 1) & (tgp_all_sites_info[tgp_sites_idx, 5] == 1))
    tgp_den_vin = np.count_nonzero((tgp_all_sites_info[tgp_sites_idx, 4] == 1) & (tgp_all_sites_info[tgp_sites_idx, 5] == 1))
    tgp_alt_cha = np.count_nonzero((tgp_all_sites_info[tgp_sites_idx, 3] == 1) & (tgp_all_sites_info[tgp_sites_idx, 2] == 1))
    tgp_alt_vin = np.count_nonzero((tgp_all_sites_info[tgp_sites_idx, 4] == 1) & (tgp_all_sites_info[tgp_sites_idx, 2] == 1))
    tgp_cha_vin = np.count_nonzero((tgp_all_sites_info[tgp_sites_idx, 4] == 1) & (tgp_all_sites_info[tgp_sites_idx, 3] == 1))
    tgp_pw_list = [
        tgp_den_alt, tgp_den_cha, tgp_den_vin,
        tgp_alt_cha, tgp_alt_vin, tgp_cha_vin,
    ]
    # Intialize a matrix.
    esl_mat = np.empty((2, 11), dtype=int)
    # Append the matrix.
    esl_mat[0, :] = np.array(arc_coord_eff_seq_len.tolist()+arc_pw_list+[0], dtype=int) # First row is the archaic data set.
    esl_mat[1, :] = np.array(tgp_coord_eff_seq_len.tolist()+tgp_pw_list+[tgp_sites_idx.size], dtype=int) # Second row is the tgp+archaic data set.
    # Export the matrix.
    np.savetxt(
        f'../meta_data/arc_and_tgp_{label}_eff_seq_len.csv',
        esl_mat, fmt='%d', delimiter=',', newline='\n',
    )
    return

# Determine the effective sequence length.
region_esl(
    chrom=int(sys.argv[1]),
    start=int(sys.argv[2]), 
    end=int(sys.argv[3]),
    label=str(sys.argv[4]),
)
