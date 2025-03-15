# Import packages.
import numpy as np
import sys

### sys.argv[1] = chromosome ###
### sys.argv[2] = region start ###
### sys.argv[3] = region end ###
### sys.argv[4] = label for output ###


# Define a function to find the effective sequence length for a set of coordinates.
def region_esl(chrom, start, end, label):
    # Intialize the paths to the all sites information.
    single_arcs = {
        'den_masked_no_aa': f'../vcf_data/bookkeeping/den_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
        'alt_masked_no_aa': f'../vcf_data/bookkeeping/alt_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
        'cha_masked_no_aa': f'../vcf_data/bookkeeping/cha_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
        'vin_masked_no_aa': f'../vcf_data/bookkeeping/vin_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
    }
    all_arcs = {
        'arcs_masked_aa': f'../vcf_data/bookkeeping/all_archaics_masked_aa_calls_all_sites_chr{chrom}.txt.gz',
        'arcs_masked_no_aa': f'../vcf_data/bookkeeping/all_archaics_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
    }
    only_hum = {
        'tgp_mod_no_aa': f'../vcf_data/bookkeeping/tgp_no_aa_calls_failed_qc_chr{chrom}.txt.gz',
        'sgdp_no_aa': f'../vcf_data/bookkeeping/sgdp_no_aa_calls_failed_qc_chr{chrom}.txt.gz',
    }
    hum_single_arcs = {
        'tgp_den_masked_aa': f'../vcf_data/bookkeeping/tgp_den_masked_aa_calls_all_sites_chr{chrom}.txt.gz',
        'tgp_den_masked_no_aa': f'../vcf_data/bookkeeping/tgp_den_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
        'tgp_alt_masked_aa': f'../vcf_data/bookkeeping/tgp_alt_masked_aa_calls_all_sites_chr{chrom}.txt.gz',
        'tgp_alt_masked_no_aa': f'../vcf_data/bookkeeping/tgp_alt_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
        'tgp_cha_masked_aa': f'../vcf_data/bookkeeping/tgp_cha_masked_aa_calls_all_sites_chr{chrom}.txt.gz',
        'tgp_cha_masked_no_aa': f'../vcf_data/bookkeeping/tgp_cha_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
        'tgp_vin_masked_aa': f'../vcf_data/bookkeeping/tgp_vin_masked_aa_calls_all_sites_chr{chrom}.txt.gz',
        'tgp_vin_masked_no_aa': f'../vcf_data/bookkeeping/tgp_vin_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
        'sgdp_den_masked_aa': f'../vcf_data/bookkeeping/sgdp_den_masked_aa_calls_all_sites_chr{chrom}.txt.gz',
        'sgdp_den_masked_no_aa': f'../vcf_data/bookkeeping/sgdp_den_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
        'sgdp_alt_masked_aa': f'../vcf_data/bookkeeping/sgdp_alt_masked_aa_calls_all_sites_chr{chrom}.txt.gz',
        'sgdp_alt_masked_no_aa': f'../vcf_data/bookkeeping/sgdp_alt_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
        'sgdp_cha_masked_aa': f'../vcf_data/bookkeeping/sgdp_cha_masked_aa_calls_all_sites_chr{chrom}.txt.gz',
        'sgdp_cha_masked_no_aa': f'../vcf_data/bookkeeping/sgdp_cha_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
        'sgdp_vin_masked_aa': f'../vcf_data/bookkeeping/sgdp_vin_masked_aa_calls_all_sites_chr{chrom}.txt.gz',
        'sgdp_vin_masked_no_aa': f'../vcf_data/bookkeeping/sgdp_vin_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
    }
    hum_all_arcs = {
        'tgp_arcs_masked_no_aa': f'../vcf_data/bookkeeping/tgp_all_archaics_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
        'sgdp_arcs_masked_no_aa': f'../vcf_data/bookkeeping/sgdp_all_archaics_masked_no_aa_calls_all_sites_chr{chrom}.txt.gz',
        'tgp_arcs_masked_aa': f'../vcf_data/bookkeeping/tgp_all_archaics_masked_aa_calls_all_sites_chr{chrom}.txt.gz',
        'sgdp_arcs_masked_aa': f'../vcf_data/bookkeeping/sgdp_all_archaics_masked_aa_calls_all_sites_chr{chrom}.txt.gz',
    }
    # For every single arachaic dataset.
    for prefix in single_arcs:
        # Load the all sites information.
        all_sites_array = np.loadtxt(
            single_arcs[prefix],
            usecols=(1), dtype=int,
        )
        # Determine how many sites fall in this region.
        region_mask = (start <= all_sites_array) & (all_sites_array <= end)
        # Export the results.
        np.savetxt(
            f'./{prefix}/{label}_eff_seq_len.txt.gz',
            [[region_mask.sum()]], fmt='%d',
        )
    # For every all archaics dataset.
    for prefix in all_arcs:
        # Load the all sites information.
        all_sites_mat = np.loadtxt(
            all_arcs[prefix],
            usecols=(0, 1, 2, 3, 4, 5), dtype=int,
        )
        # Determine how the indicies of the sites that fall in this window.
        region_idx = np.where((start <= all_sites_mat[:, 1]) & (all_sites_mat[:, 1] <= end))[0]
        # Determine the effective sequence length for each sample.
        arcs_esl = all_sites_mat[region_idx, 2:].sum(axis=0) # [ALT, CHA, VIN, DEN]
        # Determine the effective sequence length per pairwise comparison.
        den_alt = np.count_nonzero((all_sites_mat[region_idx, 2] == 1) & (all_sites_mat[region_idx, 5] == 1))
        den_cha = np.count_nonzero((all_sites_mat[region_idx, 3] == 1) & (all_sites_mat[region_idx, 5] == 1))
        den_vin = np.count_nonzero((all_sites_mat[region_idx, 4] == 1) & (all_sites_mat[region_idx, 5] == 1))
        alt_cha = np.count_nonzero((all_sites_mat[region_idx, 3] == 1) & (all_sites_mat[region_idx, 2] == 1))
        alt_vin = np.count_nonzero((all_sites_mat[region_idx, 4] == 1) & (all_sites_mat[region_idx, 2] == 1))
        cha_vin = np.count_nonzero((all_sites_mat[region_idx, 4] == 1) & (all_sites_mat[region_idx, 3] == 1))
        # Update the list.
        all_arcs_esl = arcs_esl.tolist()
        all_arcs_esl.extend([
            den_alt, den_cha, den_vin,
            alt_cha, alt_vin, cha_vin,
        ])
        # Export the results.
        np.savetxt(
            f'./{prefix}/{label}_eff_seq_len.txt.gz',
            [all_arcs_esl], fmt='%d',
        )
    # For the only only human dataset.
    for prefix in only_hum:
        # Load the qc site information.
        qc_sites_array = np.loadtxt(
            only_hum[prefix],
            usecols=(1), dtype=int,
        )
        # Determine the unique positions to avoid counting errors.
        qced_sites = np.unique(qc_sites_array)
        # Determine how many sites fall in this region.
        region_mask = (start <= qced_sites) & (qced_sites <= end)
        # Determine the effective sequence length.
        hum_esl = ((end + 1) - start) - region_mask.sum()
        # Export the results.
        np.savetxt(
            f'./{prefix}/{label}_eff_seq_len.txt.gz',
            [[hum_esl]], fmt='%d',
        )
    # For every tgp + single archaic dataset.
    for prefix in hum_single_arcs:
        # Load the all sites information.
        all_sites_array = np.loadtxt(
            hum_single_arcs[prefix],
            usecols=(1), dtype=int,
        )
        # Determine how many sites fall in this region.
        region_mask = (start <= all_sites_array) & (all_sites_array <= end)
        # Export the results.
        np.savetxt(
            f'./{prefix}/{label}_eff_seq_len.txt.gz',
            [[region_mask.sum()]], fmt='%d',
        )
    # For every tgp + all archaics dataset.
    for prefix in hum_all_arcs:
        # Load the all site information.
        all_sites_mat = np.loadtxt(
            hum_all_arcs[prefix],
            usecols=(0, 1, 2, 3, 4, 5), dtype=int,
        )
        # Determine how the indicies of the sites that fall in this window.
        region_idx = np.where((start <= all_sites_mat[:, 1]) & (all_sites_mat[:, 1] <= end))[0]
        # Determine the effective sequence length for each sample.
        arcs_esl = all_sites_mat[region_idx, 2:].sum(axis=0) # [ALT, CHA, VIN, DEN]
        # Determine the effective sequence length per pairwise comparison.
        den_alt = np.count_nonzero((all_sites_mat[region_idx, 2] == 1) & (all_sites_mat[region_idx, 5] == 1))
        den_cha = np.count_nonzero((all_sites_mat[region_idx, 3] == 1) & (all_sites_mat[region_idx, 5] == 1))
        den_vin = np.count_nonzero((all_sites_mat[region_idx, 4] == 1) & (all_sites_mat[region_idx, 5] == 1))
        alt_cha = np.count_nonzero((all_sites_mat[region_idx, 3] == 1) & (all_sites_mat[region_idx, 2] == 1))
        alt_vin = np.count_nonzero((all_sites_mat[region_idx, 4] == 1) & (all_sites_mat[region_idx, 2] == 1))
        cha_vin = np.count_nonzero((all_sites_mat[region_idx, 4] == 1) & (all_sites_mat[region_idx, 3] == 1))
        # Update the list.
        all_arcs_esl = arcs_esl.tolist()
        all_arcs_esl.extend([
            den_alt, den_cha, den_vin,
            alt_cha, alt_vin, cha_vin,
        ])
        # Export the results.
        np.savetxt(
            f'./{prefix}/{label}_eff_seq_len.txt.gz',
            [all_arcs_esl], fmt='%d',
        )
    return

# Determine the effective sequence length.
region_esl(
    chrom=int(sys.argv[1]),
    start=int(sys.argv[2]), 
    end=int(sys.argv[3]),
    label=str(sys.argv[4]),
)