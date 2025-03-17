# Import packages.
import allel
from itertools import product
import numpy as np
import pandas as pd
import sys

### sys.argv[1] = selection model ###

# Define a  function to compute Hudson's Fst estimator for a region.
def calc_fst_per_region(ac1, ac2):
    # Compute.
    num, den = allel.hudson_fst(ac1, ac2)
    # Account for the denominator being zero.
    fst = np.nansum(num) / np.nansum(den) if np.nansum(den) != 0 else 0
    # Correct for negative Fst values.
    return max(fst, 0)

# Define a function to calculate PBS per region.
def calc_pbs_per_region(pop_a_gt, pop_b_gt, pop_c_gt):
    # Determine allele counts.
    a_ac = pop_a_gt.count_alleles()
    b_ac = pop_b_gt.count_alleles()
    c_ac = pop_c_gt.count_alleles()
    # Calculate Fst.
    a_b_fst = calc_fst_per_region(a_ac, b_ac)
    a_c_fst = calc_fst_per_region(a_ac, c_ac)
    c_b_fst = calc_fst_per_region(c_ac, b_ac)
    # Correct for Fst values of 1 that will lead to inf.
    a_b_fst = min(a_b_fst, 0.99999)
    a_c_fst = min(a_c_fst, 0.99999)
    c_b_fst = min(c_b_fst, 0.99999)
    # Calculate PBS.
    pbs = (
        ((np.log(1.0 - a_b_fst) * -1.0) +\
         (np.log(1.0 - a_c_fst) * -1.0) -\
         (np.log(1.0 - c_b_fst) * -1.0)) / 2.0
    )
    return max(pbs, 0)

# Define a function to compute Hudson's estimator of Fst per site
def calc_fst_per_site(ac1, ac2):
    # Compute.
    num, den = allel.hudson_fst(ac1, ac2)
    raw_fst = num / den
    # Correct for when the denominator is zero.
    raw_fst = np.where(den == 0, 0, raw_fst)
    # Correct for negative Fst values.
    raw_fst = np.where(raw_fst < 0, 0, raw_fst)
    # Correct for Fst values of 1 that will lead to inf.
    raw_fst = np.where(raw_fst == 1, 0.99999, raw_fst)
    return raw_fst

# Define a function to calculate PBS.
def calc_pbs_per_site(pop_a_gt, pop_b_gt, pop_c_gt):
    # Determine allele counts.
    a_ac = pop_a_gt.count_alleles()
    b_ac = pop_b_gt.count_alleles()
    c_ac = pop_c_gt.count_alleles()
    # Calculate Fst with corrections.
    a_b_fst = calc_fst_per_site(a_ac, b_ac)
    a_c_fst = calc_fst_per_site(a_ac, c_ac)
    c_b_fst = calc_fst_per_site(c_ac, b_ac)
    # Calculate PBS.
    pbs = (
        ((np.log(1.0 - a_b_fst) * -1.0) +\
         (np.log(1.0 - a_c_fst) * -1.0) -\
         (np.log(1.0 - c_b_fst) * -1.0)) / 2.0
    )
    return np.where(pbs < 0, 0, pbs)


# Load the qc'ed dataframe.
qced_df = pd.read_csv(f'../data/{sys.argv[1]}_recomb_map_10k_qced_reps.csv.gz')
# Extract the qc'ed vcf files and the site info files.
qced_vcfs = qced_df['vcf_file'].values
info_files = qced_df['site_info_file'].values

# Intialize dictionaries to store the results.
rep_dicc = {}
snp_dicc = {}
# Intialize the different data partitions.
regions = ['742kb', '72kb']
snp_types = ['all', 'arc']
# For every snp type at every region.
for region, snp_type in product(regions, snp_types):
    # Intialize the dictionary.
    rep_dicc[f'geq_{snp_type}_snps_{region}'] = []
    rep_dicc[f'pbs_{snp_type}_snps_{region}'] = []
    rep_dicc[f'prop_out_pbs_{snp_type}_snps_{region}'] = []
    snp_dicc[f'{snp_type}_snp_freqs_{region}']  = []
    snp_dicc[f'pbs_{snp_type}_snps_{region}'] = []


# For every qc'ed vcf.
for path_to_qced_vcf, path_to_site_info in zip(qced_vcfs, info_files):
    # Load the qc'ed vcf.
    qced_vcf = allel.read_vcf(path_to_qced_vcf, fields='*')
    # Load the genotype matrix.
    qced_gt = allel.GenotypeArray(qced_vcf['calldata/GT'])

    # Intialize the conditions.
    isMulti = qced_vcf['variants/MULTIALLELIC']
    popOrigin = qced_vcf['variants/PO'][~isMulti]
    positions = qced_vcf['variants/POS'][~isMulti]
    mutIds = qced_vcf['variants/MT'][~isMulti]
    is72kb = (positions >= 487000) & (positions <= 559000)
    isArc = (popOrigin == 3) | (popOrigin == 2)
    
    # Intialize the genotype matricies.
    mxl_gt = qced_gt[:, 202:].compress((~isMulti), axis=0)
    eas_gt = qced_gt[:, 99:202].compress((~isMulti), axis=0)
    eur_gt = qced_gt[:, 0:99].compress((~isMulti), axis=0)
    
    # Count the number of derived alleles in mxl and compute the derived allele frequencies.
    mxl_dac = mxl_gt.count_alleles()
    mxl_daf = mxl_dac.to_frequencies()
    mxl_daf = mxl_daf[:, 0] - 1 if mxl_dac.shape[1] == 1 else mxl_daf[:, 1]
    # Determine the sites with at least one derived allele.
    has_der_mut = mxl_daf > 0
    # Determine the sites with the derived allele is at a frequency of 30% or higher.
    has_geq_daf = mxl_daf >= 0.3
    # Compute the per site pbs values.
    per_site_pbs = calc_pbs_per_site(mxl_gt, eas_gt, eur_gt)
    # Determine the simulated PBS values greater than the observed 99.95th percentile.
    is_pbs_out = per_site_pbs > 0.25886600381183444
    
    # Intialize a mask dictionary.
    mask_dicc = {
        '742kb': {
            'all': np.full(positions.size, True),
            'arc': isArc & has_der_mut,
        },
        '72kb': {
            'all': is72kb,
            'arc': isArc & is72kb & has_der_mut,
        }
    }
    
    
    # For every snp type at every region.
    for region, snp_type in product(regions, snp_types):
            # Update the dictionaries
            rep_dicc[f'geq_{snp_type}_snps_{region}'].append(
                (has_geq_daf & mask_dicc[region][snp_type]).sum()
            )
            rep_dicc[f'pbs_{snp_type}_snps_{region}'].append(calc_pbs_per_region(
                mxl_gt.compress(mask_dicc[region][snp_type], axis=0),
                eas_gt.compress(mask_dicc[region][snp_type], axis=0),
                eur_gt.compress(mask_dicc[region][snp_type], axis=0),
            ))
            rep_dicc[f'prop_out_pbs_{snp_type}_snps_{region}'].append(
                (is_pbs_out & mask_dicc[region][snp_type]).sum()
            )
            snp_dicc[f'{snp_type}_snp_freqs_{region}'].extend(mxl_daf[mask_dicc[region][snp_type]])
            snp_dicc[f'pbs_{snp_type}_snps_{region}'].extend(per_site_pbs[mask_dicc[region][snp_type]])

# Export the per replicate results.
results_path = '/users/dpeede/data/dpeede/08_muc19/muc19_results/mxl_slimulations'
pd.DataFrame(rep_dicc).to_csv(f'{results_path}/{sys.argv[1]}_per_10k_replicates.csv.gz', index=False)
# Export the per snp results.
for key, values in snp_dicc.items():
    np.savetxt(
        f'{results_path}/{sys.argv[1]}_{key}_per_snp.txt.gz',
        [values], fmt='%1.15f',
    )