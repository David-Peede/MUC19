# Import packages.
import allel
from itertools import product
import numpy as np
import pandas as pd
import sys

### sys.argv[1] = selection coefficient ###


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
    return max(pbs, 0), a_b_fst, a_c_fst, c_b_fst

# Define a function to QC the replicates.
def extract_positive_sim_info(sel_coeff):
    # Intialize the path prefix.
    sim_path = '../smodel/positive_neutral'
    # Load the qc summary.
    qc_df = pd.read_csv(f'{sim_path}_{sel_coeff}/slim_reps_qc_summary.csv.gz')
    # Convert the qc dataframe to a dictionary.
    qc_info = {
        key: np.array(value) for key, value in qc_df.to_dict(orient='list').items()
    }
    # Intialize masks.
    has_den_snps_72kb = qc_info['n_den_72kb'] > 0
    has_nea_snps_72kb = qc_info['n_nea_72kb'] > 0
    has_both_snps_72kb = has_den_snps_72kb & has_nea_snps_72kb
    # Subset the qc'ed simulatuions.
    qced_df = qc_df[has_both_snps_72kb].reset_index(drop=True)
    # Extract the qc'ed replicate ids.
    qced_reps = qced_df['org_rep_id'].values
    # Extract teh qc'ed vcf files.
    qced_vcfs = qced_df['vcf_file'].values
    
    # Intialize a dictionary.
    sel_freqs = {
        'MXB_before': [], 'EAS_before': [], 'EUR_before': [],
        'MXB_after': [], 'MXL_after': [], 'EAS_after': [], 'EUR_after': [],
    }
    # For every replicate.
    for rep in qced_reps:
        # Open the slimulation output.
        with open(f'{sim_path}_{sel_coeff}/slimulation_{rep}.txt', 'r') as slim_file:
            # For every line.
            for line in slim_file:
                # If the line starts with the identifier.
                if line.startswith('%'):
                    # Extract the selected allele frequency.
                    if 'FREQmTag_EUR_before:' in line:
                        sel_freqs['EUR_before'].append(float(line.split(':')[1].strip()))
                    elif 'FREQmTag_EAS_before:' in line:
                        sel_freqs['EAS_before'].append(float(line.split(':')[1].strip()))
                    elif 'FREQmTag_MXB_before:' in line:
                        sel_freqs['MXB_before'].append(float(line.split(':')[1].strip()))
                    elif 'FREQmTag_EUR_after:' in line:
                        sel_freqs['EUR_after'].append(float(line.split(':')[1].strip()))
                    elif 'FREQmTag_EAS_after:' in line:
                        sel_freqs['EAS_after'].append(float(line.split(':')[1].strip()))
                    elif 'FREQmTag_MXB_after:' in line:
                        sel_freqs['MXB_after'].append(float(line.split(':')[1].strip()))
                    elif 'FREQmTag_MXL_after:' in line:
                        sel_freqs['MXL_after'].append(float(line.split(':')[1].strip()))
   # Update the dataframe with the results
    for col_key, col_val in sel_freqs.items():
        qced_df[col_key] = col_val 
    
    # Intialize a dictionary to store the results.
    region_dicc = {}
    # Intialize the different data partitions.
    regions = ['742kb', '72kb']
    snp_types = ['all', 'arc']
    # For every snp type at every region.
    for region, snp_type in product(regions, snp_types):
        # Intialize the dictionary.
        region_dicc[f'n_seg_{snp_type}_snps_{region}'] = []
        region_dicc[f'pbs_{snp_type}_snps_{region}'] = []
        region_dicc[f'fst_mxl_eas_{snp_type}_snps_{region}'] = []
        region_dicc[f'fst_mxl_eur_{snp_type}_snps_{region}'] = []
        region_dicc[f'fst_eas_eur_{snp_type}_snps_{region}'] = []

    # For every qc'ed vcf file.
    for vcf_file in qced_vcfs:
        # Load the qc'ed vcf.
        qced_vcf = allel.read_vcf(f'{sim_path}_{sel_coeff}/{vcf_file}', fields='*')
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
                # Compute the number of segregating sites.
                n_seg = mask_dicc[region][snp_type].sum()
                # Compute pbs and fst.
                mxl_pbs, mxl_eas_fst, mxl_eur_fst, eas_eur_fst = calc_pbs_per_region(
                    mxl_gt.compress(mask_dicc[region][snp_type], axis=0),
                    eas_gt.compress(mask_dicc[region][snp_type], axis=0),
                    eur_gt.compress(mask_dicc[region][snp_type], axis=0),
                )
                # Update the dictionary.
                region_dicc[f'n_seg_{snp_type}_snps_{region}'].append(n_seg)
                region_dicc[f'pbs_{snp_type}_snps_{region}'].append(mxl_pbs)
                region_dicc[f'fst_mxl_eas_{snp_type}_snps_{region}'].append(mxl_eas_fst)
                region_dicc[f'fst_mxl_eur_{snp_type}_snps_{region}'].append(mxl_eur_fst)
                region_dicc[f'fst_eas_eur_{snp_type}_snps_{region}'].append(eas_eur_fst)

    # Update the dataframe with the results
    for col_key, col_val in region_dicc.items():
        qced_df[col_key] = col_val
        
    # Export the results.
    qced_df.to_csv(f'{sim_path}_{sel_coeff}/qced_slim_reps_info.csv.gz', index=False)
    return

# Compile the slimulation information.
extract_positive_sim_info(str(sys.argv[1]))