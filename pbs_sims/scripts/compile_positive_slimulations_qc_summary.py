# Import the packages.
import pandas as pd
import numpy as np
import allel
import os
import re
import sys


### sys.argv[1] = selection coefficient ###


# Define a function to qc the slimulations.
def qc_positive_sims(sel_coeff):
    # Intialize the path prefix.
    sim_path = '../smodel/positive_neutral'
    # Intialize a suffix to selection coeffiecient map.
    suffix_map = {'1': 0.1, '01': 0.01, '0015': 0.0015, '001': 0.001, '0005': 0.0005}
    # Load the vcf files.
    file_list = [
        file for file in os.listdir(f'{sim_path}_{sel_coeff}')
        if os.path.isfile(os.path.join(f'{sim_path}_{sel_coeff}', file)) and file.startswith('VCF')
    ]
    
    # Initialize lists to store results.
    vcf_files, rng_seeds, org_rep_ids = [], [], []
    n_den_742kb, n_nea_742kb, n_arc_742kb = [], [], []
    n_den_72kb, n_nea_72kb, n_arc_72kb = [], [], []
    sel_pos, sel_origin, sel_mxl_daf, sel_eas_daf, sel_eur_daf = [], [], [], [], []
    # Intialize a list of bad replicates.
    bad_reps = []
    # Intialize a dictionary.
    arc_muts = {2: 'NEA', 3: 'DEN'}
    
    # For every vcf file.
    for vcf_file in file_list:
        # Unpack the simulation information.
        seed, rep_txt = vcf_file.split('_')[-2:]
        # Load the vcf.
        sim_vcf = allel.read_vcf(f'{sim_path}_{sel_coeff}/{vcf_file}', fields='*')
        # Load the genotype matrix.
        sim_gt = allel.GenotypeArray(sim_vcf['calldata/GT'])
        # Intialize the conditions.
        isMulti = sim_vcf['variants/MULTIALLELIC']
        popOrigin = sim_vcf['variants/PO'][~isMulti]
        positions = sim_vcf['variants/POS'][~isMulti]
        sel_coeffs = sim_vcf['variants/S'][~isMulti]
        # Count the number of derived alleles in each population.
        mxl_dac = sim_gt[:, 202:].compress((~isMulti), axis=0).count_alleles()
        eas_dac = sim_gt[:, 99:202].compress((~isMulti), axis=0).count_alleles()
        eur_dac = sim_gt[:, 0:99].compress((~isMulti), axis=0).count_alleles()
        # Compute the derived allele frequencies.
        mxl_daf = mxl_dac.to_frequencies()
        mxl_daf = mxl_daf[:, 0] - 1 if mxl_dac.shape[1] == 1 else mxl_daf[:, 1]
        eas_daf = eas_dac.to_frequencies()
        eas_daf = eas_daf[:, 0] - 1 if eas_dac.shape[1] == 1 else eas_daf[:, 1]
        eur_daf = eur_dac.to_frequencies()
        eur_daf = eur_daf[:, 0] - 1 if eur_dac.shape[1] == 1 else eur_daf[:, 1]
        # Determine the sites with at least one derived allele in mxl.
        has_der_mut = mxl_daf > 0

        # Intialize the snp conditions.
        is72kb = (positions >= 487000) & (positions <= 559000)
        isDen742kb = (popOrigin == 3) & has_der_mut
        isNea742kb = (popOrigin == 2) & has_der_mut
        isArc742kb = isDen742kb | isNea742kb
        isDen72kb = isDen742kb & is72kb
        isNea72kb = isNea742kb & is72kb
        isArc72kb = isDen72kb | isNea72kb
        is_selected = sel_coeffs == suffix_map[sel_coeff]

        # If there is no or more than one selected site.
        if is_selected.sum() == 0 or is_selected.sum() > 1:
            # Update the list.
            bad_reps.append(vcf_file)
        # Else.
        else:
            # Update the lists.
            vcf_files.append(vcf_file)
            rng_seeds.append(int(seed))
            org_rep_ids.append(int(rep_txt[:-4]))
            n_den_742kb.append(isDen742kb.sum())
            n_nea_742kb.append(isNea742kb.sum())
            n_arc_742kb.append(isArc742kb.sum())
            n_den_72kb.append(isDen72kb.sum())
            n_nea_72kb.append(isNea72kb.sum())
            n_arc_72kb.append(isArc72kb.sum())
            sel_pos.append(positions[is_selected][0])
            sel_origin.append(arc_muts[popOrigin[is_selected][0]])
            sel_mxl_daf.append(mxl_daf[is_selected][0])
            sel_eas_daf.append(eas_daf[is_selected][0])
            sel_eur_daf.append(eur_daf[is_selected][0])

    # Create a dataframe from the lists.
    qc_df = pd.DataFrame({
        'vcf_file': vcf_files,
        'rng_seed': rng_seeds,
        'org_rep_id': org_rep_ids,
        'n_den_742kb': n_den_742kb,
        'n_nea_742kb': n_nea_742kb,
        'n_arc_742kb': n_arc_742kb,
        'n_den_72kb': n_den_72kb,
        'n_nea_72kb': n_nea_72kb,
        'n_arc_72kb': n_arc_72kb,
        'sel_pos': sel_pos,
        'sel_origin': sel_origin,
        'sel_mxl_daf': sel_mxl_daf,
        'sel_eas_daf': sel_eas_daf,
        'sel_eur_daf': sel_eur_daf,
    })
    # Export the results.
    qc_df.to_csv(f'{sim_path}_{sel_coeff}/slim_reps_qc_summary.csv.gz', index=False)
    np.savetxt(
        f'{sim_path}_{sel_coeff}/bad_reps.txt.gz',
        np.array([bad_reps]), fmt='%s',
    )
    return

# QC the simulations.
qc_positive_sims(str(sys.argv[1]))