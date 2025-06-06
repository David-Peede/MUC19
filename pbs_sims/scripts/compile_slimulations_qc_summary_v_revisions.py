# Import the packages.
import pandas as pd
import numpy as np
import allel
import os
import re
import sys


### sys.argv[1] = neutral or negative ###
### sys.argv[2] = uniform or recombination map ###
### sys.argv[3] = slimulation vcf prefix ###


# Define a function to qc the slimulations.
def qc_slimulations(sim_path, fprefix):
    # Initialize lists to store results.
    vcf_files, rng_seeds, org_rep_ids = [], [], []
    n_den_742kb, n_nea_742kb, n_arc_742kb = [], [], []
    n_den_72kb, n_nea_72kb, n_arc_72kb = [], [], []

    # Compile the regular expression pattern for matching filenames.
    fpattern = re.compile(fr'{fprefix}_.*_\d+\.txt$')
    # Grab all the files in the simulation directory.
    all_files = os.listdir(sim_path)
    # Find the simulation results.
    fname_list = [file for file in all_files if fpattern.match(file)]

    # For every vcf file.
    for vcf_file in fname_list:
        # Unpack the simulation information.
        seed, rep_txt = vcf_file.split('_')[-2:]
        # Load the vcf.
        sim_vcf = allel.read_vcf(f'{sim_path}/{vcf_file}', fields='*')
        # Load the genotype matrix.
        sim_gt = allel.GenotypeArray(sim_vcf['calldata/GT'])
        # Intialize the conditions.
        isMulti = sim_vcf['variants/MULTIALLELIC']
        popOrigin = sim_vcf['variants/PO'][~isMulti]
        positions = sim_vcf['variants/POS'][~isMulti]
    
        # Count the number of derived alleles in mxl.
        mxl_dac = sim_gt[:, 202:].compress((~isMulti), axis=0).count_alleles()
        # Compute the derived allele frequencies.
        mxl_daf = mxl_dac.to_frequencies()
        mxl_daf = mxl_daf[:, 0] - 1 if mxl_dac.shape[1] == 1 else mxl_daf[:, 1]
        # Determine the sites with at least one derived allele.
        has_der_mut = mxl_daf > 0
    
        # Intialize the snp conditions.
        is72kb = (positions >= 487000) & (positions <= 559000)
        isDen742kb = (popOrigin == 3) & has_der_mut
        isNea742kb = (popOrigin == 2) & has_der_mut
        isArc742kb = isDen742kb | isNea742kb
        isDen72kb = isDen742kb & is72kb
        isNea72kb = isNea742kb & is72kb
        isArc72kb = isDen72kb | isNea72kb
        
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
    })
    # Export the results.
    qc_df.to_csv(f'{sim_path}/slim_reps_qc_summary.csv.gz', index=False)
    return


# Generate the qc files for the slimulations.
qc_slimulations(f'../smodel/{str(sys.argv[1])}/{str(sys.argv[2])}', str(sys.argv[3]))