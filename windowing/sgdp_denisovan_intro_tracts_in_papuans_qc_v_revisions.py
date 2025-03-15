# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = papuan prefix ###
### sys.argv[2] = chromosome ###
### sys.argv[3] = dataset prefix ###


# Define a function to calculate the number of segregating sites.
def count_seg_sites(gt):
    # Count the number of segregsting sites.
    seg_sites = gt.count_alleles().count_segregating()
    return seg_sites

# Define a function to load genotyope and positions arrays.
def load_callset_pos(prefix, chrom):
    # Intialize the file path.
    path = f'../zarr_data/{prefix}_chr{chrom}.zarr'
    # Load the zarr array.
    zarr_array = zarr.open_group(path, mode='r')
    # Extract the genotype callset.
    callset = zarr_array[f'{chrom}/calldata/GT']
    # Load the positions.
    pos = allel.SortedIndex(zarr_array[f'{chrom}/variants/POS'])
    return callset, pos

# Define a function to qc the papuan introgressed tracts.
def pap_tracts_qc(pap, chromosome, prefix):
    # Intialize a dictionary to store the all sites path information.
    path_dicc = {
        'sgdp_den_masked_aa': '../vcf_data/bookkeeping/sgdp_den_masked_aa_calls_all_sites',
        'sgdp_den_masked_no_aa': '../vcf_data/bookkeeping/sgdp_den_masked_no_aa_calls_all_sites',
    }
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/sgdp.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Identify all the modern human indicies.
    hum_idx = meta_df.index.values
    # Intitialize the all sites path, window path, and column variable.
    hum_all_sites_path = f'{path_dicc[prefix]}_chr{chromosome}.txt.gz'
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos(prefix, chromosome)
    # Load the all site information.
    all_sites_mat = np.loadtxt(
        hum_all_sites_path,
        usecols=(0, 1),
        dtype=int,
    )
    # For each haplotype.
    for hap in ['hap1', 'hap2']:
        # Load the tracts for the papuan individual.
        pap_tracts_df = pd.read_csv(f'../hmmix_tracts/den_tracts_in_pap/den_intro_tracts_in_{pap}_{hap}.csv.gz')
        # Subset the tracts for this chromosome.
        chr_tracts_df = pap_tracts_df[pap_tracts_df['chrom'] == chromosome]
        # Extract the start and stop positions.
        starts = chr_tracts_df['start'].values
        ends = chr_tracts_df['end'].values
        # Determine the number of tracts.
        n_tracts = chr_tracts_df.shape[0]
        # Intialize the output path.
        hum_tract_path = f'./{prefix}/{pap}_{hap}_den_intro_tracts_summary_chr{chromosome}.csv.gz'
        # Intialize a results matrix to store the results.
        results_mat = np.empty((n_tracts, 6))
        # For every tracts...
        for idx in range(n_tracts):
            # Extract the start and stop positions.
            start = starts[idx]
            stop = ends[idx] # Already adjusted for inclusive indexing in the tracts file.
            # Extract the variable sites information.
            var_idx = np.where(((start <= all_pos) & (all_pos <= stop)))[0]
            # Extract the effective sequence length information.
            esl_idx = np.where(((start <= all_sites_mat[:, 1]) & (all_sites_mat[:, 1] <= stop)))[0]
            # Determine the effective sequence length.
            hum_arc_esl = esl_idx.size
            # If the window has variant sites to do calculations on...
            if (var_idx.size > 0):
                # Locate the window.
                tract_loc = all_pos.locate_range(start, stop)
                # Determine the number of segregating sites.
                s = count_seg_sites(gt=allel.GenotypeArray(callset[tract_loc]))
                # Append the results.
                results_mat[idx, :] = np.array(
                    [idx, start, stop, hum_arc_esl, s, 1], dtype=object,
                )
            # Else-if there are sites that passed qc for this window.
            elif (esl_idx.size > 0):
                # Append the results.
                results_mat[idx, :] = np.array(
                    [idx, start, stop, hum_arc_esl, 0, 1], dtype=object,
                )
            # Else.
            else:
                # Append the results.
                results_mat[idx, :] = np.array(
                    [idx, start, stop, 0, 0, 1], dtype=object,
                )
        # Export the the results matrix as a pandas dataframe.
        tract_df = pd.DataFrame(
            results_mat,
            columns=[
                'IDX',
                'START', 'STOP',
                'DEN', 'S', 'QC',
            ],
        )
        tract_df = tract_df.astype({
            'IDX': 'int',
            'START': 'int', 'STOP': 'int',
            'DEN': 'int', 'S': 'int',
            'QC': 'int',
        })
        tract_df.to_csv(hum_tract_path, index=False)
    return

# Conduct the QC.
pap_tracts_qc(pap=str(sys.argv[1]), chromosome=int(sys.argv[2]), prefix=str(sys.argv[3]))