# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###
### sys.argv[2] = dataset prefix ###


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

# Define a function to qc hg19 genes windows.
def genes_qc(chromosome, prefix):
    # Intialize a dictionary to store the all sites path information.
    path_dicc = {
        'tgp_den_masked_no_aa': '../vcf_data/bookkeeping/tgp_den_masked_no_aa_calls_all_sites',
    }
    # Intitialize the all sites path, window path, and column variable.
    hum_all_sites_path = f'{path_dicc[prefix]}_chr{chromosome}.txt.gz'
    hum_gene_path = f'./{prefix}/ncbi_refseq_genes_summary_chr{chromosome}.csv.gz'
    # Extract the archaic.
    archaic = prefix.split('_')[1]
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos(prefix, chromosome)
    # Load the all site information.
    all_sites_mat = np.loadtxt(
        hum_all_sites_path,
        usecols=(0, 1),
        dtype=int,
    )
    # Load the hg19 gene information.
    hg19_gene_df = pd.read_csv(f'./hg19_genes/ncbi_refseq_genes_chr{chromosome}.csv.gz')
    # Extract the start and stops.
    starts = hg19_gene_df.START.values
    stops = hg19_gene_df.STOP.values
    # Extract the ids.
    genes = hg19_gene_df.GENE_ID.values
    transcripts = hg19_gene_df.TRANSCRIPT_ID.values
    # Intialize a results matrix to store the results.
    results_mat = np.empty((hg19_gene_df.shape[0], 6))
    # For every gene...
    for gene in range(hg19_gene_df.shape[0]):
        # Extract the start and stop positions.
        start = starts[gene]
        stop = stops[gene]
        # Extract the variable sites information.
        var_idx = np.where(((start <= all_pos) & (all_pos <= stop)))[0]
        # Extract the effective sequence length information.
        esl_idx = np.where(((start <= all_sites_mat[:, 1]) & (all_sites_mat[:, 1] <= stop)))[0]
        # Determine the effective sequence length.
        hum_arc_esl = esl_idx.size
        # If the window has variant sites to do calculations on...
        if (var_idx.size > 0):
            # Locate the gene.
            gene_loc = all_pos.locate_range(start, stop)
            # Determine the number of segregating sites.
            s = count_seg_sites(gt=allel.GenotypeArray(callset[gene_loc]))
            # Append the results.
            results_mat[gene, :] = np.array(
                [gene, start, stop, hum_arc_esl, s, 1], dtype=object,
            )
        # Else-if there are sites that passed qc for this window.
        elif (esl_idx.size > 0):
            # Append the results.
            results_mat[gene, :] = np.array(
                [gene, start, stop, hum_arc_esl, 0, 1], dtype=object,
            )
        # Else.
        else:
            # Append the results.
            results_mat[gene, :] = np.array(
                [gene, start, stop, 0, 0, 0], dtype=object,
            )
    # Export the the results matrix as a pandas dataframe.
    gene_df = pd.DataFrame(
        results_mat,
        columns=[
            'IDX',
            'START', 'STOP',
            f'{archaic.upper()}', 'S', 'QC',
        ],
    )
    gene_df = gene_df.astype({
        'IDX': 'int',
        'START': 'int', 'STOP': 'int',
        f'{archaic.upper()}': 'int', 'S': 'int',
        'QC': 'int',
    })
    # Insert the id columns.
    gene_df.insert(loc=1, column='GENE_ID', value=genes)
    gene_df.insert(loc=2, column='TRANSCRIPT_ID', value=transcripts)
    gene_df.to_csv(hum_gene_path, index=False)
    return

# Conduct the QC.
genes_qc(chromosome=int(sys.argv[1]), prefix=str(sys.argv[2]))