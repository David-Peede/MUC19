# Import packages.
import numpy as np
import pandas as pd
import sys

### sys.argv[1] = dataset prefix ###


# Define a function to consolidate the ncbi refseq genes.
def consolidate_genes(prefix):
    # Intialzie a list to store dataframes.
    df_list = []
    # For every chromsome.
    for chrom in range(1, 23):
        # Loadt the dataframe.
        df = pd.read_csv(f'./{prefix}/ncbi_refseq_genes_summary_chr{chrom}.csv.gz')
        # Set the chromsome column.
        df['CHR'] = chrom
        # Append the data frame.
        df_list.append(df)
    # Concatenate the dataframes.
    gene_df = pd.concat(df_list)
    # Extract the archaic
    archaic = prefix.split('_')[1]
    # Intialize a list of column orders.
    col_order = [
        'IDX', 'GENE_ID', 'TRANSCRIPT_ID', 'CHR', 'START', 'STOP',
        f'{archaic.upper()}', 'S', 'QC',
    ]
    # Reorder the columns.
    gene_df = gene_df[col_order]
    # Define conditions.
    invar = (gene_df['S'] == 0)
    var = (gene_df['S'] > 0)
    qc = (gene_df['QC'] == 1)
    # Generate the invariant and variant dataframe.
    invar_df = gene_df[invar & qc]
    var_df = gene_df[var & qc]
    # Output the dataframes.
    invar_df.to_csv(f'./{prefix}/ncbi_refseq_invariant_genes.csv.gz', index=False)
    var_df.to_csv(f'./{prefix}/ncbi_refseq_variant_genes.csv.gz', index=False)
    return

# Consolidate genes.
consolidate_genes(prefix=str(sys.argv[1]))