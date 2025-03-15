# Import packages.
import numpy as np
import pandas as pd
import sys

### sys.argv[1] = papuan prefix ###
### sys.argv[2] = dataset prefix ###


# Define a function to consolidate the nonoverlapping windows.
def consolidate_windows(pap, prefix):
    # For each haplotype.
    for hap in ['hap1', 'hap2']:
        # Intialzie a list to store dataframes.
        df_list = []
        # For every chromsome.
        for chrom in range(1, 23):
            # Loadt the dataframe.
            df = pd.read_csv(f'./{prefix}/{pap}_{hap}_den_intro_tracts_summary_chr{chrom}.csv.gz')
            # Set the chromsome column.
            df['CHR'] = chrom
            # Append the data frame.
            df_list.append(df)
        # Concatenate the dataframes.
        wind_df = pd.concat(df_list)
        # Intialize a list of column orders.
        col_order = [
            'IDX', 'CHR', 'START', 'STOP',
            'DEN', 'S', 'QC',
        ]
        # Reorder the columns.
        wind_df = wind_df[col_order]
        # Define conditions.
        invar = (wind_df['S'] == 0)
        var = (wind_df['S'] > 0)
        qc = (wind_df['QC'] == 1)
        # Generate the invariant and variant dataframe.
        invar_df = wind_df[invar & qc]
        var_df = wind_df[var & qc]
        # Output the dataframes.
        invar_df.to_csv(f'./{prefix}/{pap}_{hap}_den_intro_invariant_tracts.csv.gz', index=False)
        var_df.to_csv(f'./{prefix}/{pap}_{hap}_den_intro_variant_tracts.csv.gz', index=False)
    return

# Consolidate windows.
consolidate_windows(pap=str(sys.argv[1]), prefix=str(sys.argv[2]))