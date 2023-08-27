# Import packages.
import numpy as np
import pandas as pd
import sys

### sys.argv[1] = window size ###
### sys.argv[2] = tgp or arc ###


# Define a function to consolidate the nonoverlapping windows.
def consolidate_windows(window_size, tgp_or_arc):
    # Set the file path.
    path = f'./{tgp_or_arc}/{window_size}kb_window_summary_chr'
    # Intialzie a list to store dataframes.
    df_list = []
    # For every chromsome.
    for chrom in range(1, 23):
        # Loadt the dataframe.
        df = pd.read_csv(path+f'{chrom}.csv')
        # Set the chromsome column.
        df['CHR'] = chrom
        # Append the data frame.
        df_list.append(df)
    # Concatenate the dataframes.
    wind_df = pd.concat(df_list)
    # If we are consolidating the tgp data...
    if tgp_or_arc == 'tgp':
        # Intialize a list of column orders.
        col_order = [
            'IDX', 'CHR', 'START', 'STOP',
            'ALT', 'CHA', 'VIN', 'DEN',
            'DEN-ALT', 'DEN-CHA', 'DEN-VIN', 'ALT-CHA',
            'ALT-VIN', 'CHA-VIN', 'TGP', 'S', 'QC',
        ]
    # Else...
    else:
        # Intialize a list of column orders.
        col_order = [
            'IDX', 'CHR', 'START', 'STOP',
            'ALT', 'CHA', 'VIN', 'DEN',
            'DEN-ALT', 'DEN-CHA', 'DEN-VIN', 'ALT-CHA',
            'ALT-VIN', 'CHA-VIN', 'S', 'QC',
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
    invar_df.to_csv(f'./{tgp_or_arc}/{window_size}kb_nonoverlapping_invariant_windows.csv.gz', index=False)
    var_df.to_csv(f'./{tgp_or_arc}/{window_size}kb_nonoverlapping_variant_windows.csv.gz', index=False)
    return

# Consolidate windows.
consolidate_windows(window_size=int(sys.argv[1]), tgp_or_arc=str(sys.argv[2]))
