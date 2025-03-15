# Import packages.
import numpy as np
import pandas as pd
import sys

### sys.argv[1] = window size ###
### sys.argv[2] = dataset prefix ###


# Define a function to consolidate the nonoverlapping windows.
def consolidate_windows(window_size, prefix):
    # Intialzie a list to store dataframes.
    df_list = []
    # For every chromsome.
    for chrom in range(1, 23):
        # Loadt the dataframe.
        df = pd.read_csv(f'./{prefix}/{window_size}kb_window_summary_chr{chrom}.csv.gz')
        # Set the chromsome column.
        df['CHR'] = chrom
        # Append the data frame.
        df_list.append(df)
    # Concatenate the dataframes.
    wind_df = pd.concat(df_list)
    # If we are consolidating the combined human and archaic data...
    if ('tgp_arcs' in prefix) | ('sgdp_arcs' in prefix):
        # Intialize a list of column orders.
        col_order = [
            'IDX', 'CHR', 'START', 'STOP',
            'ALT', 'CHA', 'VIN', 'DEN',
            'DEN-ALT', 'DEN-CHA', 'DEN-VIN', 'ALT-CHA',
            'ALT-VIN', 'CHA-VIN', 'HUM', 'S', 'QC',
        ]
    # Else-if we are consolidating the combined archaic data.
    elif 'arcs' in prefix:
        # Intialize a list of column orders.
        col_order = [
            'IDX', 'CHR', 'START', 'STOP',
            'ALT', 'CHA', 'VIN', 'DEN',
            'DEN-ALT', 'DEN-CHA', 'DEN-VIN', 'ALT-CHA',
            'ALT-VIN', 'CHA-VIN', 'S', 'QC',
        ]
    # Else-if we are consolidating the human data.
    elif (prefix == 'tgp_mod_no_aa') | (prefix == 'tgp_mod_aa') | (prefix == 'sgdp_no_aa') | (prefix == 'sgdp_aa'):
        # Intialize a list of column orders.
        col_order = [
            'IDX', 'CHR', 'START', 'STOP',
            'HUM', 'S', 'QC',
        ]
    # Else-if we are consolidating the human and single data.
    elif ('tgp' in prefix) | ('sgdp' in prefix):
        # Extract the archaic.
        archaic = prefix.split('_')[1]
        # Intialize a list of column orders.
        col_order = [
            'IDX', 'CHR', 'START', 'STOP',
            f'{archaic.upper()}', 'S', 'QC',
        ]
    # Else we are consolidating the archaic data.
    else:
        # Extract the archaic.
        archaic = prefix.split('_')[0]
        # Intialize a list of column orders.
        col_order = [
            'IDX', 'CHR', 'START', 'STOP',
            f'{archaic.upper()}', 'S', 'QC',
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
    invar_df.to_csv(f'./{prefix}/{window_size}kb_nonoverlapping_invariant_windows.csv.gz', index=False)
    var_df.to_csv(f'./{prefix}/{window_size}kb_nonoverlapping_variant_windows.csv.gz', index=False)
    return

# Consolidate windows.
consolidate_windows(window_size=int(sys.argv[1]), prefix=str(sys.argv[2]))