# Import packages.
import numpy as np
import pandas as pd
import sys

### sys.argv[1] = population ###


# Define a function to compute adjusted chromosome lengths.
def chr_seq_len(window_size):
    # Intialize the original assembly size for hg19.
    chr_dicc = {
        1: 249250621, 2: 243199373, 3: 198022430,
        4: 191154276, 5: 180915260, 6: 171115067,
        7: 159138663, 8: 146364022, 9: 141213431,
        10: 135534747, 11: 135006516, 12: 133851895,
        13: 115169878, 14: 107349540, 15: 102531392,
        16: 90354753, 17: 81195210, 18: 78077248,
        19: 59128983, 20: 63025520, 21: 48129895,
        22: 51304566,
    }
    # Initialize new empty dictionary to store the output.
    new_chr_dicc = {}
    # Iterate through every chromosome.
    for key in chr_dicc :
        # Floor divide every chromosome length by the window size and multiply
        # to find the new chromosome length.
        chr_len = chr_dicc[key]
        new_chr_len = (chr_len//window_size)*window_size
        # Refill dictionary with the new chromosome lengths.
        new_chr_dicc[key] = new_chr_len
    return new_chr_dicc

# Define a function to load the ihs results.
def load_ihs_tables(pop):
    # Intialize the snp-type dictionary.
    snp_type_dicc = {
        'DEN': 'denisovan_specific_snps',
        'NEA': 'neanderthal_specific_snps',
        'SHR': 'shared_archaic_snps',
        'ARC': 'archaic_specific_snps',
        'HOM': 'shared_hominin_snps',
    }
    # Intialize a list to store ihs values.
    ihs_list = []
    # For every chromosome.
    for chrom in range(1, 23):
        # Load the dataframe.
        df = pd.read_csv(
            f'../muc19_results/tgp_mod_aa/{pop.lower()}_chr{chrom}.ihs.out.100bins.norm', sep='\t',
            names=['ID', 'POS', 'DAF', 'IHH1', 'IHH0', 'U_IHS', 'N_IHS', 'CRIT'],
        )
        # Intialize the chromosome column.
        chr_col = np.full(df.shape[0], chrom)
        # Insert the new column.
        df.insert(loc=1, column='CHR', value=chr_col)
        # If this is an African population.
        if pop in ['YRI', 'LWK', 'GWD', 'MSL', 'ESN']:
            # For every snp type.
            for snp_type in snp_type_dicc:
                # Load the union of all non-AFR snp positions.
                arc_pos = np.loadtxt(f'../muc19_results/tgp_arcs_masked_no_aa/tgp_{snp_type_dicc[snp_type]}_chr{chrom}.txt.gz')
                # Create the snp mask.
                snp_mask = np.isin(df['POS'].values, arc_pos)
                # Update the dataframe.
                df[snp_type] = snp_mask
        # Else, this is one of the non-AFR populations.
        else:
            # For every snp type.
            for snp_type in snp_type_dicc:
                # Load the snp positions.
                arc_pos = np.loadtxt(f'../muc19_results/tgp_arcs_masked_no_aa/{pop.lower()}_{snp_type_dicc[snp_type]}_chr{chrom}.txt.gz')
                # Create the snp mask.
                snp_mask = np.isin(df['POS'].values, arc_pos)
                # Update the dataframe.
                df[snp_type] = snp_mask
        # Determine the human-specific snps.
        df['HUM'] = ~(df['ARC'].values | df['HOM'].values)
        # Append the list.
        ihs_list.append(df)
    # Concatenate the dataframes.
    ihs_df = pd.concat(ihs_list)
    # Extract the normalized iHS values.
    norm_ihs = ihs_df['N_IHS'].values
    # Generate the critical values mask.
    all_crit_mask = np.abs(norm_ihs) > 2
    # Update the dataframe.
    ihs_df['ALL_CRIT'] = all_crit_mask
    # For every snp partition.
    for snp_type in ['DEN', 'NEA', 'SHR', 'ARC', 'HOM', 'HUM']:
        # Extract the archaic mask.
        focal_snp_mask = ihs_df[snp_type].values
        # Update the dataframe.
        ihs_df[f'{snp_type}_CRIT'] = all_crit_mask & focal_snp_mask
    return ihs_df

# Define a function to identify iHS clusters.
def identify_ihs_clusters(ihs_df, window_size):
    # Intialize a list of snp types.
    snp_type_list = ['DEN', 'NEA', 'SHR', 'ARC', 'HOM', 'HUM']
    # Intialize a dictionary.
    ihs_clusters = {
        'CHR': [], 'START': [], 'STOP': [],
        'N_ALL_SNPS': [], 'N_DEN_SNPS': [],
        'N_NEA_SNPS': [], 'N_SHR_SNPS': [],
        'N_ARC_SNPS': [], 'N_HOM_SNPS': [], 'N_HUM_SNPS': [],
        'N_ALL_CRIT': [], 'PROP_ALL_CRIT': [],
        'N_DEN_CRIT': [], 'PROP_DEN_CRIT': [],
        'N_NEA_CRIT': [], 'PROP_NEA_CRIT': [],
        'N_SHR_CRIT': [], 'PROP_SHR_CRIT': [],
        'N_ARC_CRIT': [], 'PROP_ARC_CRIT': [],
        'N_HOM_CRIT': [], 'PROP_HOM_CRIT': [],
        'N_HUM_CRIT': [], 'PROP_HUM_CRIT': [],
    }
    # Calculate the adjusted chromosome length.
    adj_chr_len = chr_seq_len(window_size)
    # For every chromosome.
    for chrom in range(1, 23):
        # Subset the dataframe.
        chrom_df = ihs_df[ihs_df['CHR'] == chrom]
        # Extract the positions array and normalize iHS values.
        chrom_pos = chrom_df['POS'].values
        chrom_ihs = np.abs(chrom_df['N_IHS'].values)
        # Intialize a dcitionary for snp type masks.
        snp_masks = {
            snp_type: chrom_df[snp_type].values for snp_type in snp_type_list
        }
        # For every non-overlapping window.
        for window_start in range(1, int(adj_chr_len[chrom]), int(window_size)):
            # Determine the positions that fall within the window.
            window_mask = ((window_start <= chrom_pos) & (chrom_pos < (window_start + window_size)))
            # If there are no snps in this window.
            if window_mask.sum() == 0:
                # Update the dictionary.
                ihs_clusters['CHR'].append(chrom)
                ihs_clusters['START'].append(window_start)
                ihs_clusters['STOP'].append(window_start + window_size)
                ihs_clusters['N_ALL_SNPS'].append(0)
                ihs_clusters['N_ALL_CRIT'].append(0)
                ihs_clusters['PROP_ALL_CRIT'].append(np.nan)
                # For every archaic snp type.
                for snp_type in snp_type_list:
                    # Update the dictionary.
                    ihs_clusters[f'N_{snp_type}_SNPS'].append(0)
                    ihs_clusters[f'N_{snp_type}_CRIT'].append(0)
                    ihs_clusters[f'PROP_{snp_type}_CRIT'].append(np.nan)
            # Else there are snps in this window.
            else:
                # Update the dictionary.
                ihs_clusters['CHR'].append(chrom)
                ihs_clusters['START'].append(window_start)
                ihs_clusters['STOP'].append(window_start + window_size)
                ihs_clusters['N_ALL_SNPS'].append(window_mask.sum())
                ihs_clusters['N_ALL_CRIT'].append((chrom_ihs[window_mask] > 2).sum())
                ihs_clusters['PROP_ALL_CRIT'].append((chrom_ihs[window_mask] > 2).sum() / window_mask.sum())
                # For every snp type.
                for snp_type in snp_type_list:
                    # Intialize the mask.
                    snp_window_mask = snp_masks[snp_type] & window_mask
                    snp_window_ihs = chrom_ihs[snp_window_mask]
                    # Update the dictionary.
                    ihs_clusters[f'N_{snp_type}_SNPS'].append(snp_window_mask.sum())
                    ihs_clusters[f'N_{snp_type}_CRIT'].append((snp_window_ihs > 2).sum())
                    ihs_clusters[f'PROP_{snp_type}_CRIT'].append((snp_window_ihs > 2).sum() / snp_window_mask.sum())
    return pd.DataFrame(ihs_clusters)

# Load the command line argument.
focal_pop = str(sys.argv[1])
# Load the ihs results.
pop_ihs_df = load_ihs_tables(focal_pop)
# Generate the ihs clusters.
pop_ihs_clusters_742kb_df = identify_ihs_clusters(pop_ihs_df, 742_000)
pop_ihs_clusters_72kb_df = identify_ihs_clusters(pop_ihs_df, 72_000)
# Export the dataframes.
pop_ihs_df.to_csv(f'../muc19_results/tgp_mod_aa/{focal_pop.lower()}_ihs_genome_wide.csv.gz', index=False)
pop_ihs_clusters_742kb_df.to_csv(f'../muc19_results/tgp_mod_aa/{focal_pop.lower()}_ihs_windows_742kb.csv.gz', index=False)
pop_ihs_clusters_72kb_df.to_csv(f'../muc19_results/tgp_mod_aa/{focal_pop.lower()}_ihs_windows_72kb.csv.gz', index=False)