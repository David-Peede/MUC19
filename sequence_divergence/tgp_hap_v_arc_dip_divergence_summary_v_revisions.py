# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import scipy.stats as stats
import sys
import zarr

### sys.argv[1] = window size ###


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

# Define a function to exctract a genotype matrix for a specific haplotype region.
def load_hap_region(prefix, chrom, hap_start, hap_end):
    # Load the callset for the chromosome and
    # extract the genotype matrix and positions array.
    chr_callset, chr_pos = load_callset_pos(prefix, chrom)
    hap_loc = chr_pos.locate_range(hap_start, hap_end)
    hap_idx = np.where(((hap_start <= chr_pos) & (chr_pos <= hap_end)))[0]
    hap_gt = allel.GenotypeArray(chr_callset[hap_loc])
    hap_pos = chr_pos[hap_idx]
    return hap_gt, hap_pos

# Define a function to load the QC'ed windows.
def load_windows(prefix, var_or_invar, window_size=72):
    # Load the data frame.
    wind_df = pd.read_csv(f'../windowing/{prefix}/{window_size}kb_nonoverlapping_{var_or_invar}_windows.csv.gz')
    return wind_df

# Define a function to load the window indicies that passed effective sequence length QC.
def load_esl_qc_windows_idx(prefix, window_size=72):
    # Load the QC'ed windows indicies.
    qc_winds_idx = np.loadtxt(
        f'../windowing/{prefix}/{window_size}kb_esl_qced_nonoverlapping_variant_windows.txt.gz',
        dtype=int,
    )
    return qc_winds_idx

# Define a function to load the per region effective sequence lengths.
def load_region_esl(prefix, window_size=72):
    # Load the effective sequence length matrix.
    region_esl = np.loadtxt(f'../windowing/{prefix}/{window_size}kb_eff_seq_len.txt.gz', dtype=int)
    return region_esl

# Define a function to calculate alternative allele frequencies.
def calc_alt_freqs(gt):
    # If there are no altenative alleles...
    if (gt.count_alleles().shape[1] == 1):
        # Calculate alternative allele frequencies.
        alt_freqs = gt.count_alleles().to_frequencies()[:, 0] - 1
    # Else...
    else:
        # Calculate alternative allele frequencies.
        alt_freqs = gt.count_alleles().to_frequencies()[:, 1]
    return alt_freqs

# Define a function to calculate alternative allele frequencies for a single individual.
def calc_ind_alt_freqs(gt):
    # Compute the frequency for each site.
    raw_freqs = np.nansum(gt, axis=2).flatten() / 2
    # Set missing data to Nan
    alt_freqs = np.where(raw_freqs == -1, np.nan, raw_freqs)
    return alt_freqs

# Define a function to compute pairwise differences.
def pwd_per_site(px, py):
    return ((px * (1 - py)) + (py * (1 - px)))

# Define a function to calculate sequence divergence between archaic diplotypes and modern human haplotypes.
def calc_hum_hap_v_arc_dip_diffs(gt):
    # Load the meta data.
    meta_data = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store the results.
    pwd_dicc = {'hap_1': np.array([]), 'hap_2': np.array([])}
    # Intialize the archaic's allele frequency.
    arc_freq = calc_ind_alt_freqs(gt.take([2347], axis=1))
    # For every human sample.
    for samp in range(meta_data.shape[0]):
        # Extract the sample's haplotypes.
        hap_1 = gt[:, samp, 0]
        hap_2 = gt[:, samp, 1]
        # Compute the distance from the archaic diplotype.
        pwd_1 = pwd_per_site(hap_1, arc_freq)
        pwd_2 = pwd_per_site(hap_2, arc_freq)
        # Update dictionaries.
        pwd_dicc['hap_1'] = np.append(pwd_dicc['hap_1'], np.nansum(pwd_1))
        pwd_dicc['hap_2'] = np.append(pwd_dicc['hap_2'], np.nansum(pwd_2))
    return pwd_dicc

# Define a function to load the windowed haplotype divergence results.
def load_hum_hap_v_arc_dip_pwd_windows(window_size):
    # Intialize a list of archaics.
    arc_list = ['DEN', 'ALT', 'CHA', 'VIN']
    # Intialize dictionaries to store the results.
    arc_dicc = {
        'hap_1': {arc: {'pwd': [], 'div': []} for arc in arc_list},
        'hap_2': {arc: {'pwd': [], 'div': []} for arc in arc_list},
    }
    # For every archaic.
    for arc in arc_list:
        # Load the effective sequence lengths.
        esl_df = load_windows(prefix=f'tgp_{arc.lower()}_masked_no_aa', var_or_invar='variant', window_size=window_size)
        # For every haplotype.
        for hap in ['hap_1', 'hap_2']:
            # For all chromosomes.
            for chrom in range(1, 23):
                # Load the results.
                pw_mat = np.loadtxt(
                    f'../muc19_results/tgp_{arc.lower()}_masked_no_aa/{arc.lower()}_{hap}_pw_diffs_chr{chrom}_{window_size}kb.txt.gz',
                )
                # Extract the effective sequence lengths.
                esl = esl_df[esl_df['CHR'] == chrom][arc].values
                # Append the results.
                arc_dicc[hap][arc]['pwd'].append(pw_mat)
                arc_dicc[hap][arc]['div'].append(pw_mat / esl[:, np.newaxis])
            # Concatenate all the windows.
            arc_dicc[hap][arc]['pwd'] = np.concatenate(arc_dicc[hap][arc]['pwd'], axis=0)
            arc_dicc[hap][arc]['div'] = np.concatenate(arc_dicc[hap][arc]['div'], axis=0)
    return arc_dicc

# Define a function to calculate the 95% CI around the mean.
def sem_ci_of_mean(data, confidence=0.95):
    # Convert the input data to a numpy array of floats.
    data = 1.0 * data
    # Mask the np.nan's.
    data = data[~np.isnan(data)]
    # Calculate the number of observations in the data.
    n = data.size
    # Compute the SEM of the data.
    sem = stats.sem(data)
    # Compute the confidence interval around the mean.
    ci = sem * stats.t.ppf((1 + confidence) / 2., n-1)
    return sem, ci

# Define a function to compile the haplotype divergence summary.
def compile_hum_hap_v_arc_dip_div_summary(window_size):
    # Intialize the left and right haplotype positions.
    if window_size == 72:
        left = 40759001
        right = 40831000
    elif window_size == 742:
        left = 40272001
        right = 41014000
    # Intialize a list of archaics.
    arc_list = ['DEN', 'ALT', 'CHA', 'VIN']
    # Intialize dictionaries to store the results.
    obs_dicc = {
        'hap_1': {arc: {} for arc in arc_list},
        'hap_2': {arc: {} for arc in arc_list},
    }
    # For every archaic.
    for arc in arc_list:
        # Load the genotypes for the focal region.
        obs_gt, _ = load_hap_region(f'tgp_{arc.lower()}_masked_no_aa', 12, left, right)
        # Load the effective sequence length.
        arc_esl = int(load_region_esl(prefix=f'tgp_{arc.lower()}_masked_no_aa', window_size=window_size))
        # Compute the number of pairwise differences.
        obs_pwd = calc_hum_hap_v_arc_dip_diffs(obs_gt)
        # For each haplotype.
        for hap in ['hap_1', 'hap_2']:
            # Update the dictionary.
            obs_dicc[hap][arc]['pwd'] = obs_pwd[hap]
            obs_dicc[hap][arc]['div'] = obs_pwd[hap] / arc_esl
    # Load the window results.
    wind_dicc = load_hum_hap_v_arc_dip_pwd_windows(window_size=window_size)
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the sample and population arrays.
    inds = meta_df['IND'].values
    spops = meta_df['SUPERPOP'].values
    pops = meta_df['POP'].values
    # Intialize archaic labels.
    arc_dicc = {
        'DEN': 'Denisovan', 'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.', 'VIN': 'Vindija Nean.',
    }
    # Intialize a dictionary to store the results.
    df_dicc = {
        'ind': [], 'spop': [],
        'pop': [], 'arc': [],
        'obs_pwd1': [], 'obs_pwd2': [],
        'obs_div1': [], 'obs_div2': [],
        'wind_m': [], 'wind_s': [],
        'wind_se': [], 'wind_ci': [],
        'wind_p1': [], 'wind_p2': [],
    }
    # For every archaic.
    for arc in arc_list:
        # Load the window indicies of comparable effective sequence length.
        wind_idx = load_esl_qc_windows_idx(prefix=f'tgp_{arc.lower()}_masked_no_aa', window_size=window_size)
        # For every individual.
        for ind_idx, ind_tup in enumerate(zip(inds, spops, pops)):
            # Unpack.
            ind, spop, pop, = ind_tup
            # Generate the individual's distribution.
            ind_dist = np.concatenate(
                (wind_dicc['hap_1'][arc]['div'][:, ind_idx][wind_idx], wind_dicc['hap_2'][arc]['div'][:, ind_idx][wind_idx]),
            )
            # Determine the mean, standard deviation, sem, and 95% CIs for the non-overlapping window distribution.
            wind_m = np.nanmean(ind_dist)
            wind_s = np.nanstd(ind_dist)
            wind_se, wind_ci = sem_ci_of_mean(ind_dist)
            # Compute the p-values.
            wind_p1 = (np.count_nonzero(obs_dicc['hap_1'][arc]['div'][ind_idx] >= ind_dist) / np.sum(~np.isnan(ind_dist)))
            wind_p2 = (np.count_nonzero(obs_dicc['hap_2'][arc]['div'][ind_idx] >= ind_dist) / np.sum(~np.isnan(ind_dist)))
            # Fill the dictionaries.
            df_dicc['ind'].append(ind)
            df_dicc['spop'].append(spop)
            df_dicc['pop'].append(pop)
            df_dicc['arc'].append(arc_dicc[arc])
            df_dicc['obs_pwd1'].append(obs_dicc['hap_1'][arc]['pwd'][ind_idx])
            df_dicc['obs_pwd2'].append(obs_dicc['hap_2'][arc]['pwd'][ind_idx])
            df_dicc['obs_div1'].append(obs_dicc['hap_1'][arc]['div'][ind_idx])
            df_dicc['obs_div2'].append(obs_dicc['hap_2'][arc]['div'][ind_idx])
            df_dicc['wind_m'].append(wind_m)
            df_dicc['wind_s'].append(wind_s)
            df_dicc['wind_se'].append(wind_se)
            df_dicc['wind_ci'].append(wind_ci)
            df_dicc['wind_p1'].append(wind_p1)
            df_dicc['wind_p2'].append(wind_p2)
    # Convert the dictionary to a dataframe.
    div_df = pd.DataFrame(df_dicc)
    # Rename all the columns to look pretty.
    div_df.rename(
        columns={
            'ind': 'Individual', 'spop': 'Super Population',
            'pop': 'Population', 'arc': 'Archaic',
            'obs_pwd1': f'Focal {window_size}kb Region (Pairwise Diffs. Hap. 1)',
            'obs_pwd2': f'Focal {window_size}kb Region (Pairwise Diffs. Hap. 2)',
            'obs_div1': f'Focal {window_size}kb Region (Seq. Div. Hap. 1)',
            'obs_div2': f'Focal {window_size}kb Region (Seq. Div. Hap. 2)',
            'wind_m': fr'{window_size}kb Nonoverlapping Windows $\left( \mu \right)$',
            'wind_s': fr'{window_size}kb Nonoverlapping Windows $\left( \sigma \right)$',
            'wind_se': fr'{window_size}kb Nonoverlapping Windows $\left( SEM \right)$',
            'wind_ci': fr'{window_size}kb Nonoverlapping Windows $\left( \pm CI_{{95\%}} \right)$',
            'wind_p1': r'$P-value$ (Hap. 1)',
            'wind_p2': r'$P-value$ (Hap. 2)',
        }, inplace=True,
    )
    # Export the dataframe as a csv.
    div_df.to_csv(f'../science_reviews_round_01_analyses/dataframes/tgp_haplotype_archaic_diplotype_divergence_{window_size}kb.csv.gz', index=False)
    return

# Compile the haplotype divergence summary.
compile_hum_hap_v_arc_dip_div_summary(window_size=int(sys.argv[1]))