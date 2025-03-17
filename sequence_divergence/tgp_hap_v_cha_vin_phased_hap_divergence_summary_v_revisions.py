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

# Define a function to compute pairwise differences.
def pwd_per_site(px, py):
    return ((px * (1 - py)) + (py * (1 - px)))

# Define a function to calculate sequence divergence between haplotypes.
def calc_hum_hap_v_cha_vin_hap_diffs(gt, arc_hap):
    # Load the meta data.
    meta_data = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store the results.
    pwd_dicc = {'hap_1': np.array([]), 'hap_2': np.array([])}
    # For every human sample.
    for samp in range(meta_data.shape[0]):
        # Extract the sample's haplotypes.
        hap_1 = gt[:, samp, 0]
        hap_2 = gt[:, samp, 1]
        # Compute the distance.
        pwd_1 = pwd_per_site(hap_1, arc_hap)
        pwd_2 = pwd_per_site(hap_2, arc_hap)
        # Update dictionaries.
        pwd_dicc['hap_1'] = np.append(pwd_dicc['hap_1'], np.nansum(pwd_1))
        pwd_dicc['hap_2'] = np.append(pwd_dicc['hap_2'], np.nansum(pwd_2))
    return pwd_dicc

# Define a function to load the windowed psuedo-haplotype divergence results.
def load_hum_hap_v_cha_vin_phap_pwd_windows(window_size):
    # Intialize a list of archaics.
    arc_list = ['CHA', 'VIN']
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
                    f'../muc19_results/tgp_{arc.lower()}_masked_no_aa/{arc.lower()}_{hap}_psuedo_hap_diffs_chr{chrom}_{window_size}kb.txt.gz',
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
def compile_hum_hap_v_cha_vin_hap_div_summary():
    # Intialize a list of archaic keys.
    arc_keys = [
        ('CHA', 'Chagyrskaya Nean. Hap. 1'), ('CHA', 'Chagyrskaya Nean. Hap. 2'), 
        ('VIN', 'Vindija Nean. Hap. 1'), ('VIN', 'Vindija Nean. Hap. 2'),
    ]
    # Load the phased haplotypes.
    phased_df = pd.read_csv('../meta_data/altai_nean_phased_late_neanderthals_denisovan_genos_72kb.csv.gz')
    # Intialize a dictionary to store the archaic haplotypes.
    arc_haps = {
        'Chagyrskaya Nean. Hap. 1': phased_df['Chagyrskaya Nean. Hap. 1'].values,
        'Chagyrskaya Nean. Hap. 2': phased_df['Chagyrskaya Nean. Hap. 2'].values,
        'Vindija Nean. Hap. 1': phased_df['Vindija Nean. Hap. 1'].values,
        'Vindija Nean. Hap. 2': phased_df['Vindija Nean. Hap. 2'].values,
    }
    # Load the effective sequence lengths.
    arc_esl = load_region_esl('tgp_arcs_masked_no_aa', 72)
    # Adjust the effective sequence lengths to account for the sites that could not be resolved.
    cha_esl = arc_esl[1] - 2
    vin_esl = arc_esl[2] - 3
    # Load the window indicies of comparable effective sequence length.
    cha_idx = load_esl_qc_windows_idx('tgp_cha_masked_no_aa', 72)
    vin_idx = load_esl_qc_windows_idx('tgp_vin_masked_no_aa', 72)
    # Intialize the archaic info.
    arc_info = {
        'CHA': {'esl': cha_esl, 'wind_idx': cha_idx},
        'VIN': {'esl': vin_esl, 'wind_idx': vin_idx},
    }
    # Load the window results.
    wind_dicc = load_hum_hap_v_cha_vin_phap_pwd_windows(72)
    # Intialize dictionaries to store the results.
    obs_dicc = {
        'hap_1': {hap_key: {} for hap_key in arc_haps},
        'hap_2': {hap_key: {} for hap_key in arc_haps},
    }
    # Load the genotypes for the focal region.
    tgp_arcs_no_aa_72kb_gt, _ = load_hap_region('tgp_arcs_masked_no_aa', 12, 40759001, 40831000)
    # For every archaic haplotype.
    for arc_key, hap_key in arc_keys:
        # Compute the number of pairwise differences.
        obs_pwd = calc_hum_hap_v_cha_vin_hap_diffs(tgp_arcs_no_aa_72kb_gt, arc_haps[hap_key])
        # For each haplotype.
        for hap in ['hap_1', 'hap_2']:
            # Update the dictionary.
            obs_dicc[hap][hap_key]['pwd'] = obs_pwd[hap]
            obs_dicc[hap][hap_key]['div'] = obs_pwd[hap] / arc_info[arc_key]['esl']
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the sample and population arrays.
    inds = meta_df['IND'].values
    spops = meta_df['SUPERPOP'].values
    pops = meta_df['POP'].values
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
    # For every archaic haplotype.
    for arc_key, hap_key in arc_keys:
        # Load the window indicies of comparable effective sequence length.
        wind_idx = arc_info[arc_key]['wind_idx']
        # For every individual.
        for ind_idx, (ind, spop, pop) in enumerate(zip(inds, spops, pops)):
            # Generate the individual's distribution.
            ind_dist = np.concatenate(
                (wind_dicc['hap_1'][arc_key]['div'][:, ind_idx][wind_idx], wind_dicc['hap_2'][arc_key]['div'][:, ind_idx][wind_idx]),
            )
            # Determine the mean, standard deviation, sem, and 95% CIs for the non-overlapping window distribution.
            wind_m = np.nanmean(ind_dist)
            wind_s = np.nanstd(ind_dist)
            wind_se, wind_ci = sem_ci_of_mean(ind_dist)
            # Compute the p-values.
            wind_p1 = (np.count_nonzero(obs_dicc['hap_1'][hap_key]['div'][ind_idx] >= ind_dist) / np.sum(~np.isnan(ind_dist)))
            wind_p2 = (np.count_nonzero(obs_dicc['hap_2'][hap_key]['div'][ind_idx] >= ind_dist) / np.sum(~np.isnan(ind_dist)))
            # Fill the dictionaries.
            df_dicc['ind'].append(ind)
            df_dicc['spop'].append(spop)
            df_dicc['pop'].append(pop)
            df_dicc['arc'].append(hap_key)
            df_dicc['obs_pwd1'].append(obs_dicc['hap_1'][hap_key]['pwd'][ind_idx])
            df_dicc['obs_pwd2'].append(obs_dicc['hap_2'][hap_key]['pwd'][ind_idx])
            df_dicc['obs_div1'].append(obs_dicc['hap_1'][hap_key]['div'][ind_idx])
            df_dicc['obs_div2'].append(obs_dicc['hap_2'][hap_key]['div'][ind_idx])
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
            'pop': 'Population', 'arc': 'Archaic Hap.',
            'obs_pwd1': 'Focal 72kb Region (Pairwise Diffs. Hap. 1)',
            'obs_pwd2': 'Focal 72kb Region (Pairwise Diffs. Hap. 2)',
            'obs_div1': 'Focal 72kb Region (Seq. Div. Hap. 1)',
            'obs_div2': 'Focal 72kb Region (Seq. Div. Hap. 2)',
            'wind_m': r'72kb Nonoverlapping Windows $\left( \mu \right)$',
            'wind_s': r'72kb Nonoverlapping Windows $\left( \sigma \right)$',
            'wind_se': r'72kb Nonoverlapping Windows $\left( SEM \right)$',
            'wind_ci': r'72kb Nonoverlapping Windows $\left( \pm CI_{{95\%}} \right)$',
            'wind_p1': r'$P-value$ (Hap. 1)',
            'wind_p2': r'$P-value$ (Hap. 2)',
        }, inplace=True,
    )
    # Export the dataframe as a csv.
    div_df.to_csv('../analyses_nbs/dataframes/tgp_haplotype_late_neanderthal_phased_haplotype_divergence_72kb.csv.gz', index=False)
    return

# Compile the psuedo-haplotype divergence summary.
compile_hum_hap_v_cha_vin_hap_div_summary()