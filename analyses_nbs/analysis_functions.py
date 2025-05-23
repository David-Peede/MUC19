# Import packages.
from cycler import cycler
import allel
import copy
from collections import Counter
import itertools
from itertools import product
from functools import reduce
import gzip
import matplotlib
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LinearSegmentedColormap, hex2color
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
from matplotlib.patches import Patch, Ellipse, Circle
import matplotlib.lines as mlines
from matplotlib.ticker import FormatStrFormatter
import matplotlib.font_manager as fm
import numpy as np
import numcodecs
import pandas as pd
import re
import random
import zarr
from scipy.interpolate import interp1d
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.nonparametric.smoothers_lowess import lowess
import warnings

# Define the Okabe-Ito color palette
okabe_ito_palette = [
    '#E69F00',  # orange
    '#56B4E9',  # sky blue
    '#009E73',  # bluish green
    '#F0E442',  # yellow
    '#0072B2',  # blue
    '#D55E00',  # vermillion
    '#CC79A7',  # reddish purple
    '#000000',  # black
]

# Set the color cycler after updating other rcParams.
plt.rc('axes', prop_cycle=(cycler(color=okabe_ito_palette)))

# Ignore divide by 0/np.nan error and encode as np.nan's.
np.seterr(divide='ignore', invalid='ignore')
warnings.filterwarnings('ignore', category=UserWarning, message='.*loadtxt: Empty input file.*')

# Intialize my pandas preferences.
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)


#########################
### GENERAL FUNCTIONS ###
#########################

# Define a function to load genotyope and positions arrays.
def load_gt_pos(prefix, chrom):
    # Intialize the file path.
    path = f'../zarr_data/{prefix}_chr{chrom}.zarr'
    # Load the zarr array.
    zarr_array = zarr.open_group(path, mode='r')
    # Extract the genotype callset.
    callset = zarr_array[f'{chrom}/calldata/GT']
    # Convert the genotype callset to an array.
    gt = allel.GenotypeArray(callset)
    # Load the positions.
    pos = allel.SortedIndex(zarr_array[f'{chrom}/variants/POS'])
    return gt, pos

# Define a function to load the phased arachic genotyopes and positions arrays.
def load_phased_gt_pos(prefix):
    # Intialize the file path.
    path = f'../zarr_data/{prefix}.zarr'
    # Load the zarr array.
    zarr_array = zarr.open_group(path, mode='r')
    # Extract the genotype callset.
    callset = zarr_array['12/calldata/GT']
    # Convert the genotype callset to an array.
    gt = allel.GenotypeArray(callset)
    # Load the positions.
    pos = allel.SortedIndex(zarr_array['12/variants/POS'])
    return gt, pos

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
def load_windows(prefix, var_or_invar, window_size):
    # Load the data frame.
    wind_df = pd.read_csv(f'../windowing/{prefix}/{window_size}kb_nonoverlapping_{var_or_invar}_windows.csv.gz')
    return wind_df

# Define a function to load the window indicies that passed effective sequence length QC.
def load_esl_qc_windows_idx(prefix, window_size):
    # If the dataset is a single archaic dataset.
    if prefix.split('_')[0] in ['den', 'alt', 'cha', 'vin']:
        # Load the QC'ed windows indicies.
        qc_var_winds_idx = np.loadtxt(
            f'../windowing/{prefix}/{window_size}kb_esl_qced_nonoverlapping_variant_windows.txt.gz',
            dtype=int,
        )
        qc_invar_winds_idx = np.loadtxt(
            f'../windowing/{prefix}/{window_size}kb_esl_qced_nonoverlapping_invariant_windows.txt.gz',
            dtype=int,
        )
        return qc_var_winds_idx, qc_invar_winds_idx
    # Else this is any other dataset.
    else:
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

# Define a function to compile a region effective sequence lengths for the tgp + single archaic datasets.
def load_tgp_single_arc_region_esl(window_size):
    return {arc: int(load_region_esl(f'tgp_{arc.lower()}_masked_no_aa', window_size)) for arc in ['DEN', 'ALT', 'CHA', 'VIN']}

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



###########################
### ARCHAIC SNP DENSITY ###
###########################

# Define a function to count the number of archaic specific snps.
def count_arc_alleles(gt, pop_dicc, ooa_pop, ooa_freq):
    # Calculate alternative allele frequencies.
    afr_alt_freq = calc_alt_freqs(gt.take(pop_dicc['AFR'], axis=1))
    ooa_alt_freq = calc_alt_freqs(gt.take(pop_dicc[ooa_pop], axis=1))
    alt_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['ALT'], axis=1))
    cha_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['CHA'], axis=1))
    vin_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['VIN'], axis=1))
    den_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['DEN'], axis=1))
    # Intialize the conditions.
    c_ref_hum = ((1 - afr_alt_freq) < 0.01) & ((1 - ooa_alt_freq) > ooa_freq)
    c_ref_den = (den_alt_freq == 0)
    c_ref_not_den = (den_alt_freq != 0)
    c_ref_alt = (alt_alt_freq == 0)
    c_ref_not_alt = (alt_alt_freq != 0)
    c_ref_cha = (cha_alt_freq == 0)
    c_ref_not_cha = (cha_alt_freq != 0)
    c_ref_vin = (vin_alt_freq == 0)
    c_ref_not_vin = (vin_alt_freq != 0)
    c_ref_nea = c_ref_alt | c_ref_cha | c_ref_vin
    c_ref_not_nea = c_ref_not_alt & c_ref_not_cha & c_ref_not_vin
    c_ref_shared = c_ref_nea & c_ref_den
    c_alt_hum = (afr_alt_freq < 0.01) & (ooa_alt_freq > ooa_freq)
    c_alt_den = (den_alt_freq == 1)
    c_alt_not_den = (den_alt_freq != 1)
    c_alt_alt = (alt_alt_freq == 1)
    c_alt_not_alt = (alt_alt_freq != 1)
    c_alt_cha = (cha_alt_freq == 1)
    c_alt_not_cha = (cha_alt_freq != 1)
    c_alt_vin = (vin_alt_freq == 1)
    c_alt_not_vin = (vin_alt_freq != 1)
    c_alt_nea = c_alt_alt | c_alt_cha | c_alt_vin
    c_alt_not_nea = c_alt_not_alt & c_alt_not_cha & c_alt_not_vin
    c_alt_shared = c_alt_nea & c_alt_den
    # Determine the archaic specific masks.
    den_only_ref_mask = c_ref_hum & c_ref_not_nea & c_ref_den
    nea_only_ref_mask = c_ref_hum & c_ref_not_den & c_ref_nea
    shared_ref_mask = c_ref_hum & c_ref_shared
    arc_ref_mask = den_only_ref_mask | nea_only_ref_mask | shared_ref_mask
    den_only_alt_mask = c_alt_hum & c_alt_not_nea & c_alt_den
    nea_only_alt_mask = c_alt_hum & c_alt_not_den & c_alt_nea
    shared_alt_mask = c_alt_hum & c_alt_shared
    arc_alt_mask = den_only_alt_mask | nea_only_alt_mask | shared_alt_mask
    # Construct the results.
    arc_sites = np.array([
        (den_only_ref_mask | den_only_alt_mask).sum(),
        (nea_only_ref_mask | nea_only_alt_mask).sum(),
        (shared_ref_mask | shared_alt_mask).sum(),
        (arc_ref_mask | arc_alt_mask).sum(),
    ])
    return arc_sites

# Define a function to find the number of archaic snps.
def find_arc_alleles(gt, pop_dicc, ooa_pop, ooa_freq):
    # Calculate alternative allele frequencies.
    afr_alt_freq = calc_alt_freqs(gt.take(pop_dicc['AFR'], axis=1))
    ooa_alt_freq = calc_alt_freqs(gt.take(pop_dicc[ooa_pop], axis=1))
    alt_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['ALT'], axis=1))
    cha_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['CHA'], axis=1))
    vin_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['VIN'], axis=1))
    den_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['DEN'], axis=1))
    # Intialize the conditions.
    c_ref_hum = ((1 - afr_alt_freq) < 0.01) & ((1 - ooa_alt_freq) > ooa_freq)
    c_ref_den = (den_alt_freq == 0)
    c_ref_not_den = (den_alt_freq != 0)
    c_ref_alt = (alt_alt_freq == 0)
    c_ref_not_alt = (alt_alt_freq != 0)
    c_ref_cha = (cha_alt_freq == 0)
    c_ref_not_cha = (cha_alt_freq != 0)
    c_ref_vin = (vin_alt_freq == 0)
    c_ref_not_vin = (vin_alt_freq != 0)
    c_ref_nea = c_ref_alt | c_ref_cha | c_ref_vin
    c_ref_not_nea = c_ref_not_alt & c_ref_not_cha & c_ref_not_vin
    c_ref_shared = c_ref_nea & c_ref_den
    c_alt_hum = (afr_alt_freq < 0.01) & (ooa_alt_freq > ooa_freq)
    c_alt_den = (den_alt_freq == 1)
    c_alt_not_den = (den_alt_freq != 1)
    c_alt_alt = (alt_alt_freq == 1)
    c_alt_not_alt = (alt_alt_freq != 1)
    c_alt_cha = (cha_alt_freq == 1)
    c_alt_not_cha = (cha_alt_freq != 1)
    c_alt_vin = (vin_alt_freq == 1)
    c_alt_not_vin = (vin_alt_freq != 1)
    c_alt_nea = c_alt_alt | c_alt_cha | c_alt_vin
    c_alt_not_nea = c_alt_not_alt & c_alt_not_cha & c_alt_not_vin
    c_alt_shared = c_alt_nea & c_alt_den
    # Determine the archaic specific masks.
    den_only_ref_mask = c_ref_hum & c_ref_not_nea & c_ref_den
    nea_only_ref_mask = c_ref_hum & c_ref_not_den & c_ref_nea
    shared_ref_mask = c_ref_hum & c_ref_shared
    arc_ref_mask = den_only_ref_mask | nea_only_ref_mask | shared_ref_mask
    den_only_alt_mask = c_alt_hum & c_alt_not_nea & c_alt_den
    nea_only_alt_mask = c_alt_hum & c_alt_not_den & c_alt_nea
    shared_alt_mask = c_alt_hum & c_alt_shared
    arc_alt_mask = den_only_alt_mask | nea_only_alt_mask | shared_alt_mask
    # Construct a dictionary.
    arc_dicc = {
        'DEN': (den_only_ref_mask | den_only_alt_mask),
        'NEA': (nea_only_ref_mask | nea_only_alt_mask),
        'SHR': (shared_ref_mask | shared_alt_mask),
        'ARC': (arc_ref_mask | arc_alt_mask),
    }
    return arc_dicc

# Define a function to classify the snp types.
def tgp_classify_snps(gt, ooa_pop, ooa_freq):
     # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    pop_dicc = {
        'ALT': np.array([2347]), 'CHA': np.array([2348]),
        'VIN': np.array([2349]), 'DEN': np.array([2350]),
        'AFR': tgp_df[tgp_df['SUPERPOP'] == 'AFR'].index.values,
        'OOA': tgp_df[tgp_df['SUPERPOP'] != 'AFR'].index.values
    }
    # For every population...
    for pop in tgp_df['POP'].unique():
        # Append the dictionary with sample indicies.
        pop_dicc[pop] = tgp_df[tgp_df['POP'] == pop].index.values
    # Calculate alternative allele frequencies.
    afr_alt_freq = calc_alt_freqs(gt.take(pop_dicc['AFR'], axis=1))
    ooa_alt_freq = calc_alt_freqs(gt.take(pop_dicc[ooa_pop], axis=1))
    alt_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['ALT'], axis=1))
    cha_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['CHA'], axis=1))
    vin_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['VIN'], axis=1))
    den_alt_freq = calc_ind_alt_freqs(gt.take(pop_dicc['DEN'], axis=1))
    # Compute allele frequencies.
    arc_freqs = calc_alt_freqs(gt.take([2347, 2348, 2349, 2350], axis=1))
    tgp_freqs = calc_alt_freqs(gt.take(np.arange(0, 2347), axis=1))
    # Determine the sites that are invariant invariant.
    arc_invar_mask = (arc_freqs == 0) | np.isnan(arc_freqs)
    tgp_invar_mask = (tgp_freqs == 0) | np.isnan(tgp_freqs)
    # Intialize the conditions.
    c_ref_hum = ((1 - afr_alt_freq) < 0.01) & ((1 - ooa_alt_freq) > ooa_freq)
    c_ref_den = (den_alt_freq == 0)
    c_ref_not_den = (den_alt_freq != 0)
    c_ref_alt = (alt_alt_freq == 0)
    c_ref_not_alt = (alt_alt_freq != 0)
    c_ref_cha = (cha_alt_freq == 0)
    c_ref_not_cha = (cha_alt_freq != 0)
    c_ref_vin = (vin_alt_freq == 0)
    c_ref_not_vin = (vin_alt_freq != 0)
    c_ref_nea = c_ref_alt | c_ref_cha | c_ref_vin
    c_ref_not_nea = c_ref_not_alt & c_ref_not_cha & c_ref_not_vin
    c_ref_shared = c_ref_nea & c_ref_den
    c_alt_hum = (afr_alt_freq < 0.01) & (ooa_alt_freq > ooa_freq)
    c_alt_den = (den_alt_freq == 1)
    c_alt_not_den = (den_alt_freq != 1)
    c_alt_alt = (alt_alt_freq == 1)
    c_alt_not_alt = (alt_alt_freq != 1)
    c_alt_cha = (cha_alt_freq == 1)
    c_alt_not_cha = (cha_alt_freq != 1)
    c_alt_vin = (vin_alt_freq == 1)
    c_alt_not_vin = (vin_alt_freq != 1)
    c_alt_nea = c_alt_alt | c_alt_cha | c_alt_vin
    c_alt_not_nea = c_alt_not_alt & c_alt_not_cha & c_alt_not_vin
    c_alt_shared = c_alt_nea & c_alt_den
    # Determine the archaic specific masks.
    den_only_ref_mask = c_ref_hum & c_ref_not_nea & c_ref_den
    nea_only_ref_mask = c_ref_hum & c_ref_not_den & c_ref_nea
    shared_ref_mask = c_ref_hum & c_ref_shared
    arc_ref_mask = den_only_ref_mask | nea_only_ref_mask | shared_ref_mask
    den_only_alt_mask = c_alt_hum & c_alt_not_nea & c_alt_den
    nea_only_alt_mask = c_alt_hum & c_alt_not_den & c_alt_nea
    shared_alt_mask = c_alt_hum & c_alt_shared
    arc_alt_mask = den_only_alt_mask | nea_only_alt_mask | shared_alt_mask
    # Construct the archaic snps dictionary.
    arc_dicc = {
        'DEN': (den_only_ref_mask | den_only_alt_mask),
        'NEA': (nea_only_ref_mask | nea_only_alt_mask),
        'SHR': (shared_ref_mask | shared_alt_mask),
        'ARC': (arc_ref_mask | arc_alt_mask),
    }
    # Construct the human snp dictionary.
    hum_dicc = {
        'HUM': arc_invar_mask & ~tgp_invar_mask,
        'HOM': (~arc_invar_mask & ~arc_dicc['ARC']) | (arc_invar_mask & tgp_invar_mask),
        ## HOM accounts for the fact that some sites are encoded in the TGP as variant but are actually only,
        ## variant when you consider the related individuals who are in a seperate vcf.
    }
    # Construct the snp dictionary.
    snp_dicc = {**arc_dicc, **hum_dicc}
    return arc_dicc, hum_dicc, snp_dicc

# Define a function to compute the observed archaic snp density.
def tgp_arc_snp_density(gt):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    idx_dicc = {
        'ALT': np.array([2347]), 'CHA': np.array([2348]),
        'VIN': np.array([2349]), 'DEN': np.array([2350]),
        'AFR': tgp_meta_df[tgp_meta_df['SUPERPOP'] == 'AFR'].index.values,
    }
    # Intialize a list of focal populations.
    ooa_list = [
        'MXL', 'PEL', 'CLM', 'PUR', # AMR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
    ]
    # Intialize a dictionary to store the results.
    snp_dicc = {}
    # For every OOA population.
    for pop in ooa_list:
        # Update dictionaries.
        idx_dicc[pop] = tgp_meta_df[tgp_meta_df['POP'] == pop].index.values
        # Fill the dictionary.
        snp_dicc[pop] = count_arc_alleles(
            gt=gt, pop_dicc=idx_dicc,
            ooa_pop=pop, ooa_freq=0.01,
        )
    return snp_dicc

# Define a function to load the arc snp denisty windows.
def load_tgp_arc_snp_density_windows(window_size):
    # Intialize a focal population list.
    focal_pops = [
        'MXL', 'PEL', 'CLM', 'PUR', # AMR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
    ]
    # Intialize a dictionary to store the results.
    snp_dicc = {}
    # For ever focal population.
    for focal_pop in focal_pops:
        # Load the results.
        snp_mats = [
            np.loadtxt(
                f'../muc19_results/tgp_arcs_masked_no_aa/{focal_pop.lower()}_arc_snp_denisty_chr{chrom}_{window_size}kb.txt.gz',
            ) for chrom in range(1, 23)
        ]
        # Fill the dictionary.
        snp_dicc[focal_pop] = np.concatenate(snp_mats)
    return snp_dicc

# Define a function to compile the archaic snp density table.
def compile_tgp_arc_snp_denisty_summary(obs_dicc, wind_dicc, snp_type, window_size):
    # Intialize a dictionary of super populations.
    tgp_dicc = {
        'AMR': ['MXL', 'PEL', 'CLM', 'PUR'],
        'SAS': ['BEB', 'STU', 'ITU', 'PJL', 'GIH'],
        'EAS': ['CHB', 'KHV', 'CHS', 'JPT', 'CDX'],
        'EUR': ['TSI', 'CEU', 'IBS', 'GBR', 'FIN'],
    }
    # Intialize the snp-type dictionary.
    snp_type_dicc = {
        'DEN': {'lab': 'Denisovan-specific SNPs', 'idx': 0, 'export': 'denisovan_specific_snp'},
        'NEA': {'lab': 'Neanderthal-specific SNPs', 'idx': 1, 'export': 'neanderthal_specific_snp'},
        'SHR': {'lab': 'Shared Archaic SNPs', 'idx': 2, 'export': 'shared_archaic_snp'},
        'ARC': {'lab': 'All Archaic SNPs', 'idx': 3, 'export': 'archaic_specific_snp'},
    }
    # Intialize the p-value correction dictionary.
    pval_dicc = {
        742: {0: '<3.164e-04', 1: '>0.99968'},
        72: {0: '<3.389e-05', 1: '>0.9999966'},
    }
    # Intialize dictionaries to store the results.
    df_dicc = {
        'spop': [], 'pop': [], 'obs': [],
        'wind_m': [], 'wind_s': [],
        'wind_se': [], 'wind_ci': [],
        'wind_p': [],
    }
    # Load the good window indicies.
    wind_idx = load_esl_qc_windows_idx(prefix='tgp_arcs_masked_no_aa', window_size=window_size)
    # For every super population.
    for spop in tgp_dicc:
        # For every population.
        for focal_pop in tgp_dicc[spop]:
            # Grab the superpopulation results.
            winds = wind_dicc[focal_pop][:, snp_type_dicc[snp_type]['idx']][wind_idx]
            obs = obs_dicc[focal_pop][snp_type_dicc[snp_type]['idx']]
            # Determine the mean and standard deviation for non-overlapping windows.
            wind_m = np.nanmean(winds)
            wind_s = np.nanstd(winds)
            # Determine the sem and 95% ci.
            wind_se, wind_ci = sem_ci_of_mean(winds)
            # Compute the p-value.
            wind_p = (np.count_nonzero(obs <= winds) / np.sum(~np.isnan(winds)))
            # Update the dataframe dictionary.
            df_dicc['spop'].append(spop)
            df_dicc['pop'].append(focal_pop)
            df_dicc['obs'].append(int(obs))
            df_dicc['wind_m'].append(wind_m)
            df_dicc['wind_s'].append(wind_s)
            df_dicc['wind_se'].append(wind_se)
            df_dicc['wind_ci'].append(wind_ci)
            df_dicc['wind_p'].append(wind_p)
    # Convert the dictionaries to dataframes.
    snp_df = pd.DataFrame(df_dicc)
    # Adjust p-values of 0 and 1 by the number of permutations.
    snp_df['wind_p'] = np.where(snp_df['wind_p'] == 0, pval_dicc[window_size][0], snp_df['wind_p'])
    snp_df['wind_p'] = np.where(snp_df['wind_p'] == '1.0', pval_dicc[window_size][1], snp_df['wind_p'])
    # Rename all the columns to look pretty.
    snp_df.rename(
        columns={
            'spop': 'Super Population',
            'pop': 'Population',
            'obs': f'Focal {window_size}kb Region ({snp_type_dicc[snp_type]["lab"]})',
            'wind_m': fr'{window_size}kb Non-overlapping Windows $\left( \mu \right)$',
            'wind_s': fr'{window_size}kb Non-overlapping Windows $\left( \sigma \right)$',
            'wind_se': fr'{window_size}kb Non-overlapping Windows $\left( SEM \right)$',
            'wind_ci': fr'{window_size}kb Non-overlapping Windows $\left( \pm CI_{{95\%}} \right)$',
            'wind_p': fr'$P-value$',
        }, inplace=True
    )
    # Export the dataframe as a csv.
    snp_df.to_csv(f'./dataframes/tgp_{snp_type_dicc[snp_type]["export"]}_denisty_{window_size}kb.csv.gz', index=False)
    return snp_df

# Define a function to plot the archaic snp denisty windows for mxl.
def plot_mxl_arc_snp_denisty(snp_type, obs_742kb, winds_742kb, obs_72kb, winds_72kb):
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the snp-type dictionary.
    snp_type_dicc = {
        'DEN': {'lab': 'Denisovan-specific', 'idx': 0, 'export': 'denisovan_specific_snp'},
        'NEA': {'lab': 'Neanderthal-specific', 'idx': 1, 'export': 'neanderthal_specific_snp'},
    }
    # Intialize a bin dictionary.
    bin_dicc = {
        'DEN': {742: np.arange(0, 150, 10), 72: np.arange(0, 150, 10)},
        'NEA': {742: np.arange(0, 205, 5), 72: np.arange(0, 31, 1)},
    }
    # Load window indicies of, comparable effective sequence length.
    wind_idx_742kb = load_esl_qc_windows_idx(prefix='tgp_arcs_masked_no_aa', window_size=742)
    wind_idx_72kb = load_esl_qc_windows_idx(prefix='tgp_arcs_masked_no_aa', window_size=72)
    # Determine the significance threshold.
    sig_thresh = 0.05 
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 2, figsize=(7, 3.5), dpi=300,
        facecolor='white',
        sharex=False, sharey=False,
        constrained_layout=True,
    )
    # For each subplot.
    for i, tup in enumerate([
        (obs_742kb, winds_742kb, wind_idx_742kb, 742), (obs_72kb, winds_72kb, wind_idx_72kb, 72),
    ]):
        # Unpack.
        obs_dicc, wind_dicc, wind_idx, window_size = tup
        # Extract the observed value and windowed distribution.
        obs = obs_dicc['MXL'][snp_type_dicc[snp_type]['idx']]
        dist = wind_dicc['MXL'][:, snp_type_dicc[snp_type]['idx']][wind_idx]
        # Compute the p-value.
        pval = (np.count_nonzero(obs <= dist) / np.sum(~np.isnan(dist)))
        # If the p-value is significant.
        if pval < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (obs, '*', 25, 0.49)
        # Else, it is not significant.
        else:
            # Update the dictionary.
            sig_dicc[i] = (obs, r'$ns$', 10, 0.5225)
        # Plot the distribution.
        axes[i].hist(
            dist, bins=bin_dicc[snp_type][window_size],
            histtype='stepfilled', color='gray',
        )
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlabel(f'Number of {snp_type_dicc[snp_type]["lab"]}'+'\n'+f'SNPs in MXL per {window_size}kb Window')
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/mxl_{snp_type_dicc[snp_type]["export"]}_denisty_windows.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/mxl_{snp_type_dicc[snp_type]["export"]}_denisty_windows.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/mxl_{snp_type_dicc[snp_type]["export"]}_denisty_windows.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot!
    plt.show()
    return



###########
### iHS ###
###########

# Define a function to load the genome wide results.
def load_ihs_gw(pop):
    # Load the genome-wide results.
    gw_df = pd.read_csv(f'../muc19_results/tgp_mod_aa/{pop.lower()}_ihs_genome_wide.csv.gz')
    return gw_df

# Define a function load the ihs results for a chromosome.
def load_ihs_chrom(pop, chrom):
    # Load the genome-wide results.
    gw_df = load_ihs_gw(pop)
    return gw_df[gw_df['CHR'] == chrom].reset_index(drop=True)

# Define a function to compute the proportion of snps with cristical values.
def load_focal_region_ihs(window_size):
    # Intialize a dictionary to store populations.
    tgp_dicc = {
        'AMR': ['MXL', 'PEL', 'CLM', 'PUR'],
        'AFR': ['LWK', 'GWD', 'MSL', 'ESN', 'YRI'],
        'SAS': ['BEB', 'STU', 'ITU', 'PJL', 'GIH'],
        'EAS': ['CHB', 'KHV', 'CHS', 'JPT', 'CDX'],
        'EUR': ['TSI', 'CEU', 'IBS', 'GBR', 'FIN'],
    }
    # Intialize a regions dictionary.
    regions_dicc = {
        742: (40272001, 41014000),
        72: (40759001, 40831000),
    }
    # Unpack the start and end coordinates.
    start, end = regions_dicc[window_size]
    # Intialize a dictionary to store the results.
    ihs_dicc = {}
    # For every super population.
    for spop in tgp_dicc:
        # Intialize the subdictionary.
        ihs_dicc[spop] = {}
        # For every population.
        for pop in tgp_dicc[spop]:
            # Intialize a subdictionary.
            ihs_dicc[spop][pop] = {}
            # Load the chromosome 12 results.
            chr12_df = load_ihs_chrom(pop, 12)
            # Intialize masks.
            region_mask = (start <= chr12_df['POS'].values) & (chr12_df['POS'].values <= end)
            crit_mask = chr12_df['ALL_CRIT'].values
            # Determine the observed snps, observed critical snps, and the critical snp proportion.
            obs_snps = region_mask.sum()
            obs_crit = (crit_mask & region_mask).sum()
            obs_prop = obs_crit / obs_snps
            # Update the dictionary with the all snp results.
            ihs_dicc[spop][pop]['N_ALL_SNPS'] = obs_snps
            ihs_dicc[spop][pop]['N_ALL_CRIT'] = obs_crit
            ihs_dicc[spop][pop]['PROP_ALL_CRIT'] = obs_prop
            # For every snp partition.
            for snp_type in ['DEN', 'NEA', 'SHR', 'ARC', 'HOM', 'HUM']:
                # Grab the snp mask.
                snp_type_mask = chr12_df[snp_type].values
                # Create the archaic region mask.
                snp_type_region_mask = region_mask & snp_type_mask
                # Grab the absolute critical threshold mask.
                snp_type_crit_mask = chr12_df[f'{snp_type}_CRIT'].values
                # Determine the observed snps, observed critical snps, and the critical snp proportion.
                snp_type_snps = snp_type_region_mask.sum()
                snp_type_crit = (snp_type_crit_mask & snp_type_region_mask).sum()
                snp_type_prop = snp_type_crit / snp_type_snps
                # Update the dictionary with the snp set results.
                ihs_dicc[spop][pop][f'N_{snp_type}_SNPS'] = snp_type_snps
                ihs_dicc[spop][pop][f'N_{snp_type}_CRIT'] = snp_type_crit
                ihs_dicc[spop][pop][f'PROP_{snp_type}_CRIT'] = snp_type_prop
    return ihs_dicc

# Define a function to load the ihs window results.
def load_ihs_windows(window_size):
    # Intialize a dictionary to store populations.
    tgp_dicc = {
        'AMR': ['MXL', 'PEL', 'CLM', 'PUR'],
        'AFR': ['LWK', 'GWD', 'MSL', 'ESN', 'YRI'],
        'SAS': ['BEB', 'STU', 'ITU', 'PJL', 'GIH'],
        'EAS': ['CHB', 'KHV', 'CHS', 'JPT', 'CDX'],
        'EUR': ['TSI', 'CEU', 'IBS', 'GBR', 'FIN'],
    }
    # Intialize a dictionary to store the results.
    ihs_dicc = {}
    # For every super population.
    for spop in tgp_dicc:
        # Intialize the subdictionary.
        ihs_dicc[spop] = {}
        # For every population.
        for pop in tgp_dicc[spop]:
            # Update the dictionary.
            ihs_dicc[spop][pop] = pd.read_csv(f'../muc19_results/tgp_mod_aa/{pop.lower()}_ihs_windows_{window_size}kb.csv.gz')
    return ihs_dicc

# Define a function to summarize the iHS window results.
def compile_ihs_window_summary(obs_dicc, wind_dicc, window_size):
    # Intialize a dictionary to store populations.
    tgp_dicc = {
        'AMR': ['MXL', 'PEL', 'CLM', 'PUR'],
        'SAS': ['BEB', 'STU', 'ITU', 'PJL', 'GIH'],
        'EAS': ['CHB', 'KHV', 'CHS', 'JPT', 'CDX'],
        'EUR': ['TSI', 'CEU', 'IBS', 'GBR', 'FIN'],
        'AFR': ['LWK', 'GWD', 'MSL', 'ESN', 'YRI'],
    }
    # Intialize a snp set dictionary.
    snp_labs = {
        'ALL': 'All SNPs',
        'HUM': 'Human-specific SNPs',
        'HOM': 'Shared Hominin SNPs',
        'ARC': 'Archaic SNPs',
        'DEN': 'Denisovan-specific SNPs',
        'NEA': 'Neanderthal-specific SNPs',
        'SHR': 'Shared Archaic SNPs',
    }
    # Intialize a dictionary.
    ihs_window_summary = {
        'spop': [], 'pop': [], 'snp_type': [],
        'obs_snps': [], 'obs_crits': [], 'obs_prop': [], 'winds_prop_thresh': [],
        'obs_v_winds_denom': [], 'obs_v_winds_numer': [], 'obs_prop_percentile_rank': [],
        'mean_prop': [], 'std_prop': [], 'sem_prop': [], 'ci_prop': [],
        'obs_prop_is_out': [], 
    }
    # For every super population.
    for spop in tgp_dicc:
        # For every population.
        for pop in tgp_dicc[spop]:
            # For every snp partition.
            for snp_type in ['ALL', 'ARC', 'DEN', 'NEA']:
                # Grab the genomewide results.
                gw_snps = wind_dicc[spop][pop][f'N_{snp_type}_SNPS'].values
                gw_crits = wind_dicc[spop][pop][f'N_{snp_type}_CRIT'].values
                gw_props = wind_dicc[spop][pop][f'PROP_{snp_type}_CRIT'].values
                # Extract the observed values.
                obs_snps = obs_dicc[spop][pop][f'N_{snp_type}_SNPS']
                obs_crits = obs_dicc[spop][pop][f'N_{snp_type}_CRIT']
                obs_prop = obs_dicc[spop][pop][f'PROP_{snp_type}_CRIT']
                # Create the > 10 snp window threshold.
                winds_snp_mask = gw_snps > 10
                # If there are no windows that passed the snp threshold.
                if winds_snp_mask.sum() == 0:
                    # Update the dictionary.
                    ihs_window_summary['spop'].append(spop)
                    ihs_window_summary['pop'].append(pop)
                    ihs_window_summary['snp_type'].append(snp_labs[snp_type])
                    ihs_window_summary['obs_snps'].append(obs_snps)
                    ihs_window_summary['obs_crits'].append(obs_crits)
                    ihs_window_summary['obs_prop'].append(obs_prop)
                    ihs_window_summary['mean_prop'].append(np.nan)
                    ihs_window_summary['std_prop'].append(np.nan)
                    ihs_window_summary['sem_prop'].append(np.nan)
                    ihs_window_summary['ci_prop'].append(np.nan)
                    ihs_window_summary['winds_prop_thresh'].append(np.nan)
                    ihs_window_summary['obs_prop_percentile_rank'].append(np.nan)
                    ihs_window_summary['obs_v_winds_numer'].append(np.nan)
                    ihs_window_summary['obs_v_winds_denom'].append(np.nan)
                    ihs_window_summary['obs_prop_is_out'].append(False)
                # Else, there are windows to perform computations on.
                else:
                    # Apply the > 10 snp window threshold.
                    winds_snps = gw_snps[winds_snp_mask]
                    winds_crits = gw_crits[winds_snp_mask]
                    winds_props = gw_props[winds_snp_mask]
                    # Determine the mean, standard deviation, SEM and 95% CIs for the window distribution.
                    mean_prop, std_prop = np.nanmean(winds_props), np.nanstd(winds_props)
                    sem_prop, ci_prop = sem_ci_of_mean(winds_props)
                    # Determine the 1% of windows threshold.
                    winds_prop_thresh = np.percentile(winds_props, 99)
                    # Compute the percentile rank of the observed proportion.
                    obs_prop_pr = stats.percentileofscore(winds_props, obs_prop, kind='strict')
                    # Update the dictionary.
                    ihs_window_summary['spop'].append(spop)
                    ihs_window_summary['pop'].append(pop)
                    ihs_window_summary['snp_type'].append(snp_labs[snp_type])
                    ihs_window_summary['obs_snps'].append(obs_snps)
                    ihs_window_summary['obs_crits'].append(obs_crits)
                    ihs_window_summary['obs_prop'].append(obs_prop)
                    ihs_window_summary['mean_prop'].append(mean_prop if (obs_snps > 10) else np.nan)
                    ihs_window_summary['std_prop'].append(std_prop if (obs_snps > 10) else np.nan)
                    ihs_window_summary['sem_prop'].append(sem_prop if (obs_snps > 10) else np.nan)
                    ihs_window_summary['ci_prop'].append(ci_prop if (obs_snps > 10) else np.nan)
                    ihs_window_summary['winds_prop_thresh'].append(winds_prop_thresh if (obs_snps > 10) else np.nan)
                    ihs_window_summary['obs_prop_percentile_rank'].append(obs_prop_pr if (obs_snps > 10) else np.nan)
                    ihs_window_summary['obs_v_winds_numer'].append((obs_prop > winds_props).sum() if (obs_snps > 10) else np.nan)
                    ihs_window_summary['obs_v_winds_denom'].append(winds_props.size if (obs_snps > 10) else np.nan)
                    ihs_window_summary['obs_prop_is_out'].append((obs_prop > winds_prop_thresh) & (obs_snps > 10))
    # Convert to a dataframe.
    ihs_window_summary_df = pd.DataFrame(ihs_window_summary)
    # Rename all the columns to look pretty.
    ihs_window_summary_df.rename(
        columns={
            'spop': 'Super Population', 'pop': 'Population',
            'snp_type': 'SNP Set',
            'obs_snps': f'Focal {window_size}kb Region (Total SNPs)',
            'obs_crits': fr'Focal {window_size}kb Region (SNPs with $\mid iHS \mid > 2$)',
            'obs_prop': fr'Focal {window_size}kb Region (Prop. of SNPs with $\mid iHS \mid > 2$)',
            'obs_prop_percentile_rank': fr'Focal {window_size}kb Region (Percentile Rank)',
            'obs_v_winds_numer': fr'Focal {window_size}kb Region $>$ {window_size}kb Non-overlapping Windows',
            'obs_v_winds_denom': fr'{window_size}kb Non-overlapping Windows (Total SNPs $> 10$)',
            'mean_prop': fr'{window_size}kb Non-overlapping Windows $\left( \mu \right)$',
            'std_prop': fr'{window_size}kb Non-overlapping Windows $\left( \sigma \right)$',
            'sem_prop': fr'{window_size}kb Non-overlapping Windows $\left( SEM \right)$',
            'ci_prop': fr'{window_size}kb Non-overlapping Windows $\left( \pm CI_{{95\%}} \right)$',
            'winds_prop_thresh': fr'{window_size}kb Non-overlapping Windows ($99^{{th}}$ Percentile)',
            'obs_prop_is_out': f'Focal {window_size}kb Region is in Top 1%?'
        }, inplace=True
    )
    # Export the dataframe as a csv.
    ihs_window_summary_df.to_csv(f'./dataframes/tgp_ihs_critical_scores_proportions_{window_size}kb.csv.gz', index=False)
    # Subset to show the results reported in the paper.
    return ihs_window_summary_df[np.isin(ihs_window_summary_df['SNP Set'].values, ['All SNPs', 'Archaic SNPs'])].reset_index(drop=True)

# Define a function to plot the ihs results for a region of chromosome 12.
def plot_chr12_ihs_region(df, region_or_window):
    # Intialize a dictionary for plotting.
    plot_dicc = {
        'region': {
            'x_lab': 'Longest MXL Introgressed Tract (Chr12: 40272001 - 41014000)',
            'gene_dicc': {'SLC2A13': {'lab': r'$SLC2A13$'}, 'LRRK2': {'lab': r'$LRRK2$'}, 'MUC19': {'lab': r'$MUC19$'}},
            'export': 'mxl_ihs_focal_region_742kb',
        },
    }
    # Load and intialize the hg19 gene information.
    hg19_gene_df = pd.read_csv(f'../annotations/hg19_genes/ncbi_refseq_genes_chr12.csv.gz')
    gene_dicc = plot_dicc[region_or_window]['gene_dicc']
    # For every gene.
    for gene in gene_dicc:
        # Fill the dictionary.
        gene_dicc[gene]['start'] = hg19_gene_df[hg19_gene_df['GENE_ID'] == gene].START.values[0]
        gene_dicc[gene]['stop'] = hg19_gene_df[hg19_gene_df['GENE_ID'] == gene].STOP.values[0]
    # Intialize a dictionary of snp type partitions.
    snp_dicc = {
        'HUM': {'l': 'Human-specific', 'c': 'gray', 'm': df['HUM'].values},
        'HOM': {'l': 'Shared Hominin', 'c': '#000000', 'm': df['HOM'].values},
        'DEN': {'l': 'Denisovan-specific', 'c': '#E69F00', 'm': df['DEN'].values},
        'NEA': {'l': 'Neanderthal-specific', 'c': '#56B4E9', 'm': df['NEA'].values},
        'SHR': {'l': 'Shared Archaic', 'c': '#CC79A7', 'm': df['SHR'].values},
    }
    # Intialize the legend.
    legend_handles = [
        Line2D(
            [0], [0], linestyle='none', marker='o', markersize=7.5, markeredgewidth=0.1,
            color=snp_dicc[snp_type]['c'], markeredgecolor='white', label=f'{snp_dicc[snp_type]["l"]} (n = {snp_dicc[snp_type]["m"].sum()})'
        ) for snp_type in ['DEN', 'NEA', 'SHR', 'HUM', 'HOM']
    ]
    legend_handles.append(Line2D([0], [0], color='black', linestyle='dashed', label=r'|$iHS$| = 2'))
    # Extract the normalized iHS values and positions.
    norm_ihs = np.abs(df['N_IHS'].values)
    pos = df['POS'].values
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the figure and axes.
    fig = plt.figure(figsize=(7, 3.5), dpi=300, facecolor='white', constrained_layout=True)
    ax = fig.add_subplot(111)
    # Plot the 72kb haplotype region.
    ax.axvspan(max(40759001, pos[0]), min(40831000, pos[-1]), alpha=0.05, facecolor='black')
    # For all snp types.
    for snp_type in snp_dicc:
        # Plot the results.
        ax.scatter(
            pos[snp_dicc[snp_type]['m']], norm_ihs[snp_dicc[snp_type]['m']], zorder=5,
            color=snp_dicc[snp_type]['c'], marker='o', edgecolor='white', linewidth=0.1, alpha=0.75, s=20,
        )
    # Plot the threshold.
    ax.axhline(2, 0, 1, color='black', linestyle='dashed', lw=1)
    # For every gene.
    for gene in gene_dicc:
        # Plot the gene segments within the region.
        ax.plot(
            [max(gene_dicc[gene]['start'], pos[0]), min(gene_dicc[gene]['stop'], pos[-1])],
            [-0.5, -0.5], color='black', marker='|', ms=10, lw=1,
        )
        # Annotate the gene.
        ax.text(
            ((max(gene_dicc[gene]['start'], pos[0]) + min(gene_dicc[gene]['stop'], pos[-1])) / 2),
            -0.5, gene_dicc[gene]['lab'], fontsize=10,
            horizontalalignment='center',verticalalignment='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black'),
        )
    # Adjust the x-axis ticks and labels.
    if region_or_window == 'region':
        x_ticks = np.arange(40_250_000, 41_100_000, 100_000)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([f'{round(x_tick / 1e6, 3)} Mb' for x_tick in x_ticks])
    else:
        ax.set_xticks(ax.get_xticks())
        ax.set_xticklabels([f'{round(x_tick / 1e6, 3)} Mb' for x_tick in ax.get_xticks()])
    # Set the y-axis label.
    ax.set_ylim((-1, None))
    # Get the current y-ticks.
    yticks = ax.get_yticks()
    # Filter y-ticks.
    filtered_yticks = [tick for tick in yticks if tick >= 0]
    # Set the y-ticks to the filtered y-ticks.
    ax.set_yticks(filtered_yticks)
    # Set the y-axis label
    ax.set_xlabel(plot_dicc[region_or_window]['x_lab'])
    ax.set_ylabel(r'|$iHS$|')
    # Add a legend.
    ax.legend(
        handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, 1.1),
        ncol=3, frameon=True, fancybox=True, shadow=True, fontsize=8,
    )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/{plot_dicc[region_or_window]["export"]}.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/{plot_dicc[region_or_window]["export"]}.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/{plot_dicc[region_or_window]["export"]}.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot!
    plt.show()
    return




###########################
### INTROGRESSED TRACTS ###
###########################

# Define a function to find tracts with any overlap of a region.
def tract_any_overlap(tracts_df, start, end):
    # Intialize a set to store tract indicies.
    tract_idx = set()
    # For every row in the qc'ed tracts.
    for _, row in tracts_df.iterrows():
        # If this tract has any overlap with the gene.
        if ((row['START'] <= start < row['END']) or         # t_s---s---t_e---e
            (row['START'] <= end < row['END']) or           # s---t_s---e---t_e
            (start >= row['START'] and end < row['END']) or # t_s---s---e---t_e
            (start <= row['START'] and end > row['END'])):  # s---t_s---t_e---e:
            # Update the set.
            tract_idx.add(row.name)
    return tracts_df.iloc[np.sort(np.array(list(tract_idx))), :].reset_index(drop=True)

# Define a function to find genes with any overlap of a region.
def genes_any_overlap(genes_df, start, end):
    # Intialize a set to store genes indicies.
    genes_idx = set()
    # For every row in the qc'ed genes.
    for _, row in genes_df.iterrows():
        # If this gene has any overlap with the region.
        if ((row['START'] <= start < row['STOP']) or         # g_s---s---g_e---e
            (row['START'] <= end < row['STOP']) or           # s---g_s---e---g_e
            (start >= row['START'] and end < row['STOP']) or # g_s---s---e---g_e
            (start <= row['START'] and end > row['STOP'])):  # s---g_s---g_e---e:
            # Update the set.
            genes_idx.add(row.name)
    return genes_df.iloc[np.sort(np.array(list(genes_idx))), :].reset_index(drop=True)

# Define a function to summarize the introgressed tract frequency.
def compile_tgp_introgressed_tract_freqs(tracts_df, midfix):
    # Load the meta data file for the TGP.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    tgp_info = {key: np.array(value) for key, value in tgp_df.to_dict(orient='list').items()}
    # Intialize an ordered dictionary.
    tgp_dicc = {
        'AMR': ['MXL', 'PEL', 'CLM', 'PUR'],
        'SAS': ['BEB', 'STU', 'ITU', 'PJL', 'GIH'],
        'EAS': ['CHB', 'KHV', 'CHS', 'JPT', 'CDX'],
        'EUR': ['TSI', 'CEU', 'IBS', 'GBR', 'FIN'],
        'AFR': ['LWK', 'GWD', 'MSL', 'ESN', 'YRI'],
    }
    # Intialize a dictionaries.
    pop_dicc = {
        'spop': [], 'pop': [], 'n_tot': [],
        'n_tracts': [],  'freq': [], 'len': []
    }
    spop_dicc = {
        'spop': [], 'n_tot': [],
        'n_tracts': [],  'freq': [], 'len': []
    }
    # For every super population.
    for spop in tgp_dicc:
        # Determine the total number of chromosomes.
        spop_tot = (tgp_info['SUPERPOP'] == spop).sum() * 2
        # Determine the number of tracts and mean length.
        spop_tracts = (tracts_df['SUPERPOP'].values == spop).sum()
        spop_avg_len = tracts_df[tracts_df['SUPERPOP'].values == spop]['LENGTH'].mean()
        # Update the dictionary.
        spop_dicc['spop'].append(spop)
        spop_dicc['n_tot'].append(spop_tot)
        spop_dicc['n_tracts'].append(spop_tracts)
        spop_dicc['freq'].append(spop_tracts / spop_tot)
        spop_dicc['len'].append(spop_avg_len)
        # For every population.
        for pop in tgp_dicc[spop]:
            # Determine the total number of chromosomes.
            pop_tot = (tgp_info['POP'] == pop).sum() * 2
            # Determine the number of tracts and mean length.
            pop_tracts = (tracts_df['POP'].values == pop).sum()
            pop_avg_len = tracts_df[tracts_df['POP'].values == pop]['LENGTH'].mean()
            # Update the dictionary.
            pop_dicc['spop'].append(spop)
            pop_dicc['pop'].append(pop)
            pop_dicc['n_tot'].append(pop_tot)
            pop_dicc['n_tracts'].append(pop_tracts)
            pop_dicc['freq'].append(pop_tracts / pop_tot)
            pop_dicc['len'].append(pop_avg_len)
    # Convert to dataframes.
    pop_df = pd.DataFrame(pop_dicc)
    spop_df = pd.DataFrame(spop_dicc)
    # Cleanup the column names.
    spop_df.rename(
        columns={
            'spop': 'Super Population',
            'n_tot': 'Total Number of Chromosomes',
            'n_tracts': 'Number of Introgressed Tracts',
            'freq': 'Introgressed Tract Frequency',
            'len': 'Mean Tract Length',
        }, inplace=True
    )
    pop_df.rename(
        columns={
            'spop': 'Super Population', 'pop': 'Population',
            'n_tot': 'Total Number of Chromosomes',
            'n_tracts': 'Number of Introgressed Tracts',
            'freq': 'Introgressed Tract Frequency',
            'len': 'Mean Tract Length',
        }, inplace=True
    )
    # Export.
    spop_df.to_csv(f'./dataframes/tgp_{midfix}_introgressed_tract_frequency_per_super_population.csv.gz', index=False)
    pop_df.to_csv(f'./dataframes/tgp_{midfix}_introgressed_tract_frequency_per_population.csv.gz', index=False)
    # Determine the total number of chromosomes.
    amr_chroms = spop_df[
        spop_df['Super Population'] == 'AMR'
    ]['Total Number of Chromosomes'].values[0]
    non_amr_chroms = spop_df[
        (spop_df['Super Population'] != 'AMR') & (spop_df['Super Population'] != 'AFR')
    ]['Total Number of Chromosomes'].sum()
    amr_v_non_amr_chroms = np.array([amr_chroms, non_amr_chroms])
    # Determine the total number of tracts.
    amr_tracts = spop_df[
        spop_df['Super Population'] == 'AMR'
    ]['Number of Introgressed Tracts'].values[0]
    non_amr_tracts = spop_df[
        (spop_df['Super Population'] != 'AMR') & (spop_df['Super Population'] != 'AFR')
    ]['Number of Introgressed Tracts'].sum()
    amr_v_non_amr_tracts = np.array([amr_tracts, non_amr_tracts])
    # Perform a z-proprtions test.
    z_stat, z_p = proportions_ztest(
        amr_v_non_amr_tracts, amr_v_non_amr_chroms, alternative='larger',
    )
    # Build contigency tables.
    fet_table = [
        [amr_tracts, amr_chroms - amr_tracts],
        [non_amr_tracts, non_amr_chroms - non_amr_tracts],
    ]
    # Run a FET.
    odds_ratio, f_p = stats.fisher_exact(fet_table, alternative='greater')
    # Print the Results.
    amr_v_non_amr_freqs = amr_v_non_amr_tracts / amr_v_non_amr_chroms
    tract_summary = f"""
    AMR vs Non-AMR Introgressed Tracts Summary
    ===========================================


    Proportions Z-Test
    ------------------
    Population   Chromosomes   Tracts   Frequency
    ----------   -----------   ------   ---------
    AMR          {amr_chroms}           {amr_tracts}      {amr_v_non_amr_freqs[0]}
    Non-AMR      {non_amr_chroms}          {non_amr_tracts}      {amr_v_non_amr_freqs[1]}

    Z-statistic: {z_stat}
    P-value:     {z_p}


    Fisher's Exact Test
    -------------------
    Contingency Table:
                        Introgressed   Non-Introgressed
                        Tracts         Tracts
                        -----------   ---------------
                   AMR| {fet_table[0][0]}           {fet_table[0][1]}
               Non-AMR| {fet_table[1][0]}           {fet_table[1][1]}

    Odds Ratio: {odds_ratio}
    P-value:    {f_p}
    """
    print(tract_summary)
    return spop_df, pop_df

# Define a function to find longest contiguous segment conditioned on the minimum numbers of tracts.
def find_longest_contiguous_segments(positions, tract_counts, min_tract_count):
    # Determine the tract_counts that meet the minimum tract count.
    valid_tract_counts = tract_counts >= min_tract_count
    # Edge case for when no tract_counts meet the minimum tract count.
    if not np.any(valid_tract_counts):
        return np.nan, np.nan, np.nan
    
    # Find contiguous segements of true values.
    contiguous_segments = np.split(                   # Split the contiguous segements into subarrays.
        np.where(valid_tract_counts)[0],              # Finds the indicies of all true values.
        np.where(                                     # Finds where the difference isn't 1, ie not contiguos.
            np.diff(                                  # Computes the difference between adjacent indicies.
                np.where(valid_tract_counts)[0]) != 1 # Finds where the difference isn't 1, ie not contiguos.
        )[0] + 1)                                     # Point to the start of the next contiguos run, ie + 1.
    # Find the longest contiguous segment.
    longest_contiguous_segment = max(contiguous_segments, key=len)
    # Determine the coordinates and segement length.
    start = positions[longest_contiguous_segment[0]]
    end = positions[longest_contiguous_segment[-1]]
    segment_length = end - start + 1          # + 1 since both coordinates are inclusive.
    return start, end, segment_length

# Define a function to compute the probability of ILS.
def calc_prob_ils(r, m):
    # Intialize parameters.
    b = 550_000 # Archaic - Modern human split time from Prüfer et al., 2017; Peyrégne et al., 2022.
    g = 29 # Human generation time from Langergraber et al., 2012; Zeberg & Pääbo, 2021.
    L = 1 / (r * (b / g))
    # Compute the probability of ILS.
    pr_ils = 1 - stats.gamma.cdf(m, a=2, scale=(1 / L))
    # If the probability reported by scipy is 0.
    if pr_ils == 0:
        # Determine the smallest value scipy would report.
        pr_ils = np.nextafter(pr_ils, 1)
    # Print a summary.
    print(f'E[L | r] = {L}')
    print(f'Pr(ILS | m) = {pr_ils}')
    return

# Define a function load and order introgressed tracts by lengths for a given region.
def load_region_tracts(tracts_df, pop_list):
    # Sort the dataframe.
    sorted_tracts_df = tracts_df.sort_values(by=['LENGTH'], ascending=False, ignore_index=True)
    # Extract the population column.
    pop_tracts = sorted_tracts_df['POP'].values
    # Determine the population masks.
    pop_df_mask = np.in1d(pop_tracts, pop_list)
    # Subset the dataframe.
    pop_tracts_df = sorted_tracts_df[pop_df_mask].reset_index(drop=True)
    # Extract the individuals.
    pop_inds = pop_tracts_df['IND'].values
    # Determine the unique individual counts.
    pop_u_inds, pop_ind_counts = np.unique(pop_inds, return_counts=True)
    # Determine the inividuals with multiple and single tracts.
    pop_multi_tracts_inds = pop_u_inds[(pop_ind_counts > 1)]
    pop_single_tracts_inds = pop_u_inds[~(pop_ind_counts > 1)]
    # Intialize a data frame dictionary.
    pop_df_dicc = {
        'IND': [], 'SUPERPOP': [], 'POP': [],
        'START': [], 'END': [],
    }
    # For every unique individual...
    for ind in pop_u_inds:
        # Subset the dataframe.
        subset_df = pop_tracts_df[pop_tracts_df['IND'] == ind]
        # If the individual only has one tract...
        if ind in pop_single_tracts_inds:
            # Fill the dictionary.
            pop_df_dicc['IND'].append(ind)
            pop_df_dicc['SUPERPOP'].append(subset_df['SUPERPOP'].values[0])
            pop_df_dicc['POP'].append(subset_df['POP'].values[0])
            pop_df_dicc['START'].append(subset_df['START'].values[0])
            pop_df_dicc['END'].append(subset_df['END'].values[0])
        # Else...
        else:
            # Fill the dictionary.
            pop_df_dicc['IND'].append(ind)
            pop_df_dicc['SUPERPOP'].append(subset_df['SUPERPOP'].values[0])
            pop_df_dicc['POP'].append(subset_df['POP'].values[0])
            pop_df_dicc['START'].append(np.min(subset_df['START'].values))
            pop_df_dicc['END'].append(np.max(subset_df['END'].values))
    # Convert the dictionary into a dataframe.
    pop_df = pd.DataFrame(pop_df_dicc)
    # Determine the sequence lengths.
    pop_df['LENGTH'] = pop_df['END'] - pop_df['START']
    # Sort the dataframe.
    pop_collasped_tracts = pop_df.sort_values(by=['LENGTH'], ascending=False, ignore_index=True)
    # Extract all the introgressed tracts.
    pop_starts = pop_collasped_tracts['START'].values
    pop_ends = pop_collasped_tracts['END'].values
    # Find the minimum start and maximum end.
    min_pos = pop_starts.min()
    max_pos = pop_ends.max()
    # Intialize the bins and bin counts.
    tract_bins = np.arange(min_pos, max_pos+1)
    tract_counts = np.zeros(tract_bins.size)
    # Iterate through each interval.
    for start, end in zip(pop_starts, pop_ends):
        # Update the positions counts.
        tract_counts[start-min_pos:end-min_pos+1] += 1
    return pop_collasped_tracts, pop_starts, pop_ends, tract_bins, tract_counts

# Define a function to plot mxl tracts for muc19.
def plot_mxl_muc19_tracts(tracts_df):
    # Subset the mxl tracts.
    mxl_df = tracts_df[tracts_df['POP'] == 'MXL'].sort_values(by=['LENGTH'], ascending=False, ignore_index=True)
    # Extract the start and end positions.
    starts, ends = mxl_df['START'].values, mxl_df['END'].values
    # Load the gene information.
    hg19_gene_df = pd.read_csv(f'../annotations/hg19_genes/ncbi_refseq_genes_chr12.csv.gz')
    genes_df = genes_any_overlap(hg19_gene_df, starts.min(), ends.max())
    gene_dicc = {
        gene_id: [genes_df['START'].values[i], genes_df['STOP'].values[i]] 
        for i, gene_id in enumerate(genes_df['GENE_ID'].values)
    }
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the figure.
    fig = plt.figure(
        figsize=(6, 4), dpi=300,
        facecolor='white', constrained_layout=True,
    )
    ax = fig.add_subplot(111)
    # Plot the 72kb focal region.
    ax.axvspan(40759001, 40831000, alpha=0.05, facecolor='black')
    # For every introgressed tract.
    for i in range(mxl_df.shape[0]):
        # Plot the tract.
        ax.plot([starts[i], ends[i]],[i+1, i+1], color='#009E73', lw=1.25)
    # For every gene.
    for gene in gene_dicc:
        # Plot the gene segments within the region.
        ax.plot(
            [max(gene_dicc[gene][0], starts.min()), min(gene_dicc[gene][-1], ends.max())],
            [-2, -2], color='black', marker='|', ms=10, lw=1, zorder=5,
        )
        # Annotate the gene.
        ax.text(
            ((max(gene_dicc[gene][0], starts.min()) + min(gene_dicc[gene][-1], ends.max())) / 2), -2,
            fr'${gene}$', fontsize=8, horizontalalignment='center',verticalalignment='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black'), zorder=5,
        )
    # Label the axes.
    ax.set_xlabel(f'Chromosome 12: {int(starts.min())} - {int(ends.max() - 1)}')
    ax.set_ylabel('MXL Introgressed Tracts')
    # Turn off the y-ticks for the first subplot.
    ax.set_yticks([])
    # Rescale x-axis to Mb.
    x_ticks = np.arange(40_250_000, 41_100_000, 100_000)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{round(x_tick / 1e6, 3)} Mb' for x_tick in x_ticks])
    # Export the plot.
    plt.savefig(
        './supp_figures/png/mxl_introgressed_tracts_muc19.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/mxl_introgressed_tracts_muc19.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/mxl_introgressed_tracts_muc19.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to plot the introgressed tract denisty.
def plot_tgp_muc19_tract_density(tracts_df, tgp_tract_pos, tgp_tract_counts):
    # Intialize the distribution data for the tgp.
    tgp_spops = tracts_df['SUPERPOP'].values
    tgp_min, tgp_max = tracts_df['START'].min(), tracts_df['END'].max()
    tgp_bins = np.arange(tgp_min, tgp_max+1)
    tgp_counts = {
        'EUR': np.zeros(tgp_bins.size),
        'SAS': np.zeros(tgp_bins.size),
        'EAS': np.zeros(tgp_bins.size),
        'AMR': np.zeros(tgp_bins.size),
    }
    # For every superpopulation.
    for spop in tgp_counts:
        # Determine the superpopulation mask.
        spop_df_mask = np.in1d(tgp_spops, [spop])
        # Subset the dataframe.
        spop_tracts_df = tracts_df[spop_df_mask].reset_index(drop=True)
        # Extract all the introgressed tracts.
        spop_starts = spop_tracts_df['START'].values
        spop_ends = spop_tracts_df['END'].values
        # Iterate through each interval.
        for start, end in zip(spop_starts, spop_ends):
            # Update the positions counts.
            tgp_counts[spop][start-tgp_min:end-tgp_min+1] += 1

    # Intialize dictionaries for plotting.
    color_dicc = {
        'AMR': '#009E73', 'SAS': '#CC79A7', 'EAS': '#0072B2', 'EUR': '#D55E00',
    }
    plot_dicc = {
        'min':tgp_min, 'max': tgp_max, 'gene_loc': -20,
        'bins': np.arange(tgp_bins.min(), tgp_bins.max()+1.5),
        'pos': [tgp_bins.astype(int), tgp_bins.astype(int), tgp_bins.astype(int), tgp_bins.astype(int)],
        'weights': [tgp_counts[spop].astype(int) for spop in tgp_counts],
        'colors': [color_dicc[spop] for spop in tgp_counts],
        'legend': [Patch(color=color_dicc[spop], label=f'{spop}') for spop in color_dicc], 
    }
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
   # Intialize the figure and axes.
    fig = plt.figure(
        figsize=(6, 4), dpi=300,
        facecolor='white', constrained_layout=True,
    )
    ax = fig.add_subplot(111)
    # Plot the 72kb focal region.
    ax.axvspan(40759001, 40831000, alpha=0.05, facecolor='black', zorder=5)
    # Plot the distribution.
    ax.hist(
        plot_dicc['pos'], bins=plot_dicc['bins'], weights=plot_dicc['weights'],
        histtype='stepfilled', color=plot_dicc['colors'], stacked=True,
    )
    ax.hist(
        tgp_tract_pos.astype(int), bins=plot_dicc['bins'], weights=tgp_tract_counts.astype(int),
        histtype='step', color='black', linewidth=0.75,
    )
    # Add axes labels.
    ax.set_ylabel('Introgressed Tract Count')
    ax.set_xlabel(f'Chromosome 12: {int(tgp_min)} - {int(tgp_max - 1)}')
    # Rescale x-axis to Mb.
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels([f'{round(x_tick / 1e6, 3)} Mb' for x_tick in ax.get_xticks()])
    # Add a figure lgend.
    legend = ax.legend(
        handles=plot_dicc['legend'],
        loc='upper left', bbox_to_anchor=(0, 1),
        frameon=True, fancybox=True, shadow=True, fontsize=8,
    )
    # Export the plot.
    plt.savefig(
        './supp_figures/png/tgp_introgressed_tract_density_muc19.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/tgp_introgressed_tract_density_muc19.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/tgp_introgressed_tract_density_muc19.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to plot the super population tracts and tract denisty.
def plot_tract_info_per_super_population(tracts_df, spop):
    # Intialize a dictionary to store populations.
    spop_dicc = {
        'AMR': ['MXL', 'PEL', 'CLM', 'PUR'],
        'SAS': ['BEB', 'STU', 'ITU', 'PJL', 'GIH'],
        'EAS': ['CHB', 'KHV', 'CHS', 'JPT', 'CDX'],
        'EUR': ['TSI', 'CEU', 'IBS', 'GBR', 'FIN'],
    }
    # Intialize dictionaries for plotting.
    color_dicc = {
        'AMR': '#009E73', 'SAS': '#CC79A7', 'EAS': '#0072B2', 'EUR': '#D55E00',
    }
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        2, len(spop_dicc[spop]), figsize=(12, 6), dpi=300,
        facecolor='white',
        sharex=True, sharey=False,
        constrained_layout=True,
    )
    # For every population.
    for i, pop in enumerate(spop_dicc[spop]):
        # Subset the tracts.
        pop_df = tracts_df[tracts_df['POP'] == pop].sort_values(by=['LENGTH'], ascending=False, ignore_index=True)
        # Extract the start and end positions.
        starts, ends = pop_df['START'].values, pop_df['END'].values
        pop_min, pop_max = starts.min(), ends.max()
        pop_bins = np.arange(pop_min, pop_max+1)
        pop_counts = np.zeros(pop_bins.size)
        # Iterate through each interval.
        for start, end in zip(starts, ends):
            # Update the positions counts.
            pop_counts[start-pop_min:end-pop_min+1] += 1
        # Plot the tract denisty.
        axes[0, i].hist(
            pop_bins.astype(int), bins=np.arange(pop_min, pop_max+1.5),
            weights=pop_counts.astype(int), histtype='step', color=color_dicc[spop],
        )
        # For every introgressed tract.
        for j in range(pop_df.shape[0]):
            # Plot the tract.
            axes[1, i].plot([starts[j], ends[j]],[j+1, j+1], color=color_dicc[spop], lw=1.25)
        # Add axes labels.
        axes[0, i].set_ylabel(f'{pop} Introgressed Tract Count')
        axes[0, i].set_xlabel('Chromosome 12')
        axes[1, i].set_xlabel('Chromosome 12')
        axes[1, i].set_ylabel(f'{pop} Introgressed Tracts')
        axes[1, i].set_yticks([])
    # For every subplot.
    for i, ax in enumerate(axes.flatten()):
        # Plot the 72kb focal region.
        ax.axvspan(40759001, 40831000, alpha=0.05, facecolor='black')
        # Plot the panel label.
        ax.set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Force labels to appear on all subplots.
        ax.xaxis.set_tick_params(which='both', labelbottom=True)
        # Rescale x-axis to Mb.
        ax.set_xticks(ax.get_xticks())
        ax.set_xticklabels(
            [f'{round(x_tick / 1e6, 1)} Mb' for x_tick in ax.get_xticks()],
            rotation=45, ha='right', rotation_mode='anchor',
        )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/{spop.lower()}_introgressed_tract_info_muc19.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/{spop.lower()}_introgressed_tract_info_muc19.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/{spop.lower()}_introgressed_tract_info_muc19.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to plot the introgressed tract denisty.
def plot_tgp_sr_rr_tract_density(tracts_df, tgp_tract_pos, tgp_tract_counts):
    # Intialize the distribution data for the tgp.
    tgp_spops = tracts_df['SUPERPOP'].values
    tgp_min, tgp_max = tracts_df['START'].min(), tracts_df['END'].max()
    tgp_bins = np.arange(tgp_min, tgp_max+1)
    tgp_counts = {
        'EUR': np.zeros(tgp_bins.size),
        'SAS': np.zeros(tgp_bins.size),
        'EAS': np.zeros(tgp_bins.size),
        'AMR': np.zeros(tgp_bins.size),
    }
    # For every superpopulation.
    for spop in tgp_counts:
        # Determine the superpopulation mask.
        spop_df_mask = np.in1d(tgp_spops, [spop])
        # Subset the dataframe.
        spop_tracts_df = tracts_df[spop_df_mask].reset_index(drop=True)
        # Extract all the introgressed tracts.
        spop_starts = spop_tracts_df['START'].values
        spop_ends = spop_tracts_df['END'].values
        # Iterate through each interval.
        for start, end in zip(spop_starts, spop_ends):
            # Update the positions counts.
            tgp_counts[spop][start-tgp_min:end-tgp_min+1] += 1

    # Intialize dictionaries for plotting.
    color_dicc = {
        'AMR': '#009E73', 'SAS': '#CC79A7', 'EAS': '#0072B2', 'EUR': '#D55E00',
    }
    plot_dicc = {
        'min':tgp_min, 'max': tgp_max, 'gene_loc': -20,
        'bins': np.arange(tgp_bins.min(), tgp_bins.max()+1.5),
        'pos': [tgp_bins.astype(int), tgp_bins.astype(int), tgp_bins.astype(int), tgp_bins.astype(int)],
        'weights': [tgp_counts[spop].astype(int) for spop in tgp_counts],
        'colors': [color_dicc[spop] for spop in tgp_counts],
        'legend': [Patch(color=color_dicc[spop], label=f'{spop}') for spop in color_dicc], 
    }
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
   # Intialize the figure and axes.
    fig = plt.figure(
        figsize=(6, 4), dpi=300,
        facecolor='white', constrained_layout=True,
    )
    ax = fig.add_subplot(111)
    # Plot the short read repeat region.
    ax.axvspan(40875941, 40885367, alpha=0.25, facecolor='black', zorder=5)
    # Plot the distribution.
    ax.hist(
        plot_dicc['pos'], bins=plot_dicc['bins'], weights=plot_dicc['weights'],
        histtype='stepfilled', color=plot_dicc['colors'], stacked=True,
    )
    ax.hist(
        tgp_tract_pos.astype(int), bins=plot_dicc['bins'], weights=tgp_tract_counts.astype(int),
        histtype='step', color='black', linewidth=0.75,
    )
    # Add axes labels.
    ax.set_ylabel('Introgressed Tract Count')
    ax.set_xlabel(f'Chromosome 12: {int(tgp_min)} - {int(tgp_max - 1)}')
    # Rescale x-axis to Mb.
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels([f'{round(x_tick / 1e6, 3)} Mb' for x_tick in ax.get_xticks()])
    # Add a figure lgend.
    legend = ax.legend(
        handles=plot_dicc['legend'],
        loc='upper left', bbox_to_anchor=(0, 1),
        frameon=True, fancybox=True, shadow=True, fontsize=8,
    )
    # Export the plot.
    plt.savefig(
        './supp_figures/png/tgp_introgressed_tract_density_short_read_repeat_region.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/tgp_introgressed_tract_density_short_read_repeat_region.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/tgp_introgressed_tract_density_short_read_repeat_region.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return



############
### VNTR ###
############

# Define a function to load the short read vntr data.
def load_short_read_vntr_data():
    # Load the meta data file for the TGP.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a list of non-Afr populations.
    non_afr_pops = [
        'MXL', 'PEL', 'CLM', 'PUR',
        'BEB', 'STU', 'ITU', 'PJL', 'GIH',
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX',
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN',
    ]
    # Load the introgressed tracts for chromosome 12.
    chr12_tracts_df = pd.read_csv('../hmmix_tracts/tgp_hmmix_haplotype_tracts_chr12.csv.gz')
    # Subset and clean the introgressed tracts overlapping the short read repeat region.
    short_read_repeat_region_tracts_df, _, _, _, _ = load_region_tracts(
        tract_any_overlap(chr12_tracts_df, 40876395, 40885001), non_afr_pops,
    )
    # Load in the repeats dataframe.
    short_read_repeats_df = pd.read_csv('../vntr/hg38_muc19_short_read_repeats.csv')
    # Filter out ACB and ASW.
    not_acb_or_asw = (short_read_repeats_df['POP'].values != 'ACB') & (short_read_repeats_df['POP'].values != 'ASW') 
    short_read_repeats_df = short_read_repeats_df[not_acb_or_asw].reset_index(drop=True)
    # Extract the individuals with an introgressed tract at the repeat region.
    short_read_repeat_region_tract_haps = short_read_repeat_region_tracts_df['IND'].values
    short_read_repeat_region_tract_inds = np.array([hap.split('-')[0] for hap in short_read_repeat_region_tract_haps])
    # Intialize a list.
    n_short_read_reapeat_tracts = []
    # For every individual.
    for ind in short_read_repeats_df['IND'].values:
        # If this individual doesn't have any repeat tracts.
        if ind not in np.unique(short_read_repeat_region_tract_inds):
            # Set the number of tracts to 0.
            n_short_read_reapeat_tracts.append(0)
        # Else.
        else:
            # Create a mask for the tract dataframe.
            tract_mask = short_read_repeat_region_tract_inds == ind
            # Extract the tract lengths.
            ind_tracts = short_read_repeat_region_tracts_df[tract_mask]['LENGTH'].values
            # If there are two tracts.
            if ind_tracts.size == 2:
                # Update the lists.
                n_short_read_reapeat_tracts.append(2)
            # Else there is only one tract.
            else:
                # Update the lists.
                n_short_read_reapeat_tracts.append(1)
    # Update the dataframe.
    short_read_repeats_df['N_SR_TRACTS'] = n_short_read_reapeat_tracts
    # Determine the elevated threshold.
    short_read_thresh = np.percentile(short_read_repeats_df['REPEAT_COPIES'].values, 95)
    # Add the outlier column.
    short_read_repeats_df['OUT'] = short_read_repeats_df['REPEAT_COPIES'].values > short_read_thresh
    # Subset the dataframe to only include outliers.
    short_read_repeats_outlier_df = short_read_repeats_df[short_read_repeats_df['OUT'].values].reset_index(drop=True)
    # Print the outlier threshold.
    print(f'Repeat Copies Outlier Threshold: {short_read_thresh}')
    # Rename all the columns to look pretty.
    short_read_repeats_export_df = short_read_repeats_df.rename(
        columns={
            'IND': 'Individual', 'POP': 'Population', 'SUPERPOP': 'Super Population',
            'REPEAT_COPIES': 'Repeat Copies',
            'N_SR_TRACTS': 'Number of Tracts Overlapping the Repeat Region',
            'OUT': 'Repeat Copies > 95th Percentile',
        }
    )
    # Export the dataframe as a csv.
    short_read_repeats_export_df.to_csv('./dataframes/tgp_muc19_short_read_vntr_info.csv.gz', index=False)
    return short_read_repeats_df, short_read_repeats_outlier_df

# Define a function to generate the super population, population, and haplotype group short read repeat count summary.
def compile_short_read_vntr_summary(short_read_repeats_df):
    # Intialize dictionaries.
    short_read_hap_summary = {'hap_grp': [], 'n': [], 'mean': [], 'std':[], 'sem': [], '95_ci': [], 'n_outs': [], 'p_outs': []}
    short_read_spop_summary = {'spop': [], 'n': [], 'mean': [], 'std':[], 'sem': [], '95_ci': [], 'n_outs': [], 'p_outs': []}
    short_read_pop_summary = {'spop': [], 'pop': [], 'n': [],  'mean': [], 'std':[], 'sem': [],  '95_ci': [], 'n_outs': [], 'p_outs': []}
    # For every super population.
    for spop in short_read_repeats_df['SUPERPOP'].unique():
        # Subset the dataframe.
        spop_df = short_read_repeats_df[short_read_repeats_df['SUPERPOP'] == spop]
        # Extract the repeat counts.
        repeat_counts = spop_df['REPEAT_COPIES'].values
        # Determine the sem and 95% ci.
        repeat_se, repeat_ci = sem_ci_of_mean(repeat_counts)
        # Fill the dictionary.
        short_read_spop_summary['spop'].append(spop)
        short_read_spop_summary['n'].append(spop_df.shape[0])
        short_read_spop_summary['mean'].append(np.mean(repeat_counts))
        short_read_spop_summary['std'].append(np.std(repeat_counts))
        short_read_spop_summary['sem'].append(repeat_se)
        short_read_spop_summary['95_ci'].append(repeat_ci)
        short_read_spop_summary['n_outs'].append(spop_df['OUT'].sum())
        short_read_spop_summary['p_outs'].append(spop_df['OUT'].sum() / spop_df.shape[0])
        # For every population.
        for pop in spop_df['POP'].unique():
            # Subset the dataframe.
            pop_df = spop_df[spop_df['POP'] == pop]
            # Extract the repeat counts.
            repeat_counts = pop_df['REPEAT_COPIES'].values
            # Determine the sem and 95% ci.
            repeat_se, repeat_ci = sem_ci_of_mean(repeat_counts)
            # Fill the dictionary.
            short_read_pop_summary['spop'].append(spop)
            short_read_pop_summary['pop'].append(pop)
            short_read_pop_summary['n'].append(pop_df.shape[0])
            short_read_pop_summary['mean'].append(np.mean(repeat_counts))
            short_read_pop_summary['std'].append(np.std(repeat_counts))
            short_read_pop_summary['sem'].append(repeat_se)
            short_read_pop_summary['95_ci'].append(repeat_ci)
            short_read_pop_summary['n_outs'].append(pop_df['OUT'].sum())
            short_read_pop_summary['p_outs'].append(pop_df['OUT'].sum() / pop_df.shape[0])
    # For every haplotype group.
    for col, hap_grp, grp_lab in [
        ('N_SR_TRACTS', 0, 'Inds. with No Tracts Overlapping the Repeat Region'),
        ('N_SR_TRACTS', 1, 'Inds. with One Tract Overlapping the Repeat Region'),
        ('N_SR_TRACTS', 2, 'Inds. with Two Tracts Overlapping the Repeat Region'),
    ]:
        # Subset the dataframe.
        sub_df = short_read_repeats_df[short_read_repeats_df[col] == hap_grp]
        # Extract the repeat counts.
        repeat_counts = sub_df['REPEAT_COPIES'].values
        # Determine the sem and 95% ci.
        repeat_se, repeat_ci = sem_ci_of_mean(repeat_counts)
        # Fill the dictionary.
        short_read_hap_summary['hap_grp'].append(grp_lab)
        short_read_hap_summary['n'].append(repeat_counts.size)
        short_read_hap_summary['mean'].append(np.mean(repeat_counts))
        short_read_hap_summary['std'].append(np.std(repeat_counts))
        short_read_hap_summary['sem'].append(repeat_se)
        short_read_hap_summary['95_ci'].append(repeat_ci)
        short_read_hap_summary['n_outs'].append(sub_df['OUT'].sum())
        short_read_hap_summary['p_outs'].append(sub_df['OUT'].sum() / sub_df.shape[0])
    # Conver the dictionaries to dataframes.
    short_read_hap_summary_df = pd.DataFrame(short_read_hap_summary)
    short_read_spop_summary_df = pd.DataFrame(short_read_spop_summary).sort_values('mean', ascending=False).reset_index(drop=True)
    short_read_pop_summary_df = pd.DataFrame(short_read_pop_summary).sort_values('mean', ascending=False).reset_index(drop=True)
    # Rename all the columns to look pretty.
    short_read_hap_summary_df.rename(
        columns={
            'hap_grp': 'Group', 'n': 'Number of Inds.',
            'mean': r'Repeat Copies $\left( \mu \right)$',
            'std': r'Repeat Copies $\left( \sigma \right)$',
            'sem': r'Repeat Copies $\left( SEM \right)$',
            '95_ci': r'Repeat Copies $\left( \pm CI_{95\%} \right)$',
            'n_outs': 'Number of Inds. with Elevated Repeats',
            'p_outs': 'Proportion of Inds. with Elevated Repeats',
        }, inplace=True
    )
    short_read_spop_summary_df.rename(
        columns={
            'spop': 'Super Population', 'n': 'Number of Inds.',
            'mean': r'Repeat Copies $\left( \mu \right)$',
            'std': r'Repeat Copies $\left( \sigma \right)$',
            'sem': r'Repeat Copies $\left( SEM \right)$',
            '95_ci': r'Repeat Copies $\left( \pm CI_{95\%} \right)$',
            'n_outs': 'Number of Inds. with Elevated Repeats',
            'p_outs': 'Proportion of Inds. with Elevated Repeats',
        }, inplace=True
    )
    short_read_pop_summary_df.rename(
        columns={
            'spop': 'Super Population', 'pop': 'Population', 'n': 'Number of Inds.',
            'mean': r'Repeat Copies $\left( \mu \right)$',
            'std': r'Repeat Copies $\left( \sigma \right)$',
            'sem': r'Repeat Copies $\left( SEM \right)$',
            '95_ci': r'Repeat Copies $\left( \pm CI_{95\%} \right)$',
            'n_outs': 'Number of Inds. with Elevated Repeats',
            'p_outs': 'Proportion of Inds. with Elevated Repeats',
        }, inplace=True
    )
    # Export the dataframe as a csv.
    short_read_hap_summary_df.to_csv('./dataframes/tgp_short_read_vntr_introgressed_tract_summary.csv.gz', index=False)
    short_read_spop_summary_df.to_csv('./dataframes/tgp_short_read_vntr_superpopulation_summary.csv.gz', index=False)
    short_read_pop_summary_df.to_csv('./dataframes/tgp_short_read_vntr_population_summary.csv.gz', index=False)
    return short_read_spop_summary_df, short_read_pop_summary_df, short_read_hap_summary_df

# Define a function to generate the elevated repeats summary for all individuals.
def compile_short_read_vntr_all_repeats_summary(short_read_repeats_df):
    # Intialize dictionaries to store the results.
    short_read_all_z_and_u_tests = {
        'grp1': [], 'grp2': [],
        'tot_grp1': [], 'tot_grp2': [],
        'out_grp1': [], 'out_grp2': [],
        'p_grp1': [], 'p_grp2': [],
        'm_grp1': [], 'm_grp2': [],
        's_grp1': [], 's_grp2': [],
        'se_grp1': [], 'se_grp2': [],
        'ci_grp1': [], 'ci_grp2': [],
        'z': [], 'z_p': [], 'u': [], 'u_p': [], 
    }
    short_read_all_z_and_u_plots = {
        'TRACTS': {}, 'AMR': {}
    }
    # Intialize a list of comparisons.
    short_read_all_z_and_u_configs = [
        (short_read_repeats_df, 'N_SR_TRACTS', 0, r'All Inds.: $\geq 1$ Overlapping Tract', 'All Inds.: No Overlapping Tract', 'TRACTS'),
        (short_read_repeats_df, 'SUPERPOP', 'AMR', 'All AMR Inds.', 'All non-AMR Inds.', 'AMR'),
    ]
    # For every configuration.
    for df, col, val, grp1, grp2, plot_key in short_read_all_z_and_u_configs:
        # Extract the count data.
        if col == 'SUPERPOP':
            grp1_counts = df[df[col] == val]['REPEAT_COPIES'].values
            grp2_counts = df[df[col] != val]['REPEAT_COPIES'].values
            grp1_outs = df[(df[col] == val) & (df['OUT'].values)]['REPEAT_COPIES'].values
            grp2_outs = df[(df[col] != val) & (df['OUT'].values)]['REPEAT_COPIES'].values
        else:
            grp1_counts = df[df[col] != val]['REPEAT_COPIES'].values
            grp2_counts = df[df[col] == val]['REPEAT_COPIES'].values
            grp1_outs = df[(df[col] != val) & (df['OUT'].values)]['REPEAT_COPIES'].values
            grp2_outs = df[(df[col] == val) & (df['OUT'].values)]['REPEAT_COPIES'].values
        # Run tests.
        u, u_p = stats.mannwhitneyu(grp1_counts, grp2_counts, alternative='greater')
        z, z_p = proportions_ztest(
            [grp1_outs.size, grp2_outs.size],
            [grp1_counts.size, grp2_counts.size],
            alternative='larger',
        )
        # Determine the sem and 95% ci.
        se_grp1, ci_grp1 = sem_ci_of_mean(grp1_counts)
        se_grp2, ci_grp2 = sem_ci_of_mean(grp2_counts)
        # Fill the dictionary.
        short_read_all_z_and_u_tests['grp1'].append(grp1)
        short_read_all_z_and_u_tests['grp2'].append(grp2)
        short_read_all_z_and_u_tests['tot_grp1'].append(grp1_counts.size)
        short_read_all_z_and_u_tests['tot_grp2'].append(grp2_counts.size)
        short_read_all_z_and_u_tests['out_grp1'].append(grp1_outs.size)
        short_read_all_z_and_u_tests['out_grp2'].append(grp2_outs.size)
        short_read_all_z_and_u_tests['p_grp1'].append(grp1_outs.size / grp1_counts.size)
        short_read_all_z_and_u_tests['p_grp2'].append(grp2_outs.size / grp2_counts.size)
        short_read_all_z_and_u_tests['m_grp1'].append(np.mean(grp1_counts))
        short_read_all_z_and_u_tests['m_grp2'].append(np.mean(grp2_counts))
        short_read_all_z_and_u_tests['s_grp1'].append(np.std(grp1_counts))
        short_read_all_z_and_u_tests['s_grp2'].append(np.std(grp2_counts))
        short_read_all_z_and_u_tests['se_grp1'].append(se_grp1)
        short_read_all_z_and_u_tests['se_grp2'].append(se_grp2)
        short_read_all_z_and_u_tests['ci_grp1'].append(ci_grp1)
        short_read_all_z_and_u_tests['ci_grp2'].append(ci_grp2)
        short_read_all_z_and_u_tests['z'].append(z)
        short_read_all_z_and_u_tests['z_p'].append(z_p)
        short_read_all_z_and_u_tests['u'].append(u)
        short_read_all_z_and_u_tests['u_p'].append(u_p)
        # Fill the plotting dictionary.
        short_read_all_z_and_u_plots[plot_key]['counts'] = [grp1_counts, grp2_counts]
        short_read_all_z_and_u_plots[plot_key]['means'] = [np.mean(grp1_counts), np.mean(grp2_counts)]
        short_read_all_z_and_u_plots[plot_key]['cis'] = [ci_grp1, ci_grp2]
        short_read_all_z_and_u_plots[plot_key]['pvals'] = [u_p, z_p]
    # Convert to a dataframe.
    short_read_all_z_and_u_tests_df = pd.DataFrame(short_read_all_z_and_u_tests)
    short_read_all_z_and_u_tests_df.rename(
        columns={
            'grp1': 'Group 1', 'grp2': 'Group 2',
            'tot_grp1': 'Number of Inds. (Group 1)', 'tot_grp2': 'Number of Inds. (Group 2)',
            'out_grp1': 'Number of Inds. with Elevated Repeats (Group 1)', 'out_grp2': 'Number of Inds. with Elevated Repeats (Group 2)',
            'p_grp1': 'Prop. of Inds. with Elevated Repeats (Group 1)', 'p_grp2': 'Prop. of Inds. with Elevated Repeats (Group 2)',
            'm_grp1': r'Group 1 Repeat Copies $\left( \mu \right)$',
            'm_grp2': r'Group 2 Repeat Copies $\left( \mu \right)$',
            's_grp1': r'Group 1 Repeat Copies $\left( \sigma \right)$',
            's_grp2': r'Group 2 Repeat Copies $\left( \sigma \right)$',
            'se_grp1': r'Group 1 Repeat Copies $\left( SEM \right)$',
            'se_grp2': r'Group 2 Repeat Copies $\left( SEM \right)$',
            'ci_grp1': r'Group 1 Repeat Copies $\left( \pm CI_{95\%} \right)$',
            'ci_grp2': r'Group 2 Repeat Copies $\left( \pm CI_{95\%} \right)$',
            'z': r'$Z-Statistic$', 'z_p': r'$P-value$ (Prop. $Z$-Test)',
            'u': r'$U-Statistic$', 'u_p': r'$P-value$ (Mann-Whitney $U$-Test)', 
        }, inplace=True
    )
    short_read_all_z_and_u_tests_df.to_csv('./dataframes/tgp_short_read_vntr_group_comparisons_all_individuals.csv.gz', index=False)
    return short_read_all_z_and_u_tests_df, short_read_all_z_and_u_plots

# Define a function to generate the elevated repeats summary for individuals with elevated repeats.
def compile_short_read_vntr_elevated_repeats_summary(short_read_repeats_outlier_df):
    # Intialize dictionaries to store the results.
    short_read_outlier_z_and_u_tests = {
        'grp1': [], 'grp2': [], 'tot': [],
        'n_grp1': [], 'n_grp2': [],
        'p_grp1': [], 'p_grp2': [],
        'm_grp1': [], 'm_grp2': [],
        's_grp1': [], 's_grp2': [],
        'se_grp1': [], 'se_grp2': [],
        'ci_grp1': [], 'ci_grp2': [],
        'z': [], 'z_p': [], 'u': [], 'u_p': [], 
    }
    short_read_outlier_z_and_u_plots = {
        'TRACTS': {}, 'AMR': {}
    }
    # Intialize a list of comparisons.
    short_read_outlier_z_and_u_configs = [
        (short_read_repeats_outlier_df, 'N_SR_TRACTS', 0, r'Outlier Inds.: $\geq 1$ Overlapping Tract', 'Outlier Inds.: No Overlapping Tract', 'TRACTS'),
        (short_read_repeats_outlier_df, 'SUPERPOP', 'AMR', 'Outlier AMR Inds.', 'Outlier non-AMR Inds.', 'AMR'),
    ]
    # For every configuration.
    for df, col, val, grp1, grp2, plot_key in short_read_outlier_z_and_u_configs:
        # Extract the count data.
        if col == 'SUPERPOP':
            grp1_counts = df[df[col] == val]['REPEAT_COPIES'].values
            grp2_counts = df[df[col] != val]['REPEAT_COPIES'].values
        else:
            grp1_counts = df[df[col] != val]['REPEAT_COPIES'].values
            grp2_counts = df[df[col] == val]['REPEAT_COPIES'].values
        # Run tests.
        u, u_p = stats.mannwhitneyu(grp1_counts, grp2_counts, alternative='greater')
        z, z_p = proportions_ztest([grp1_counts.size, grp2_counts.size], df.shape[0], alternative='larger')
        # Determine the sem and 95% ci.
        se_grp1, ci_grp1 = sem_ci_of_mean(grp1_counts)
        se_grp2, ci_grp2 = sem_ci_of_mean(grp2_counts)
        # Fill the dictionary.
        short_read_outlier_z_and_u_tests['grp1'].append(grp1)
        short_read_outlier_z_and_u_tests['grp2'].append(grp2)
        short_read_outlier_z_and_u_tests['tot'].append(df.shape[0])
        short_read_outlier_z_and_u_tests['n_grp1'].append(grp1_counts.size)
        short_read_outlier_z_and_u_tests['n_grp2'].append(grp2_counts.size)
        short_read_outlier_z_and_u_tests['p_grp1'].append(grp1_counts.size / df.shape[0])
        short_read_outlier_z_and_u_tests['p_grp2'].append(grp2_counts.size / df.shape[0])
        short_read_outlier_z_and_u_tests['m_grp1'].append(np.mean(grp1_counts))
        short_read_outlier_z_and_u_tests['m_grp2'].append(np.mean(grp2_counts))
        short_read_outlier_z_and_u_tests['s_grp1'].append(np.std(grp1_counts))
        short_read_outlier_z_and_u_tests['s_grp2'].append(np.std(grp2_counts))
        short_read_outlier_z_and_u_tests['se_grp1'].append(se_grp1)
        short_read_outlier_z_and_u_tests['se_grp2'].append(se_grp2)
        short_read_outlier_z_and_u_tests['ci_grp1'].append(ci_grp1)
        short_read_outlier_z_and_u_tests['ci_grp2'].append(ci_grp2)
        short_read_outlier_z_and_u_tests['z'].append(z)
        short_read_outlier_z_and_u_tests['z_p'].append(z_p)
        short_read_outlier_z_and_u_tests['u'].append(u)
        short_read_outlier_z_and_u_tests['u_p'].append(u_p)
        # Fill the plotting dictionary.
        short_read_outlier_z_and_u_plots[plot_key]['counts'] = [grp1_counts, grp2_counts]
        short_read_outlier_z_and_u_plots[plot_key]['means'] = [np.mean(grp1_counts), np.mean(grp2_counts)]
        short_read_outlier_z_and_u_plots[plot_key]['cis'] = [ci_grp1, ci_grp2]
        short_read_outlier_z_and_u_plots[plot_key]['pvals'] = [u_p, z_p]
    # Convert to a dataframe.
    short_read_outlier_z_and_u_tests_df = pd.DataFrame(short_read_outlier_z_and_u_tests)
    short_read_outlier_z_and_u_tests_df.rename(
        columns={
            'grp1': 'Group 1', 'grp2': 'Group 2', 'tot': 'Total Number of Inds. with Elevated Repeats',
            'n_grp1': 'Number of Inds. (Group 1)', 'n_grp2': 'Number of Inds. (Group 2)',
            'p_grp1': 'Prop. of Inds. (Group 1)', 'p_grp2': 'Prop. of Inds. (Group 2)',
            'm_grp1': r'Group 1 Repeat Copies $\left( \mu \right)$',
            'm_grp2': r'Group 2 Repeat Copies $\left( \mu \right)$',
            's_grp1': r'Group 1 Repeat Copies $\left( \sigma \right)$',
            's_grp2': r'Group 2 Repeat Copies $\left( \sigma \right)$',
            'se_grp1': r'Group 1 Repeat Copies $\left( SEM \right)$',
            'se_grp2': r'Group 2 Repeat Copies $\left( SEM \right)$',
            'ci_grp1': r'Group 1 Repeat Copies $\left( \pm CI_{95\%} \right)$',
            'ci_grp2': r'Group 2 Repeat Copies $\left( \pm CI_{95\%} \right)$',
            'z': r'$Z-Statistic$', 'z_p': r'$P-value$ (Prop. $Z$-Test)',
            'u': r'$U-Statistic$', 'u_p': r'$P-value$ (Mann-Whitney $U$-Test)', 
        }, inplace=True
    )
    short_read_outlier_z_and_u_tests_df.to_csv('./dataframes/tgp_short_read_vntr_group_comparisons_outlier_individuals.csv.gz', index=False)
    return short_read_outlier_z_and_u_tests_df, short_read_outlier_z_and_u_plots

# Define a function to plot repeat count group comparisons.
def plot_vntr_group_comparisons(plot_key, all_inds, out_inds):
    # Intialize the x-axis positions.
    x_pos = [1, 3]
    # Intialize a plot info dictionary.
    plot_info = {
        'TRACTS': {
            0: {
                'label': 'All Individuals',
                'ticks': [r'$\geq 1$ Overlapping Tract', 'No Overlapping Tract'],
            },
            1: {
                'label': 'Outlier Individuals',
                'ticks': [r'$\geq 1$ Overlapping Tract', 'No Overlapping Tract'],
            },
            'export': 'tgp_short_read_vntr_introgressed_tract_group_comparison',
        },
        'AMR': {
            0: {
                'label': 'All Individuals',
                'ticks': ['AMR', 'Non-AMR'],
            },
            1: {
                'label': 'Outlier Individuals',
                'ticks': ['AMR', 'Non-AMR'],
            },
            'export': 'tgp_short_read_vntr_superpopulation_group_comparison',
        }
    }
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 2, figsize=(7, 3.5), dpi=300,
        facecolor='white',
        sharex=False, sharey=True,
        constrained_layout=True,
    )
    # For each subplot.
    for i, plot_dicc in enumerate([all_inds, out_inds]):
        # Plot the distributions.
        vp = axes[i].violinplot(plot_dicc[plot_key]['counts'], positions=x_pos, widths=0.9, showextrema=False)
        # Adjust the violin plot color.
        for pc in vp['bodies']:
            pc.set_facecolor('grey')
            pc.set_edgecolor('grey')
            pc.set_alpha(0.75)
        # For every group.
        for j, grp in enumerate(plot_dicc[plot_key]['counts']):
            # Plot the mean and 95% ci.
            axes[i].errorbar(
                x_pos[j], plot_dicc[plot_key]['means'][j], yerr=plot_dicc[plot_key]['cis'][j],
                color='white', marker='o', markersize=5,
                elinewidth=1.5, capthick=1.5, capsize=5, zorder=5,
            )
            # Generate some jitter to the x-axis.
            jitter = np.random.normal(x_pos[j], 0.1, size=grp.size)
            # Plot the points!
            axes[i].scatter(
                jitter, grp, color='black',
                marker='x', s=10, linewidths=0.75,
            )
        # Add a figure lgend.
        mw_p, zp_p = plot_dicc[plot_key]['pvals']
        axes[i].legend(
            handles=[
                Line2D([0], [0], color='none', linestyle='none', label=fr'Proportions $Z$-Test$ = {zp_p:.3e}$'),
                Line2D([0], [0], color='none', linestyle='none', label=fr'Mann–Whitney $U$-Test$ = {mw_p:.3e}$'),
            ], 
            loc='upper right', bbox_to_anchor=(1, 1), handlelength=0, handletextpad=0,
            frameon=True, fancybox=True, shadow=True, fontsize=8,
        )
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Set the axes labels.
        axes[i].set_ylabel('Repeat Copies')
        axes[i].set_xlabel(plot_info[plot_key][i]['label'])
        # Set the x-axis tick positions and labels.
        axes[i].set_xticks(x_pos)
        axes[i].set_xticklabels(plot_info[plot_key][i]['ticks'])
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/{plot_info[plot_key]["export"]}.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/{plot_info[plot_key]["export"]}.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/{plot_info[plot_key]["export"]}.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to summarize the correlations in repeat counts.
def compile_short_read_vntr_correlation_summary(short_read_repeats_df):
    # Intialize a dictionaries.
    short_read_corr_dicc = {
        'grp': [], 'rho': [], 'pval': [],
    }
    mxl_dicc = {}
    corr_info = {'df': short_read_repeats_df, 'col': 'N_SR_TRACTS'}
    tgp_dicc = {
        'AFR': ['LWK', 'GWD', 'MSL', 'ESN', 'YRI'],
        'AMR': ['MXL', 'PEL', 'CLM', 'PUR'],
        'SAS': ['BEB', 'STU', 'ITU', 'PJL', 'GIH'],
        'EAS': ['CHB', 'KHV', 'CHS', 'JPT', 'CDX'],
        'EUR': ['TSI', 'CEU', 'IBS', 'GBR', 'FIN'],
    }
    # Compile the results for the TGP.
    short_read_corr_dicc['grp'].append('TGP')
    # Unpack.
    df, col = corr_info['df'], corr_info['col']
    # Extract the data.
    n_haps, n_repeats = df[col].values, df['REPEAT_COPIES'].values
    # Run a correlation test.
    rho, pval = stats.spearmanr(n_haps, n_repeats)
    # If the p-value reported by scipy is 0.
    if pval == 0:
        # Determine the smallest p-value scipy would report.
        pval = np.nextafter(pval, 1)
    # Update the dictionary.
    short_read_corr_dicc['rho'].append(rho)
    short_read_corr_dicc['pval'].append(pval)
    # For every superpopulation.
    for spop in tgp_dicc:
        # Update the dictionary.
        short_read_corr_dicc['grp'].append(spop)
        # Unpack.
        df, col = corr_info['df'], corr_info['col']
        df = df[df['SUPERPOP'] == spop]
        # Extract the data.
        n_haps, n_repeats = df[col].values, df['REPEAT_COPIES'].values
        # If there is only one type of haplotype in the group.
        if np.unique(n_haps).size == 1:
            # Update the dictionary.
            short_read_corr_dicc['rho'].append(np.nan)
            short_read_corr_dicc['pval'].append(np.nan)
        # Else, there are at least two haplotype groups.
        else:
            # Run a correlation test.
            rho, pval = stats.spearmanr(n_haps, n_repeats)
            # If the p-value reported by scipy is 0.
            if pval == 0:
                # Determine the smallest p-value scipy would report.
                pval = np.nextafter(pval, 1)
            # Update the dictionary.
            short_read_corr_dicc['rho'].append(rho)
            short_read_corr_dicc['pval'].append(pval)
    # For every superpopulation.
    for spop in tgp_dicc:
        # For every population.
        for pop in tgp_dicc[spop]:
            # Update the dictionary.
            short_read_corr_dicc['grp'].append(pop)
            # Unpack.
            df, col = corr_info['df'], corr_info['col']
            df = df[df['POP'] == pop]
            # Extract the data.
            n_haps, n_repeats = df[col].values, df['REPEAT_COPIES'].values
            # If there is only one type of haplotype in the group.
            if np.unique(n_haps).size == 1:
                # Update the dictionary.
                short_read_corr_dicc['rho'].append(np.nan)
                short_read_corr_dicc['pval'].append(np.nan)
            # Else, there are at least two haplotype groups.
            else:
                # Run a correlation test.
                rho, pval = stats.spearmanr(n_haps, n_repeats)
                # If the p-value reported by scipy is 0.
                if pval == 0:
                    # Determine the smallest p-value scipy would report.
                    pval = np.nextafter(pval, 1)
                # Update the dictionary.
                short_read_corr_dicc['rho'].append(rho)
                short_read_corr_dicc['pval'].append(pval)
            # If the current population is mxl.
            if pop == 'MXL':
                # Fill the dictionary.
                mxl_dicc['data'] = [
                    df[df[col] == 0]['REPEAT_COPIES'].values,
                    df[df[col] == 1]['REPEAT_COPIES'].values,
                    df[df[col] == 2]['REPEAT_COPIES'].values,
                ]
                mxl_dicc['corr'] = [rho, pval]
    # Convert the dictionary to a dataframe.
    short_read_corr_df = pd.DataFrame(short_read_corr_dicc)
    short_read_corr_df.rename(
        columns={
            'grp': 'Group',
            'rho': r"Spearman's $\rho$",
            'pval': r'$P-value$',
        }, inplace=True
    )
    short_read_corr_df.to_csv('./dataframes/tgp_short_read_vntr_introgressed_tract_correlations.csv.gz', index=False)
    return short_read_corr_df, mxl_dicc

# Define a function to plot repeat count correlations for mxl.
def plot_vntr_mxl_correlations(plot_dicc):
    # Intialize the x-axis positions and labels.
    x_pos = np.array([0, 1, 2])
    x_ticklabs = x_pos.astype(str)
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the figure and axes.
    fig = plt.figure(
        figsize=(6, 4), dpi=300,
        facecolor='white', constrained_layout=True,
    )
    ax = fig.add_subplot(111)
    # Unpack.
    repeat_counts = plot_dicc['data']
    rho, pval = plot_dicc['corr']
    # Plot the alternating background.
    ax.axvspan(x_pos[1] - 0.5, x_pos[1] + 0.5, facecolor='black', alpha=0.05, zorder=1)
    # Plot the distributions.
    vp = ax.violinplot(repeat_counts, positions=x_pos, widths=0.9, showextrema=False)
    # Adjust the violin plot color.
    for pc in vp['bodies']:
        pc.set_facecolor('grey')
        pc.set_edgecolor('grey')
        pc.set_alpha(0.75)
    # For every group.
    for i, grp in enumerate(repeat_counts):
        # Compute 95% ci.
        _, grp_ci = sem_ci_of_mean(grp)
        # Plot the mean and 95% ci.
        ax.errorbar(
            x_pos[i], np.mean(grp), yerr=grp_ci,
            color='black', marker='o',markersize=5,
            elinewidth=1.5, capthick=1.5, capsize=5, zorder=5,
        )
        # Generate some jitter to the x-axis.
        jitter = np.random.normal(x_pos[i], 0.1, size=grp.size)
        # Plot the points!
        ax.scatter(
            jitter, grp, facecolor='white', marker='X', s=25,
            edgecolor='black', linewidths=0.25, zorder=4,
        )
    # Add a figure lgend.
    ax.legend(
        handles=[
            Line2D([0], [0], color='none', linestyle='none', label=f"Spearman's $\\rho = {rho:.3f}$"),
            Line2D([0], [0], color='none', linestyle='none', label=fr'$P-value = {pval:.3e}$')
        ], 
        loc='upper left', bbox_to_anchor=(0, 1), handlelength=0, handletextpad=0,
        frameon=True, fancybox=True, shadow=True, fontsize=8,
    )
    # Adjust the axes limits.
    ax.set_xlim(-0.5, 2.5)
    # Set the axes labels.
    ax.set_ylabel('Repeat Copies')
    ax.set_xlabel('Tracts Overlapping the Repeat Region in MXL Inds.')
    # Set the x-axis tick positions and labels.
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_ticklabs)
    # Export the plot.
    plt.savefig(
        './supp_figures/png/mxl_short_read_vntr_introgressed_tract_correlations.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/mxl_short_read_vntr_introgressed_tract_correlations.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/mxl_short_read_vntr_introgressed_tract_correlations.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to summarize the correlations between ancestry proportion and repeat counts.
def compile_short_read_vntr_amr_ancestry_summary(short_read_repeats_df):
    # Subset the dataframe for amr populations.
    amr_short_read_repeats_df = short_read_repeats_df[short_read_repeats_df['SUPERPOP'] == 'AMR'].reset_index(drop=True)
    # Load the ancestry proportions for the short read repeat region.
    short_read_anc_props_df = pd.read_csv('../amr_lai/anc_props/amr_short_read_repeat_region_props.csv.gz')
    # Intialize a dictionaries.
    anc_labs = {
        'NAT': 'Indigenous American',
        'EUR': 'European',
        'AFR': 'African',
    }
    anc_props_col = {
        'NAT': [], 'EUR': [],
        'AFR': [], 'UNK': [],
    }
    anc_prop_corr_dicc = {
        'pop': [], 'anc': [], 'rho': [], 'pval': [],
    }
    mxl_dicc = {anc: {} for anc in anc_labs}
    # Intialize a dictionary.
    anc_props_summary = {
        'pop': [], 'anc': [],
        'mean': [], 'std': [],
        'sem': [], '95_ci': [],
    }
    # For every amr population.
    for pop in ['MXL', 'PEL', 'CLM', 'PUR']:
        # Subset the dataframe.
        pop_df = short_read_anc_props_df[short_read_anc_props_df['POP'] == pop]
        # For every ancestry componenet.
        for anc in anc_labs:
            # Extract the ancestry proportions.
            anc_props = pop_df[anc].values
            # Determine the sem and 95% ci.
            anc_se, anc_ci = sem_ci_of_mean(anc_props)
            # Fill the dictionary.
            anc_props_summary['pop'].append(pop)
            anc_props_summary['anc'].append(anc_labs[anc])
            anc_props_summary['mean'].append(np.mean(anc_props))
            anc_props_summary['std'].append(np.std(anc_props))
            anc_props_summary['sem'].append(anc_se)
            anc_props_summary['95_ci'].append(anc_ci)
    # Extract the individuals.
    cnv_inds = amr_short_read_repeats_df['IND'].values
    anc_inds = short_read_anc_props_df['IND'].values
    # For every amr individual.
    for ind in cnv_inds:
        # Create a mask.
        ind_mask = anc_inds == ind
        # For every ancestry componenet.
        for anc in anc_props_col:
            # Update the dictionary.
            anc_props_col[anc].append(short_read_anc_props_df[ind_mask][anc].values[0])
    # For every ancestry proportion.
    for anc in anc_props_col:
        # Update the dataframe.
        amr_short_read_repeats_df[anc] = anc_props_col[anc]
    # For every amr population.
    for pop in ['MXL', 'PEL', 'CLM', 'PUR']:
        # Subset the dataframe.
        df = amr_short_read_repeats_df[amr_short_read_repeats_df['POP'] == pop]
        # For every ancestry component.
        for anc in anc_labs:
            # Extract the data.
            anc_props, n_repeats = df[anc].values, df['REPEAT_COPIES'].values
            # Run a correlation test.
            rho, pval = stats.spearmanr(anc_props, n_repeats)
            # If the p-value reported by scipy is 0.
            if pval == 0:
                # Determine the smallest p-value scipy would report.
                pval = np.nextafter(pval, 1)
            # Update the dictionary.
            anc_prop_corr_dicc['pop'].append(pop)
            anc_prop_corr_dicc['anc'].append(anc_labs[anc])
            anc_prop_corr_dicc['rho'].append(rho)
            anc_prop_corr_dicc['pval'].append(pval)
            # If the current population is mxl.
            if pop == 'MXL':
                # If only two ancestry proprtions are present.
                if np.unique(anc_props).size == 2:
                    # Intialize data for plotting.
                    anc_pos = [0, 1]
                    data = [
                        df[df[anc] == 0]['REPEAT_COPIES'].values,
                        df[df[anc] == 0.5]['REPEAT_COPIES'].values,
                    ]
                # Else all three are present.
                else:
                    # Intialize data for plotting.
                    anc_pos = [0, 1, 2]
                    data = [
                        df[df[anc] == 0]['REPEAT_COPIES'].values,
                        df[df[anc] == 0.5]['REPEAT_COPIES'].values,
                        df[df[anc] == 1]['REPEAT_COPIES'].values,
                    ]
                # Fill the dictionary.
                mxl_dicc[anc]['data'] = data
                mxl_dicc[anc]['anc_pos'] = anc_pos
                mxl_dicc[anc]['corr'] = [rho, pval]
    # Convert the dictionaries to a dataframes.
    anc_prop_corr_df = pd.DataFrame(anc_prop_corr_dicc)
    anc_prop_corr_df.rename(
        columns={
            'pop': 'Population', 'anc': 'Ancestry Componenet',
            'rho': "Spearman's $\\rho$", 'pval': '$P-value$',
        }, inplace=True
    )
    anc_props_summary_df = pd.DataFrame(anc_props_summary)
    # Rename all the columns to look pretty.
    amr_short_read_repeats_export_df = amr_short_read_repeats_df.rename(
        columns={
            'IND': 'Individual', 'POP': 'Population', 'SUPERPOP': 'Super Population',
            'REPEAT_COPIES': 'Repeat Copies',
            'N_SR_TRACTS': 'Number of Tracts Overlapping the Repeat Region',
            'OUT': 'Repeat Copies > 95th Percentile',
            'NAT': 'Indigenous American Anc. Prop.', 'EUR': 'European Anc. Prop.',
            'AFR': 'African Anc. Prop.', 'UNK': 'Unknown Anc. Prop.',
        }
    )
    anc_props_summary_df.rename(
        columns={
            'pop': 'Population', 'anc': 'Ancestry Componenet',
            'mean': r'Ancestry Proportion $\left( \mu \right)$',
            'std': r'Ancestry Proportion $\left( \sigma \right)$',
            'sem': r'Ancestry Proportion $\left( SEM \right)$',
            '95_ci': r'Ancestry Proportion $\left( \pm CI_{95\%} \right)$',
        }, inplace=True
    )
    # Export the dataframe as a csv.
    amr_short_read_repeats_export_df.to_csv('./dataframes/tgp_muc19_short_read_vntr_info_with_amr_ancestry_proportions.csv.gz', index=False)
    anc_prop_corr_df.to_csv('./dataframes/tgp_short_read_vntr_amr_ancestry_proportions_correlations.csv.gz', index=False)
    anc_props_summary_df.to_csv('./dataframes/tgp_short_read_vntr_amr_ancestry_proportions_summary.csv.gz', index=False)
    return anc_props_summary_df, anc_prop_corr_df, mxl_dicc

# Define a function to plot repeat count correlations for mxl.
def plot_vntr_mxl_ancestry_proportions(plot_dicc):
    # Intialize the x-axis positions, labels, and legend.
    x_pos = np.array([0, 1, 2])
    x_ticklabs = ['0%', '50%', '100%']
    x_labs = [
        'Indigenous American Anc. Props. in MXL Inds.',
        'European Anc. Props. in MXL Inds.',
        'African Anc. Props. in MXL Inds.',
    ]
    legend_info = [('upper left', (0, 1)), ('upper right', (1, 1)), ('upper right', (1, 1))]
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 3, figsize=(9, 3.5), dpi=300,
        facecolor='white',
        sharex=False, sharey=True,
        constrained_layout=True,
    )
    # For each subplot.
    for i, plot_key in enumerate(plot_dicc):
        # Unpack.
        repeat_counts = plot_dicc[plot_key]['data']
        anc_pos = plot_dicc[plot_key]['anc_pos']
        rho, pval = plot_dicc[plot_key]['corr']
        # Plot the alternating background.
        axes[i].axvspan(x_pos[1] - 0.5, x_pos[1] + 0.5, facecolor='black', alpha=0.05, zorder=1)
        # Plot the distributions.
        vp = axes[i].violinplot(repeat_counts, positions=anc_pos, widths=0.9, showextrema=False)
        # Adjust the violin plot color.
        for pc in vp['bodies']:
            pc.set_facecolor('grey')
            pc.set_edgecolor('grey')
            pc.set_alpha(0.75)
        # For every group.
        for j, grp in enumerate(repeat_counts):
            # Compute 95% ci.
            _, grp_ci = sem_ci_of_mean(grp)
            # Plot the mean and 95% ci.
            axes[i].errorbar(
                x_pos[j], np.mean(grp), yerr=grp_ci,
                color='black', marker='o',markersize=5,
                elinewidth=1.5, capthick=1.5, capsize=5, zorder=5,
            )
            # Generate some jitter to the x-axis.
            jitter = np.random.normal(x_pos[j], 0.1, size=grp.size)
            # Plot the points!
            axes[i].scatter(
                jitter, grp, facecolor='white', marker='X', s=25,
                edgecolor='black', linewidths=0.25, zorder=4,
            )
        # Add a figure lgend.
        axes[i].legend(
            handles=[
                Line2D([0], [0], color='none', linestyle='none', label=f"Spearman's $\\rho = {rho:.3f}$"),
                Line2D([0], [0], color='none', linestyle='none', label=fr'$P-value = {pval:.3e}$')
            ], 
            loc=legend_info[i][0], bbox_to_anchor=legend_info[i][1], handlelength=0, handletextpad=0,
            frameon=True, fancybox=True, shadow=True, fontsize=8,
        )
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Adjust the axes limits.
        axes[i].set_xlim(-0.5, 2.5)
        # Set the axes labels.
        axes[i].set_ylabel('Repeat Copies')
        axes[i].set_xlabel(x_labs[i])
        # Set the x-axis tick positions and labels.
        axes[i].set_xticks(x_pos)
        axes[i].set_xticklabels(x_ticklabs)
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
    # Export the plot.
    plt.savefig(
        './supp_figures/png/mxl_short_read_vntr_ancestry_proportions_correlations.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/mxl_short_read_vntr_ancestry_proportions_correlations.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/mxl_short_read_vntr_ancestry_proportions_correlations.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to plot repeat count correlations between sequencing technologies.
def plot_vntr_dataset_correlation():
    # Load in the repeats dataframe.
    df = pd.read_csv('../vntr/hg38_muc19_long_vs_short_read_repeats.csv')
    # Extract the data.
    sr_copies = df['SR_REPEAT_COPIES'].values
    lr_copies = df['LR_REPEAT_COPIES'].values
    # Perform a least squares regression.
    slope, intercept, rho, pval, se = stats.linregress(lr_copies, sr_copies)
    # Compute the line of best fir for the linear regression model.
    line_of_best_fit = (slope * lr_copies) + intercept
    # Sort the values based on the the average long read repeat length.
    sorted_indices = np.argsort(lr_copies)
    sorted_lr_copies = lr_copies[sorted_indices]
    sorted_line_of_best_fit = line_of_best_fit[sorted_indices]
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the figure and axes.
    fig = plt.figure(
        figsize=(6, 4), dpi=300,
        facecolor='white', constrained_layout=True,
    )
    ax = fig.add_subplot(111)
    # Plot the points!
    ax.scatter(
        lr_copies, sr_copies, facecolor='none', marker='o', s=55,
        edgecolor='grey',
    )
    # Plot the line of best fit.
    ax.plot(sorted_lr_copies, sorted_line_of_best_fit, 'k--', lw=1)
    # Add a figure lgend.
    ax.legend(
        handles=[
            Line2D([0], [0], color='none', linestyle='none', label=f"Pearson's $\\rho = {rho:.3f}$"),
            Line2D([0], [0], color='none', linestyle='none', label=fr'$P-value = {pval:.3e}$')
        ], 
        loc='upper left', bbox_to_anchor=(0, 1), handlelength=0, handletextpad=0,
        frameon=True, fancybox=True, shadow=True, fontsize=8,
    )
    # Set the axes labels.
    ax.set_xlabel('Average Repeat Copies (Long-Read Data)')
    ax.set_ylabel('Repeat Copies (Short-Read Data)')
    # Export the plot.
    plt.savefig(
        './supp_figures/png/vntr_short_read_v_long_read_correlation.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/vntr_short_read_v_long_read_correlation.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/vntr_short_read_v_long_read_correlation.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return



#####################
### INTROGRESSION ###
#####################

# Define a function to compute U_ABC(w, x, y).
def U_ABC(gt, A, B, C, w, x, y):
    # Calculate alternative allele frequencies.
    A_alt_freq = calc_alt_freqs(gt.take(A, axis=1))
    B_alt_freq = calc_alt_freqs(gt.take(B, axis=1))
    C_alt_freq = calc_ind_alt_freqs(gt.take(C, axis=1))
    # Intialize conditions.
    A_ref = (1 - A_alt_freq) < w
    A_alt = A_alt_freq < w
    B_ref = (1 - B_alt_freq) > x
    B_alt = B_alt_freq > x
    C_ref = (1 - C_alt_freq) == y
    C_alt = C_alt_freq == y
    # Determine the U sites.
    U_ref = A_ref & B_ref & C_ref
    U_alt = A_alt & B_alt & C_alt
    U_sites = U_ref | U_alt
    return U_sites.sum()

# Define a function to compute observed U30 - statistics.
def tgp_U30(gt):
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a list of focal populations.
    ooa_list = [
        'MXL', 'PEL', 'CLM', 'PUR', # AMR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
    ]
    # Intialize a dictionary to store all sample indicies.
    idx_dicc = {
        'ARC': np.array([2347]),
        'AFR': meta_df[meta_df['SUPERPOP'] == 'AFR'].index.values,
    }
    # For every OOA population.
    for pop in ooa_list:
        # Update the dictionary.
        idx_dicc[pop] = meta_df[meta_df['POP'] == pop].index.values
    # Intialize a list to store the reseults.
    U30_list = []
    # For every OOA population.
    for pop in ooa_list:
        # Compute the U30 statistic.
        U30_list.append(U_ABC(
            gt=gt, A=idx_dicc['AFR'], B=idx_dicc[pop],
            C=idx_dicc['ARC'], w=0.01, x=0.3, y=1,
        ))
    return np.array(U30_list)

# Define a function to load the windowed U30 results.
def load_U30_windows(window_size):
    # Load the window results.
    U30_wind_list = [np.loadtxt(f'../muc19_results/tgp_den_masked_no_aa/u30_afr_b_den_chr{chrom}_{window_size}kb.txt.gz') for chrom in range(1, 23)]
    # Concatenate all the chromosomes.
    U30_wind_mat = np.concatenate(U30_wind_list)
    return U30_wind_mat

# Deifne a function to compile the U30 results summary.
def compile_U30_summary(obs_vals, wind_mat, window_size):
    # Intialize a list of focal populations.
    ooa_list = [
        'MXL', 'PEL', 'CLM', 'PUR', # AMR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
    ]
    spop_list = [
        'AMR', 'AMR', 'AMR', 'AMR',
        'SAS', 'SAS', 'SAS', 'SAS', 'SAS',
        'EAS', 'EAS', 'EAS', 'EAS', 'EAS',   
        'EUR', 'EUR', 'EUR', 'EUR', 'EUR',
    ]
    # Intialize the p-value correction dictionary.
    pval_dicc = {
        742: {0: '<3.139e-04', 1: '>0.999686'},
        72: {0: '<3.284e-05', 1: '>0.999967'},
    }
    # Intialize a dictionary to store the results.
    df_dicc = {
        'spop': [], 'pop': [], 'obs': [],
        'wind_m': [], 'wind_s': [],
        'wind_se': [], 'wind_ci': [],
        'wind_p': [],
    }
    # Load the good window indicies.
    wind_idx = load_esl_qc_windows_idx(prefix='tgp_den_masked_no_aa', window_size=window_size)
    # For every population.
    for i, pop in enumerate(ooa_list):
        # Grab the windowed and observed results.
        winds = wind_mat[:, i][wind_idx]
        obs = obs_vals[i]
        # Determine the mean and standard deviation for non-overlapping windows.
        wind_m = np.nanmean(winds)
        wind_s = np.nanstd(winds)
        # Determine the sem and 95% ci.
        wind_se, wind_ci = sem_ci_of_mean(winds)
        # Compute the p-value.
        wind_p = (np.count_nonzero(obs <= winds) / np.sum(~np.isnan(winds)))
        # Update the dataframe dictionary.
        df_dicc['spop'].append(spop_list[i])
        df_dicc['pop'].append(pop)
        df_dicc['obs'].append(int(obs))
        df_dicc['wind_m'].append(wind_m)
        df_dicc['wind_s'].append(wind_s)
        df_dicc['wind_se'].append(wind_se)
        df_dicc['wind_ci'].append(wind_ci)
        df_dicc['wind_p'].append(wind_p)
    # Convert the dictionary to a dataframe.
    U_df = pd.DataFrame(df_dicc)
    # Adjust p-values of 0 and 1 by the number of permutations.
    U_df['wind_p'] = np.where(U_df['wind_p'] == 0, pval_dicc[window_size][0], U_df['wind_p'])
    U_df['wind_p'] = np.where(U_df['wind_p'] == '1.0', pval_dicc[window_size][1], U_df['wind_p'])
    # Rename all the columns to look pretty.
    U_df.rename(
        columns={
            'spop': 'Super Population',
            'pop': r'Population $\left( B \right)$',
            'obs': f'Focal {window_size}kb Region '+r'$\left( U_{AFR:B:DEN}(1\%, 30\%, 100\%) \right)$',
            'wind_m': fr'{window_size}kb Non-overlapping Windows $\left( \mu \right)$',
            'wind_s': fr'{window_size}kb Non-overlapping Windows $\left( \sigma \right)$',
            'wind_se': fr'{window_size}kb Non-overlapping Windows $\left( SEM \right)$',
            'wind_ci': fr'{window_size}kb Non-overlapping Windows $\left( \pm CI_{{95\%}} \right)$',
            'wind_p': r'$P-value$',
        }, inplace=True
    )
    # Export the dataframe as a csv.
    U_df.to_csv(f'./dataframes/u30_afr_b_den_{window_size}kb.csv.gz', index=False)
    return U_df

# Define a function to plot the U30 windows for mxl.
def plot_U30_afr_mxl_den(obs_742kb, winds_742kb, obs_72kb, winds_72kb):
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Load window indicies of, comparable effective sequence length.
    wind_idx_742kb = load_esl_qc_windows_idx(prefix='tgp_den_masked_no_aa', window_size=742)
    wind_idx_72kb = load_esl_qc_windows_idx(prefix='tgp_den_masked_no_aa', window_size=72)
    # Determine the significance threshold.
    sig_thresh = 0.05 
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 2, figsize=(7, 3.5), dpi=300,
        facecolor='white',
        sharex=False, sharey=False,
        constrained_layout=True,
    )
    # For each subplot.
    for i, tup in enumerate([
        (obs_742kb, winds_742kb, wind_idx_742kb, 742), (obs_72kb, winds_72kb, wind_idx_72kb, 72),
    ]):
        # Unpack.
        obs_vals, wind_mat, wind_idx, window_size = tup
        # Extract the observed value and windowed distribution.
        obs = obs_vals[0]
        dist = wind_mat[:, 0][wind_idx]
        # Compute the p-value.
        pval = (np.count_nonzero(obs <= dist) / np.sum(~np.isnan(dist)))
        # If the p-value is significant.
        if pval < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (obs, '*', 25, 0.49)
        # Else, it is not significant.
        else:
            # Update the dictionary.
            sig_dicc[i] = (obs, r'$ns$', 10, 0.5225)
        # Plot the distribution.
        axes[i].hist(
            dist, bins=np.arange(0, 145, 5),
            histtype='stepfilled', color='gray',
        )
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlabel(fr'$U_{{AFR:MXL:DEN}}(1\%, 30\%, 100\%)$ per {window_size}kb Window')
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        './supp_figures/png/u30_afr_mxl_den_windows.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/u30_afr_mxl_den_windows.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/u30_afr_mxl_den_windows.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot!
    plt.show()
    return

# Define a funct to compute Q95_ABC(w, y).
def Q95_ABC(gt, A, B, C, w, y):
    # Calculate alternative allele frequencies.
    A_alt_freq = calc_alt_freqs(gt.take(A, axis=1))
    B_alt_freq = calc_alt_freqs(gt.take(B, axis=1))
    C_alt_freq = calc_ind_alt_freqs(gt.take(C, axis=1))
    # Intialize conditions.
    A_ref = (1 - A_alt_freq) < w
    A_alt = A_alt_freq < w
    C_ref = (1 - C_alt_freq) == y
    C_alt = C_alt_freq == y
    Q_ref = A_ref & C_ref
    Q_alt = A_alt & C_alt
    # If there are no archaic sites.
    if (Q_ref | Q_alt).sum() == 0:
        return 0
    # Else there are archaic sites.
    else:
        # Subset the archaic alleles.
        Q_freqs = np.concatenate([B_alt_freq[Q_alt], (1 - B_alt_freq)[Q_ref]])
        # Find the 95% percentile.
        Q95_thresh = np.percentile(Q_freqs, 95)
        return Q95_thresh

# Define a function to compute the observed Q95 and U20/U30 statistics.
def tgp_Q95_U30(gt):
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a list of focal populations.
    ooa_list = [
        'MXL', 'PEL', 'CLM', 'PUR', # AMR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
    ]
    # Intialize a dictionary to store all sample indicies.
    idx_dicc = {
        'ARC': np.array([2347]),
        'AFR': meta_df[meta_df['SUPERPOP'] == 'AFR'].index.values,
    }
    # For every OOA population.
    for pop in ooa_list:
        # Update the dictionary.
        idx_dicc[pop] = meta_df[meta_df['POP'] == pop].index.values
    # Intilize a dictionary to store the results.
    Q95_U_dicc = {}
    # For every population.
    for pop in ooa_list:
        # Compute Q95.
        Q95 = Q95_ABC(
            gt, idx_dicc['AFR'], idx_dicc[pop], idx_dicc['ARC'], 0.01, 1,
        )
        # Compute the U_ABC stats
        U30 = U_ABC(gt, idx_dicc['AFR'], idx_dicc[pop], idx_dicc['ARC'], 0.01, 0.3, 1)
        Q95_U_dicc[pop] = np.array([Q95, U30])
    return Q95_U_dicc

# Define a function to load the windowed Q95 and U20/U30 results.
def load_Q95_U30_windows(window_size):
    # Intialize a list of focal populations.
    ooa_list = [
        'MXL', 'PEL', 'CLM', 'PUR', # AMR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
    ]
    # Intilize a dictionary to store the results.
    Q95_U_dicc = {}
    # Load the good window indicies.
    wind_idx = load_esl_qc_windows_idx(prefix=f'tgp_den_masked_no_aa', window_size=window_size)
    # Mask the region partially overlapping the focal region for clarity when plotting.
    if window_size == 72:
        wind_idx = wind_idx[wind_idx != 26235]
    elif window_size == 742:
        wind_idx = wind_idx[wind_idx != 2581]
    # For every population.
    for pop in ooa_list:
        # Load the windowed results.
        Q95_U_winds = [
            np.loadtxt(f'../muc19_results/tgp_den_masked_no_aa/q95_u20_u30_afr_{pop.lower()}_den_chr{chrom}_{window_size}kb.txt.gz') for chrom in range(1, 23)
        ]
        # Update the dictionary.
        Q95_U_dicc[pop] = np.concatenate(Q95_U_winds)[wind_idx]
    return Q95_U_dicc

# Define a function to plot the Q95 results.
def plot_Q95_x_U30_by_spop(obs_dicc, wind_dicc, window_size):
    # Intialize dictionaries for plotting.
    ooa_dicc = {
        'AMR': ['MXL', 'PEL', 'CLM', 'PUR'],
        'SAS': ['BEB', 'STU', 'ITU', 'PJL', 'GIH'],
        'EAS': ['CHB', 'KHV', 'CHS', 'JPT', 'CDX'],   
        'EUR': ['TSI', 'CEU', 'IBS', 'GBR', 'FIN'],
    }
    spop_colors = {
        'AMR': '#009E73', 'SAS': '#CC79A7', 'EAS': '#0072B2', 'EUR': '#D55E00',
    }
    pop_markers = {
        'MXL': 'o', 'PEL': '^', 'CLM': 's', 'PUR': 'p',
        'BEB': 'o', 'STU': '^', 'ITU': 's', 'PJL': 'p', 'GIH': 'h',
        'CHB': 'o', 'KHV': '^', 'CHS': 's', 'JPT': 'p', 'CDX': 'h',
        'TSI': 'o', 'CEU': '^', 'IBS': 's', 'GBR': 'p', 'FIN': 'h',
    }
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'ytick.labelsize': 10,
        'xtick.labelsize': 10,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        2, 2, figsize=(11, 11), dpi=300,
        facecolor='white',
        sharex=True, sharey=True,
        constrained_layout=True,
    )
    # Flatten the axes.
    axes = axes.flatten()
    # For every superpopulation.
    for i, spop in enumerate(spop_colors):
        # For every population.
        for pop in ooa_dicc[spop][::-1]:
            # Plot the windowed results.
            axes[i].scatter(
                wind_dicc[pop][:, 0], wind_dicc[pop][:, -1], color=spop_colors[spop],
                marker=pop_markers[pop], edgecolor='grey', linewidth=0.5, alpha=0.25, s=125,
            )
        # For every population.
        for pop in ooa_dicc[spop][::-1]:
            # Plot the obsereved results.
            axes[i].scatter(
                obs_dicc[pop][0], obs_dicc[pop][-1], color='#F0E442',
                marker=pop_markers[pop], edgecolor='grey', linewidth=0.5, s=100,
            )
    # For every superpopulation.
    for i, spop in enumerate(spop_colors):
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        axes[i].set_xlabel(fr"$Q95_{{AFR:{spop}:Denisovan}}(1\%, 100\%)$ per {window_size}kb Window")
        axes[i].set_ylabel(fr"$U_{{AFR:{spop}:Denisovan}}(1\%, 30\%, 100\%)$ per {window_size}kb Window")
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
        # Generate the legend.
        legend_handles = [
            Line2D(
                [0], [0], linestyle='none', marker=pop_markers[pop], markersize=10, markeredgewidth=0.25,
                color=spop_colors[spop], markeredgecolor='grey', label=f'{pop}'
            ) for pop in ooa_dicc[spop]
        ]
        # Add the legend.
        legend = axes[i].legend(
            handles=legend_handles, loc='upper right', bbox_to_anchor=(1, 1),
            frameon=True, fancybox=True, shadow=True, fontsize=10,
        )
        legend.set_title(spop, prop={'size': 12, 'weight': 'bold'})
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/u30_x_q95_afr_mxl_den_windows_{window_size}kb.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/u30_x_q95_afr_mxl_den_windows_{window_size}kb.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/u30_x_q95_afr_mxl_den_windows_{window_size}kb.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to summarize the Q95 results.
def compile_Q95_summary(obs_72kb, obs_742kb):
    # Intialize a list of focal populations.
    ooa_list = [
        'MXL', 'PEL', 'CLM', 'PUR', # AMR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
    ]
    spop_list = [
        'AMR', 'AMR', 'AMR', 'AMR',
        'SAS', 'SAS', 'SAS', 'SAS', 'SAS',
        'EAS', 'EAS', 'EAS', 'EAS', 'EAS',   
        'EUR', 'EUR', 'EUR', 'EUR', 'EUR',
    ]
    # Intialize a dictionary to store the results.
    df_dicc = {
        'spop': [], 'pop': [], 'obs_72kb': [], 'obs_742kb': [],
    }
    # For every population.
    for i, pop in enumerate(ooa_list):
        # Update the dataframe dictionary.
        df_dicc['spop'].append(spop_list[i])
        df_dicc['pop'].append(pop)
        df_dicc['obs_72kb'].append(obs_72kb[pop][0])
        df_dicc['obs_742kb'].append(obs_742kb[pop][0])
    # Convert the dictionary to a dataframe.
    Q95_df = pd.DataFrame(df_dicc)
    # Rename all the columns to look pretty.
    Q95_df.rename(
        columns={
            'spop': 'Super Population',
            'pop': r'Population $\left( B \right)$',
            'obs_72kb': 'Focal 72kb Region '+r'$\left( Q95_{AFR:B:DEN}(1\%, 100\%) \right)$',
            'obs_742kb': 'Focal 742kb Region '+r'$\left( Q95_{AFR:B:DEN}(1\%, 100\%) \right)$',
        }, inplace=True
    )
    # Export the dataframe as a csv.
    Q95_df.to_csv(f'./dataframes/q95_u30_afr_b_den_72kb_and_742kb.csv.gz', index=False)
    return Q95_df

# Define a site pattern function.
def site_patterns(p1, p2, p3):
    # Calculate site pattern counts.
    abba = np.sum((1 - p1) * (p2) * (p3))
    baba = np.sum((p1) * (1 - p2) * (p3))
    bbaa = np.sum((p1) * (p2) * (1 - p3))
    baaa = np.sum((p1) * (1 - p2) * (1 - p3))
    abaa = np.sum((1 - p1) * (p2) * (1 - p3))
    aaba = np.sum((1 - p1) * (1 - p2) * (p3))
    return abba, baba, bbaa, baaa, abaa, aaba

# Define a function to calculate site patterns between the mxl and archaics.
def mxl_arc_dplus(gt, p2_ind):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Determine the sample indicies.
    p1_idx = tgp_meta_df[tgp_meta_df['POP'] == 'YRI'].index.values
    p2_idx = tgp_meta_df[tgp_meta_df['IND'] == p2_ind].index.values
    p3_idx = np.array([2347])
    # Determine the sample list.
    samp_list = np.concatenate([p1_idx, p2_idx, p3_idx])
    # Determine the indicies where all samples are called.
    called_mask = (gt.take(samp_list, axis=1).is_called() == True).all(axis=1)
    # If there are no sites called between all three samples...
    if (called_mask.sum() == 0):
        # Set the results to 0 since we are iterating over QC'ed regions.
        results = np.zeros(6)
    # Else...
    else:
        # Determine the indicies where we have varibale sites.
        var_mask = gt.take(samp_list, axis=1).compress(called_mask, axis=0).count_alleles().is_segregating()
        # If there are no variable sites...
        if (var_mask.sum() == 0):
            # Set the results to 0 since we are iterating over QC'ed regions.
            results = np.zeros(6)
        # Else...
        else:
            # Calculate the alternative allele frequencies.
            p1_alt_freqs = calc_alt_freqs(gt.take(p1_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p2_alt_freqs = calc_alt_freqs(gt.take(p2_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p3_alt_freqs = calc_ind_alt_freqs(gt.take(p3_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            anc_freqs = calc_ind_alt_freqs(gt.take([-1], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            # Polarize the samples.
            p1_der_freqs = np.where(anc_freqs == 1, np.abs(p1_alt_freqs - 1), p1_alt_freqs)
            p2_der_freqs = np.where(anc_freqs == 1, np.abs(p2_alt_freqs - 1), p2_alt_freqs)
            p3_der_freqs = np.where(anc_freqs == 1, np.abs(p3_alt_freqs - 1), p3_alt_freqs)
            # Calculate the site pattern counts.
            abba, baba, bbaa, baaa, abaa, aaba = site_patterns(
                p1_der_freqs, p2_der_freqs, p3_der_freqs,
            )
            # Calculate D+.
            dplus = (((abba - baba) + (baaa - abaa)) / ((abba + baba) + (baaa + abaa)))
            # Intialize a results array.
            results = np.array([abba, baba, bbaa, baaa, abaa, aaba, dplus])
    return results

# Define a funtion to load the mxl observed reuslts.
def build_mxl_arc_scenario_obs_dicc(gt_dicc):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a list of focal individuals.
    p2_list = ['NA19664']
    # Intialize a dictionary to store all of the results.
    obs_dicc = {}
    # For every p2.
    for p2 in p2_list:
        # For every p3.
        for p3 in ['DEN', 'ALT', 'CHA', 'VIN']:
            # Intialize the configuration.
            config = f'{p2}-{p3}'
            # Load the window results.
            obs_dicc[config] = mxl_arc_dplus(gt_dicc[p3], p2)
    return obs_dicc


# Define a function to load the tgp site pattern windowed results.
def load_mxl_dplus_windows(p2_ind, p3_arc, window_size):
    # Intialize a list to store the results.
    dplus_list = []
    # Intialize the file path prefix.
    path_prefix = f'../muc19_results/tgp_{p3_arc.lower()}_masked_aa/yri_{p2_ind.lower()}_{p3_arc.lower()}'
    # For all chromosomes.
    for chrom in range(1, 23):
        # Load the site pattern results.
        site_patterns = np.loadtxt(f'{path_prefix}_chr{chrom}_{window_size}kb.txt.gz')
        # Extract the site pattern arrays.
        abba = site_patterns[:, 0]
        baba = site_patterns[:, 1]
        bbaa = site_patterns[:, 2]
        baaa = site_patterns[:, 3]
        abaa = site_patterns[:, 4]
        aaba = site_patterns[:, 5]
        # Calculate D+.
        dplus = (((abba - baba) + (baaa - abaa)) / ((abba + baba) + (baaa + abaa)))
        # Append the list
        dplus_list.append(dplus)
    return np.concatenate(dplus_list)

# Define a funtion to load the tgp window reuslts.
def build_mxl_arc_scenario_window_dicc(window_size):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a list of focal individuals.
    p2_list = ['NA19664']
    # Intialize a dictionary to store all of the results.
    wind_dicc = {}
    # For every p2.
    for p2 in p2_list:
        # For every p3.
        for p3 in ['DEN', 'ALT', 'CHA', 'VIN']:
            # Intialize the configuration.
            config = f'{p2}-{p3}'
            # Load the window results.
            wind_dicc[config] = load_mxl_dplus_windows(p2_ind=p2, p3_arc=p3, window_size=window_size)
    return wind_dicc

# Define a function to compile the archaic snp density table.
def compile_mxl_arc_dplus_scenario_summary(obs_dicc, wind_dicc, window_size):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a list of focal individuals.
    p2_list = ['NA19664']
    # Intialize an archaic information dictionary.
    arc_dicc = {
        'DEN': {'lab': 'Denisovan', 'export': 'dplus_yri_mxl_den'},
        'ALT': {'lab': 'Altai Nean.', 'export': 'dplus_yri_mxl_alt'},
        'CHA': {'lab': 'Chagyrskaya Nean.', 'export': 'dplus_yri_mxl_cha'},
        'VIN': {'lab': 'Vindija Nean.', 'export': 'dplus_yri_mxl_vin'},
    }
    # Intialize dictionaries to store the results.
    df_dicc = {
        'DEN': {'p1': [], 'p2': [], 'p3': [],
                'abba': [], 'baba': [], 'bbaa': [],
                'baaa': [], 'abaa': [], 'aaba': [],
                'dplus': [], 'wind_m': [], 'wind_s': [],
                'wind_se': [], 'wind_ci': [], 'wind_p': []},
        'ALT': {'p1': [], 'p2': [], 'p3': [],
                'abba': [], 'baba': [], 'bbaa': [],
                'baaa': [], 'abaa': [], 'aaba': [],
                'dplus': [], 'wind_m': [], 'wind_s': [],
                'wind_se': [], 'wind_ci': [], 'wind_p': []},
        'CHA': {'p1': [], 'p2': [], 'p3': [],
                'abba': [], 'baba': [], 'bbaa': [],
                'baaa': [], 'abaa': [], 'aaba': [],
                'dplus': [], 'wind_m': [], 'wind_s': [],
                'wind_se': [], 'wind_ci': [], 'wind_p': []},
        'VIN': {'p1': [], 'p2': [], 'p3': [],
                'abba': [], 'baba': [], 'bbaa': [],
                'baaa': [], 'abaa': [], 'aaba': [],
                'dplus': [], 'wind_m': [], 'wind_s': [],
                'wind_se': [], 'wind_ci': [], 'wind_p': []},
        
    }
    # For every p3.
    for p3 in ['DEN', 'ALT', 'CHA', 'VIN']:
        # Load the good window indicies.
        wind_idx = load_esl_qc_windows_idx(prefix=f'tgp_{p3.lower()}_masked_aa', window_size=window_size)
        # For every p2.
        for p2 in p2_list:
            # Intialize the configuration.
            config = f'{p2}-{p3}'
            # Grab the observed and windowed results.
            winds = wind_dicc[config][wind_idx]
            abba, baba, bbaa, baaa, abaa, aaba, dplus = obs_dicc[config]
            # Determine the mean and standard deviation for non-overlapping windows.
            wind_m = np.nanmean(winds)
            wind_s = np.nanstd(winds)
            # Determine the sem and 95% ci.
            wind_se = wind_s / np.sqrt(np.sum(~np.isnan(winds)))
            wind_ci = wind_s * stats.norm.ppf(0.95)
            # Compute the p-value.
            wind_p = stats.norm.sf(x=dplus, loc=0, scale=wind_s)
            # Update the dataframe dictionary.
            df_dicc[p3]['p1'].append('YRI')
            df_dicc[p3]['p2'].append(p2)
            df_dicc[p3]['p3'].append(arc_dicc[p3]['lab'])
            df_dicc[p3]['abba'].append(abba)
            df_dicc[p3]['baba'].append(baba)
            df_dicc[p3]['bbaa'].append(bbaa)
            df_dicc[p3]['baaa'].append(baaa)
            df_dicc[p3]['abaa'].append(abaa)
            df_dicc[p3]['aaba'].append(aaba)
            df_dicc[p3]['dplus'].append(dplus)
            df_dicc[p3]['wind_m'].append(wind_m)
            df_dicc[p3]['wind_s'].append(wind_s)
            df_dicc[p3]['wind_se'].append(wind_se)
            df_dicc[p3]['wind_ci'].append(wind_ci)
            df_dicc[p3]['wind_p'].append(wind_p)
    # Intialize a list to store the dataframes.
    df_list = []
    # For every arcahic.
    for arc in ['DEN', 'ALT', 'CHA', 'VIN']:
        # Convert the dictionary to a dataframe.
        dplus_df = pd.DataFrame(df_dicc[arc]).sort_values('wind_p').reset_index(drop=True)
        # Rename all the columns to look pretty.
        dplus_df.rename(
            columns={
                'p1': r'$P1$', 'p2': r'$P2$', 'p3': r'$P3$',
                'abba': r'$ABBA$', 'baba': r'$BABA$', 'bbaa': r'$BBAA$',
                'baaa': r'$BAAA$', 'abaa': r'$ABAA$', 'aaba': r'$AABA$',
                'dplus': f'Focal {window_size}kb Region '+r'$\left( D+ \right)$',
                'wind_m': fr'{window_size}kb Non-overlapping Windows $\left( \mu \right)$',
                'wind_s': fr'{window_size}kb Non-overlapping Windows $\left( \sigma \right)$',
                'wind_se': fr'{window_size}kb Non-overlapping Windows $\left( SEM \right)$',
                'wind_ci': fr'{window_size}kb Non-overlapping Windows $\left( \pm CI_{{95\%}} \right)$',
                'wind_p': r'$P-value$',
            }, inplace=True
        )
        # Update the list.
        df_list.append(dplus_df)
    # Concatenate.
    all_dplus_df = pd.concat(df_list)
    # Export the dataframe as a csv.
    all_dplus_df.to_csv(f'./dataframes/dplus_yri_mxl_arcs_{window_size}kb.csv.gz', index=False)
    return pd.concat(df_list)

# Define a function to plot the focal mxl sequence divergence results.
def plot_mxl_arc_dplus(obs_dicc, wind_dicc, window_size):
    # Intialize a dictionary of archaic labels.
    x_labs = {
        'DEN': fr'$D+\left(\left(YRI, NA19664 \right), Denisovan \right)$ per {window_size}kb Window',
        'ALT': fr'$D+\left(\left(YRI, NA19664 \right), Altai \; Nean. \right)$ per {window_size}kb Window',
        'CHA': fr'$D+\left(\left(YRI, NA19664 \right), Chagyrskaya \; Nean. \right)$ per {window_size}kb Window',
        'VIN': fr'$D+\left(\left(YRI, NA19664 \right), Vindija \; Nean. \right)$ per {window_size}kb Window',
    }
    # Determine the significance threshold.
    sig_thresh = 0.05
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        2, 2, figsize=(10, 5), dpi=300,
        facecolor='white',
        sharex=True, sharey=True,
        constrained_layout=True,
    )
    # Flatten the axes.
    axes = axes.flatten()
    # For every archaic.
    for i, p3 in enumerate(['DEN', 'ALT', 'CHA', 'VIN']):
        # Load the good window indicies.
        wind_idx = load_esl_qc_windows_idx(prefix=f'tgp_{p3.lower()}_masked_aa', window_size=window_size)
        # Grab the observed and windowed results.
        dist = wind_dicc[f'NA19664-{p3}'][wind_idx]
        obs = obs_dicc[f'NA19664-{p3}'][-1]
        # Compute the p-value.
        pval = stats.norm.sf(x=obs, loc=0, scale=np.nanstd(dist))
        # If the p-value is significant.
        if pval < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (obs, '*', 25, 0.49)
        # Else, it is not significant
        else:
            # Update the dictionary.
            sig_dicc[i] = (obs, r'$ns$', 10, 0.5225)
        # Plot the distribution.
        axes[i].hist(
            dist, bins=np.linspace(-1, 1, 100),
            histtype='stepfilled', color='gray',
        )
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlabel(x_labs[p3])
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/dplus_mxl_arcs_{window_size}kb.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/dplus_mxl_arcs_{window_size}kb.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/dplus_mxl_arcs_{window_size}kb.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot!
    plt.show()
    return

# Define a function to filter and calculate site pattern counts for the archaics.
def arc_dplus(gt, p1, p2, p3):
    # Intialize a dictionary of indicies.
    idx_dicc = {
        'ALT': 0, 'CHA': 1,
        'VIN': 2, 'DEN': 3,
    }
    # Determine the sample list.
    samp_list = [idx_dicc[p1], idx_dicc[p2], idx_dicc[p3]]
    # Determine the indicies where all samples are called.
    called_mask = (gt.take(samp_list, axis=1).is_called() == True).all(axis=1)
    # If there are no sites called between all three samples...
    if (called_mask.sum() == 0):
        # Set the results to 0 since we are iterating over QC'ed regions.
        results = np.zeros(6)
    # Else...
    else:
        # Determine the indicies where we have varibale sites.
        var_mask = gt.take(samp_list, axis=1).compress(called_mask, axis=0).count_alleles().is_segregating()
        # If there are no variable sites...
        if (var_mask.sum() == 0):
            # Set the results to 0 since we are iterating over QC'ed regions.
            results = np.zeros(6)
        # Else...
        else:
            # Calculate the alternative allele frequencies.
            p1_alt_freqs = calc_ind_alt_freqs(gt.take([idx_dicc[p1]], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p2_alt_freqs = calc_ind_alt_freqs(gt.take([idx_dicc[p2]], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p3_alt_freqs = calc_ind_alt_freqs(gt.take([idx_dicc[p3]], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            anc_freqs = calc_ind_alt_freqs(gt.take([-1], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            # Polarize the samples.
            p1_der_freqs = np.where(anc_freqs == 1, np.abs(p1_alt_freqs - 1), p1_alt_freqs)
            p2_der_freqs = np.where(anc_freqs == 1, np.abs(p2_alt_freqs - 1), p2_alt_freqs)
            p3_der_freqs = np.where(anc_freqs == 1, np.abs(p3_alt_freqs - 1), p3_alt_freqs)
            # Calculate the site pattern counts.
            abba, baba, bbaa, baaa, abaa, aaba = site_patterns(
                p1_der_freqs, p2_der_freqs, p3_der_freqs,
            )
            # Calculate D+.
            dplus = (((abba - baba) + (baaa - abaa)) / ((abba + baba) + (baaa + abaa)))
            # Intialize a results array.
            results = np.array([abba, baba, bbaa, baaa, abaa, aaba, dplus])
    return results

# Define a funtion to load the archaic observed reuslts.
def build_arc_scenario_obs_dicc(gt):
    # Intialize a list of site pattern configurations.
    config_list = [
        ('ALT', 'CHA', 'DEN'),
        ('ALT', 'VIN', 'DEN'),
#         ('CHA', 'VIN', 'DEN'),
#         ('CHA', 'VIN', 'ALT'),
    ]
    # Intialize a dictionary to store all of the results.
    obs_dicc = {}
    # For every site pattern configuration
    for config in config_list:
        # Unpack.
        p1, p2, p3 = config
        # Load the window results.
        obs_dicc[config] = arc_dplus(gt=gt, p1=p1, p2=p2, p3=p3)
    return obs_dicc

# Define a function to load the archaic site pattern windowed results.
def load_arc_dplus_windows(p1_arc, p2_arc, p3_arc, window_size):
    # Intialize a list to store the results.
    dplus_list = []
    # Intialize the file path prefix.
    path_prefix = f'../muc19_results/arcs_masked_aa/{p1_arc.lower()}_{p2_arc.lower()}_{p3_arc.lower()}'
    # For all chromosomes.
    for chrom in range(1, 23):
        # Load the site pattern results.
        site_patterns = np.loadtxt(f'{path_prefix}_chr{chrom}_{window_size}kb.txt.gz')
        # Extract the site pattern arrays.
        abba = site_patterns[:, 0]
        baba = site_patterns[:, 1]
        bbaa = site_patterns[:, 2]
        baaa = site_patterns[:, 3]
        abaa = site_patterns[:, 4]
        aaba = site_patterns[:, 5]
        # Calculate D+.
        dplus = (((abba - baba) + (baaa - abaa)) / ((abba + baba) + (baaa + abaa)))
        # Append the list
        dplus_list.append(dplus)
    return np.concatenate(dplus_list)

# Define a funtion to load the archaic window reuslts.
def build_arc_scenario_window_dicc(window_size):
    # Intialize a list of site pattern configurations.
    config_list = [
        ('ALT', 'CHA', 'DEN'),
        ('ALT', 'VIN', 'DEN'),
    ]
    # Intialize a dictionary to store all of the results.
    wind_dicc = {}
    # For every site pattern configuration
    for config in config_list:
        # Unpack.
        p1, p2, p3 = config
        # Load the window results.
        wind_dicc[config] = load_arc_dplus_windows(p1_arc=p1, p2_arc=p2, p3_arc=p3, window_size=window_size)
    return wind_dicc

# Define a function to compile the archaic snp density table.
def compile_arc_dplus_scenario_summary(obs_dicc, wind_dicc, window_size):
    # Intialize a list of site pattern configurations.
    config_list = [
        ('ALT', 'CHA', 'DEN'),
        ('ALT', 'VIN', 'DEN'),
    ]
    # Intialize a dictionary of archaic labels.
    arc_dicc = {
        'DEN': 'Denisovan', 'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.', 'VIN': 'Vindija Nean.',
    }
    # Intialize dictionaries to store the results.
    df_dicc = {
        'p1': [], 'p2': [], 'p3': [],
        'abba': [], 'baba': [], 'bbaa': [],
        'baaa': [], 'abaa': [], 'aaba': [],
        'dplus': [], 'wind_m': [], 'wind_s': [],
        'wind_se': [], 'wind_ci': [], 'wind_p': [],
    }
    # Load the good window indicies.
    wind_idx = load_esl_qc_windows_idx(prefix='arcs_masked_aa', window_size=window_size)
    # For every configuration.
    for config in config_list:
        # Unpack.
        p1, p2, p3 = config
        # Grab the observed and windowed results.
        winds = wind_dicc[config][wind_idx]
        abba, baba, bbaa, baaa, abaa, aaba, dplus = obs_dicc[config]
        # Determine the mean and standard deviation for non-overlapping windows.
        wind_m = np.nanmean(winds)
        wind_s = np.nanstd(winds)
        # Determine the sem and 95% ci.
        wind_se = wind_s / np.sqrt(np.sum(~np.isnan(winds)))
        wind_ci = wind_s * stats.norm.ppf(0.95)
        # Compute the p-value.
        wind_p = stats.norm.sf(x=dplus, loc=0, scale=wind_s)
        # Update the dataframe dictionary.
        df_dicc['p1'].append(arc_dicc[p1])
        df_dicc['p2'].append(arc_dicc[p2])
        df_dicc['p3'].append(arc_dicc[p3])
        df_dicc['abba'].append(abba)
        df_dicc['baba'].append(baba)
        df_dicc['bbaa'].append(bbaa)
        df_dicc['baaa'].append(baaa)
        df_dicc['abaa'].append(abaa)
        df_dicc['aaba'].append(aaba)
        df_dicc['dplus'].append(dplus)
        df_dicc['wind_m'].append(wind_m)
        df_dicc['wind_s'].append(wind_s)
        df_dicc['wind_se'].append(wind_se)
        df_dicc['wind_ci'].append(wind_ci)
        df_dicc['wind_p'].append(wind_p)
    # Convert the dictionaries to dataframes.
    dplus_df = pd.DataFrame(df_dicc)
    # Rename all the columns to look pretty.
    dplus_df.rename(
        columns={
            'p1': r'$P1$', 'p2': r'$P2$', 'p3': r'$P3$',
            'abba': r'$ABBA$', 'baba': r'$BABA$', 'bbaa': r'$BBAA$',
            'baaa': r'$BAAA$', 'abaa': r'$ABAA$', 'aaba': r'$AABA$',
            'dplus': f'Focal {window_size}kb Region '+r'$\left( D+ \right)$',
            'wind_m': fr'{window_size}kb Non-overlapping Windows $\left( \mu \right)$',
            'wind_s': fr'{window_size}kb Non-overlapping Windows $\left( \sigma \right)$',
            'wind_se': fr'{window_size}kb Non-overlapping Windows $\left( SEM \right)$',
            'wind_ci': fr'{window_size}kb Non-overlapping Windows $\left( \pm CI_{{95\%}} \right)$',
            'wind_p': r'$P-value$',
        }, inplace=True
    )
    # Export the dataframe as a csv.
    dplus_df.to_csv(f'./dataframes/dplus_altai_neanderthal_p2_denisovan_{window_size}kb.csv.gz', index=False)
    return dplus_df

# Define a function to plot the archaic D+ results.
def plot_arc_dplus(obs_dicc, wind_dicc, window_size):
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize a list of site pattern configurations.
    config_list = [
        ('ALT', 'CHA', 'DEN'),
        ('ALT', 'VIN', 'DEN'),
    ]
    # Intialize a dictionary of archaics x-axis labels.
    arc_x_labs = {
        ('ALT', 'CHA', 'DEN'): fr'$D+\left(\left(Altai \; Nean., Chagyrskaya \; Nean. \right), Denisovan \right)$ per {window_size}kb Window',
        ('ALT', 'VIN', 'DEN'): fr'$D+\left(\left(Altai \; Nean., Vindija \; Nean. \right), Denisovan \right)$ per {window_size}kb Window',
    }
    # Intialize a bin dictionary.
    bin_dicc = {
        742: np.arange(-1, 1.02, 0.01),
        72: np.arange(-1, 1.03, 0.075),
    }
    # Load window indicies of, comparable effective sequence length.
    wind_idx = load_esl_qc_windows_idx(prefix='arcs_masked_aa', window_size=window_size)
    # Determine the significance threshold.
    sig_thresh = 0.05
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 2, figsize=(9, 3.5), dpi=300,
        facecolor='white',
        sharex=True, sharey=True,
        constrained_layout=True,
    )
    # For every site pattern configuration.
    for i, config in enumerate(config_list):
        # Grab the observed and windowed distribution.
        dist = wind_dicc[config][wind_idx]
        obs = obs_dicc[config][-1]
        # Compute the p-value.
        pval = stats.norm.sf(x=obs, loc=0, scale=np.nanstd(dist))
        # If the p-value is significant.
        if pval < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (obs, '*', 25, 0.49)
        # Else, it is not significant
        else:
            # Update the dictionary.
            sig_dicc[i] = (obs, r'$ns$', 10, 0.5225)
        # Plot the distribution.
        axes[i].hist(
            dist, bins=bin_dicc[window_size],
            histtype='stepfilled', color='gray',
        )
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlabel(arc_x_labs[config])
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/dplus_altai_neanderthal_p2_denisovan_{window_size}kb.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/dplus_altai_neanderthal_p2_denisovan_{window_size}kb.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/dplus_altai_neanderthal_p2_denisovan_{window_size}kb.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot!
    plt.show()
    return



###########
### PAP ###
###########

# Define a function to calculate pap scores.
def psuedo_ancestry_painting(gt, target, source_1, source_2):
    # Determine the alternative allele frequencies for the target individual.
    t_aaf = calc_ind_alt_freqs(gt.take(target, axis=1))
    # Create a mask for the heterozygous sites.
    het_mask = t_aaf == 0.5
    # If there are no sites heterozygous sites/
    if het_mask.sum() == 0:
        pap = np.array([np.nan, 0])
    # Determine the alternative allele frequencies for the sites of interest.
    t_aaf = t_aaf[het_mask]
    s1_aaf = calc_ind_alt_freqs(gt.take(source_1, axis=1).compress(het_mask, axis=0))
    s2_aaf = calc_ind_alt_freqs(gt.take(source_2, axis=1).compress(het_mask, axis=0))
    # Determine the sites that can be painted.
    s1_ref_s2_alt = (s1_aaf == 0) & (s2_aaf == 1) & (t_aaf == 0.5)
    s1_alt_s2_ref = (s1_aaf == 1) & (s2_aaf == 0) & (t_aaf == 0.5)
    # Compute the pap sites and total possible sites.
    pap = np.array([(s1_ref_s2_alt | s1_alt_s2_ref).sum(), het_mask.sum()])
    return pap

# Define a function to determine pap masks for plotting 
def psuedo_ancestry_painting_mask(gt, target, source_1, source_2):
    # Determine the alternative allele frequencies for the target individual.
    t_aaf = calc_ind_alt_freqs(gt.take(target, axis=1))
    # Create a mask for the heterozygous sites.
    het_mask = t_aaf == 0.5
    # Determine the alternative allele frequencies for the sites of interest.
    t_aaf = t_aaf[het_mask]
    s1_aaf = calc_ind_alt_freqs(gt.take(source_1, axis=1).compress(het_mask, axis=0))
    s2_aaf = calc_ind_alt_freqs(gt.take(source_2, axis=1).compress(het_mask, axis=0))
    # Determine the sites that can be painted.
    s1_ref_s2_alt = (s1_aaf == 0) & (s2_aaf == 1) & (t_aaf == 0.5)
    s1_alt_s2_ref = (s1_aaf == 1) & (s2_aaf == 0) & (t_aaf == 0.5)
    # Create the pap mask and pap matrix.
    pap_mask = s1_ref_s2_alt | s1_alt_s2_ref
    pap_mat = np.array([s1_aaf, t_aaf, s2_aaf])
    return pap_mask, pap_mat

# Define a function to load the observed pap scores.
def calc_pap_scores(arc_gt, tgp_gt):
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionaries of focal individuals.
    arc_idx = {
        'ALT': [0], 'CHA': [1],
        'VIN': [2], 'DEN': [3],
    }
    tgp_idx = {
        'ALT': np.array([2347]), 'CHA': np.array([2348]),
        'VIN': np.array([2349]), 'DEN': np.array([2350]),
        'MXL': tgp_df[tgp_df['IND'] == 'NA19664'].index.values,
        'YRI': tgp_df[tgp_df['IND'] == 'NA19190'].index.values,
    }
    # Intialize the pap configurations.
    arc_configs = [('CHA', 'DEN', 'ALT'), ('VIN', 'DEN', 'ALT')]
    tgp_configs = [
        ('CHA', 'MXL', 'YRI'), ('VIN', 'MXL', 'YRI'),
        ('DEN', 'MXL', 'YRI'), ('ALT', 'MXL', 'YRI'),
    ]
    # Intialize a dictionaries to store the results.
    pap_dicc, plot_dicc = {}, {}
    # For every pap configuration.
    for config in arc_configs:
        # Unpack.
        t, s1, s2 = config
        # Compute the results.
        pap_sites, het_sites = psuedo_ancestry_painting(arc_gt, arc_idx[t], arc_idx[s1], arc_idx[s2])
        pap_mask, pap_mat = psuedo_ancestry_painting_mask(arc_gt, arc_idx[t], arc_idx[s1], arc_idx[s2])
        # Fill the the results.
        pap_dicc[config] = pap_sites / het_sites
        plot_dicc[config] = {'mask': pap_mask, 'mat': pap_mat}
    # For every pap configuration.
    for config in tgp_configs:
        # Unpack.
        t, s1, s2 = config
        # Compute the results.
        pap_sites, het_sites = psuedo_ancestry_painting(tgp_gt, tgp_idx[t], tgp_idx[s1], tgp_idx[s2])
        pap_mask, pap_mat = psuedo_ancestry_painting_mask(tgp_gt, tgp_idx[t], tgp_idx[s1], tgp_idx[s2])
        # Fill the the results.
        pap_dicc[config] = pap_sites / het_sites
        plot_dicc[config] = {'mask': pap_mask, 'mat': pap_mat}
    return pap_dicc, plot_dicc

# Define a function to load pap scores in windows.
def load_pap_scores_windows(window_size):
    # Intialize a dictionary to store the results.
    pap_dicc = {}
    # Intialize the pap configurations.
    pap_configs = [
        ('arcs_masked_no_aa', 'CHA', 'DEN', 'ALT'), ('arcs_masked_no_aa', 'VIN', 'DEN', 'ALT'),
        ('tgp_arcs_masked_no_aa', 'CHA', 'MXL', 'YRI'), ('tgp_arcs_masked_no_aa', 'VIN', 'MXL', 'YRI'),
        ('tgp_arcs_masked_no_aa', 'DEN', 'MXL', 'YRI'), ('tgp_arcs_masked_no_aa', 'ALT', 'MXL', 'YRI'),
    ]
    # For every pap configuration.
    for config in pap_configs:
        # Unpack.
        data_dir, t, s1, s2 = config
        # Intialize a list.
        all_pap = []
        # For all chromosomes...
        for chrom in range(1, 23):
            # Load and append the pap data.
            all_pap.append(
                np.loadtxt(
                    f'../muc19_results/{data_dir}/{t.lower()}_{s1.lower()}_{s2.lower()}_pap_counts_chr{chrom}_{window_size}kb.txt.gz',
            ))
        # Concatenate the results.
        all_pap = np.concatenate(all_pap)
        # Fill the dictionary.
        pap_dicc[(t, s1, s2)] = all_pap[:, 0] / all_pap[:, 1]
    return pap_dicc

# Define a function to compile the pap score results per window.
def compile_pap_score_summary(obs_dicc, wind_dicc, window_size):
    # Intialize labels.
    label_dicc = {
        'DEN': 'Denisovan',
        'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.',
        'VIN': 'Vindija Nean.',
        'MXL': 'MXL (NA19664)',
        'YRI': 'YRI (NA19190)',
    }
    # Intialize dictionaries to store the results.
    df_dicc = {
        's1': [], 't': [], 's2': [], 'obs': [],
        'wind_m': [], 'wind_s': [],
        'wind_se': [], 'wind_ci': [],
        'wind_p': [],
    }
    # Load the good window indicies.
    arc_idx = load_esl_qc_windows_idx(prefix='arcs_masked_no_aa', window_size=window_size)
    tgp_idx = load_esl_qc_windows_idx(prefix='tgp_arcs_masked_no_aa', window_size=window_size)
    # Intialize the pap configurations.
    pap_configs = [
        (arc_idx, 'CHA', 'DEN', 'ALT'), (arc_idx, 'VIN', 'DEN', 'ALT'),
        (tgp_idx, 'CHA', 'MXL', 'YRI'), (tgp_idx, 'VIN', 'MXL', 'YRI'),
        (tgp_idx, 'DEN', 'MXL', 'YRI'), (tgp_idx, 'ALT', 'MXL', 'YRI'),
    ]
    # For every pap configuration.
    for config in pap_configs:
        # Unpack.
        wind_idx, t, s1, s2 = config
        # Grab the results.
        winds = wind_dicc[(t, s1, s2)][wind_idx]
        obs = obs_dicc[(t, s1, s2)]
        # Determine the mean and standard deviation for non-overlapping windows.
        wind_m = np.nanmean(winds)
        wind_s = np.nanstd(winds)
        # Determine the sem and 95% ci.
        wind_se, wind_ci = sem_ci_of_mean(winds)
        # Compute the p-value.
        wind_p = (np.count_nonzero(obs <= winds) / np.sum(~np.isnan(winds)))
        # Update the dataframe dictionary.
        df_dicc['s1'].append(label_dicc[s1])
        df_dicc['t'].append(label_dicc[t])
        df_dicc['s2'].append(label_dicc[s2])
        df_dicc['obs'].append(obs)
        df_dicc['wind_m'].append(wind_m)
        df_dicc['wind_s'].append(wind_s)
        df_dicc['wind_se'].append(wind_se)
        df_dicc['wind_ci'].append(wind_ci)
        df_dicc['wind_p'].append(wind_p)
    # Convert the dictionary to a dataframe.
    pap_df = pd.DataFrame(df_dicc)
    # Rename all the columns to look pretty.
    pap_df.rename(
        columns={
            's1': r'$Source^{1}$',
            't': r'$Target$',
            's2': r'$Source^{2}$',
            'obs': fr'Focal {window_size}kb Region $\left( PAP \text{{ Score}} \right)$',
            'wind_m': fr'{window_size}kb Non-overlapping Windows $\left( \mu \right)$',
            'wind_s': fr'{window_size}kb Non-overlapping Windows $\left( \sigma \right)$',
            'wind_se': fr'{window_size}kb Non-overlapping Windows $\left( SEM \right)$',
            'wind_ci': fr'{window_size}kb Non-overlapping Windows $\left( \pm CI_{{95\%}} \right)$',
            'wind_p': r'$P-value$',
        }, inplace=True
    )
    # Export the dataframe as a csv.
    pap_df.to_csv(f'./dataframes/pseudo_ancestry_painting_scores_{window_size}kb.csv.gz', index=False)
    return pap_df

# Define a function to plot the pap scores distribution.
def plot_pap_score_dist(obs_dicc, wind_dicc, window_size):
    # Intialize labels.
    label_dicc = {
        'DEN': 'Denisovan',
        'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.',
        'VIN': 'Vindija Nean.',
        'MXL': 'MXL',
        'YRI': 'YRI',
    }
    # Load the good window indicies.
    arc_idx = load_esl_qc_windows_idx(prefix='arcs_masked_no_aa', window_size=window_size)
    tgp_idx = load_esl_qc_windows_idx(prefix='tgp_arcs_masked_no_aa', window_size=window_size)
    # Intialize the pap configurations.
    pap_configs = [
        (arc_idx, 'CHA', 'DEN', 'ALT'), (arc_idx, 'VIN', 'DEN', 'ALT'),
        (tgp_idx, 'CHA', 'MXL', 'YRI'), (tgp_idx, 'VIN', 'MXL', 'YRI'),
        (tgp_idx, 'DEN', 'MXL', 'YRI'), (tgp_idx, 'ALT', 'MXL', 'YRI'),
    ]
    # Determine the significance threshold.
    sig_thresh = 0.05
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        3, 2, figsize=(9, 6), dpi=300,
        facecolor='white',
        sharex=True, sharey=True,
        constrained_layout=True,
    )
    # Flatten the axes.
    axes = axes.flatten()
    # For every subplot
    for i, (wind_idx, t, s1, s2) in enumerate(pap_configs):
        # Intialize the configuration.
        config = (t, s1, s2)
        # Extract the observed value and windowed distribution.
        obs = obs_dicc[config]
        dist = wind_dicc[config][wind_idx]
        # Compute the p-value.
        pval = (np.count_nonzero(obs <= dist) / np.sum(~np.isnan(dist)))
        # If the p-value is significant.
        if pval < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (obs, '*', 25, 0.49)
        # Else, it is not significant.
        else:
            # Update the dictionary.
            sig_dicc[i] = (obs, r'$ns$', 10, 0.5225)
        # Plot the distribution.
        axes[i].hist(
            dist, bins=np.arange(0, 1.1, 0.1),
            histtype='stepfilled', color='gray',
        )
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlabel(r'$PAP$'+f'({label_dicc[s1]}, {label_dicc[t]}, {label_dicc[s2]}) Score per {window_size}kb Window')
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/psuedo_ancestry_painting_scores_{window_size}kb.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/psuedo_ancestry_painting_scores_{window_size}kb.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/psuedo_ancestry_painting_scores_{window_size}kb.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return



###############
### PHASING ###
###############

# Define a function to plot the synthetic neanderthal phasing results.
def plot_syn_nea_phased_df(df):
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
    })
    # Intialize the figure and axes.
    fig = plt.figure(
        figsize=(9, 3), dpi=300,
        facecolor='white', constrained_layout=True,
    )
    ax = fig.add_subplot(111)
    # Intialize color map.
    cmap = ListedColormap(['#0072B2', '#CC79A7'])
    # Plot the weird haplotypes.
    im = ax.imshow(df.iloc[:, 1:].T.to_numpy(), aspect='auto', cmap=cmap)
    # Label the rows.
    ax.set_yticks(np.arange(df.iloc[:, 1:].columns.size))
    ax.set_yticklabels(df.iloc[:, 1:].columns.values, size=10)
    # Seperate each box and add a grid.
    ax.set_xticks(np.arange(0, df.iloc[:, 1:].T.to_numpy().shape[1], 1))
    ax.set_yticks(np.arange(0, df.iloc[:, 1:].T.to_numpy().shape[0], 1))
    ax.set_xticks(np.arange(-0.5, df.iloc[:, 1:].T.to_numpy().shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, df.iloc[:, 1:].T.to_numpy().shape[0], 1), minor=True)
    ax.grid(which='minor', color='black', linestyle='-', linewidth=0.5)
    # Remove the ticks.
    ax.tick_params(bottom=False, labelbottom=False)
    ax.tick_params(which='minor', left=False, bottom=False, labelbottom=False)
    # Export the plot.
    plt.savefig(
        './supp_figures/png/synthetic_neanderthal_phasing.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/synthetic_neanderthal_phasing.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/synthetic_neanderthal_phasing.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to computed the number of pairwise differences for the phased archaic dataset.
def compile_phased_archaic_pwd_summary(df):
    # Intialize a dictionary.
    df_dicc = {'arc1': [], 'arc2': [], 'pwd': []}
    # Intialize a list of phased archaics.
    phased_arcs = [
        'Denisovan', 'Altai Nean.',
        'Chagyrskaya Nean. Hap. 1', 'Chagyrskaya Nean. Hap. 2', 
        'Vindija Nean. Hap. 1', 'Vindija Nean. Hap. 2',
    ]
    # Generate unique combinations of pairs (n=2)
    u_pairs = list(itertools.combinations(phased_arcs, 2))
    # For all pairwise comparisons of archaics.
    for (arc1, arc2) in u_pairs:
        # Compute the updated number of pairwise differences after phasing.
        pwd = np.nansum(pwd_per_site(df[arc1].values, df[arc2].values))
        # Update the dictionary.
        df_dicc['arc1'].append(arc1)
        df_dicc['arc2'].append(arc2)
        df_dicc['pwd'].append(pwd)
    # Convert to a dataframe.
    pwd_df = pd.DataFrame(df_dicc)
    # Rename all the columns to look pretty.
    pwd_df.rename(
        columns={
            'arc1': 'Archaic 1',
            'arc2': 'Archaic 2',
            'pwd': 'Focal 72kb Region (Pairwise Diffs.)',
        }, inplace=True
    )
    # Export the dataframe as a csv.
    pwd_df.to_csv('./dataframes/phased_late_neanderthals_archaic_pairwise_diffs_72kb.csv.gz', index=False)
    return pwd_df



########################
### CODING MUTATIONS ###
########################

# Define a function to compile the frequency of denisovan-specific missense mutations in modern humans.
def compile_tgp_sgdp_den_mis_mut_info(tgp_or_sgdp):
    # Intialize the meta information.
    tgp_or_sgdp_dicc = {
        'tgp': {
            'meta_path': '../meta_data/tgp_mod.txt',
            'spops': ['AFR', 'AMR', 'SAS', 'EAS', 'EUR'],
        },
        'sgdp': {
            'meta_path': '../meta_data/sgdp.txt',
            'spops': ['AFR', 'AMR', 'OCN', 'SAS', 'EAS', 'CAS', 'WER'],
        },
    }
    # Load the meta data.
    meta_data = pd.read_csv(
        tgp_or_sgdp_dicc[tgp_or_sgdp]['meta_path'], sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize the missense mutations dataframe.
    arc_coding_muts_df = pd.read_csv('../meta_data/muc19_72kb_denisovan_coding_muts_info.csv')
    arc_mis_muts_df = arc_coding_muts_df[arc_coding_muts_df['Mut. Type'] == 'Missense'].reset_index(drop=True)
    mis_mut_pos = arc_mis_muts_df['Chr12 Position (Hg19)'].values
    # Load the genotypes and the positions.
    gt, pos = load_hap_region(f'{tgp_or_sgdp}_arcs_masked_no_aa', 12, 40759001, 40831000)
    # Subset for the missenses muts.
    mis_mut_gt = gt.compress(np.isin(pos, mis_mut_pos), axis=0)
    # For every superpopulation.
    for spop in tgp_or_sgdp_dicc[tgp_or_sgdp]['spops']:
        # Update the dataframe.
        arc_mis_muts_df[spop] = calc_alt_freqs(
            mis_mut_gt.take(meta_data[meta_data['SUPERPOP'] == spop].index.values, axis=1),
        )
    # For every superpopulation.
    for spop in tgp_or_sgdp_dicc[tgp_or_sgdp]['spops']:
        # For every population.
        for pop in meta_data[meta_data['SUPERPOP'] == spop]['POP'].unique():
            # Update the dataframe.
            arc_mis_muts_df[pop] = calc_alt_freqs(
                mis_mut_gt.take(meta_data[meta_data['POP'] == pop].index.values, axis=1),
            )
    # Export the dataframe as a csv.    
    arc_mis_muts_df.to_csv(f'./dataframes/{tgp_or_sgdp}_denisovan_specific_missense_mutation_summary.csv.gz', index=False)
    return arc_mis_muts_df

# Define a function to compile the frequency of denisovan-specific missense mutations per ancient sample.
def compile_anc_amr_mis_mut_info():
    # Intialize the missense mutations dataframe and information.
    arc_coding_muts_df = pd.read_csv('../meta_data/muc19_72kb_denisovan_coding_muts_info.csv')
    arc_mis_muts_df = arc_coding_muts_df[arc_coding_muts_df['Mut. Type'] == 'Missense'].reset_index(drop=True)
    mis_mut_pos = arc_mis_muts_df['Chr12 Position (Hg19)'].values
    mis_mut_refs = arc_mis_muts_df['Ref. Allele'].values
    mis_mut_alts = arc_mis_muts_df['Denisovan Allele'].values
    # Load the ancient sample ids.
    anc_samps = np.array([
        'USR1', 'CK-13', 'Anzick',
        'CR-01', 'CT-01', 'NC_C',
        'PS-06', 'SC-05', 'SM-02',
        'SN-13', 'SN-17', 'SN-44',
        '2417J', '2417Q', '333B',
        'SRR8144644', 'SRR8144646',
        'SRR8144647', 'SRR8144650',
        'ERR2270782', 'ERR2270783',
        'ERR2270784', 'ERR2270785',

    ])
    # Intialize a dictionary.
    anc_ind_mis_muts_ac = {'POS': mis_mut_pos}
    # Intialize a dictionary to determine the allele from angsd.
    angsd_allele_dicc = {
        0: 'A', 1: 'C',
        2: 'G', 3: 'T',
    }
    # For all ancient samples.
    for samp in anc_samps:
        # Load the read data.
        anc_df = pd.read_csv(f'../ancient_americans/{samp}_ac.txt.gz', sep='\t')
        # Extract the positions array.
        anc_pos = anc_df['pos'].values
        # Intialize a column for the dataframes.
        df_col = []
        anc_ind_mis_muts_ac[samp] = []
        # For every missense mutation.
        for i, pos in enumerate(mis_mut_pos):
            # Extract the missense info.
            ref = mis_mut_refs[i]
            alt = mis_mut_alts[i]
            # Determine the index in the ancient sample.
            pos_idx = np.where(anc_pos == pos)[0]
            # If the site passed QC in the ancient sample...
            if pos_idx.size > 0:
                # Extract read counts.
                read_counts = anc_df.iloc[pos_idx].to_numpy()[:, -4:].flatten()
                # Determine which allele indicies have more than two reads.
                allele_idx = np.where(read_counts >= 2)[0]
                # If this site is not mono- or bi-allelic.
                if allele_idx.size > 2:
                    # Append the list and update the dictionary.
                    df_col.append('.')
                    anc_ind_mis_muts_ac[samp].append(np.nan)
                # Else-if this is a mono-allelic site...
                elif allele_idx.size == 1:
                    # Determine the allele identity.
                    anc_allele_call = angsd_allele_dicc[allele_idx[0]]
                    # If the ancient allele is the refernce allele...
                    if anc_allele_call == ref:
                        # Append the list and update the dictionary.
                        df_col.append(ref+'/'+ref)
                        anc_ind_mis_muts_ac[samp].append(0)
                    # Else-if the ancient allele is the alternative allele...
                    elif anc_allele_call == alt:
                        # Append the list and update the dictionary.
                        df_col.append(alt+'/'+alt)
                        anc_ind_mis_muts_ac[samp].append(2)
                    # Else...
                    else:
                        # Append the list and update the dictionary.
                        df_col.append('.')
                        anc_ind_mis_muts_ac[samp].append(np.nan)
                # Else-if this is a bi-allelic site...
                elif allele_idx.size == 2:
                    # Determine the allele identities.
                    anc_allele_calls = np.array([angsd_allele_dicc[j] for j in allele_idx])
                    # If the ancient alleles are the same as the hg19 alleles...
                    if np.array_equal(np.sort(anc_allele_calls), np.sort(np.array([ref, alt]))):
                        # Append the list and update the dictionary.
                        df_col.append(ref+'/'+alt)
                        anc_ind_mis_muts_ac[samp].append(1)
                    # Else...
                    else:
                        # Append the list and update the dictionary.
                        df_col.append('.')
                        anc_ind_mis_muts_ac[samp].append(np.nan)
                # Else...
                else:
                    # Append the list and update the dictionary.
                    df_col.append('.')
                    anc_ind_mis_muts_ac[samp].append(np.nan)
            # Else...
            else:
                # Append the list and update the dictionary.
                df_col.append('.')
                anc_ind_mis_muts_ac[samp].append(np.nan)
        # Update the the dataframe.
        arc_mis_muts_df[samp] = df_col
    # Export the dataframe as a csv.    
    arc_mis_muts_df.to_csv('./dataframes/ancient_americans_denisovan_specific_missense_mutation_summary.csv.gz', index=False)
    # Convert the dictionary to a dataframe.
    arc_mis_muts_aac_df = pd.DataFrame(anc_ind_mis_muts_ac)
    # Export the dataframe as a csv.    
    arc_mis_muts_aac_df.to_csv('./dataframes/ancient_americans_denisovan_specific_missense_mutation_aac.csv.gz', index=False)
    return arc_mis_muts_df, arc_mis_muts_aac_df

# Define a function to assess if missense mutations are elevated in AMR populations.
def compile_amr_v_non_amr_mis_mut_summary():
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    tgp_info = {key: np.array(value) for key, value in tgp_df.to_dict(orient='list').items()}
    # Extract the focal masks.
    is_amr = tgp_info['SUPERPOP'] == 'AMR'
    is_non_amr = (tgp_info['SUPERPOP'] != 'AMR') & (tgp_info['SUPERPOP'] != 'AFR')
    # Intialize the missense mutations dataframe and information.
    arc_coding_muts_df = pd.read_csv('../meta_data/muc19_72kb_denisovan_coding_muts_info.csv')
    arc_mis_muts_df = arc_coding_muts_df[arc_coding_muts_df['Mut. Type'] == 'Missense'].reset_index(drop=True)
    mis_mut_pos = arc_mis_muts_df['Chr12 Position (Hg19)'].values
    # Load the genotypes and the positions.
    gt, pos = load_hap_region('tgp_arcs_masked_no_aa', 12, 40759001, 40831000)
    # Subset the genotype matricies.
    amr_gt = gt.take(tgp_df.index.values[is_amr], axis=1)
    non_amr_gt = gt.take(tgp_df.index.values[is_non_amr], axis=1)
    # Intialize a dictionary to store the results.
    df_dicc = {
        'mis_pos': [],
        'grp1_chroms': [], 'grp1_aac': [], 'grp1_prop': [],
        'grp2_chroms': [], 'grp2_aac': [], 'grp2_prop': [],
        'z_stat': [], 'z_p': [], 'odds_ratio': [], 'f_p': [],
    }
    # For every missense mutation.
    for mut_pos in mis_mut_pos:
        # Extract the missesnse mutation info per group.
        amr_ind_aac = amr_gt.compress(pos == mut_pos, axis=0).to_n_alt().flatten()
        amr_aac = np.nansum(amr_ind_aac)
        amr_chroms = amr_ind_aac.size * 2
        non_amr_ind_aac = non_amr_gt.compress(pos == mut_pos, axis=0).to_n_alt().flatten()
        non_amr_aac = np.nansum(non_amr_ind_aac)
        non_amr_chroms = non_amr_ind_aac.size * 2
        # Perform a proportions Z-test.
        z_stat, z_p = proportions_ztest(
            np.array([amr_aac, non_amr_aac]),
            np.array([amr_chroms, non_amr_chroms]),
            alternative='larger',
        )
        # Perform a Fisher's exact test.
        odds_ratio, f_p = stats.fisher_exact(
            [[amr_aac, amr_chroms - amr_aac], [non_amr_aac, non_amr_chroms - non_amr_aac]], 
            alternative='greater',
        )
        # Update the dictionary.
        df_dicc['mis_pos'].append(mut_pos)
        df_dicc['grp1_chroms'].append(amr_chroms)
        df_dicc['grp1_aac'].append(amr_aac)
        df_dicc['grp1_prop'].append(amr_aac / amr_chroms)
        df_dicc['grp2_chroms'].append(non_amr_chroms)
        df_dicc['grp2_aac'].append(non_amr_aac)
        df_dicc['grp2_prop'].append(non_amr_aac / non_amr_chroms)
        df_dicc['z_stat'].append(z_stat)
        df_dicc['z_p'].append(z_p)
        df_dicc['odds_ratio'].append(odds_ratio)
        df_dicc['f_p'].append(f_p)
    # Convert to a dataframe.
    summary_df = pd.DataFrame(df_dicc)
    # Rename the columns to look pretty.
    summary_df.rename(
        columns={
            'mis_pos': 'Denisovan-specific Mis. Mut. Pos.',
            'grp1_chroms': 'Number of Chroms. (AMR)',
            'grp1_aac': 'Denisovan-specific Mis. Mut. Count (AMR)',
            'grp1_prop': 'Denisovan-specific Mis. Freq. (AMR)',
            'grp2_chroms': 'Number of Chroms. (Non-AMR)',
            'grp2_aac': 'Denisovan-specific Mis. Mut. Count (Non-AMR)',
            'grp2_prop': 'Denisovan-specific Mis. Freq. (Non-AMR)',
            'z_stat': r'$Z-Statistic$',
            'z_p': r'$P-value$ (Prop. $Z$-Test)',
            'odds_ratio': 'Odds Ratio',
            'f_p': r"$P-value$ (Fisher's Exact Test)",
        }, inplace=True
    )
    # Export.
    summary_df.to_csv('./dataframes/amr_v_non_amr_denisovan_specific_missense_mutation_proportions_ztest_fet.csv.gz', index=False)
    return summary_df

# Define a function to summarize the introgressed tract frequency.
def compile_introgressed_tract_freqs(tracts_df):
    # Load the meta data file for the TGP.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    tgp_info = {key: np.array(value) for key, value in tgp_df.to_dict(orient='list').items()}
    # Intialize an ordered dictionary.
    tgp_dicc = {
        'AMR': ['MXL', 'PEL', 'CLM', 'PUR'],
        'SAS': ['BEB', 'STU', 'ITU', 'PJL', 'GIH'],
        'EAS': ['CHB', 'KHV', 'CHS', 'JPT', 'CDX'],
        'EUR': ['TSI', 'CEU', 'IBS', 'GBR', 'FIN'],
        'AFR': ['LWK', 'GWD', 'MSL', 'ESN', 'YRI'],
    }
    # Intialize a list of non-Afr populations.
    non_afr_pops = [
        'MXL', 'PEL', 'CLM', 'PUR',
        'BEB', 'STU', 'ITU', 'PJL', 'GIH',
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX',
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN',
    ]
    # Clean the introgressed tracts.
    cleaned_tracts_df, _, _, _, _ = load_region_tracts(tracts_df, non_afr_pops)
    # Intialize a dictionary.
    tract_summary = {
        'spop': [], 'pop': [], 'n_tot': [],
        'n_tracts': [],  'freq': [],
    }
    # For every super population.
    for spop in tgp_dicc:
        # For every population.
        for pop in tgp_dicc[spop]:
            # Determine the total number of chromosomes.
            n_tot = (tgp_info['POP'] == pop).sum() * 2
            # Determine the number of tracts.
            n_tracts = (cleaned_tracts_df['POP'].values == pop).sum()
            # Update the dictionary.
            tract_summary['spop'].append(spop)
            tract_summary['pop'].append(pop)
            tract_summary['n_tot'].append(n_tot)
            tract_summary['n_tracts'].append(n_tracts)
            tract_summary['freq'].append(n_tracts / n_tot)
    return pd.DataFrame(tract_summary)

# Define a fucntion to plot the relationship between denisovan-like haplotype frequency and the mean frequency of the missense mutations.
def plot_tract_freq_x_mean_mis_mut_freq():
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize the missense mutations dataframe.
    arc_coding_muts_df = pd.read_csv('../meta_data/muc19_72kb_denisovan_coding_muts_info.csv')
    arc_mis_muts_df = arc_coding_muts_df[arc_coding_muts_df['Mut. Type'] == 'Missense'].reset_index(drop=True)
    mis_mut_pos = arc_mis_muts_df['Chr12 Position (Hg19)'].values
    # Load the genotypes and the positions.
    gt, pos = load_hap_region(f'tgp_arcs_masked_no_aa', 12, 40759001, 40831000)
    # Subset for the missenses muts.
    mis_mut_gt = gt.compress(np.isin(pos, mis_mut_pos), axis=0)
    # Load the introgressed tract frequency dataframe and extract the information.
    muc19_tracts_df = pd.read_csv('../hmmix_tracts/tgp_hmmix_haplotype_tracts_muc19.csv.gz')
    tract_freq_df = compile_introgressed_tract_freqs(muc19_tracts_df)
    tract_freqs = tract_freq_df['freq'].values
    spops = tract_freq_df['spop'].values
    pops = tract_freq_df['pop'].values
    # Intialize dictionaries for plotting.
    spop_colors = {
        'AFR': 'grey', 'AMR': '#009E73', 'SAS': '#CC79A7', 'EAS': '#0072B2', 'EUR': '#D55E00',
    }
    marker_dicc = {0: '^', 1: 'v'}
    # Intialize lists.
    mis_mut_freqs, mis_mut_cis, colors = [], [], []
    # For every population.
    for spop, pop in zip(spops, pops):
        # Extract the allele frequency per individual per site.
        ind_freqs = mis_mut_gt.take(tgp_df[tgp_df['POP'] == pop].index.values, axis=1).to_n_alt().flatten() / 2
        # Compute the 95% confidence intervals.
        _, freqs_ci = sem_ci_of_mean(ind_freqs)
        # Update the lists.
        mis_mut_freqs.append(ind_freqs.mean())
        mis_mut_cis.append(freqs_ci)
        colors.append(spop_colors[spop])
    # Perform a least squares regression.
    slope, intercept, rho, pval, se = stats.linregress(tract_freqs, mis_mut_freqs)
    # Compute the line of best fit for the linear regression model.
    line_of_best_fit = (slope * tract_freqs) + intercept
    # Sort the values based on the tract frequencies.
    sorted_indices = np.argsort(tract_freqs)
    sorted_tract_freqs = tract_freqs[sorted_indices]
    sorted_line_of_best_fit = line_of_best_fit[sorted_indices]
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the figure and axes.
    fig = plt.figure(
        figsize=(6, 4), dpi=300,
        facecolor='white', constrained_layout=True,
    )
    ax = fig.add_subplot(111)
    # Plot the regression results.
    ax.plot(sorted_tract_freqs, sorted_line_of_best_fit,  'k--')
    # Intialize the current marker.
    c_marker = 0
    # For every population sorted by haplotype frequency.
    for idx in sorted_indices:
        # Plot the results.
        ax.errorbar(
            tract_freqs[idx], mis_mut_freqs[idx], yerr=mis_mut_cis[idx],
            color=colors[idx], marker=marker_dicc[c_marker], markersize=5,
            markeredgewidth=1.25, elinewidth=1.25, capthick=1.25, capsize=5, zorder=5,
        )
        # Move the current marker forward.
        c_marker = abs(c_marker - 1)
    # Label the axes.
    ax.set_xlabel(r'$MUC19$ Introgressed Tract Frequency')
    ax.set_ylabel('Denisovan-specific Missense Mut. Frequency')
    # Add the legend.
    legend = ax.legend(
        handles=[Patch(color=spop_colors[spop], label=spop) for spop in spop_colors], 
        loc='upper left', bbox_to_anchor=(0, 1),
        frameon=True, fancybox=True, shadow=True, fontsize=8,
    )
    legend.set_title(f"Pearson's $\\rho = {rho:.3f}$\n$P-value = {pval:.3e}$", prop={'size': 8})
    legend._legend_box.align = 'left'
    # Export the plot.
    plt.savefig(
        './supp_figures/png/tgp_introgressed_tract_frequency_x_missense_mutation_frequency.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/tgp_introgressed_tract_frequency_x_missense_mutation_frequency.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/tgp_introgressed_tract_frequency_x_missense_mutation_frequency.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to compile the frequency of denisovan-specific missense mutations per individual.
def compute_mxl_anc_amr_sgdp_mis_mut_freq_per_ind():
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    sgdp_df = pd.read_csv(
        '../meta_data/sgdp.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Load the ancient sample ids.
    anc_samps = np.array([
        'USR1', 'CK-13', 'Anzick',
        'CR-01', 'CT-01', 'NC_C',
        'PS-06', 'SC-05', 'SM-02',
        'SN-13', 'SN-17', 'SN-44',
        '2417J', '2417Q', '333B',
        'SRR8144644', 'SRR8144646',
        'SRR8144647', 'SRR8144650',
        'ERR2270782', 'ERR2270783',
        'ERR2270784', 'ERR2270785',

    ])
    # Intialize a dictionary to determine the allele from angsd.
    angsd_allele_dicc = {
        0: 'A', 1: 'C',
        2: 'G', 3: 'T',
    }
    # Intialize the missense mutations dataframe and information.
    arc_coding_muts_df = pd.read_csv('../meta_data/muc19_72kb_denisovan_coding_muts_info.csv')
    arc_mis_muts_df = arc_coding_muts_df[arc_coding_muts_df['Mut. Type'] == 'Missense'].reset_index(drop=True)
    mis_mut_pos = arc_mis_muts_df['Chr12 Position (Hg19)'].values
    mis_mut_refs = arc_mis_muts_df['Ref. Allele'].values
    mis_mut_alts = arc_mis_muts_df['Denisovan Allele'].values
    # Load the genotypes and the positions.
    tgp_gt, tgp_pos = load_hap_region('tgp_arcs_masked_no_aa', 12, 40759001, 40831000)
    sgdp_gt, sgdp_pos = load_hap_region('sgdp_arcs_masked_no_aa', 12, 40759001, 40831000)
    # Subset for the focal populations.
    mxl_gt = tgp_gt.take(tgp_df[tgp_df['POP'] == 'MXL'].index.values, axis=1)
    amr_gt = sgdp_gt.take(sgdp_df[sgdp_df['SUPERPOP'] == 'AMR'].index.values, axis=1)
    # Intialize a dictionary to store the allele frequencies per individual.
    mis_mut_ind_freqs = {pos: {} for pos in mis_mut_pos}
    # For evry missense mutation.
    for pos in mis_mut_pos:
        # Update the dictionary.
        mis_mut_ind_freqs[pos]['MXL'] = mxl_gt.compress(tgp_pos == pos, axis=0).to_n_alt().flatten() / 2
        mis_mut_ind_freqs[pos]['AMR'] = amr_gt.compress(sgdp_pos == pos, axis=0).to_n_alt().flatten() / 2
        mis_mut_ind_freqs[pos]['ANC'] = np.array([])
    # For all ancient samples.
    for samp in anc_samps:
        # Load the read data.
        anc_df = pd.read_csv(f'../ancient_americans/{samp}_ac.txt.gz', sep='\t')
        # Extract the positions array.
        anc_pos = anc_df['pos'].values
        # For every missense mutation.
        for i, pos in enumerate(mis_mut_pos):
            # Extract the missense info.
            ref = mis_mut_refs[i]
            alt = mis_mut_alts[i]
            # Determine the index in the ancient sample.
            pos_idx = np.where(anc_pos == pos)[0]
            # If the site passed QC in the ancient sample...
            if pos_idx.size > 0:
                # Extract read counts.
                read_counts = anc_df.iloc[pos_idx].to_numpy()[:, -4:].flatten()
                # Determine which allele indicies have more than two reads.
                allele_idx = np.where(read_counts >= 2)[0]
                # If this site is not mono- or bi-allelic.
                if allele_idx.size > 2:
                    # Update the dictionary.
                    mis_mut_ind_freqs[pos]['ANC'] = np.append(mis_mut_ind_freqs[pos]['ANC'], np.nan)
                # Else-if this is a mono-allelic site...
                elif allele_idx.size == 1:
                    # Determine the allele identity.
                    anc_allele_call = angsd_allele_dicc[allele_idx[0]]
                    # If the ancient allele is the refernce allele...
                    if anc_allele_call == ref:
                        # Update the dictionary.
                        mis_mut_ind_freqs[pos]['ANC'] = np.append(mis_mut_ind_freqs[pos]['ANC'], 0)
                    # Else-if the ancient allele is the alternative allele...
                    elif anc_allele_call == alt:
                        # Update the dictionary.
                        mis_mut_ind_freqs[pos]['ANC'] = np.append(mis_mut_ind_freqs[pos]['ANC'], 1)
                    # Else...
                    else:
                        # Update the dictionary.
                        mis_mut_ind_freqs[pos]['ANC'] = np.append(mis_mut_ind_freqs[pos]['ANC'], np.nan)
                # Else-if this is a bi-allelic site...
                elif allele_idx.size == 2:
                    # Determine the allele identities.
                    anc_allele_calls = np.array([angsd_allele_dicc[j] for j in allele_idx])
                    # If the ancient alleles are the same as the hg19 alleles...
                    if np.array_equal(np.sort(anc_allele_calls), np.sort(np.array([ref, alt]))):
                        # Update the dictionary.
                        mis_mut_ind_freqs[pos]['ANC'] = np.append(mis_mut_ind_freqs[pos]['ANC'], 0.5)
                    # Else...
                    else:
                        # Update the dictionary.
                        mis_mut_ind_freqs[pos]['ANC'] = np.append(mis_mut_ind_freqs[pos]['ANC'], np.nan)
                # Else...
                else:
                    # Update the dictionary.
                    mis_mut_ind_freqs[pos]['ANC'] = np.append(mis_mut_ind_freqs[pos]['ANC'], np.nan)
            # Else...
            else:
                # Update the dictionary.
                mis_mut_ind_freqs[pos]['ANC'] = np.append(mis_mut_ind_freqs[pos]['ANC'], np.nan)
    return mis_mut_ind_freqs

# Define a function to compile the a summary for the frequency of denisovan-specific missense mutations per individual.
def compile_mxl_anc_amr_sgdp_mis_mut_freq_per_ind_summary():
    # Compile the per indiviual frequencies.
    ind_freqs = compute_mxl_anc_amr_sgdp_mis_mut_freq_per_ind()
    # Intialize a label dictionary.
    pop_labs = {
        'MXL': 'MXL', 'AMR': 'SDGP AMR', 'ANC': 'Ancient Americans'
    }
    # Intialize a list of pairwise comparisons.
    pw_combos = [('MXL', 'ANC'), ('MXL', 'AMR'), ('AMR', 'ANC')]
    # Intialize dictionaries.
    mis_muts_info = {
        'mis_pos': [], 'pop': [], 'n_chroms': [],
        'freq': [], 'std': [], 'sem': [], '95_ci': [],
    }
    df_dicc = {
        'mis_pos': [], 'grp1': [], 'grp2': [],
        'grp1_chroms': [], 'grp1_aac': [], 'grp1_prop': [],
        'grp2_chroms': [], 'grp2_aac': [], 'grp2_prop': [],
        'z_stat': [], 'z_p': [], 'odds_ratio': [], 'f_p': [],
    }
    # For every missense mutation.
    for mis_pos in ind_freqs:
        # For every population.
        for pop in ind_freqs[mis_pos]:
            # Compute the allele frequency for the population, std, sem, and 95% CIs.
            n_chroms = (~np.isnan(ind_freqs[mis_pos][pop])).sum() * 2
            pop_freq = np.nanmean(ind_freqs[mis_pos][pop])
            pop_std = np.nanstd(ind_freqs[mis_pos][pop])
            pop_sem, pop_ci = sem_ci_of_mean(ind_freqs[mis_pos][pop])
            # Fill the dictionary.
            mis_muts_info['mis_pos'].append(mis_pos)
            mis_muts_info['pop'].append(pop_labs[pop])
            mis_muts_info['n_chroms'].append(n_chroms)
            mis_muts_info['freq'].append(pop_freq)
            mis_muts_info['std'].append(pop_std)
            mis_muts_info['sem'].append(pop_sem)
            mis_muts_info['95_ci'].append(pop_ci)
    # For every pairwise combination.
    for grp1, grp2 in pw_combos:
        # For every missense mutation.
        for mis_pos in ind_freqs:
            # Extract the information for group 1.
            grp1_freqs = ind_freqs[mis_pos][grp1][~np.isnan(ind_freqs[mis_pos][grp1])]
            grp1_chroms = grp1_freqs.size * 2
            grp1_aac = (grp1_freqs * 2).sum()
            # Extract the information for group 2.
            grp2_freqs = ind_freqs[mis_pos][grp2][~np.isnan(ind_freqs[mis_pos][grp2])]
            grp2_chroms = grp2_freqs.size * 2
            grp2_aac = (grp2_freqs * 2).sum()
            # Perform a proportions Z-test.
            z_stat, z_p = proportions_ztest(
                np.array([grp1_aac, grp2_aac]),
                np.array([grp1_chroms, grp2_chroms]),
            )
            # Perform a Fisher's exact test.
            odds_ratio, f_p = stats.fisher_exact(
                [[grp1_aac, grp1_chroms - grp1_aac],
                 [grp2_aac, grp2_chroms - grp2_aac]],
            )
            # Fill the dictionary.
            df_dicc['mis_pos'].append(mis_pos)
            df_dicc['grp1'].append(pop_labs[grp1])
            df_dicc['grp2'].append(pop_labs[grp2])
            df_dicc['grp1_chroms'].append(grp1_chroms)
            df_dicc['grp1_aac'].append(grp1_aac)
            df_dicc['grp1_prop'].append(grp1_aac / grp1_chroms)
            df_dicc['grp2_chroms'].append(grp2_chroms)
            df_dicc['grp2_aac'].append(grp2_aac)
            df_dicc['grp2_prop'].append(grp2_aac / grp2_chroms)
            df_dicc['z_stat'].append(z_stat)
            df_dicc['z_p'].append(z_p)
            df_dicc['odds_ratio'].append(odds_ratio)
            df_dicc['f_p'].append(f_p)
    # Convert the dictionaries to dataframes.
    mis_muts_info_df = pd.DataFrame(mis_muts_info)
    summary_df = pd.DataFrame(df_dicc)
    # Rename the columns to look pretty.
    mis_muts_info_df.rename(
        columns={
            'mis_pos': 'Denisovan-specific Mis. Mut. Pos.',
            'pop': 'Group',
            'n_chroms': 'Number of Chroms.',
            'freq': r'Denisovan-specific Mis. Mut. Freq.',
            'std': r'Denisovan-specific Mis. Mut. $\left( \sigma \right)$',
            'sem': r'Denisovan-specific Mis. Mut. $\left( SEM \right)$',
            '95_ci': r'Denisovan-specific Mis. Mut. $\left( \pm CI_{95\%} \right)$',
        }, inplace=True
    )
    summary_df.rename(
        columns={
            'mis_pos': 'Denisovan-specific Mis. Mut. Pos.',
            'grp1': 'Group 1', 'grp2': 'Group 2', 
            'grp1_chroms': 'Number of Chroms. (Group 1)',
            'grp1_aac': 'Denisovan-specific Mis. Mut. Count (Group 1)',
            'grp1_prop': 'Denisovan-specific Mis. Freq. (Group 1)',
            'grp2_chroms': 'Number of Chroms. (Group 2)',
            'grp2_aac': 'Denisovan-specific Mis. Mut. Count (Group 2)',
            'grp2_prop': 'Denisovan-specific Mis. Freq. (Group 2)',
            'z_stat': r'$Z-Statistic$',
            'z_p': r'$P-value$ (Prop. $Z$-Test)',
            'odds_ratio': 'Odds Ratio',
            'f_p': r"$P-value$ (Fisher's Exact Test)"
        }, inplace=True
    )
    # Export the dataframes as csvs.
    mis_muts_info_df.to_csv('./dataframes/mxl_sgdp_amr_anc_amr_denisovan_specific_missense_mutation_summary.csv.gz', index=False)
    summary_df.to_csv('./dataframes/mxl_sgdp_amr_anc_amr_denisovan_specific_missense_mutation_proportions_ztest_fet.csv.gz', index=False)
    return mis_muts_info_df, summary_df

# Define a function to compile tgp amr individual missense mutation and indigenous american ancestry proportions.
def compile_tgp_amr_mis_mut_region_info():
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize the missense mutations dataframe.
    den_coding_muts_df = pd.read_csv('../meta_data/muc19_72kb_denisovan_coding_muts_info.csv')
    den_mis_muts_df = den_coding_muts_df[den_coding_muts_df['Mut. Type'] == 'Missense'].reset_index(drop=True)
    all_mis_mut_pos = den_mis_muts_df['Chr12 Position (Hg19)'].values
    is_focal_pos = all_mis_mut_pos == 40808726
    mis_mut_pos = all_mis_mut_pos[is_focal_pos]
    mis_mut_refs = den_mis_muts_df['Ref. Allele'].values[is_focal_pos]
    mis_mut_alts = den_mis_muts_df['Denisovan Allele'].values[is_focal_pos]
    # Load the ancestry proportions.
    amr_anc_props_df = pd.read_csv('../amr_lai/anc_props/amr_72kb_region_props.csv.gz')
    # Load the genotypes and the positions.
    gt, pos = load_hap_region('tgp_arcs_masked_no_aa', 12, 40759001, 40831000)
    # Subset for the missenses muts.
    mis_mut_gt = gt.compress(np.isin(pos, mis_mut_pos), axis=0)
    # Intialize a dictionary to store the results.
    ind_info = {'POP': [], 'IND': [], 'FREQ': [], 'PROP': []}
    # For every amr population.
    for pop in ['MXL', 'PEL', 'CLM', 'PUR']:
        # Subset the the dataframes.
        pop_df = tgp_df[tgp_df['POP'] == pop]
        anc_df = amr_anc_props_df[amr_anc_props_df['POP'] == pop]
        # Extract the individuals, indicies, and ancestry proportions.
        pop_inds, pop_idx = pop_df['IND'].values, pop_df.index.values
        anc_inds, nat_props = anc_df['IND'].values, anc_df['NAT'].values
        # Subset the genotypes.
        pop_gt = mis_mut_gt.take(pop_idx, axis=1)
        # For every individual.
        for i, ind in enumerate(pop_inds):
            # Subset the genotypes.
            ind_gt = pop_gt[:, i, :]
            # Compute the mean frequency across all nine mutations.
            ind_freq = ind_gt.sum() / ind_gt.size
            # Create a mask.
            is_ind = anc_inds == ind
            # Extract the indigenous american ancestry proportion.
            ind_prop = nat_props[is_ind][0]
            # Update the dictionary.
            ind_info['POP'].append(pop)
            ind_info['IND'].append(ind)
            ind_info['FREQ'].append(ind_freq)
            ind_info['PROP'].append(ind_prop)
    # Load the ancient sample ids.
    anc_samps = np.array([
        'USR1', 'CK-13', 'Anzick',
        'CR-01', 'CT-01', 'NC_C',
        'PS-06', 'SC-05', 'SM-02',
        'SN-13', 'SN-17', 'SN-44',
        '2417J', '2417Q', '333B',
        'SRR8144644', 'SRR8144646',
        'SRR8144647', 'SRR8144650',
        'ERR2270782', 'ERR2270783',
        'ERR2270784', 'ERR2270785',

    ])
    # Intialize a dictionary to determine the allele from angsd.
    angsd_allele_dicc = {
        0: 'A', 1: 'C',
        2: 'G', 3: 'T',
    }
    # For all ancient samples.
    for samp in anc_samps:
        # Load the read data.
        anc_df = pd.read_csv(f'../ancient_americans/{samp}_ac.txt.gz', sep='\t')
        # Extract the positions array.
        anc_pos = anc_df['pos'].values
        # Intialize the numerator and denominator.
        anc_numer, anc_denom = 0, 0
        # For every missense mutation.
        for i, pos in enumerate(mis_mut_pos):
            # Extract the missense info.
            ref = mis_mut_refs[i]
            alt = mis_mut_alts[i]
            # Determine the index in the ancient sample.
            pos_idx = np.where(anc_pos == pos)[0]
            # If the site passed QC in the ancient sample...
            if pos_idx.size > 0:
                # Extract read counts.
                read_counts = anc_df.iloc[pos_idx].to_numpy()[:, -4:].flatten()
                # Determine which allele indicies have more than two reads.
                allele_idx = np.where(read_counts >= 2)[0]
                # If this site is not mono- or bi-allelic.
                if allele_idx.size > 2:
                    # Continue to the next site.
                    continue
                # Else-if this is a mono-allelic site...
                elif allele_idx.size == 1:
                    # Determine the allele identity.
                    anc_allele_call = angsd_allele_dicc[allele_idx[0]]
                    # If the ancient allele is the refernce allele...
                    if anc_allele_call == ref:
                        # Update the denominator.
                        anc_denom += 2 
                    # Else-if the ancient allele is the alternative allele...
                    elif anc_allele_call == alt:
                        # Update the numertaor and denominator.
                        anc_numer += 2 
                        anc_denom += 2 
                    # Else...
                    else:
                        # Continue to the next site.
                        continue
                # Else-if this is a bi-allelic site...
                elif allele_idx.size == 2:
                    # Determine the allele identities.
                    anc_allele_calls = np.array([angsd_allele_dicc[j] for j in allele_idx])
                    # If the ancient alleles are the same as the hg19 alleles...
                    if np.array_equal(np.sort(anc_allele_calls), np.sort(np.array([ref, alt]))):
                        # Update the numertaor and denominator.
                        anc_numer += 1 
                        anc_denom += 2 
                    # Else...
                    else:
                        # Continue to the next site.
                        continue
                # Else...
                else:
                    # Continue to the next site.
                    continue
            # Else...
            else:
                # Continue to the next site.
                continue
        # If there is information.
        if anc_denom > 0:
            # Update the dictionary.
            ind_info['POP'].append('ANC')
            ind_info['IND'].append(samp)
            ind_info['FREQ'].append(anc_numer / anc_denom)
            ind_info['PROP'].append(1)
    # Convert the dictionary to a dataframe.
    return pd.DataFrame(ind_info)

# Define a function to plot missense mutation frequency x indigenous american ancestry proportion.
def plot_mis_mut_freq_x_iaa_prop(ind_info_df):
    # Convert the dataframe to a dictionary.
    ind_info = {key: np.array(value) for key, value in ind_info_df.to_dict(orient='list').items()}
    # Intialize a dictionary for plotting.
    c_dicc = {
        'MXL': '#009E73',
        'PEL': '#E69F00',
        'CLM': '#0072B2',
        'PUR': '#D55E00',
        'ANC': '#CC79A7',
    }
    # Perform a least squares regression.
    slope, intercept, rho, pval, se = stats.linregress(ind_info['PROP'], ind_info['FREQ'])
    # Compute the line of best fit for the linear regression model.
    line_of_best_fit = (slope * ind_info['PROP']) + intercept
    # Sort the values based on the IAA ancestry proportions.
    sorted_indices = np.argsort(ind_info['PROP'])
    sorted_anc_props = ind_info['PROP'][sorted_indices]
    sorted_line_of_best_fit = line_of_best_fit[sorted_indices]
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the figure and axes.
    fig = plt.figure(
        figsize=(6, 4), dpi=300,
        facecolor='white', constrained_layout=True,
    )
    ax = fig.add_subplot(111)
    # Plot the regression results.
    ax.plot(sorted_anc_props, sorted_line_of_best_fit,  'k--')
    # Intialize the jitter threshold.
    jitter_thresh = 0.005
    # Add jitter.
    x_jittered = ind_info['PROP'] + np.random.uniform(-jitter_thresh, jitter_thresh, size=ind_info['IND'].size)
    y_jittered = ind_info['FREQ'] + np.random.uniform(-jitter_thresh, jitter_thresh, size=ind_info['IND'].size)
    # Shuffle the indicies.
    shuffled_idx = np.random.choice(np.arange(ind_info['IND'].size), size=ind_info['IND'].size, replace=False)
    # For every individual.
    for i in range(shuffled_idx.size):
        # Determine the population.
        pop = ind_info['POP'][shuffled_idx[i]]
        # Plot the results.
        ax.scatter(
            x_jittered[shuffled_idx[i]], y_jittered[shuffled_idx[i]],
            facecolor='none', edgecolor=c_dicc[pop], marker='X', s=75, linewidths=0.75, alpha=0.5,
        )
    # Label the axes.
    ax.set_xlabel('Indigenous American Ancestry Proportion')
    ax.set_ylabel('Denisovan-specific Missense Mut. Frequency')
    # Add the legend.
    legend = ax.legend(
        handles=[Patch(color=c_dicc[pop], label=pop) for pop in c_dicc], 
        loc='upper left', bbox_to_anchor=(0, 1),
        frameon=True, fancybox=True, shadow=True, fontsize=8,
    )
    legend.set_title(f"Pearson's $\\rho = {rho:.3f}$\n$P-value = {pval:.3e}$", prop={'size': 8})
    legend._legend_box.align = 'left'
    # Export the plot.
    plt.savefig(
        './supp_figures/png/amr_ind_iaa_prop_x_missense_mutation_frequency.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/amr_ind_iaa_prop_x_missense_mutation_frequency.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/amr_ind_iaa_prop_x_missense_mutation_frequency.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return



######################
### HETEROZYGOSITY ###
######################

# Define a function to load the genome wide heterozygosity results for the archaics.
def load_arc_gw_het_info():
    # Intialize a list of archaics.
    arc_list = ['DEN', 'ALT', 'CHA', 'VIN']
    # Load the results.
    het_dicc = {
        arc: np.vstack([np.loadtxt(f'../muc19_results/{arc.lower()}_masked_no_aa/archaic_het_sites_heterozygosity_eff_seq_len_chr{chrom}.txt.gz') for chrom in range(1, 23)]).sum(axis=0)
        for arc in arc_list
    }
    return het_dicc

# Define a function to load the genome wide heterozygosity results for the tgp.
def load_tgp_gw_het_info():
    # Load the meta data file for the tgp.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize the tgp supepopulation list.
    spop_list = ['AFR', 'SAS', 'EAS', 'EUR', 'AMR']
    # Load the results.
    het_counts = np.vstack([
        np.loadtxt(f'../muc19_results/tgp_mod_no_aa/tgp_het_sites_chr{chrom}.txt.gz', dtype=int) for chrom in range(1, 23)
    ]).sum(axis=0)
    # Intialize a dictionary.
    het_dicc = {spop: het_counts[tgp_df[tgp_df['SUPERPOP'] == spop].index.values] for spop in spop_list}
    return het_dicc

# Define a function to compile a summary of the genome-wide heterozygosity.
def compile_afr_arc_gw_het_summary():
    # Load the afr results.
    tgp_het_info = load_tgp_gw_het_info()
    afr_het_info = np.vstack([
        np.loadtxt(f'../muc19_results/tgp_mod_no_aa/afr_heterozygosity_eff_seq_len_chr{chrom}.txt.gz') for chrom in range(1, 23)
    ]).sum(axis=0)
    # Load the archaic results.
    arc_het_info = load_arc_gw_het_info()
    # Intialize dictionaries.
    arc_labs = {
        'DEN': 'Denisovan', 'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.', 'VIN': 'Vindija Nean.',
    }
    het_dicc = {'grp': [], 'obs': []}
    # Update the dictionary with the afr results.
    het_dicc['grp'].append('AFR')
    het_dicc['obs'].append(np.mean(tgp_het_info['AFR'] / afr_het_info[1]))
    # For every archaic.
    for arc in arc_labs:
        # Update the dictionary with the afr results.
        het_dicc['grp'].append(arc_labs[arc])
        het_dicc['obs'].append(arc_het_info[arc][0] / arc_het_info[arc][2])
    # Convert to a dataframe.
    het_df = pd.DataFrame(het_dicc)
    # Rename all the columns to look pretty.
    het_df.rename(
        columns={
            'grp': 'Group',
            'obs': 'Genome-Wide Heterozygosity',
        }, inplace=True,
    )
    # Export the dataframe.
    het_df.to_csv(f'./dataframes/afr_archaic_heterozygosity_genome_wide.csv.gz', index=False)
    return het_df

# Define a function to plot the genome-wide distribution of heterozygous sites.
def plot_arc_tgp_gw_het():
    # Intialize the archaic labels.
    arc_labs = {
        'DEN': 'Denisovan', 'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.', 'VIN': 'Vindija Nean.',
    }
    arc_list = [arc_labs[arc] for arc in arc_labs]
    # Intialize the tgp supepopulation list.
    spop_list = ['AFR', 'SAS', 'EAS', 'EUR', 'AMR']
    # Load the genome-wide results.
    arc_het_info = load_arc_gw_het_info()
    tgp_het_info = load_tgp_gw_het_info()
    # Intialize lists for plotting.
    arc_het_counts = [arc_het_info[arc][0] for arc in arc_labs]
    tgp_het_counts = [tgp_het_info[spop] for spop in spop_list]
    xtick_labs = arc_list + spop_list
    x_pos = np.arange(len(xtick_labs))
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the figure and axes.
    fig = plt.figure(
        figsize=(6, 4), dpi=300,
        facecolor='white', constrained_layout=True,
    )
    ax = fig.add_subplot(111)
    # For every other x-axis label.
    for loc in np.arange(1, len(xtick_labs), 2):
        # Plot an alternating background.
        ax.axvspan(loc-0.5, loc+0.5, facecolor='black', alpha=0.05, zorder=1)
    # Plot the archaic results.
    ax.scatter(x_pos[:4], arc_het_counts, color='black', marker='X', s=50)
    # Plot the tgp distributions.
    vp = ax.violinplot(tgp_het_counts, positions=x_pos[4:], widths=0.9, showextrema=False)
    # Adjust the violin plot color.
    for pc in vp['bodies']:
        pc.set_facecolor('grey')
        pc.set_edgecolor('grey')
        pc.set_alpha(0.75)
    # For every superpopulation.
    for i, spop_counts in enumerate(tgp_het_counts):
        # Generate some jitter to the x-axis.
        jitter = np.random.normal(x_pos[4:][i], 0.1, size=spop_counts.size)
        # Plot the points!
        ax.scatter(
            jitter, spop_counts, facecolor='none', alpha=0.5,
            edgecolor='black', marker='X', s=5, linewidths=0.25, zorder=5,
        )
    # Set the x-axis limits, tick positions, and labels.
    ax.set_xlim(x_pos[0]-0.5, x_pos[-1]+0.5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(xtick_labs, rotation=45, ha='right', rotation_mode='anchor')
    # Plot the y-axis limits.
    ax.set_ylabel('Number of Heterozygous Sites (Genome-Wide)')
    plt.savefig(
        f'./supp_figures/png/tgp_superpopulation_archaic_heterozygous_sites_genome_wide.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/tgp_superpopulation_archaic_heterozygous_sites_genome_wide.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/tgp_superpopulation_archaic_heterozygous_sites_genome_wide.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to load the observed number of heterozygous sites for the archaics.
def arc_obs_het(gt_dicc):
    # Intialize a dictionary with the observed number of heterozygous sites.
    het_dicc = {key: gt.count_het(axis=0)[0] for key, gt in gt_dicc.items()}
    return het_dicc

# Define a function to load the windowed heterozygous site counts for the archaics.
def load_arc_het_windows(window_size):
    het_dicc = {
        arc: np.concatenate([np.loadtxt(f'../muc19_results/{arc.lower()}_masked_no_aa/archaic_het_sites_heterozygosity_chr{chrom}_{window_size}kb.txt.gz') for chrom in range(1, 23)])[:, 0]
        for arc in ['DEN', 'ALT', 'CHA', 'VIN']
    }
    return het_dicc

# Define a function to compile the archaic heterozygous per window summary.
def compile_arc_het_summary(obs_dicc, wind_dicc, window_size):
    # Intialize the archaic labels.
    arc_labs = {
        'DEN': 'Denisovan', 'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.', 'VIN': 'Vindija Nean.',
    }
    # Intialize a dictionary to store the results.
    df_dicc = {
        'arc': [], 'obs': [],
        'wind_m': [], 'wind_s': [],
        'wind_se': [], 'wind_ci': [], 'wind_p': [],
    }
    # For every archaic.
    for arc in arc_labs:
        # Load window indicies of, comparable effective sequence length.
        var_idx, invar_idx = load_esl_qc_windows_idx(f'{arc.lower()}_masked_no_aa', window_size)
        # Intialize the window distribution for the sites that passed qc.
        winds = np.concatenate((wind_dicc[arc][var_idx], np.zeros(invar_idx.size)))
        # Determine the mean and standard deviation for non-overlapping windows.
        wind_m = np.nanmean(winds)
        wind_s = np.nanstd(winds)
        # Determine the sem and 95% ci.
        wind_se, wind_ci = sem_ci_of_mean(winds)
        # Compute the p-value.
        wind_p = (np.count_nonzero(obs_dicc[arc] <= winds) / np.sum(~np.isnan(winds)))
        # Append the results.
        df_dicc['arc'].append(arc_labs[arc])
        df_dicc['obs'].append(int(obs_dicc[arc]))
        df_dicc['wind_m'].append(wind_m)
        df_dicc['wind_s'].append(wind_s)
        df_dicc['wind_se'].append(wind_se)
        df_dicc['wind_ci'].append(wind_ci)
        df_dicc['wind_p'].append(wind_p)
    # Convert the dictionary to a dataframe.
    het_df = pd.DataFrame(data=df_dicc)
    # Rename all the columns to look pretty.
    het_df.rename(
        columns={
            'arc': 'Archaic',
            'obs': f'Focal {window_size}kb Region (Het. Sites)',
            'wind_m': r'Non-overlapping Windows $\left( \mu \right)$',
            'wind_s': r'Non-overlapping Windows $\left( \sigma \right)$',
            'wind_se': r'Non-overlapping Windows $\left( SEM \right)$',
            'wind_ci': r'Non-overlapping Windows $\left( \pm CI_{95\%} \right)$',
            'wind_p': r'$P-value$',
        }, inplace=True,
    )
    # Export the dataframe.
    het_df.to_csv(f'./dataframes/archaic_heterozygous_sites_{window_size}kb.csv.gz', index=False)
    return het_df


# Define a function to plot the distribution of heterozygous sites per archaic.
def plot_arc_het_dist(obs_dicc, wind_dicc, window_size):
    # Intialize the archaic labels.
    arc_labs = {
        'DEN': 'Denisovan', 'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.', 'VIN': 'Vindija Nean.',
    }
    # Intialize a bin dictionary.
    bin_dicc = {742: 'fd', 72: np.arange(0, 177, 2)}
    # Determine the significance threshold.
    sig_thresh = 0.05 
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        2, 2, figsize=(8, 4), dpi=300,
        facecolor='white',
        sharex=True, sharey=True,
        constrained_layout=True,
    )
    # Flatten the axes.
    axes = axes.flatten()
    # For every archaic.
    for i, arc in enumerate(arc_labs):
        # Load window indicies of, comparable effective sequence length.
        var_idx, invar_idx = load_esl_qc_windows_idx(f'{arc.lower()}_masked_no_aa', window_size)
        # Extract the observed value and windowed distribution.
        obs = obs_dicc[arc]
        dist = np.concatenate((wind_dicc[arc][var_idx], np.zeros(invar_idx.size)))
        # Compute the p-value.
        pval = (np.count_nonzero(obs <= dist) / np.sum(~np.isnan(dist)))
        # If the p-value is significant.
        if pval < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (obs, '*', 25, 0.49)
        # Else, it is not significant.
        else:
            # Update the dictionary.
            sig_dicc[i] = (obs, r'$ns$', 10, 0.5225)
        # Plot the distribution.
        axes[i].hist(
            dist, bins=bin_dicc[window_size],
            histtype='stepfilled', color='gray',
        )
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlabel(f'{arc_labs[arc]} Het. Sites per {window_size}kb Window')
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/archaic_heterozygous_sites_{window_size}kb.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/archaic_heterozygous_sites_{window_size}kb.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/archaic_heterozygous_sites_{window_size}kb.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot!
    plt.show()
    return

# Define a function to load the observed number of heterozygous sites for the tgp focal groups.
def tgp_obs_het(gt):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary of sample indicies.
    samp_idx_dicc = {
        'AFR': tgp_meta_df[tgp_meta_df['SUPERPOP'] == 'AFR'].index.values,
        'HET': np.loadtxt('../meta_data/tgp_het_den_like_ind_idx_72kb.txt.gz', dtype=int),
        'HOM': np.loadtxt('../meta_data/tgp_hom_den_like_ind_idx_72kb.txt.gz', dtype=int),
    }
    # Compute the observed number of heterozygous sites.
    het_dicc = {key: gt.take(idx, axis=1).count_het(axis=0) for key, idx in samp_idx_dicc.items()}
    return het_dicc 

# Define a function to load the windowed heterozygous sites for the tgp focal groups.
def load_tgp_het_windows(window_size):
    het_dicc = {
        grp: np.concatenate([np.loadtxt(f'../muc19_results/tgp_mod_no_aa/{grp.lower()}_het_sites_chr{chrom}_{window_size}kb.txt.gz') for chrom in range(1, 23)])
        for grp in ['AFR', 'HET', 'HOM']
    }
    return het_dicc

# Define a function to compile the the tgp focal groups heterozygous per window summary.
def compile_tgp_het_summary(obs_dicc, wind_dicc, window_size):
    # Intialize a dictionary of group labels.
    group_labs = {
        'AFR': 'African Inds.',
        'HET': r'Inds. with One $Denisovan-like$ Hap.',
        'HOM': r'Inds. with Two $Denisovan-like$ Haps.',
    }
    # Intialize a dictionary to store the results.
    df_dicc = {
        'grp': [], 'obs': [],
        'wind_m': [], 'wind_s': [],
        'wind_se': [], 'wind_ci': [], 'wind_p': [],
    }
    # Load window indicies of, comparable effective sequence length.
    wind_idx = load_esl_qc_windows_idx('tgp_mod_no_aa', window_size)
    # For every archaic.
    for grp in group_labs:
        # Intialize the observed value and the window distribution for the sites that passed qc.
        obs = np.nanmean(obs_dicc[grp])
        winds = np.nanmean(wind_dicc[grp][wind_idx], axis=1)
        # Determine the mean and standard deviation for non-overlapping windows.
        wind_m = np.nanmean(winds)
        wind_s = np.nanstd(winds)
        # Determine the sem and 95% ci.
        wind_se, wind_ci = sem_ci_of_mean(winds)
        # Compute the p-value.
        if grp == 'HOM':
            wind_p = (np.count_nonzero(obs >= winds) / np.sum(~np.isnan(winds)))
        else:
            wind_p = (np.count_nonzero(obs <= winds) / np.sum(~np.isnan(winds)))
        # Append the results.
        df_dicc['grp'].append(group_labs[grp])
        df_dicc['obs'].append(obs)
        df_dicc['wind_m'].append(wind_m)
        df_dicc['wind_s'].append(wind_s)
        df_dicc['wind_se'].append(wind_se)
        df_dicc['wind_ci'].append(wind_ci)
        df_dicc['wind_p'].append(wind_p)
    # Convert the dictionary to a dataframe.
    het_df = pd.DataFrame(data=df_dicc)
    # Rename all the columns to look pretty.
    het_df.rename(
        columns={
            'grp': 'Group',
            'obs': f'Focal {window_size}kb Region (Average Number of Het. Sites)',
            'wind_m': r'Non-overlapping Windows $\left( \mu \right)$',
            'wind_s': r'Non-overlapping Windows $\left( \sigma \right)$',
            'wind_se': r'Non-overlapping Windows $\left( SEM \right)$',
            'wind_ci': r'Non-overlapping Windows $\left( \pm CI_{95\%} \right)$',
            'wind_p': r'$P-value$',
        }, inplace=True,
    )
    # Export the dataframe.
    het_df.to_csv(f'./dataframes/afr_het_hom_heterozygous_sites_{window_size}kb.csv.gz', index=False)
    return het_df

# Define a function to plot the distribution of heterozygous sites per archaic.
def plot_tgp_het_dist(obs_dicc, wind_dicc, window_size):
    # Intialize a dictionary of group labels.
    group_labs = {
        'AFR': 'Avg. Number of Het. Sites Among'+'\n'+f'African Inds. per {window_size}kb Window',
        'HET': 'Avg. Number of Het. Sites Among Inds. with'+'\n'+fr'One $Denisovan-like$ Hap. per {window_size}kb Window',
        'HOM': 'Avg. Number of Het. Sites Among Inds. with'+'\n'+fr'Two $Denisovan-like$ Haps. per {window_size}kb Window',
    }
    # Intialize a bin dictionary.
    bin_dicc = {742: 'fd', 72: np.arange(0, 302, 2)}
    # Determine the significance threshold.
    sig_thresh = 0.05 
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Load window indicies of, comparable effective sequence length.
    wind_idx = load_esl_qc_windows_idx('tgp_mod_no_aa', window_size)
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 3, figsize=(9, 3.5), dpi=300,
        facecolor='white',
        sharex=True, sharey=True,
        constrained_layout=True,
    )
    # For every focal group.
    for i, grp in enumerate(group_labs):
        # Extract the observed value and windowed distribution.
        obs = np.nanmean(obs_dicc[grp])
        dist = np.nanmean(wind_dicc[grp][wind_idx], axis=1)
        # Compute the p-value.
        if grp == 'HOM':
            pval = (np.count_nonzero(obs >= dist) / np.sum(~np.isnan(dist)))
        else:
            pval = (np.count_nonzero(obs <= dist) / np.sum(~np.isnan(dist)))
        # If the p-value is significant.
        if pval < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (obs, '*', 25, 0.49)
        # Else, it is not significant.
        else:
            # Update the dictionary.
            sig_dicc[i] = (obs, r'$ns$', 10, 0.5225)
        # Plot the distribution.
        axes[i].hist(
            dist, bins=bin_dicc[window_size],
            histtype='stepfilled', color='gray',
        )
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlabel(group_labs[grp])
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/afr_het_hom_heterozygous_sites_{window_size}kb.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/afr_het_hom_heterozygous_sites_{window_size}kb.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/afr_het_hom_heterozygous_sites_{window_size}kb.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot!
    plt.show()
    return



###########
### PBS ###
###########

# Define a function to compute Hudson's estimator of Fst per site
def calc_fst_per_site(ac1, ac2):
    # Compute.
    num, den = allel.hudson_fst(ac1, ac2)
    raw_fst = num / den
    # Correct for when the denominator is zero.
    raw_fst = np.where(den == 0, 0, raw_fst)
    # Correct for negative Fst values.
    raw_fst = np.where(raw_fst < 0, 0, raw_fst)
    # Correct for Fst values of 1 that will lead to inf.
    raw_fst = np.where(raw_fst == 1, 0.99999, raw_fst)
    return raw_fst

# Define a function to calculate PBS.
def calc_pbs_per_site(gt, pop_a, pop_b, pop_c):
    # Determine the sample list.
    samp_list = np.concatenate([pop_a, pop_b, pop_c])
    # Determine if any of the sites are segregating.
    seg_mask = gt.take(samp_list, axis=1).count_alleles().is_segregating()
    # Determine allele counts.
    a_ac = gt.take(pop_a, axis=1).count_alleles()
    b_ac = gt.take(pop_b, axis=1).count_alleles()
    c_ac = gt.take(pop_c, axis=1).count_alleles()
    # Calculate Fst with corrections.
    a_b_fst = calc_fst_per_site(a_ac, b_ac)
    a_c_fst = calc_fst_per_site(a_ac, c_ac)
    c_b_fst = calc_fst_per_site(c_ac, b_ac)
    # Calculate PBS.
    pbs = (
        ((np.log(1.0 - a_b_fst) * -1.0) +\
         (np.log(1.0 - a_c_fst) * -1.0) -\
         (np.log(1.0 - c_b_fst) * -1.0)) / 2.0
    )
    # Set the sites that are not segregating to undefined.
    pbs[~seg_mask] = np.nan
    return np.where(pbs < 0, 0, pbs)

# Define a  function to compute Hudson's Fst estimator for a region.
def calc_fst_per_region(ac1, ac2):
    # Compute.
    num, den = allel.hudson_fst(ac1, ac2)
    # Account for the denominator being zero.
    fst = np.nansum(num) / np.nansum(den) if np.nansum(den) != 0 else 0
    # Correct for negative Fst values.
    return max(fst, 0)

# Define a function to calculate PBS per region.
def calc_pbs_per_region(gt, pop_a, pop_b, pop_c):
    # Determine the sample list.
    samp_list = np.concatenate([pop_a, pop_b, pop_c])
    # Determine if any of the sites are segregating.
    seg_mask = gt.take(samp_list, axis=1).count_alleles().is_segregating()
    # If no sites are segregating.
    if seg_mask.sum() == 0:
        # Return undefinded.
        return np.nan
    # Else there are segreagting sites.
    else:
        # Determine allele counts.
        a_ac = gt.take(pop_a, axis=1).count_alleles()
        b_ac = gt.take(pop_b, axis=1).count_alleles()
        c_ac = gt.take(pop_c, axis=1).count_alleles()
        # Calculate Fst.
        a_b_fst = calc_fst_per_region(a_ac, b_ac)
        a_c_fst = calc_fst_per_region(a_ac, c_ac)
        c_b_fst = calc_fst_per_region(c_ac, b_ac)
        # Correct for Fst values of 1 that will lead to inf.
        a_b_fst = min(a_b_fst, 0.99999)
        a_c_fst = min(a_c_fst, 0.99999)
        c_b_fst = min(c_b_fst, 0.99999)
        # Calculate PBS.
        pbs = (
            ((np.log(1.0 - a_b_fst) * -1.0) +\
             (np.log(1.0 - a_c_fst) * -1.0) -\
             (np.log(1.0 - c_b_fst) * -1.0)) / 2.0
        )
    return max(pbs, 0)

# Define a function to compute pbs per region results for mxl.
def compute_mxl_per_region_pbs_all_snps(gt):
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {
        'MXL_NAT': np.loadtxt('../amr_lai/anc_props/mxl_nat_idx.txt.gz', dtype=int),
        'MXL_NOT': np.loadtxt('../amr_lai/anc_props/mxl_not_idx.txt.gz', dtype=int),
    }
    # For every population.
    for pop in ['MXL', 'CEU', 'CHB']:
        # Append the dictionary with sample indicies.
        samp_idx_dicc[pop] = meta_df[meta_df['POP'] == pop].index.values
    # Intialize a dictionary.
    pbs_dicc = {}
    # For ever population.
    for pop in ['MXL', 'MXL_NAT', 'MXL_NOT']:
        # Update the dictionary.
        pbs_dicc[pop] = calc_pbs_per_region(
            gt=gt, pop_a=samp_idx_dicc[pop], pop_b=samp_idx_dicc['CHB'], pop_c=samp_idx_dicc['CEU'],
        )
    return pbs_dicc

# Define a function to load pbs per region per window results for mxl.
def load_mxl_per_region_pbs_all_snps_windows(window_size):
    # Load the window results.
    mxl_pbs = [np.loadtxt(f'../muc19_results/tgp_mod_no_aa/mxl_chb_ceu_pbs_partitions_chr{chrom}_{window_size}kb.txt.gz') for chrom in range(1, 23)]
    # Concatenate the chromosome results.
    mxl_winds = np.concatenate(mxl_pbs)
    # Initialize a dictionary.
    pbs_winds = {}
    # For every ancestry partition.
    for i, mxl in enumerate(['MXL', 'MXL_NAT', 'MXL_NOT']):
        # Fill the dictionary.
        pbs_winds[mxl] = np.where(mxl_winds[:, i] < 0, 0, mxl_winds[:, i])
    return pbs_winds

# Define a function to compile the pbs per region results for mxl.
def compile_mxl_per_region_pbs_all_snps_summary(obs_dicc, wind_dicc, window_size):
    # Intialize a dictionary for the amr labels.
    amr_labs = {
        'MXL': r'$A =$ All MXL Inds.', 'MXL_NAT': r'$A =$ MXL Inds. >50% IAA',
        'MXL_NOT': r'$A =$ MXL Inds. <50% IAA',
    }
    # Intialize dictionaries to store the results.
    df_dicc = {
        'pop': [], 'obs': [],
        'wind_m': [], 'wind_s': [],
        'wind_se': [], 'wind_ci': [],
        'wind_p': [],
    }
    # Load the good window indicies.
    wind_idx = load_esl_qc_windows_idx(prefix='tgp_mod_no_aa', window_size=window_size)
    # For every population.
    for focal_pop in amr_labs:
        # Grab the superpopulation results.
        winds = wind_dicc[focal_pop][wind_idx]
        obs = obs_dicc[focal_pop]
        # Determine the mean and standard deviation for non-overlapping windows.
        wind_m = np.nanmean(winds)
        wind_s = np.nanstd(winds)
        # Determine the sem and 95% ci.
        wind_se, wind_ci = sem_ci_of_mean(winds)
        # Compute the p-value.
        wind_p = (np.count_nonzero(obs <= winds) / np.sum(~np.isnan(winds)))
        # Update the dataframe dictionary.
        df_dicc['pop'].append(amr_labs[focal_pop])
        df_dicc['obs'].append(obs)
        df_dicc['wind_m'].append(wind_m)
        df_dicc['wind_s'].append(wind_s)
        df_dicc['wind_se'].append(wind_se)
        df_dicc['wind_ci'].append(wind_ci)
        df_dicc['wind_p'].append(wind_p)
    # Convert the dictionaries to dataframes.
    pbs_df = pd.DataFrame(df_dicc)
    # Rename all the columns to look pretty.
    pbs_df.rename(
        columns={
            'pop': r'$PBS_{A:CHB:CEU}$',
            'obs': fr'Focal {window_size}kb Region $\left( PBS \right)$',
            'wind_m': fr'{window_size}kb Non-overlapping Windows $\left( \mu \right)$',
            'wind_s': fr'{window_size}kb Non-overlapping Windows $\left( \sigma \right)$',
            'wind_se': fr'{window_size}kb Non-overlapping Windows $\left( SEM \right)$',
            'wind_ci': fr'{window_size}kb Non-overlapping Windows $\left( \pm CI_{{95\%}} \right)$',
            'wind_p': fr'$P-value$',
        }, inplace=True
    )
    # Export the dataframe as a csv.
    pbs_df.to_csv(f'./dataframes/mxl_chb_ceu_pbs_per_region_{window_size}kb.csv.gz', index=False)
    return pbs_df

# Define a function to plot the archaic snp denisty windows for mxl.
def plot_mxl_per_region_pbs_all_snps(obs_742kb, winds_742kb, obs_72kb, winds_72kb):
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Load window indicies of, comparable effective sequence length.
    wind_idx_742kb = load_esl_qc_windows_idx(prefix='tgp_mod_no_aa', window_size=742)
    wind_idx_72kb = load_esl_qc_windows_idx(prefix='tgp_mod_no_aa', window_size=72)
    # Determine the significance threshold.
    sig_thresh = 0.05 
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 2, figsize=(7, 3.5), dpi=300,
        facecolor='white',
        sharex=False, sharey=False,
        constrained_layout=True,
    )
    # For each subplot.
    for i, tup in enumerate([
        (obs_742kb, winds_742kb, wind_idx_742kb, 742), (obs_72kb, winds_72kb, wind_idx_72kb, 72),
    ]):
        # Unpack.
        obs_dicc, wind_dicc, wind_idx, window_size = tup
        # Extract the observed value and windowed distribution.
        obs = obs_dicc['MXL']
        dist = wind_dicc['MXL'][wind_idx]
        # Compute the p-value.
        pval = (np.count_nonzero(obs <= dist) / np.sum(~np.isnan(dist)))
        # If the p-value is significant.
        if pval < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (obs, '*', 25, 0.49)
        # Else, it is not significant.
        else:
            # Update the dictionary.
            sig_dicc[i] = (obs, r'$ns$', 10, 0.5225)
        # Plot the distribution.
        axes[i].hist(
            dist, bins=np.arange(0, 0.16, 0.01),
            histtype='stepfilled', color='gray',
        )
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlabel(fr'$PBS_{{MXL:CHB:CEU}}$ per {window_size}kb Window')
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        './supp_figures/png/mxl_chb_ceu_pbs_per_region_windows.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/mxl_chb_ceu_pbs_per_region_windows.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/mxl_chb_ceu_pbs_per_region_windows.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot!
    plt.show()
    return

# Define a function to compute pbs per region results for mxl.
def compute_mxl_per_region_pbs_sprime_snps(gt, pos):
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {}
    # For every population...
    for pop in ['MXL', 'CHB', 'CEU']:
        # Append the dictionary with sample indicies.
        samp_idx_dicc[pop] = meta_df[meta_df['POP'] == pop].index.values
    # Intialize a list of sprime site types.
    sprime_sites = [
        'all_arc',
    ]
    # Intilize a dictionary.
    sprime_dicc = {}
    # For every sprime site.
    for site in sprime_sites:
        # Load the information.
        sprime_pos = np.loadtxt(f'../meta_data/mxl_sprime_{site}_sites_chr12.txt.gz', dtype=int)
        # Update the dictioanry.
        sprime_dicc[site] = np.isin(pos, sprime_pos)
    # Intialize a dictionary.
    pbs_dicc = {}
    # For every combination of population b and c.
    for site, mask in sprime_dicc.items():
        # Compute PBS.
        pbs_dicc[site] = calc_pbs_per_region(
            gt=gt.compress(mask, axis=0),
            pop_a=samp_idx_dicc['MXL'], pop_b=samp_idx_dicc['CHB'], pop_c=samp_idx_dicc['CEU'],
        )
    return pbs_dicc

# Define a function to load pbs per region per window results for mxl.
def load_mxl_per_region_pbs_sprime_snps_windows(window_size):
    # Intialize the snp-type dictionary.
    sprime_sites = [
        'all_arc',
    ]
    # Initialize a dictionary.
    pbs_winds = {}
    # For every snp type.
    for site in sprime_sites:
        # Load the windowed reults.
        mxl_winds = np.concatenate([
            np.loadtxt(f'../muc19_results/tgp_mod_no_aa/mxl_chb_ceu_pbs_sprime_{site}_chr{chrom}_{window_size}kb.txt.gz') for chrom in range(1, 23)
        ])
        winds_mxl = mxl_winds[:, 0]
        # Fill the dictionary.
        pbs_winds[site] = np.where(winds_mxl < 0, 0, winds_mxl)
    return pbs_winds

# Define a function to compile the pbs per region results for mxl.
def compile_mxl_per_region_pbs_sprime_snps_summary(obs_742kb, winds_742kb, obs_72kb, winds_72kb):
    # Intialize dictionaries to store the results.
    df_dicc = {
        'region': [], 'obs': [],
        'wind_m': [], 'wind_s': [],
        'wind_se': [], 'wind_ci': [],
        'wind_p': [],
    }
    # Load window indicies of, comparable effective sequence length.
    wind_idx_742kb = load_esl_qc_windows_idx(prefix='tgp_mod_no_aa', window_size=742)
    wind_idx_72kb = load_esl_qc_windows_idx(prefix='tgp_mod_no_aa', window_size=72)
    # For every sprime site.
    for obs_dicc, wind_dicc, wind_idx, window_size in [
        (obs_742kb, winds_742kb, wind_idx_742kb, '742kb'), (obs_72kb, winds_72kb, wind_idx_72kb, '72kb'),
    ]:
        # Grab the superpopulation results.
        winds = wind_dicc['all_arc'][wind_idx]
        obs = obs_dicc['all_arc']
        # Determine the mean and standard deviation for non-overlapping windows.
        wind_m = np.nanmean(winds)
        wind_s = np.nanstd(winds)
        # Determine the sem and 95% ci.
        wind_se, wind_ci = sem_ci_of_mean(winds)
        # Compute the p-value.
        wind_p = (np.count_nonzero(obs <= winds) / np.sum(~np.isnan(winds)))
        # Update the dataframe dictionary.
        df_dicc['region'].append(window_size)
        df_dicc['obs'].append(obs)
        df_dicc['wind_m'].append(wind_m)
        df_dicc['wind_s'].append(wind_s)
        df_dicc['wind_se'].append(wind_se)
        df_dicc['wind_ci'].append(wind_ci)
        df_dicc['wind_p'].append(wind_p)
    # Convert the dictionaries to dataframes.
    pbs_df = pd.DataFrame(df_dicc)
    # Rename all the columns to look pretty.
    pbs_df.rename(
        columns={
            'region': 'Focal Region',
            'obs': r'$PBS_{MXL:CHB:CEU}$',
            'wind_m': r'Non-overlapping Windows $\left( \mu \right)$',
            'wind_s': r'Non-overlapping Windows $\left( \sigma \right)$',
            'wind_se': r'Non-overlapping Windows $\left( SEM \right)$',
            'wind_ci': r'Non-overlapping Windows $\left( \pm CI_{{95\%}} \right)$',
            'wind_p': r'$P-value$',
        }, inplace=True
    )
    # Export the dataframe as a csv.
    pbs_df.to_csv(f'./dataframes/sprime_mxl_chb_ceu_pbs_per_region_742kb_and_72kb.csv.gz', index=False)
    return pbs_df

# Define a function to plot the archaic snp denisty windows for mxl.
def plot_mxl_per_region_pbs_sprime_snps(obs_742kb, winds_742kb, obs_72kb, winds_72kb):
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Load window indicies of, comparable effective sequence length.
    wind_idx_742kb = load_esl_qc_windows_idx(prefix='tgp_mod_no_aa', window_size=742)
    wind_idx_72kb = load_esl_qc_windows_idx(prefix='tgp_mod_no_aa', window_size=72)
    # Determine the significance threshold.
    sig_thresh = 0.05 
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 2, figsize=(7, 3.5), dpi=300,
        facecolor='white',
        sharex=False, sharey=False,
        constrained_layout=True,
    )
    # For each subplot.
    for i, tup in enumerate([
        (obs_742kb, winds_742kb, wind_idx_742kb, 742), (obs_72kb, winds_72kb, wind_idx_72kb, 72),
    ]):
        # Unpack.
        obs_dicc, wind_dicc, wind_idx, window_size = tup
        # Extract the observed value and windowed distribution.
        obs = obs_dicc['all_arc']
        dist = wind_dicc['all_arc'][wind_idx]
        # Compute the p-value.
        pval = (np.count_nonzero(obs <= dist) / np.sum(~np.isnan(dist)))
        # If the p-value is significant.
        if pval < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (obs, '*', 25, 0.49)
        # Else, it is not significant.
        else:
            # Update the dictionary.
            sig_dicc[i] = (obs, r'$ns$', 10, 0.5225)
        # Plot the distribution.
        axes[i].hist(
            dist, bins=np.arange(0, 0.31, 0.01),
            histtype='stepfilled', color='gray',
        )
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlabel(fr'$PBS_{{MXL:CHB:CEU}}$ per {window_size}kb Window (SPrime Sites)')
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        './supp_figures/png/sprime_mxl_chb_ceu_pbs_per_region_windows.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/sprime_mxl_chb_ceu_pbs_per_region_windows.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/sprime_mxl_chb_ceu_pbs_per_region_windows.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot!
    plt.show()
    return

# Define a function to compute pbs per region results for amr.
def compute_amr_per_region_pbs_all_snps(gt):
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the asian, european, and admixed american populations.
    asn_pops = np.unique(meta_df[
        (meta_df['SUPERPOP'] == 'SAS') | (meta_df['SUPERPOP'] == 'EAS')
    ]['POP'].values)
    eur_pops = np.unique(meta_df[meta_df['SUPERPOP'] == 'EUR']['POP'].values)
    amr_pops = np.unique(meta_df[meta_df['SUPERPOP'] == 'AMR']['POP'].values)
    # Find all unique combinations.
    amr_asn_eur_combos = list(product(amr_pops, asn_pops, eur_pops))
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {}
    # For every population...
    for pop in np.concatenate((amr_pops, asn_pops, eur_pops)):
        # Append the dictionary with sample indicies.
        samp_idx_dicc[pop] = meta_df[meta_df['POP'] == pop].index.values
    # Intialize a dictionary.
    pbs_dicc = {}
    # For every combination of population b and c.
    for a_pop, b_pop, c_pop in amr_asn_eur_combos:
        # Compute PBS.
        pbs_dicc[(a_pop, b_pop, c_pop)] = calc_pbs_per_region(
            gt=gt, pop_a=samp_idx_dicc[a_pop], pop_b=samp_idx_dicc[b_pop], pop_c=samp_idx_dicc[c_pop],
        )
    return pbs_dicc

# Define a function to load pbs per region per window results for mxl.
def load_amr_per_region_pbs_all_snps_windows(window_size):
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the asian, european, and admixed american populations.
    asn_pops = np.unique(meta_df[
        (meta_df['SUPERPOP'] == 'SAS') | (meta_df['SUPERPOP'] == 'EAS')
    ]['POP'].values)
    eur_pops = np.unique(meta_df[meta_df['SUPERPOP'] == 'EUR']['POP'].values)
    amr_pops = np.unique(meta_df[meta_df['SUPERPOP'] == 'AMR']['POP'].values)
    # Find all unique combinations.
    amr_asn_eur_combos = list(product(amr_pops, asn_pops, eur_pops))
    # Initialize a dictionary.
    pbs_winds = {}
    # For every combination of population b and c.
    for a_pop, b_pop, c_pop in amr_asn_eur_combos:
        # Load the windowed reults.
        amr_winds = np.concatenate([
            np.loadtxt(f'../muc19_results/tgp_mod_no_aa/{a_pop.lower()}_{b_pop.lower()}_{c_pop.lower()}_pbs_chr{chrom}_{window_size}kb.txt.gz') for chrom in range(1, 23)
        ])
        # Fill the dictionary.
        pbs_winds[(a_pop, b_pop, c_pop)] = np.where(amr_winds < 0, 0, amr_winds)
    return pbs_winds

# Define a function to compile the pbs per region results for mxl.
def compile_amr_per_region_pbs_all_snps_summary(obs_dicc, wind_dicc, window_size):
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the asian, european, and admixed american populations.
    asn_pops = np.unique(meta_df[
        (meta_df['SUPERPOP'] == 'SAS') | (meta_df['SUPERPOP'] == 'EAS')
    ]['POP'].values)
    eur_pops = np.unique(meta_df[meta_df['SUPERPOP'] == 'EUR']['POP'].values)
    amr_pops = np.unique(meta_df[meta_df['SUPERPOP'] == 'AMR']['POP'].values)
    # Find all unique combinations.
    amr_asn_eur_combos = list(product(amr_pops, asn_pops, eur_pops))
    # Intialize dictionaries to store the results.
    df_dicc = {
        'pop_a': [], 'pop_b': [], 'pop_c': [],
        'obs': [], 'wind_m': [], 'wind_s': [],
        'wind_se': [], 'wind_ci': [], 'wind_p': [],
        
    }
    # Load the good window indicies.
    wind_idx = load_esl_qc_windows_idx(prefix='tgp_mod_no_aa', window_size=window_size)
    # For every combination of population b and c.
    for a_pop, b_pop, c_pop in amr_asn_eur_combos:
        # Grab the superpopulation results.
        winds = wind_dicc[(a_pop, b_pop, c_pop)][wind_idx]
        obs = obs_dicc[(a_pop, b_pop, c_pop)]
        # Determine the mean and standard deviation for non-overlapping windows.
        wind_m = np.nanmean(winds)
        wind_s = np.nanstd(winds)
        # Determine the sem and 95% ci.
        wind_se, wind_ci = sem_ci_of_mean(winds)
        # Compute the p-value.
        wind_p = (np.count_nonzero(obs <= winds) / np.sum(~np.isnan(winds)))
        # Update the dataframe dictionary.
        df_dicc['pop_a'].append(a_pop)
        df_dicc['pop_b'].append(b_pop)
        df_dicc['pop_c'].append(c_pop)
        df_dicc['obs'].append(obs)
        df_dicc['wind_m'].append(wind_m)
        df_dicc['wind_s'].append(wind_s)
        df_dicc['wind_se'].append(wind_se)
        df_dicc['wind_ci'].append(wind_ci)
        df_dicc['wind_p'].append(wind_p)
    # Convert the dictionaries to dataframes.
    pbs_df = pd.DataFrame(df_dicc)
    # Rename all the columns to look pretty.
    pbs_df.rename(
        columns={
            'pop_a': r'$A$ Pop.',
            'pop_b': r'$B$ Pop.',
            'pop_c': r'$C$ Pop.',
            'obs': fr'Focal {window_size}kb Region $\left( PBS_{{A:B:C}} \right)$',
            'wind_m': fr'{window_size}kb Non-overlapping Windows $\left( \mu \right)$',
            'wind_s': fr'{window_size}kb Non-overlapping Windows $\left( \sigma \right)$',
            'wind_se': fr'{window_size}kb Non-overlapping Windows $\left( SEM \right)$',
            'wind_ci': fr'{window_size}kb Non-overlapping Windows $\left( \pm CI_{{95\%}} \right)$',
            'wind_p': fr'$P-value$',
        }, inplace=True
    )
    # Export the dataframe as a csv.
    pbs_df.to_csv(f'./dataframes/amr_asn_eur_pbs_per_region_{window_size}kb.csv.gz', index=False)
    return pbs_df

# Define a function to load all pbs values for mxl.
def load_gw_mxl_chb_ceu_pbs_all_snps():
    # Intialize the pbs data.
    all_pbs = [np.loadtxt(f'../muc19_results/tgp_mod_no_aa/mxl_chb_ceu_pbs_chr{chrom}.txt.gz') for chrom in range(1, 23)]
    # Concatenate the chromosome results.
    pbs_chroms = np.concatenate(all_pbs)[:, 0]
    return np.where(pbs_chroms < 0, 0, pbs_chroms)

# Define a function to compile the pbs information for mxl
def init_mxl_chb_ceu_pbs_per_site_info(tgp_gt_742kb, tgp_pos_742kb, tgp_arcs_gt_742kb, tgp_arcs_pos_742kb):
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    pop_dicc = {
        'MXL': tgp_df[tgp_df['POP'] == 'MXL'].index.values,
        'CHB': tgp_df[tgp_df['POP'] == 'CHB'].index.values,
        'CEU': tgp_df[tgp_df['POP'] == 'CEU'].index.values,
        'MXL_NAT': np.loadtxt('../amr_lai/anc_props/mxl_nat_idx.txt.gz', dtype=int),
        'MXL_NOT': np.loadtxt('../amr_lai/anc_props/mxl_not_idx.txt.gz', dtype=int),
    }
    # Classify each snp type.
    _, _, arc_hum_snps = tgp_classify_snps(tgp_arcs_gt_742kb, 'MXL', 0.01)
    # Intialize a snp type mask dictionary.
    snp_dicc = {
        'ARC': np.isin(tgp_pos_742kb, tgp_arcs_pos_742kb[arc_hum_snps['ARC']]),
        'DEN': np.isin(tgp_pos_742kb, tgp_arcs_pos_742kb[arc_hum_snps['DEN']]),
        'NEA': np.isin(tgp_pos_742kb, tgp_arcs_pos_742kb[arc_hum_snps['NEA']]),
        'SHR': np.isin(tgp_pos_742kb, tgp_arcs_pos_742kb[arc_hum_snps['SHR']]),
        'HOM': np.isin(tgp_pos_742kb, tgp_arcs_pos_742kb[arc_hum_snps['HOM']]),
    }
    snp_dicc['HUM'] = ~(snp_dicc['ARC'] | snp_dicc['HOM'])
    # Load the genome wide results
    gw_pbs = load_gw_mxl_chb_ceu_pbs_all_snps()
    # Clean the pbs values.
    cleaned_gw_pbs = gw_pbs[~np.isnan(gw_pbs)]
    # Compute the outlier threshold.
    gw_out_thresh = np.percentile(cleaned_gw_pbs, 99.95)
    # Compute the per site pbs values.
    mxl_pbs_per_site = calc_pbs_per_site(tgp_gt_742kb, pop_dicc['MXL'], pop_dicc['CHB'], pop_dicc['CEU'])
    mxl_nat_pbs_per_site = calc_pbs_per_site(tgp_gt_742kb, pop_dicc['MXL_NAT'], pop_dicc['CHB'], pop_dicc['CEU'])
    mxl_not_pbs_per_site = calc_pbs_per_site(tgp_gt_742kb, pop_dicc['MXL_NOT'], pop_dicc['CHB'], pop_dicc['CEU'])
    # Clean the pbs information.
    cleaned_mxl_pbs_per_site_mask = ~np.isnan(mxl_pbs_per_site)
    cleaned_mxl_pbs_per_site = mxl_pbs_per_site[cleaned_mxl_pbs_per_site_mask]
    cleaned_mxl_nat_pbs_per_site = mxl_nat_pbs_per_site[cleaned_mxl_pbs_per_site_mask]
    cleaned_mxl_not_pbs_per_site = mxl_not_pbs_per_site[cleaned_mxl_pbs_per_site_mask]
    cleaned_pbs_masks = {
        'IS_ARC': snp_dicc['ARC'][cleaned_mxl_pbs_per_site_mask],
        'IS_DEN': snp_dicc['DEN'][cleaned_mxl_pbs_per_site_mask],
        'IS_NEA': snp_dicc['NEA'][cleaned_mxl_pbs_per_site_mask],
        'IS_SHR': snp_dicc['SHR'][cleaned_mxl_pbs_per_site_mask],
        'IS_HOM': snp_dicc['HOM'][cleaned_mxl_pbs_per_site_mask],
        'IS_HUM': snp_dicc['HUM'][cleaned_mxl_pbs_per_site_mask],
    }
    cleaned_pos = tgp_pos_742kb[cleaned_mxl_pbs_per_site_mask]
    # Compute alternative allele frequencies.
    cleaned_mxl_aaf = calc_alt_freqs(tgp_gt_742kb.take(pop_dicc['MXL'], axis=1).compress(cleaned_mxl_pbs_per_site_mask, axis=0))
    cleaned_mxl_nat_aaf = calc_alt_freqs(tgp_gt_742kb.take(pop_dicc['MXL_NAT'], axis=1).compress(cleaned_mxl_pbs_per_site_mask, axis=0))
    cleaned_mxl_not_aaf = calc_alt_freqs(tgp_gt_742kb.take(pop_dicc['MXL_NOT'], axis=1).compress(cleaned_mxl_pbs_per_site_mask, axis=0))
    cleaned_chb_aaf = calc_alt_freqs(tgp_gt_742kb.take(pop_dicc['CHB'], axis=1).compress(cleaned_mxl_pbs_per_site_mask, axis=0))
    cleaned_ceu_aaf = calc_alt_freqs(tgp_gt_742kb.take(pop_dicc['CEU'], axis=1).compress(cleaned_mxl_pbs_per_site_mask, axis=0))
    # Determine what positions are in the 72kb region.
    is_72kb = (40759001 <= cleaned_pos) & (cleaned_pos <= 40831000)
    # Intialize lists to store the percentile ranks and uncorrected p-values.
    pbs_prs, unc_pvals = [], []
    # For every observed pbs values.
    for obs_pbs in cleaned_mxl_pbs_per_site:
        # Update the lists.
        pbs_prs.append(stats.percentileofscore(cleaned_gw_pbs, obs_pbs, kind='strict'))
        unc_pvals.append((cleaned_gw_pbs >= obs_pbs).sum() / cleaned_gw_pbs.size)
    # Convert the lists to arrays.
    pbs_prs = np.array(pbs_prs)
    unc_pvals = np.array(unc_pvals)
    # Deterime what snps are outliers and significnat based on the uncorrected p-value.
    out_mask = cleaned_mxl_pbs_per_site > gw_out_thresh
    unc_mask = unc_pvals < 0.05
    # Compute adjusted p-values using the Bonferroni correction.
    bon_pvals = np.minimum((unc_pvals * unc_pvals.size), 1)
    bon_mask = bon_pvals < 0.05
    # Compute adjusted p-values using the Benjamini-Hochberg procedure.
    _, bhp_pvals, _, _ = multipletests(unc_pvals, alpha=0.01, method='fdr_bh')
    bhp_mask = bhp_pvals < 0.01
    # Intialize a dictionary.
    pbs_info = {
        'POS': cleaned_pos, 'IS_72KB': is_72kb, 'IS_ARC': cleaned_pbs_masks['IS_ARC'], 'IS_DEN': cleaned_pbs_masks['IS_DEN'],
        'IS_NEA': cleaned_pbs_masks['IS_NEA'], 'IS_SHR': cleaned_pbs_masks['IS_SHR'], 'IS_HUM': cleaned_pbs_masks['IS_HUM'],
        'IS_HOM': cleaned_pbs_masks['IS_HOM'], 'MXL_AAF': cleaned_mxl_aaf, 'PBS': cleaned_mxl_pbs_per_site, 'PR': pbs_prs,
        'UNC_PVAL': unc_pvals, 'BON_PVAL': bon_pvals, 'BHP_PVAL': bhp_pvals,
        'IS_OUT': out_mask, 'IS_UNC': unc_mask, 'IS_BON': bon_mask, 'IS_BHP': bhp_mask,
    }
    mxl_pbs_partitions = {
        'POS': cleaned_pos, 'IS_ARC': cleaned_pbs_masks['IS_ARC'], 'IS_DEN': cleaned_pbs_masks['IS_DEN'], 'IS_NEA': cleaned_pbs_masks['IS_NEA'], 
        'IS_SHR': cleaned_pbs_masks['IS_SHR'], 'IS_HUM': cleaned_pbs_masks['IS_HUM'], 'IS_HOM': cleaned_pbs_masks['IS_HOM'], 
        'MXL': cleaned_mxl_pbs_per_site, 'MXL_NAT': cleaned_mxl_nat_pbs_per_site, 'MXL_NOT': cleaned_mxl_not_pbs_per_site,
        'MXL_AAF': cleaned_mxl_aaf, 'MXL_NAT_AAF': cleaned_mxl_nat_aaf, 'MXL_NOT_AAF': cleaned_mxl_not_aaf, 'CHB_AAF': cleaned_chb_aaf, 'CEU_AAF': cleaned_ceu_aaf,
    }
    # Intialize a summary dataframe.
    pbs_summary_df = pd.DataFrame({
        'Description': [
            'Number of PBS SNPs Genome Wide',
            'Number of PBS SNPs 742kb',
            'Number of PBS SNPs 72kb',
            '99.95th Percentile',
            'Number of SNPs > 99.95th Percentile',
            'Uncorrected Signficance Level',
            'Number of Sig. SNPs Before Correcting for Multiple Comparisons',
            'Bonferroni-Corrected Significance Level',
            'Number of Sig. SNPs After Bonferroni Correction',
            'False Discovery Rate',
            'Number of Sig. SNPs After Benjamini-Hochberg Procedure',
        ],
        'Value': [
            cleaned_gw_pbs.size,
            cleaned_pos.size,
            is_72kb.sum(),
            gw_out_thresh,
            out_mask.sum(),
            0.05,
            unc_mask.sum(),
            0.05 / cleaned_pos.size,
            bon_mask.sum(),
            0.01,
            bhp_mask.sum()
        ],
    })
    pbs_summary_df['Value'] = pbs_summary_df['Value'].map(lambda x: f'{x:.17f}')
    # Convert the partioned pbs dictionary to a dataframe.
    mxl_pbs_partitions_df = pd.DataFrame(mxl_pbs_partitions)
    # Rename all the columns to look pretty.
    mxl_pbs_partitions_df.rename(
        columns={
            'POS': 'Position', 'IS_ARC': 'Archaic SNP', 'IS_DEN': 'Denisovan-specific SNP',
            'IS_NEA': 'Neanderthal-specific SNP', 'IS_SHR': 'Shared Archaic SNP',
            'IS_HUM': 'Human-specific SNP', 'IS_HOM': 'Shared Hominin SNP', 
            'MXL': r'$A =$ All MXL Inds.', 'MXL_NAT': r'$A =$ MXL Inds. >50% Indigenous American Anc.',
            'MXL_NOT': r'$A =$ MXL Inds. <50% Indigenous American Anc.',
            'MXL_AAF': 'All MXL Inds. Alt. Allele Freq.',
            'MXL_NAT_AAF': 'MXL Inds. >50% Indigenous American Anc. Alt. Allele Freq.',
            'MXL_NOT_AAF': 'MXL Inds. <50% Indigenous American Anc. Alt. Allele Freq.',
            'CHB_AAF': 'CHB Alt. Allele Freq.', 'CEU_AAF': 'CEU Alt. Allele Freq.',
        }, inplace=True
    )
    # Export the dataframe as a csv.
    mxl_pbs_partitions_df.to_csv(f'./dataframes/mxl_chb_ceu_partitioned_pbs_per_snp_742kb.csv.gz', index=False)
    return mxl_pbs_partitions, mxl_pbs_partitions_df, pbs_info, gw_out_thresh, pd.DataFrame(pbs_info), pbs_summary_df

# Define a function to plot the mxl vs partiontied pbs results.
def plot_mxl_pbs_partitions(mxl_dicc):
    # Intialize a dictionary of snp type partitions.
    snp_dicc = {
        'IS_HUM': {'l': 'Human-specific', 'c': 'gray'},
        'IS_HOM': {'l': 'Shared Hominin', 'c': '#000000'},
        'IS_DEN': {'l': 'Denisovan-specific', 'c': '#E69F00'},
        'IS_NEA': {'l': 'Neanderthal-specific', 'c': '#56B4E9'},
        'IS_SHR': {'l': 'Shared Archaic', 'c': '#CC79A7'},
    }
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 2, figsize=(7, 3.5), dpi=300,
        facecolor='white',
        sharex=True, sharey=True,
        constrained_layout=True,
    )
    # For each subplot.
    for i, tup in enumerate([
        ('MXL_NAT', r'$PBS_{A:CHB:CEU}$ ($A =$ MXL Inds. >50% IAA)', 'lower right', (1, 0)),
        ('MXL_NOT', r'$PBS_{A:CHB:CEU}$ ($A =$ MXL Inds. <50% IAA)', 'upper left', (0, 1)),
    ]):
        # Unpack.
        mxl_key, y_lab, legend_loc, legend_bbox = tup
        # For all snp types.
        for snp_type, zorder in [
            ('IS_DEN', 5), ('IS_NEA', 4), ('IS_SHR', 4), ('IS_HOM', 1), ('IS_HUM', 1), 
        ]:
            # Plot the results.
            axes[i].scatter(
                mxl_dicc['MXL'][mxl_dicc[snp_type]], mxl_dicc[mxl_key][mxl_dicc[snp_type]], zorder=zorder,
                color=snp_dicc[snp_type]['c'], marker='o', edgecolor='white', linewidth=0.1, s=20,
            )
        # Plot x = y.
        lims = [
            np.min([axes[i].get_xlim(), axes[i].get_ylim()]),
            np.max([axes[i].get_xlim(), axes[i].get_ylim()]),
        ]
        axes[i].plot(lims, lims, 'k--')
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel(y_lab)
        axes[i].set_xlabel(r'$PBS_{A:CHB:CEU}$ ($A =$ All MXL Inds.)')
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
        # Add a legend.
        axes[i].legend(
            handles=[
                Line2D(
                    [0], [0], linestyle='none', marker='o', markersize=7.5, markeredgewidth=0.1,
                    color=snp_dicc[snp_type]['c'], markeredgecolor='white', label=f'{snp_dicc[snp_type]["l"]} (n = {mxl_dicc[snp_type].sum() - np.isnan(mxl_dicc[mxl_key][mxl_dicc[snp_type]]).sum()})'
                ) for snp_type in ['IS_DEN', 'IS_NEA', 'IS_SHR', 'IS_HUM', 'IS_HOM']
            ], loc=legend_loc, bbox_to_anchor=legend_bbox,
            ncol=1, frameon=True, fancybox=True, shadow=True, fontsize=8,
        )
    # For each subplot.
    for i in range(2):
        # Adjust the axes limits.
        axes[i].set_ylim(-0.05, 0.65)
        axes[i].set_xlim(-0.05, 0.65)
    
    # Export the plot.
    plt.savefig(
        './supp_figures/png/mxl_chb_ceu_partitioned_pbs_per_snp_742kb.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/mxl_chb_ceu_partitioned_pbs_per_snp_742kb.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/mxl_chb_ceu_partitioned_pbs_per_snp_742kb.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot!
    plt.show()
    return

# Define a function to summarize the per site PBS_{MXL:CHB:CEU} results.
def summarize_mxl_chb_ceu_pbs_per_site_info(info_dicc):
    # Intialize snp type label dictionary.
    snp_type_labs = {
        'ALL': 'All SNPs',
        'ARC': 'Archaic SNPs',
        'DEN': 'Denisovan-specific SNPs',
        'NEA': 'Neanderthal-specific SNPs',
        'HUM': 'Human-specific SNPs',
        'HOM': 'Shared Hominin SNPs',
    }
    # Intialize a list for the significance columns.
    sig_cols = [
        'snp_type', 'tot', 'sig',
        'min_pbs', 'max_pbs', 'mean_pbs', 'std_pbs', 'sem_pbs', 'ci_pbs',
        'min_pr', 'max_pr', 'mean_pr', 'std_pr', 'sem_pr', 'ci_pr',
    ]
    # Intialize dictionaries.
    sig_dicc = {key: [] for key in ['snp_type', 'tot', 'out', 'unc', 'bon', 'bhp']}
    out_dicc = {key: [] for key in sig_cols}
    unc_dicc = {key: [] for key in sig_cols}
    bhp_dicc = {key: [] for key in sig_cols}
    mask_dicc = {'ALL': np.full(info_dicc['POS'].size, True)}
    for snp_type in ['ARC', 'DEN', 'NEA', 'HUM', 'HOM']:
        mask_dicc[snp_type] = info_dicc[f'IS_{snp_type}']
    for sig_type in ['OUT', 'UNC', 'BON', 'BHP']:
        mask_dicc[sig_type] = info_dicc[f'IS_{sig_type}']
    # Extract the pbs values and percentile ranks.
    pbs, prs = info_dicc['PBS'], info_dicc['PR']
    # Intialize lists for the significance dictionaries.
    sig_type_dicc_list = [('OUT', out_dicc), ('UNC', unc_dicc), ('BHP', bhp_dicc)]
    # For every snp type.
    for snp_type in snp_type_labs:
        # Intialize the masks.
        snp_mask = mask_dicc[snp_type]
        # Update the summary dictionary.
        sig_dicc['snp_type'].append(snp_type_labs[snp_type])
        sig_dicc['tot'].append(snp_mask.sum())
        # For every significance type.
        for sig_col in ['out', 'unc', 'bon', 'bhp']:
            # Update the summary dictionary.
            sig_dicc[sig_col].append((snp_mask & mask_dicc[sig_col.upper()]).sum())
        # For every significance dictionary.
        for sig_type, sig_type_dicc in sig_type_dicc_list:
            # Intialize the mask and exctarct the pbs and percentile rank values.
            snp_sig_mask = snp_mask & mask_dicc[sig_type]
            snp_sig_pbs, snp_sig_prs = pbs[snp_sig_mask], prs[snp_sig_mask]
            # Determine the sem and 95% ci.
            pbs_se, pbs_ci = sem_ci_of_mean(snp_sig_pbs)
            prs_se, prs_ci = sem_ci_of_mean(snp_sig_prs)
            # Update the dictionary.
            sig_type_dicc['snp_type'].append(snp_type_labs[snp_type])
            sig_type_dicc['tot'].append(snp_mask.sum())
            sig_type_dicc['sig'].append(snp_sig_mask.sum())
            sig_type_dicc['min_pbs'].append(snp_sig_pbs.min())
            sig_type_dicc['max_pbs'].append(snp_sig_pbs.max())
            sig_type_dicc['mean_pbs'].append(snp_sig_pbs.mean())
            sig_type_dicc['std_pbs'].append(snp_sig_pbs.std())
            sig_type_dicc['sem_pbs'].append(pbs_se)
            sig_type_dicc['ci_pbs'].append(pbs_ci)
            sig_type_dicc['min_pr'].append(snp_sig_prs.min())
            sig_type_dicc['max_pr'].append(snp_sig_prs.max())
            sig_type_dicc['mean_pr'].append(snp_sig_prs.mean())
            sig_type_dicc['std_pr'].append(snp_sig_prs.std())
            sig_type_dicc['sem_pr'].append(prs_se)
            sig_type_dicc['ci_pr'].append(prs_ci)
    # Convert the dictionaries to dataframes.
    sig_df = pd.DataFrame(sig_dicc)
    out_df = pd.DataFrame(out_dicc)
    unc_df = pd.DataFrame(unc_dicc)
    bhp_df = pd.DataFrame(bhp_dicc)
    # Rename all the columns to look pretty.
    sig_df.rename(
        columns={
            'snp_type': 'SNP Set',
            'tot': 'Total SNPs',
            'out': r'$PBS > 99.95^{th}$ Genome-Wide Percentile',
            'unc': r'Uncorrected $P-values < 0.05$',
            'bon': r'Bonferroni Adjusted $P-values < 0.05$',
            'bhp': r'Benjamini-Hochberg Adjusted $P-values < 0.01$',
        }, inplace=True
    )
    out_df.rename(
        columns={
            'snp_type': 'SNP Set',
            'tot': 'Total SNPs',
            'sig': r'$PBS > 99.95^{th}$ Genome-Wide Percentile',
            'min_pbs': r'$PBS$ (Min)',
            'max_pbs': r'$PBS$ (Max)',
            'mean_pbs': r'$PBS \; \left( \mu \right)$',
            'std_pbs': r'$PBS \; \left( \sigma \right)$',
            'sem_pbs': r'$PBS \; \left( SEM \right)$',
            'ci_pbs': r'$PBS \; \left( \pm CI_{95\%} \right)$',
            'min_pr': 'Percentile Rank (Min)',
            'max_pr': 'Percentile Rank (Max)',
            'mean_pr': r'Percentile Rank $\left( \mu \right)$',
            'std_pr': r'Percentile Rank $\left( \sigma \right)$',
            'sem_pr': r'Percentile Rank $\left( SEM \right)$',
            'ci_pr': r'Percentile Rank $\left( \pm CI_{95\%} \right)$',
        }, inplace=True
    )
    unc_df.rename(
        columns={
            'snp_type': 'SNP Set',
            'tot': 'Total SNPs',
            'sig': r'Uncorrected $P-values < 0.05$',
            'min_pbs': r'$PBS$ (Min)',
            'max_pbs': r'$PBS$ (Max)',
            'mean_pbs': r'$PBS \; \left( \mu \right)$',
            'std_pbs': r'$PBS \; \left( \sigma \right)$',
            'sem_pbs': r'$PBS \; \left( SEM \right)$',
            'ci_pbs': r'$PBS \; \left( \pm CI_{95\%} \right)$',
            'min_pr': 'Percentile Rank (Min)',
            'max_pr': 'Percentile Rank (Max)',
            'mean_pr': r'Percentile Rank $\left( \mu \right)$',
            'std_pr': r'Percentile Rank $\left( \sigma \right)$',
            'sem_pr': r'Percentile Rank $\left( SEM \right)$',
            'ci_pr': r'Percentile Rank $\left( \pm CI_{95\%} \right)$',
        }, inplace=True
    )
    bhp_df.rename(
        columns={
            'snp_type': 'SNP Set',
            'tot': 'Total SNPs',
            'sig': r'Benjamini-Hochberg Adjusted $P-values < 0.01$',
            'min_pbs': r'$PBS$ (Min)',
            'max_pbs': r'$PBS$ (Max)',
            'mean_pbs': r'$PBS \; \left( \mu \right)$',
            'std_pbs': r'$PBS \; \left( \sigma \right)$',
            'sem_pbs': r'$PBS \; \left( SEM \right)$',
            'ci_pbs': r'$PBS \; \left( \pm CI_{95\%} \right)$',
            'min_pr': 'Percentile Rank (Min)',
            'max_pr': 'Percentile Rank (Max)',
            'mean_pr': r'Percentile Rank $\left( \mu \right)$',
            'std_pr': r'Percentile Rank $\left( \sigma \right)$',
            'sem_pr': r'Percentile Rank $\left( SEM \right)$',
            'ci_pr': r'Percentile Rank $\left( \pm CI_{95\%} \right)$',
        }, inplace=True
    )
    # Export.
    sig_df.to_csv('./dataframes/pbs_empirical_significance_summary.csv.gz', index=False)
    out_df.to_csv('./dataframes/pbs_empirical_outlier_summary.csv.gz', index=False)
    bhp_df.to_csv('./dataframes/pbs_empirical_bhp_significance_summary.csv.gz', index=False)
    return sig_df, out_df, unc_df, bhp_df

# Define a function to plot the per site PBS_{MXL:CHB:CEU} results for the 742kb region.
def plot_pbs_mxl_chb_ceu_info_742kb(info_dicc, plot_key, thresh):
    # Intialize a dictionary for plotting.
    plot_dicc = {
        'PBS': {
            'y_val': info_dicc['PBS'], 'y_lab': r'$PBS_{MXL:CHB:CEU}$',
            'gene_loc': -0.05, 'y_lim': (-0.1, None), 'legend_loc': (0.5, 1.1), 
            'thresh': Line2D([0], [0], color='black', linestyle='dashed', label=r'$99.95^{th}$ Percentile (Genome-Wide)'),
            'export': 'pbs_mxl_chb_ceu_742kb',
        },
        'PR': {
            'y_val': info_dicc['PR'], 'y_lab': r'$PBS_{MXL:CHB:CEU}$ Percentile Rank',
            'gene_loc': 99.945, 'y_lim': (99.94, None), 'legend_loc': (0.5, 1.15), 
            'thresh': Line2D([0], [0], color='black', linestyle='dashed', label=r'$99.95^{th}$ Percentile (Genome-Wide)'),
            'export': 'pbs_mxl_chb_ceu_percentile_ranks_742kb',
        },
        'UNC_PVAL': {
            'y_val': -1 * np.log10(info_dicc['UNC_PVAL']), 'y_lab': r'$- \log_{10}(P-value_{Uncorrected})$',
            'gene_loc': -0.5,  'y_lim': (-1, None), 'legend_loc': (0.5, 1.15),
            'thresh': Line2D([0], [0], color='black', linestyle='dashed', label=r'$P-value_{Uncorrected} = 0.05$'),
            'export': 'pbs_mxl_chb_ceu_uncorrected_pvalues_742kb',
        },
        'BON_PVAL': {
            'y_val': -1 * np.log10(info_dicc['BON_PVAL']), 'y_lab': r'$- \log_{10}(P-value_{Bonferroni})$',
            'gene_loc': -0.15,  'y_lim': (-0.3, None), 'legend_loc': (0.5, 1.15),
            'thresh': Line2D([0], [0], color='black', linestyle='dashed', label=r'$P-value_{Bonferroni} = 0.05$'),
            'export': 'pbs_mxl_chb_ceu_bonferroni_corrected_pvalues_742kb',
        },
        'BHP_PVAL': {
            'y_val': -1 * np.log10(info_dicc['BHP_PVAL']), 'y_lab': r'$- \log_{10}(P-value_{Benjamini-Hochberg})$',
            'gene_loc': -0.25,  'y_lim': (-0.5, None), 'legend_loc': (0.5, 1.15),
            'thresh': Line2D([0], [0], color='black', linestyle='dashed', label=r'$P-value_{Benjamini-Hochberg} = 0.01$'),
            'export': 'pbs_mxl_chb_ceu_benjamini_hochberg_corrected_pvalues_742kb',
        },
    }
    # Load the hg19 gene information.
    hg19_gene_df = pd.read_csv(f'../annotations/hg19_genes/ncbi_refseq_genes_chr12.csv.gz')
    # Intialize the gene dictionary and labels.
    gene_dicc = {
        'SLC2A13': {'lab': r'$SLC2A13$'}, 'LRRK2': {'lab': r'$LRRK2$'},
        'MUC19': {'lab': r'$MUC19$'},
    }
    region = 'Longest MXL Introgressed Tract (Chr12: 40272001 - 41014000)'
    # For every gene.
    for gene in gene_dicc:
        # Fill the dictionary.
        gene_dicc[gene]['start'] = hg19_gene_df[hg19_gene_df['GENE_ID'] == gene].START.values[0]
        gene_dicc[gene]['stop'] = hg19_gene_df[hg19_gene_df['GENE_ID'] == gene].STOP.values[0]
    # Intialize a dictionary of snp type partitions.
    snp_dicc = {
        'IS_HUM': {'l': 'Human-specific', 'c': 'gray'},
        'IS_HOM': {'l': 'Shared Hominin', 'c': '#000000'},
        'IS_DEN': {'l': 'Denisovan-specific', 'c': '#E69F00'},
        'IS_NEA': {'l': 'Neanderthal-specific', 'c': '#56B4E9'},
        'IS_SHR': {'l': 'Shared Archaic', 'c': '#CC79A7'},
    }
    # Intialize the legend.
    if plot_key == 'PR':
        legend_handles = [
            Line2D(
                [0], [0], linestyle='none', marker='o', markersize=7.5, markeredgewidth=0.1,
                color=snp_dicc[snp_type]['c'], markeredgecolor='white', label=f'{snp_dicc[snp_type]["l"]} (n = {(info_dicc[snp_type] & info_dicc["IS_OUT"]).sum()})'
            ) for snp_type in ['IS_DEN', 'IS_NEA', 'IS_SHR', 'IS_HUM', 'IS_HOM']
        ]
    else:
        legend_handles = [
            Line2D(
                [0], [0], linestyle='none', marker='o', markersize=7.5, markeredgewidth=0.1,
                color=snp_dicc[snp_type]['c'], markeredgecolor='white', label=f'{snp_dicc[snp_type]["l"]} (n = {info_dicc[snp_type].sum()})'
            ) for snp_type in ['IS_DEN', 'IS_NEA', 'IS_SHR', 'IS_HUM', 'IS_HOM']
        ]
    legend_handles.append(plot_dicc[plot_key]['thresh'])
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the figure and axes.
    fig = plt.figure(figsize=(7, 3.5), dpi=300, facecolor='white', constrained_layout=True)
    ax = fig.add_subplot(111)
    # Plot the 72kb haplotype region.
    ax.axvspan(40759001, 40831000, alpha=0.05, facecolor='black')
    # For all snp types.
    for snp_type in snp_dicc:
        # Plot the results.
        if plot_key == 'PR':
            ax.scatter(
                info_dicc['POS'][info_dicc[snp_type] & info_dicc['IS_OUT']],
                plot_dicc[plot_key]['y_val'][info_dicc[snp_type] & info_dicc['IS_OUT']], zorder=5,
                color=snp_dicc[snp_type]['c'], marker='o', edgecolor='white', linewidth=0.1, alpha=0.75, s=20,
            )
        else:
            ax.scatter(
                info_dicc['POS'][info_dicc[snp_type]], plot_dicc[plot_key]['y_val'][info_dicc[snp_type]], zorder=5,
                color=snp_dicc[snp_type]['c'], marker='o', edgecolor='white', linewidth=0.1, alpha=0.75, s=20,
            )
    # Plot the threshold.
    ax.axhline(thresh, 0, 1, color='black', linestyle='dashed', lw=1)
    # For every gene.
    for gene in gene_dicc:
        # Plot the gene segments within the region.
        ax.plot(
            [max(gene_dicc[gene]['start'], info_dicc['POS'][0]), min(gene_dicc[gene]['stop'], info_dicc['POS'][-1])],
            [plot_dicc[plot_key]['gene_loc'], plot_dicc[plot_key]['gene_loc']], color='black', marker='|', ms=10, lw=1,
        )
        # Annotate the gene.
        ax.text(
            ((max(gene_dicc[gene]['start'], info_dicc['POS'][0]) + min(gene_dicc[gene]['stop'], info_dicc['POS'][-1])) / 2),
            plot_dicc[plot_key]['gene_loc'], gene_dicc[gene]['lab'], fontsize=10,
            horizontalalignment='center',verticalalignment='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black'),
        )
    # Rescale x-axis to Mb.
    x_ticks = np.arange(40_250_000, 41_100_000, 100_000)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{round(x_tick / 1e6, 3)} Mb' for x_tick in x_ticks])
    # Set the y-axis label.
    ax.set_ylim(plot_dicc[plot_key]['y_lim'])
    # Get the current y-ticks.
    yticks = ax.get_yticks()
    # Filter y-ticks.
    if plot_key == 'PR':
        filtered_yticks = [99.95, 99.96, 99.97, 99.98, 99.99, 100]
    else:
        filtered_yticks = [tick for tick in yticks if tick >= 0]
    # Set the y-ticks to the filtered y-ticks.
    ax.set_yticks(filtered_yticks)
    # Set the y-axis label
    ax.set_xlabel(region)
    ax.set_ylabel(plot_dicc[plot_key]['y_lab'])
    # Add a legend.
    ax.legend(
        handles=legend_handles, loc='upper center', bbox_to_anchor=plot_dicc[plot_key]['legend_loc'],
        ncol=3, frameon=True, fancybox=True, shadow=True, fontsize=8,
    )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/{plot_dicc[plot_key]["export"]}.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/{plot_dicc[plot_key]["export"]}.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/{plot_dicc[plot_key]["export"]}.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot!
    plt.show()
    return

# Define a function to load the simulated per snp results.
def load_mxl_slimulated_per_snp_results():
    # Intialize a list of models.
    model_list = ['neutral', 'negative']
    # Intialize a list of snp type partitions.
    snp_types = ['all', 'arc']
    # Intialize a dictionary to store the slimulated results.
    sim_dicc = {
        model: {'freqs': {}, 'pbs': {}} for model in model_list
    }
    # For every simulated model.
    for model in model_list:
        # For every snp set partition.
        for snp_type in snp_types:
            # Fill the dictionaries.
            sim_dicc[model]['freqs'][snp_type] = np.loadtxt(
                f'../muc19_results/mxl_slimulations/{model}_{snp_type}_snp_freqs_742kb_per_snp.txt.gz'
            )
            sim_dicc[model]['pbs'][snp_type] = np.loadtxt(
                f'../muc19_results/mxl_slimulations/{model}_pbs_{snp_type}_snps_742kb_per_snp.txt.gz'
            )
    return sim_dicc

# Define a function to summarize the per snp frequencies >= 30%
def summarize_snps_geq_thirty_perct(sim_dicc):
    # Intialize a list of models.
    model_list = ['neutral', 'negative']
    # Intialize a list of snp type partitions.
    snp_types = ['all', 'arc']
    # Intialize a label dictionary.
    lab_dicc = {
        'all': 'All SNPs', 'arc': 'Archaic SNPs',
        'neutral': 'Neutral',
        'negative': 'Negative',
    }
    # Intialize a dictionary to store the results.
    df_dicc = {
        'model': [], 'snp_type': [], 'tot': [], 'geq': [], 'pval': [],
    }
    # For every simulated model.
    for model in model_list:
        # For every snp set partition.
        for snp_type in snp_types:
            # Compute the total number of snps and the number of snps with an aaf >= 30%.
            n_tot = (~np.isnan(sim_dicc[model]['freqs']['all'])).sum()
            n_geq = (sim_dicc[model]['freqs'][snp_type] >= 0.3).sum()
            # Compute the p-value while adjusting for the number of permutations.
            pval = n_geq / n_tot
            if pval == 0:
                pval = '< {:.5e}'.format((1 / (n_tot + 1)))
            # Update the dictionary.
            df_dicc['model'].append(lab_dicc[model])
            df_dicc['snp_type'].append(lab_dicc[snp_type])
            df_dicc['tot'].append(int(n_tot))
            df_dicc['geq'].append(int(n_geq))
            df_dicc['pval'].append(pval)
    # Convert to a dataframe.
    geq_df = pd.DataFrame(df_dicc)
    # Rename the columns to look pretty.
    geq_df.rename(
        columns={
            'model': 'Selection Model',
            'snp_type': 'SNP Set',
            'tot': 'Total SNPs',
            'geq': r'MXL SNPs $\geq 30\%$',
            'pval': r'$P-value$',
        }, inplace=True
    )
    # Export.
    geq_df.to_csv('./dataframes/slim_per_snp_geq_freq_summary.csv.gz', index=False)
    return geq_df

# Define a function to summarize the p-values for the 417 observed outlier pbs snps.
def summarize_per_snp_pbs_pvals(info_dicc, sim_dicc):
    # Intialize a list of models.
    model_list = ['neutral', 'negative']
    # Intialize a list of snp type partitions.
    snp_types = ['all', 'arc']
    # Intialize a label dictionary.
    lab_dicc = {
        'all': 'All SNPs', 'arc': 'Archaic SNPs',
        'neutral': 'Neutral',
        'negative': 'Negative',
    }
    # Intialize a dictionary with the observed pbs values.
    obs_dicc = {
        'all': info_dicc['PBS'][info_dicc['IS_OUT']],
        'arc': info_dicc['PBS'][info_dicc['IS_ARC'] & info_dicc['IS_OUT']],
    }
    # Intialize a dictionary to store the p-values.
    df_dicc = {
        'model': [], 'snp_type': [], 'pbs': [], 'numer': [], 'denom': [],
        'unc_pval': [], 'unc_thresh': [], 'is_unc': [],
        'bon_pval': [], 'bon_thresh': [], 'is_bon': [],
        'bhp_pval': [], 'bhp_thresh': [], 'is_bhp': [],
    }
    # For every model.
    for model in model_list:
        # For every snp type.
        for snp_type in snp_types:
            # Determine the number of pbs values.
            n_tests = obs_dicc[snp_type].size
            # Intialize a list to store the numerator of p-values.
            numer = []
            # For every observered pbs value.
            for pbs in obs_dicc[snp_type]:
                # Update the list with the numerator.
                numer.append((sim_dicc[model]['pbs'][snp_type] >= pbs).sum())
            numer = np.array(numer)
            # Compute the uncorrected p-value.
            denom = (~np.isnan(sim_dicc[model]['pbs']['all'])).sum()
            unc_pvals = numer / denom
            unc_mask = unc_pvals < 0.05
            # Intialize the Bonferroni correction info.
            bon_thresh = 0.05 / n_tests
            bon_pvals = np.minimum((unc_pvals * n_tests), 1)
            bon_mask = bon_pvals < 0.05
            # Compute adjusted p-values using the Benjamini-Hochberg procedure.
            _, bhp_pvals, _, _ = multipletests(unc_pvals, alpha=0.01, method='fdr_bh')
            bhp_mask = bhp_pvals < 0.01
            # Update the dictionary.
            df_dicc['model'].extend(np.full(n_tests, lab_dicc[model]))
            df_dicc['snp_type'].extend(np.full(n_tests, lab_dicc[snp_type]))
            df_dicc['pbs'].extend(obs_dicc[snp_type])
            df_dicc['numer'].extend(numer)
            df_dicc['denom'].extend(np.full(n_tests, denom))
            df_dicc['unc_pval'].extend(unc_pvals)
            df_dicc['unc_thresh'].extend(np.full(n_tests, 0.05))
            df_dicc['is_unc'].extend(unc_mask)
            df_dicc['bon_pval'].extend(bon_pvals)
            df_dicc['bon_thresh'].extend(np.full(n_tests, bon_thresh))
            df_dicc['is_bon'].extend(bon_mask)
            df_dicc['bhp_pval'].extend(bhp_pvals)
            df_dicc['bhp_thresh'].extend(np.full(n_tests, 0.01))
            df_dicc['is_bhp'].extend(bhp_mask)
    # Convert the dictionary to a dataframe.
    pval_df = pd.DataFrame(df_dicc)
    # Rename the columns to look pretty.
    pval_df.rename(
        columns={
            'model': 'Selection Model',
            'snp_type': 'SNP Set',
            'pbs': r'$PBS_{OBS}$',
            'numer': r'$PBS_{OBS} \geq PBS_{SIM}$',
            'denom': 'Total Number of Simulated SNPs',
            'unc_pval': r'Uncorrected $P-value$',
            'unc_thresh': r'Uncorrected Significance Level',
            'is_unc': r'Uncorrected $P-value < 0.05$',
            'bon_pval': r'Bonferroni Adj. $P-value$',
            'bon_thresh': r'Bonferroni-Corrected Significance Level',
            'is_bon': r'Bonferroni Adj. $P-value < 0.05$',
            'bhp_pval': r'Benjamini-Hochberg Adj. $P-value$',
            'bhp_thresh': 'False Discovery Rate',
            'is_bhp': r'Benjamini-Hochberg Adj. $P-value < 0.01$',
        }, inplace=True
    )
    # Intialize the summary dataframe.
    summary_df = pval_df.groupby(['Selection Model', 'SNP Set']).agg({
        r'$PBS_{OBS}$': ['count'],
        r'Uncorrected $P-value < 0.05$': ['sum'],
        r'Bonferroni Adj. $P-value < 0.05$': ['sum'],
        r'Benjamini-Hochberg Adj. $P-value < 0.01$': ['sum']
    }).reset_index()
    # Reset the columns.
    summary_df.columns = [
        'Selection Model', 'SNP Set', r'Number of $PBS_{MXL:CHB:CEU}$ Values',
        r'Uncorrected $P-value < 0.05$', r'Bonferroni Adj. $P-value < 0.05$', r'Benjamini-Hochberg Adj. $P-value < 0.01$',
    ]
    # Export the results.
    summary_df.to_csv('./dataframes/slim_per_snp_pbs_significance_summary.csv.gz', index=False)
    pval_df.to_csv('./dataframes/slim_per_snp_pbs_pvalues.csv.gz', index=False)
    return summary_df, pval_df

# Define a function to summarize the p-values for the 417 observed outlier pbs snps.
def summarize_per_snp_unique_pbs_pvals(info_dicc, sim_dicc):
    # Intialize a list of models.
    model_list = ['neutral', 'negative']
    # Intialize a list of snp type partitions.
    snp_types = ['all', 'arc']
    # Intialize a label dictionary.
    lab_dicc = {
        'all': 'All SNPs', 'arc': 'Archaic SNPs',
        'neutral': 'Neutral',
        'negative': 'Negative',
    }
    # Intialize a dictionary with the observed pbs values.
    obs_dicc = {
        'all': info_dicc['PBS'][info_dicc['IS_OUT']],
        'arc': info_dicc['PBS'][info_dicc['IS_ARC'] & info_dicc['IS_OUT']],
    }
    # Intialize a dictionary with the observed pbs values.
    obs_dicc = {
        'all': info_dicc['PBS'][info_dicc['IS_OUT']],
        'arc': info_dicc['PBS'][info_dicc['IS_ARC'] & info_dicc['IS_OUT']],
    }
    # Intialize a dictionary to store the p-values.
    df_dicc = {
        'model': [], 'snp_type': [], 'pbs': [], 'numer': [], 'denom': [],
        'unc_pval': [], 'unc_thresh': [], 'is_unc': [],
        'bon_pval': [], 'bon_thresh': [], 'is_bon': [],
        'bhp_pval': [], 'bhp_thresh': [], 'is_bhp': [],
    }
    # For every model.
    for model in model_list:
        # For every snp type.
        for snp_type in snp_types:
            # Extract the number of unique pbs outliers.
            u_pbs = np.unique(obs_dicc[snp_type])
            # Determine the number of pbs values.
            n_tests = u_pbs.size
            # Intialize a list to store the numerator of p-values.
            numer = []
            # For every observered pbs value.
            for pbs in u_pbs:
                # Update the list with the numerator.
                numer.append((sim_dicc[model]['pbs'][snp_type] >= pbs).sum())
            numer = np.array(numer)
            # Compute the uncorrected p-value.
            denom = (~np.isnan(sim_dicc[model]['pbs']['all'])).sum()
            unc_pvals = numer / denom
            unc_mask = unc_pvals < 0.05
            # Intialize the Bonferroni correction info.
            bon_thresh = 0.05 / n_tests
            bon_pvals = np.minimum((unc_pvals * n_tests), 1)
            bon_mask = bon_pvals < 0.05
             # Compute adjusted p-values using the Benjamini-Hochberg procedure.
            _, bhp_pvals, _, _ = multipletests(unc_pvals, alpha=0.01, method='fdr_bh')
            bhp_mask = bhp_pvals < 0.01
            # Update the dictionary.
            df_dicc['model'].extend(np.full(n_tests, lab_dicc[model]))
            df_dicc['snp_type'].extend(np.full(n_tests, lab_dicc[snp_type]))
            df_dicc['pbs'].extend(u_pbs)
            df_dicc['numer'].extend(numer)
            df_dicc['denom'].extend(np.full(n_tests, denom))
            df_dicc['unc_pval'].extend(unc_pvals)
            df_dicc['unc_thresh'].extend(np.full(n_tests, 0.05))
            df_dicc['is_unc'].extend(unc_mask)
            df_dicc['bon_pval'].extend(bon_pvals)
            df_dicc['bon_thresh'].extend(np.full(n_tests, bon_thresh))
            df_dicc['is_bon'].extend(bon_mask)
            df_dicc['bhp_pval'].extend(bhp_pvals)
            df_dicc['bhp_thresh'].extend(np.full(n_tests, 0.01))
            df_dicc['is_bhp'].extend(bhp_mask)
    # Convert the dictionary to a dataframe.
    pval_df = pd.DataFrame(df_dicc)
    # Rename the columns to look pretty.
    pval_df.rename(
        columns={
            'model': 'Selection Model',
            'snp_type': 'SNP Set',
            'pbs': r'$PBS_{OBS}$',
            'numer': r'$PBS_{OBS} \geq PBS_{SIM}$',
            'denom': 'Total Number of Simulated SNPs',
            'unc_pval': r'Uncorrected $P-value$',
            'unc_thresh': r'Uncorrected Significance Level',
            'is_unc': r'Uncorrected $P-value < 0.05$',
            'bon_pval': r'Bonferroni Adj. $P-value$',
            'bon_thresh': r'Bonferroni-Corrected Significance Level',
            'is_bon': r'Bonferroni Adj. $P-value < 0.05$',
            'bhp_pval': r'Benjamini-Hochberg Adj. $P-value$',
            'bhp_thresh': 'False Discovery Rate',
            'is_bhp': r'Benjamini-Hochberg Adj. $P-value < 0.01$',
        }, inplace=True
    )
    # Intialize the summary dataframe.
    summary_df = pval_df.groupby(['Selection Model', 'SNP Set']).agg({
        r'$PBS_{OBS}$': ['count'],
        r'Uncorrected $P-value < 0.05$': ['sum'],
        r'Bonferroni Adj. $P-value < 0.05$': ['sum'],
        r'Benjamini-Hochberg Adj. $P-value < 0.01$': ['sum']
    }).reset_index()
    # Reset the columns.
    summary_df.columns = [
        'Selection Model', 'SNP Set', r'Unique $PBS_{MXL:CHB:CEU}$ Values',
        r'Uncorrected $P-value < 0.05$', r'Bonferroni Adj. $P-value < 0.05$', r'Benjamini-Hochberg Adj. $P-value < 0.01$',
    ]
    # Export the results.
    summary_df.to_csv('./dataframes/slim_per_snp_unique_pbs_significance_summary.csv.gz', index=False)
    pval_df.to_csv('./dataframes/slim_per_snp_unique_pbs_pvalues.csv.gz', index=False)
    return summary_df, pval_df

# Define a function to plot the per-snp distribution of mxl aaf.
def plot_snps_mxl_aaf(sim_dicc):
    # Intialize a list of models.
    model_list = ['neutral', 'negative']
    # Intialize a list of snp type partitions.
    snp_types = ['all', 'arc']
    # Intialize a list of combinations.
    model_snp_combos = []
    # For every simulated model.
    for model in model_list:
        # For every snp partiton.
        for snp_type in snp_types:
            # Update the list.
            model_snp_combos.append((model, snp_type))
    # Intialize a label dictionary.
    lab_dicc = {
        'all': 'All', 'arc': 'Archaic',
        'neutral': 'Neutral',
        'negative': 'Negative',
    }
    # Intialize the y-axis labels.
    snp_labs = [lab_dicc[snp_type] for snp_type in snp_types]
    # Intialize the positions to plot.
    snp_pos = np.arange(0, len(snp_types)*2, 2)
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        2, 2, figsize=(8, 4), dpi=300, facecolor='white',
        sharex=True, sharey=True, constrained_layout=True,
    )
    # Flatten the axes.
    axes = axes.flatten()
    # For every model and snp partition.
    for i, (model, snp_type) in enumerate(model_snp_combos):
        # Plot the simulated distribution.
        axes[i].hist(
            sim_dicc[model]['freqs'][snp_type],bins=np.arange(0, 1.05, 0.05),
            weights= 1 / sim_dicc[model]['freqs'][snp_type].size * np.ones(sim_dicc[model]['freqs'][snp_type].size),
            histtype='step', color='gray', linewidth=2,
        )
        # Plot the threshold.
        axes[i].axvline(x=0.3, color='black', linestyle='--', linewidth=1.5)
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel('Density')
        axes[i].set_xlabel(f'MXL AAF - {lab_dicc[snp_type]} SNPs ({lab_dicc[model]})')
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
    # Export the plot.
    plt.savefig(
        './supp_figures/png/slim_per_snp_aaf.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/slim_per_snp_aaf.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/slim_per_snp_aaf.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to load the observed per snp results.
def init_mxl_obs_region_results(tgp_gt_742kb, tgp_pos_742kb, tgp_arcs_gt_742kb, tgp_arcs_pos_742kb):
    # Intialize a dictionary to store the observed results.
    obs_dicc = {
        'freqs': {'742kb': {}, '72kb': {}},
        'pbs': {'742kb': {}, '72kb': {}},
        'out': {'742kb': {}, '72kb': {}}
    }
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    pop_dicc = {
        'MXL': tgp_df[tgp_df['POP'] == 'MXL'].index.values,
        'CHB': tgp_df[tgp_df['POP'] == 'CHB'].index.values,
        'CEU': tgp_df[tgp_df['POP'] == 'CEU'].index.values,
    }
    # Classify each snp type.
    _, _, arc_hum_snps = tgp_classify_snps(tgp_arcs_gt_742kb, 'MXL', 0.01)
    # Intialize a snp type mask dictionary.
    mask_dicc = {
        'all': np.full(tgp_pos_742kb.size, True),
        'arc': np.isin(tgp_pos_742kb, tgp_arcs_pos_742kb[arc_hum_snps['ARC']]),
    }
    # Determine what positions are in the 72kb region.
    is_72kb = (40759001 <= tgp_pos_742kb) & (tgp_pos_742kb <= 40831000)
    # Compute the per site pbs values.
    mxl_pbs_per_site = calc_pbs_per_site(tgp_gt_742kb, pop_dicc['MXL'], pop_dicc['CHB'], pop_dicc['CEU'])
    # Compute alternative allele frequencies.
    mxl_aaf = calc_alt_freqs(tgp_gt_742kb.take(pop_dicc['MXL'], axis=1))
    # Determine the sites with an aaf >= 30%.
    is_geq = mxl_aaf >= 0.3
    # For every snp type.
    for snp_type in mask_dicc:
        # Fill the dictionary.
        obs_dicc['freqs']['742kb'][snp_type] = (is_geq & mask_dicc[snp_type]).sum()
        obs_dicc['freqs']['72kb'][snp_type] = (is_geq & mask_dicc[snp_type] & is_72kb).sum()
        obs_dicc['pbs']['742kb'][snp_type] = calc_pbs_per_region(
            tgp_gt_742kb.compress(mask_dicc[snp_type], axis=0), 
            pop_dicc['MXL'], pop_dicc['CHB'], pop_dicc['CEU']
        )
        obs_dicc['pbs']['72kb'][snp_type] = calc_pbs_per_region(
            tgp_gt_742kb.compress((mask_dicc[snp_type] & is_72kb), axis=0), 
            pop_dicc['MXL'], pop_dicc['CHB'], pop_dicc['CEU']
        )
    return obs_dicc

# Define a function to load the simulated per region.
def load_mxl_slimulated_region_results():
    # Intialize a list of models.
    model_list = ['neutral', 'negative']
    # Intialize the different data partitions.
    regions = ['742kb', '72kb']
    snp_types = ['all', 'arc']
    # Intialize a dictionary to store the slimulated results.
    sim_dicc = {
        model: {'freqs': {'742kb': {}, '72kb': {}}, 'pbs': {'742kb': {}, '72kb': {}}}
        for model in model_list
    }
    # For every simulated model.
    for model in model_list:
        # Load the results.
        sim_df = pd.read_csv(f'../muc19_results/mxl_slimulations/{model}_per_10k_replicates.csv.gz')
        # For every region.
        for region in regions:
            # For every snp set partition.
            for snp_type in snp_types:
                # Fill the dictionaries.
                sim_dicc[model]['freqs'][region][snp_type] = sim_df[f'geq_{snp_type}_snps_{region}'].values
                sim_dicc[model]['pbs'][region][snp_type] = sim_df[f'pbs_{snp_type}_snps_{region}'].values
    return sim_dicc

# Define a function to summarize the region frequencies >= 30%
def summarize_region_geq_thirty_perct(obs_dicc, sim_dicc):
    # Intialize a list of models.
    model_list = ['neutral', 'negative']
    # Intialize the different data partitions.
    regions = ['742kb', '72kb']
    snp_types = ['all', 'arc']
    # Intialize a label dictionary.
    lab_dicc = {
        'all': 'All SNPs', 'arc': 'Archaic SNPs',
        'neutral': 'Neutral',
        'negative': 'Negative',
    }
    # Intialize a dictionary to store the results.
    df_dicc = {
        'model': [], 'snp_type': [], 'region': [], 'obs': [], 'tot': [], 'geq': [], 'pval': [],
    }
    # For every simulated model.
    for model in model_list:
        # For every region.
        for region in regions:
            # For every snp set partition.
            for snp_type in snp_types:
                # Compute the total number of replicates and p-values.
                n_tot = (~np.isnan(sim_dicc[model]['freqs'][region][snp_type])).sum()
                n_geq = (sim_dicc[model]['freqs'][region][snp_type] >= obs_dicc['freqs'][region][snp_type]).sum()
                # Compute the p-value while adjusting for the number of permutations.
                pval = n_geq / n_tot
                if pval == 0:
                    pval = '<0.0001'
                # Update the dictionary.
                df_dicc['model'].append(lab_dicc[model])
                df_dicc['region'].append(region)
                df_dicc['obs'].append(int(obs_dicc['freqs'][region][snp_type]))
                df_dicc['snp_type'].append(lab_dicc[snp_type])
                df_dicc['tot'].append(int(n_tot))
                df_dicc['geq'].append(int(n_geq))
                df_dicc['pval'].append(pval)
    # Convert to a dataframe.
    geq_df = pd.DataFrame(df_dicc)
    # Rename the columns to look pretty.
    geq_df.rename(
        columns={
            'model': 'Selection Model',
            'snp_type': 'SNP Set',
            'region': 'Region',
            'obs': r'MXL AAF $\geq 30\%$ (Observed)',
            'tot': 'Total Replicates',
            'geq': r'Simulated Replicates $\geq$ Observed',
            'pval': r'$P-value$',
        }, inplace=True
    )
    # Export the dataframe.
    geq_df.to_csv('./dataframes/slim_per_region_geq_freq_summary.csv.gz', index=False)
    return geq_df

# Define a function to summarize the region pbs values.
def summarize_region_pbs(obs_dicc, sim_dicc):
    # Intialize a list of models.
    model_list = ['neutral', 'negative']
    # Intialize the different data partitions.
    regions = ['742kb', '72kb']
    snp_types = ['all', 'arc']
    # Intialize a label dictionary.
    lab_dicc = {
        'all': 'All SNPs', 'arc': 'Archaic SNPs',
        'neutral': 'Neutral',
        'negative': 'Negative',
    }
    # Intialize a dictionary to store the results.
    df_dicc = {
        'model': [], 'snp_type': [], 'region': [], 'obs': [], 'tot': [], 'geq': [], 'pval': [],
    }
    # For every simulated model.
    for model in model_list:
        # For every region.
        for region in regions:
            # For every snp set partition.
            for snp_type in snp_types:
                # Compute the total number of replicates and p-values.
                n_tot = (~np.isnan(sim_dicc[model]['pbs'][region][snp_type])).sum()
                n_geq = (sim_dicc[model]['pbs'][region][snp_type] >= obs_dicc['pbs'][region][snp_type]).sum()
                # Compute the p-value while adjusting for the number of permutations.
                pval = n_geq / n_tot
                if pval == 0:
                    pval = '<0.0001'
                # Update the dictionary.
                df_dicc['model'].append(lab_dicc[model])
                df_dicc['region'].append(region)
                df_dicc['obs'].append(obs_dicc['pbs'][region][snp_type])
                df_dicc['snp_type'].append(lab_dicc[snp_type])
                df_dicc['tot'].append(int(n_tot))
                df_dicc['geq'].append(int(n_geq))
                df_dicc['pval'].append(pval)
    # Convert to a dataframe.
    pbs_df = pd.DataFrame(df_dicc)
    # Rename the columns to look pretty.
    pbs_df.rename(
        columns={
            'model': 'Selection Model',
            'snp_type': 'SNP Set',
            'region': 'Region',
            'obs': r'$PBS$ (Observed)',
            'tot': 'Total Replicates',
            'geq': r'Simulated Replicates $\geq$ Observed',
            'pval': r'$P-value$',
        }, inplace=True
    )
    # Export the dataframe.
    pbs_df.to_csv('./dataframes/slim_per_region_pbs_summary.csv.gz', index=False)
    return pbs_df

# Define a function to plot the per-snp distribution of mxl aaf.
def plot_region_mxl_aaf(obs_dicc, sim_dicc, region):
    # Intialize a list of models.
    model_list = ['neutral', 'negative']
    # Intialize a list of snp type partitions.
    snp_types = ['arc']
    # Intialize a list of combinations.
    model_snp_combos = []
    # For every simulated model.
    for model in model_list:
        # For every snp partiton.
        for snp_type in snp_types:
            # Update the list.
            model_snp_combos.append((model, snp_type))
    # Intialize a label dictionary.
    lab_dicc = {
        'all': 'All', 'arc': 'Archaic',
        'neutral': 'Neutral',
        'negative': 'Negative',
    }
    # Determine the significance threshold.
    sig_thresh = 0.05 
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 2, figsize=(7, 3), dpi=300, facecolor='white',
        sharex='row', sharey=True, constrained_layout=True,
    )
    # Flatten the axes.
    axes = axes.flatten()
    # For every model and snp partition.
    for i, (model, snp_type) in enumerate(model_snp_combos):
        # Plot the simulated distribution.
        if region == '72kb':
            axes[i].hist(
                sim_dicc[model]['freqs'][region][snp_type], bins=np.arange(0, 25, 2.5),
                histtype='stepfilled', color='gray',
            )
        else:
            axes[i].hist(
                sim_dicc[model]['freqs'][region][snp_type], bins=np.arange(0, 160, 10),
                histtype='stepfilled', color='gray',
            )
        # Compute the total number of replicates and p-values.
        n_tot = (~np.isnan(sim_dicc[model]['freqs'][region][snp_type])).sum()
        n_geq = (sim_dicc[model]['freqs'][region][snp_type] >= obs_dicc['freqs'][region][snp_type]).sum()
        # Compute the p-value.
        pval = n_geq / n_tot
        # If the p-value is significant.
        if pval < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (obs_dicc['freqs'][region][snp_type], '*', 25, 0.49)
        # Else, it is not significant.
        else:
            # Update the dictionary.
            sig_dicc[i] = (obs_dicc['freqs'][region][snp_type], r'$ns$', 10, 0.5225)
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Set the x-axis limits.
        if region == '72kb':
            axes[i].set_xlim(right=150)
        else:
            axes[i].set_xlim(right=220)
        # Add axes labels.
        axes[i].set_ylabel(f'Frequency (per {region} Region)')
        axes[i].set_xlabel(fr'Number of {lab_dicc[snp_type]} SNPs $\geq$ 30% ({lab_dicc[model]})')
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/slim_region_aaf_geq_30_percent_{region}.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/slim_region_aaf_geq_30_percent_{region}.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/slim_region_aaf_geq_30_percent_{region}.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to plot the per-snp distribution of mxl aaf.
def plot_region_mxl_pbs(obs_dicc, sim_dicc, region):
    # Intialize a list of models.
    model_list = ['neutral', 'negative']
    # Intialize a list of snp type partitions.
    snp_types = ['all', 'arc']
    # Intialize a list of combinations.
    model_snp_combos = []
    # For every simulated model.
    for model in model_list:
        # For every snp partiton.
        for snp_type in snp_types:
            # Update the list.
            model_snp_combos.append((model, snp_type))
    # Intialize a label dictionary.
    lab_dicc = {
        'all': 'All', 'arc': 'Archaic',
        'neutral': 'Neutral',
        'negative': 'Negative',
    }
    # Determine the significance threshold.
    sig_thresh = 0.05 
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        2, 2, figsize=(9, 4.5), dpi=300, facecolor='white',
        sharex='row', sharey=True, constrained_layout=True,
    )
    # Flatten the axes.
    axes = axes.flatten()
    # For every model and snp partition.
    for i, (model, snp_type) in enumerate(model_snp_combos):
        # Plot the simulated distribution.
        axes[i].hist(
            sim_dicc[model]['pbs'][region][snp_type], bins=np.arange(0, 0.33, 0.01),
            histtype='stepfilled', color='gray',
        )
        # Compute the total number of replicates and p-values.
        n_tot = (~np.isnan(sim_dicc[model]['pbs'][region][snp_type])).sum()
        n_geq = (sim_dicc[model]['pbs'][region][snp_type] >= obs_dicc['pbs'][region][snp_type]).sum()
        # Compute the p-value.
        pval = n_geq / n_tot
        # If the p-value is significant.
        if pval < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (obs_dicc['pbs'][region][snp_type], '*', 25, 0.49)
        # Else, it is not significant.
        else:
            # Update the dictionary.
            sig_dicc[i] = (obs_dicc['pbs'][region][snp_type], r'$ns$', 10, 0.5225)
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel(f'Frequency (per {region} Region)')
        axes[i].set_xlabel(fr'$PBS$ for {lab_dicc[snp_type]} SNPs ({lab_dicc[model]})')
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/slim_region_pbs_{region}.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/slim_region_pbs_{region}.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/slim_region_pbs_{region}.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to compute the power.
def compute_pbs_power():
    # Load the results for no selection.
    neu_df = pd.read_csv(f'../muc19_results/mxl_slimulations/neutral_per_10k_replicates.csv.gz')
    # Convert the dataframe to a dictionary.
    neu_info = {
        key: np.array(value) for key, value in neu_df.to_dict(orient='list').items()
    }
    # Find the the upper 5% quantile.
    q5_all_72kb = np.quantile(neu_info['pbs_all_snps_72kb'], 0.95)
    q5_arc_72kb = np.quantile(neu_info['pbs_arc_snps_72kb'], 0.95)
    q5_all_742kb = np.quantile(neu_info['pbs_all_snps_742kb'], 0.95)
    q5_arc_742kb = np.quantile(neu_info['pbs_arc_snps_742kb'], 0.95)
    # Load the positive selection results.
    s1_df = pd.read_csv(f'../muc19_results/mxl_slimulations/positive_s1_per_1k_replicates.csv.gz')
    s01_df = pd.read_csv(f'../muc19_results/mxl_slimulations/positive_s01_per_1k_replicates.csv.gz')
    s0015_df = pd.read_csv(f'../muc19_results/mxl_slimulations/positive_s0015_per_1k_replicates.csv.gz')
    # Intialize a dictionary to store the results.
    power_dicc = {
        's': [0.1, 0.01, 0.0015],
        'prop_all_snps_72kb': [], 'prop_arc_snps_72kb': [],
        'prop_all_snps_742kb': [], 'prop_arc_snps_742kb': [],
    }
    # For every selection coeficient.
    for sel_df in [s1_df, s01_df, s0015_df]:
        # Determine how many snps fall outside the neutral distribution.
        is_gt_q5_all_72kb = sel_df['pbs_all_snps_72kb'].values > q5_all_72kb
        is_gt_q5_arc_72kb = sel_df['pbs_arc_snps_72kb'].values > q5_arc_72kb
        is_gt_q5_all_742kb = sel_df['pbs_all_snps_742kb'].values > q5_all_742kb
        is_gt_q5_arc_742kb = sel_df['pbs_arc_snps_742kb'].values > q5_arc_742kb
        # Fill the dictionary.
        power_dicc['prop_all_snps_72kb'].append(is_gt_q5_all_72kb.sum() / sel_df.shape[0])
        power_dicc['prop_arc_snps_72kb'].append(is_gt_q5_arc_72kb.sum() / sel_df.shape[0])
        power_dicc['prop_all_snps_742kb'].append(is_gt_q5_all_742kb.sum() / sel_df.shape[0])
        power_dicc['prop_arc_snps_742kb'].append(is_gt_q5_arc_742kb.sum() / sel_df.shape[0])
    # Convert to a dataframe.
    power_df = pd.DataFrame(power_dicc)
    # Rename the columns to look pretty.
    power_df.rename(
        columns={
            's': r'$s$',
            'prop_all_snps_72kb': 'All SNPs (72kb)',
            'prop_arc_snps_72kb': 'Archaic SNPs (72kb)',
            'prop_all_snps_742kb': 'All SNPs (742kb)',
            'prop_arc_snps_742kb': 'Archaic SNPs (742kb)',
        }, inplace=True
    )
    # Export the dataframe.
    power_df.to_csv('./dataframes/slim_pbs_power_summary.csv.gz', index=False)
    return power_df



###########################
### SEQUENCE DIVERGENCE ###
###########################

# Define a function to show the results for the focal mxl individual with the 742kb introgressed tract.
def show_focal_mxl_seq_div_742kb():
    # Load the sequence divergence results.
    tgp_hap_div_742kb_df = pd.read_csv('./dataframes/tgp_haplotype_archaic_diplotype_divergence_742kb.csv.gz')
    # Subset the focal individual.
    mxl_hap_div_742kb_df = tgp_hap_div_742kb_df[tgp_hap_div_742kb_df['Individual'] == 'NA19725']
    # Export the the dataframe.
    mxl_hap_div_742kb_df.to_csv('./dataframes/focal_mxl_sequence_divergence_742kb.csv.gz', index=False)
    return mxl_hap_div_742kb_df

# Define a function to show the results for the focal mxl individual at the 72kb region.
def show_focal_mxl_seq_div_72kb():
    # Load the sequence divergence results.
    tgp_hap_div_72kb_df = pd.read_csv('./dataframes/tgp_haplotype_archaic_diplotype_divergence_72kb.csv.gz')
    # Subset the focal individual.
    mxl_hap_div_72kb_df = tgp_hap_div_72kb_df[tgp_hap_div_72kb_df['Individual'] == 'NA19664']
    # Export the the dataframe.
    mxl_hap_div_72kb_df.to_csv('./dataframes/focal_mxl_sequence_divergence_72kb.csv.gz', index=False)
    return mxl_hap_div_72kb_df

# Define a function to show the results for the focal mxl individual at the phased 72kb region.
def show_focal_mxl_hap_div_72kb():
    # Load the sequence divergence results.
    tgp_hap_div_72kb_df = pd.read_csv('./dataframes/tgp_haplotype_late_neanderthal_phased_haplotype_divergence_72kb.csv.gz')
    # Subset the focal individual.
    mxl_hap_div_72kb_df = tgp_hap_div_72kb_df[tgp_hap_div_72kb_df['Individual'] == 'NA19664']
    # Export the the dataframe.
    mxl_hap_div_72kb_df.to_csv('./dataframes/focal_mxl_haplotype_divergence_72kb.csv.gz', index=False)
    return mxl_hap_div_72kb_df

# Define a function to load the windowed sequence divergence results.
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
        esl_df = load_windows(f'tgp_{arc.lower()}_masked_no_aa', 'variant', window_size)
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

# Define a function to plot the focal mxl sequence divergence results.
def plot_focal_mxl_seq_div(window_size):
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary for the focal individuals.
    mxl_dicc = {742: 'NA19725', 72: 'NA19664'}
    # Determine the focal individual's index.
    mxl_idx = np.where(tgp_df['IND'].values == mxl_dicc[window_size])[0][0]
    # Load the sequence divergence results.
    div_df = pd.read_csv(f'./dataframes/tgp_haplotype_archaic_diplotype_divergence_{window_size}kb.csv.gz')
    # Subset the focal individual.
    mxl_div_df = div_df[div_df['Individual'] == mxl_dicc[window_size]]
    # Extract the archaic ids, sequence divergence and p-values.
    # Note for NA19725 the longest haplotype is the second haplotype and NA19664 carries two haplotypes and no heterozygous sites.
    arc_ids = mxl_div_df['Archaic'].values
    obs_divs = mxl_div_df[f'Focal {window_size}kb Region (Seq. Div. Hap. 2)'].values
    pvals = mxl_div_df[r'$P-value$ (Hap. 2)'].values
    # Load the windowed results.
    wind_dicc = load_hum_hap_v_arc_dip_pwd_windows(window_size)
    # Determine the significance threshold.
    sig_thresh = 0.05 / 2
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        2, 2, figsize=(8, 8), dpi=300,
        facecolor='white',
        sharex=True, sharey=True,
        constrained_layout=True,
    )
    # Flatten the axes.
    axes = axes.flatten()
    # For every archaic.
    for i, arc in enumerate(['DEN', 'ALT', 'CHA', 'VIN']):
        # Extract the windows of comparable effective sequence length.
        wind_idx = load_esl_qc_windows_idx(f'tgp_{arc.lower()}_masked_no_aa', window_size)
        # Extract the distribution of sequence divergence.
        dist = np.concatenate((
            wind_dicc['hap_1'][arc]['div'][:, mxl_idx][wind_idx],
            wind_dicc['hap_2'][arc]['div'][:, mxl_idx][wind_idx],
        ))
        # If the p-value is significant.
        if pvals[i] < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (obs_divs[i], '*', 25, 0.49)
        # Else, it is not significant.
        else:
            # Update the dictionary.
            sig_dicc[i] = (obs_divs[i], r'$ns$', 10, 0.5225)
        # Plot the distribution.
        axes[i].hist(
            dist, bins=np.linspace(0, 0.005, 100),
            histtype='stepfilled', color='gray',
        )
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlabel(f'{arc_ids[i]} v MXL ({mxl_dicc[window_size]})'+'\n'+f'Seq. Div. per {window_size}kb Window')
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/focal_mxl_v_archaics_seq_div_{window_size}kb.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/focal_mxl_v_archaics_seq_div_{window_size}kb.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/focal_mxl_v_archaics_seq_div_{window_size}kb.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot!
    plt.show()
    return

# Define a function to load the windowed psuedo-haplotype divergence results.
def load_hum_hap_v_cha_vin_phap_pwd_windows():
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
        esl_df = load_windows(f'tgp_{arc.lower()}_masked_no_aa', 'variant', 72)
        # For every haplotype.
        for hap in ['hap_1', 'hap_2']:
            # For all chromosomes.
            for chrom in range(1, 23):
                # Load the results.
                pw_mat = np.loadtxt(
                    f'../muc19_results/tgp_{arc.lower()}_masked_no_aa/{arc.lower()}_{hap}_psuedo_hap_diffs_chr{chrom}_72kb.txt.gz',
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

# Define a function to plot the sequence divergence results for the phased archaics and focal mxl individual.
def plot_phased_archaics_mxl_seq_div():
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Determine the focal individual's index who is homozygous for the denisovan-like haplotype.
    mxl_idx = np.where(tgp_df['IND'].values == 'NA19664')[0][0]
    # Load the windowed results.
    wind_dicc = load_hum_hap_v_cha_vin_phap_pwd_windows()
    # Load the window indicies of comparable effective sequence length.
    cha_idx = load_esl_qc_windows_idx('tgp_cha_masked_no_aa', 72)
    vin_idx = load_esl_qc_windows_idx('tgp_vin_masked_no_aa', 72)
    # Grab the individual's windowed distributions.
    cha_dist = np.concatenate((
        wind_dicc['hap_1']['CHA']['div'][:, mxl_idx][cha_idx],
        wind_dicc['hap_2']['CHA']['div'][:, mxl_idx][cha_idx],
    ))
    vin_dist = np.concatenate((
        wind_dicc['hap_1']['VIN']['div'][:, mxl_idx][vin_idx],
        wind_dicc['hap_2']['VIN']['div'][:, mxl_idx][vin_idx],
    ))
    # Load the observed phased results.
    obs_df = pd.read_csv('./dataframes/tgp_haplotype_late_neanderthal_phased_haplotype_divergence_72kb.csv.gz')
    # Intialize masks.
    mxl_mask = np.isin(obs_df['Individual'].values, 'NA19664')
    cha_mask = np.isin(obs_df['Archaic Hap.'].values, 'Chagyrskaya Nean. Hap. 2')
    vin_mask = np.isin(obs_df['Archaic Hap.'].values, 'Vindija Nean. Hap. 2')
    # Grab the results (Note: this individual has no heterozygous sites so the haplotype is arbitrary).
    cha_div = obs_df[mxl_mask & cha_mask]['Focal 72kb Region (Seq. Div. Hap. 1)'].values[0]
    cha_pval = obs_df[mxl_mask & cha_mask][r'$P-value$ (Hap. 1)'].values[0]
    vin_div = obs_df[mxl_mask & vin_mask]['Focal 72kb Region (Seq. Div. Hap. 1)'].values[0]
    vin_pval = obs_df[mxl_mask & vin_mask][r'$P-value$ (Hap. 1)'].values[0]
    # Intialize the plotting dictionary.
    plot_dicc = {
        'Chagyrskaya Nean. v MXL (NA19664)'+'\n'+'Pseudo-Hap. Div. per 72kb Window': {
            'obs': cha_div, 'dist': cha_dist, 'pval': cha_pval,
        },
        'Vindija Nean. v MXL (NA19664)'+'\n'+'Pseudo-Hap. Div. per 72kb Window': {
            'obs': vin_div, 'dist': vin_dist, 'pval': vin_pval,
        },
    }
    # Determine the significance threshold.
    sig_thresh = 0.05 / 4 
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 2, figsize=(7, 3.5), dpi=300,
        facecolor='white',
        sharex=True, sharey=True,
        constrained_layout=True,
    )
    # For each subplot.
    for i, key in enumerate(plot_dicc):
        # If the p-value is significant.
        if plot_dicc[key]['pval'] < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (plot_dicc[key]['obs'], '*', 25, 0.49)
        # Else, it is not significant.
        else:
            # Update the dictionary.
            sig_dicc[i] = (plot_dicc[key]['obs'], r'$ns$', 10, 0.5225)
        # Plot the distribution.
        axes[i].hist(
            plot_dicc[key]['dist'], bins=np.linspace(0, 0.005, 100),
            histtype='stepfilled', color='gray',
        )
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlabel(key)
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        './supp_figures/png/mxl_denisovan_like_hap_v_late_neanderthal_phased_hap_div_72kb.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/mxl_denisovan_like_hap_v_late_neanderthal_phased_hap_div_72kb.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/mxl_denisovan_like_hap_v_late_neanderthal_phased_hap_div_72kb.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to annotate haplotype identities.
def annotate_hap_identities(den_thresh, rec_thresh, tgp_or_sgdp):
    # Intialize the meta information.
    tgp_or_sgdp_dicc = {
        'tgp': '../meta_data/tgp_mod.txt', 'sgdp': '../meta_data/sgdp.txt',
    }
    meta_data = pd.read_csv(
        tgp_or_sgdp_dicc[tgp_or_sgdp], sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the inidcies.
    all_inds = meta_data['IND'].values
    # Load the sequence divergence results.
    div_72kb = pd.read_csv(f'./dataframes/{tgp_or_sgdp}_haplotype_archaic_diplotype_divergence_72kb.csv.gz')
    div_742kb = pd.read_csv(f'./dataframes/{tgp_or_sgdp}_haplotype_archaic_diplotype_divergence_742kb.csv.gz')
    # Extract the denisovan divergence results.
    den_div_72kb = div_72kb[div_72kb['Archaic'] == 'Denisovan']
    den_hap1_72kb = den_div_72kb['Focal 72kb Region (Pairwise Diffs. Hap. 1)'].values
    den_hap2_72kb = den_div_72kb['Focal 72kb Region (Pairwise Diffs. Hap. 2)'].values
    # Intialize masks for the unphased haplotype identities at the 72kb region.
    # Note: this threshold was chosen from the S-curves.
    is_den_hap1_72kb = np.isin(all_inds, den_div_72kb[den_hap1_72kb < den_thresh]['Individual'].values)
    is_den_hap2_72kb = np.isin(all_inds, den_div_72kb[den_hap2_72kb < den_thresh]['Individual'].values)
    is_rec_hap1_72kb = np.isin(
        all_inds, den_div_72kb[(den_hap1_72kb > den_thresh) & (den_hap1_72kb < rec_thresh)]['Individual'].values,
    )
    is_rec_hap2_72kb = np.isin(
        all_inds, den_div_72kb[(den_hap2_72kb > den_thresh) & (den_hap2_72kb < rec_thresh)]['Individual'].values,
    )
    is_hum_hap1_72kb = ~(is_den_hap1_72kb | is_rec_hap1_72kb)
    is_hum_hap2_72kb = ~(is_den_hap2_72kb | is_rec_hap2_72kb)
    # Update the dataframe.
    meta_data['IS_HAP1_DEN_72KB'] = is_den_hap1_72kb
    meta_data['IS_HAP2_DEN_72KB'] = is_den_hap2_72kb
    meta_data['IS_HAP1_HUM_72KB'] = is_hum_hap1_72kb
    meta_data['IS_HAP2_HUM_72KB'] = is_hum_hap2_72kb
    meta_data['IS_HAP1_REC_72KB'] = is_rec_hap1_72kb
    meta_data['IS_HAP2_REC_72KB'] = is_rec_hap2_72kb
    meta_data['N_DEN_HAPS_72KB'] = is_den_hap1_72kb.astype(int) + is_den_hap2_72kb.astype(int)
    meta_data['N_HUM_HAPS_72KB'] = is_hum_hap1_72kb.astype(int) + is_hum_hap2_72kb.astype(int)
    meta_data['N_REC_HAPS_72KB'] = is_rec_hap1_72kb.astype(int) + is_rec_hap2_72kb.astype(int)
    # Intialize masks for the unphased haplotype identities at the 742kb region.
    # Note: here we correct for two multiple comparisons per individual (ie two modern haps x one archaic genotype).
    is_cha_742kb = div_742kb['Archaic'].values == 'Chagyrskaya Nean.'
    is_vin_742kb = div_742kb['Archaic'].values == 'Vindija Nean.'
    is_hap1_sig_742kb = div_742kb['$P-value$ (Hap. 1)'].values < 0.025
    is_hap2_sig_742kb = div_742kb['$P-value$ (Hap. 2)'].values < 0.025
    is_cha_hap1_742kb = np.isin(all_inds, div_742kb[is_cha_742kb & is_hap1_sig_742kb]['Individual'].values)
    is_vin_hap1_742kb = np.isin(all_inds, div_742kb[is_vin_742kb & is_hap1_sig_742kb]['Individual'].values)
    is_arc_hap1_742kb = is_cha_hap1_742kb | is_vin_hap1_742kb
    is_cha_hap2_742kb = np.isin(all_inds, div_742kb[is_cha_742kb & is_hap2_sig_742kb]['Individual'].values)
    is_vin_hap2_742kb = np.isin(all_inds, div_742kb[is_vin_742kb & is_hap2_sig_742kb]['Individual'].values)
    is_arc_hap2_742kb = is_cha_hap2_742kb | is_vin_hap2_742kb
    # Update the dataframe.
    meta_data['IS_HAP1_ARC_742KB'] = is_arc_hap1_742kb
    meta_data['IS_HAP2_ARC_742KB'] = is_arc_hap2_742kb
    meta_data['N_ARC_HAPS_742KB'] = is_arc_hap1_742kb.astype(int) + is_arc_hap2_742kb.astype(int)
    return meta_data, {key: np.array(value) for key, value in meta_data.to_dict(orient='list').items()}

# Define a function to summarize the introgressed haplotype frequency.
def compile_tgp_introgressed_hap_freqs(window_size):
    # Intialize the hap key.
    if window_size == 72:
        hap_key = 'N_DEN_HAPS_72KB'
    else:
        hap_key = 'N_ARC_HAPS_742KB'
    # Load the haplotype identities.
    _, hap_info = annotate_hap_identities(70, 100, 'tgp')
    # Load the meta data file for the TGP.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    tgp_info = {key: np.array(value) for key, value in tgp_df.to_dict(orient='list').items()}
    # Intialize an ordered dictionary.
    tgp_dicc = {
        'AMR': ['MXL', 'PEL', 'CLM', 'PUR'],
        'SAS': ['BEB', 'STU', 'ITU', 'PJL', 'GIH'],
        'EAS': ['CHB', 'KHV', 'CHS', 'JPT', 'CDX'],
        'EUR': ['TSI', 'CEU', 'IBS', 'GBR', 'FIN'],
        'AFR': ['LWK', 'GWD', 'MSL', 'ESN', 'YRI'],
    }
    # Intialize a dictionaries.
    pop_dicc = {
        'spop': [], 'pop': [], 'n_tot': [],
        'n_haps': [],  'freq': [],
    }
    spop_dicc = {
        'spop': [], 'n_tot': [],
        'n_haps': [],  'freq': [],
    }
    # For every super population.
    for spop in tgp_dicc:
        # Determine the total number of chromosomes.
        spop_tot = (tgp_info['SUPERPOP'] == spop).sum() * 2
        # Determine the number of haplotypes.
        spop_haps = (hap_info[hap_key][hap_info['SUPERPOP'] == spop]).sum()
        # Update the dictionary.
        spop_dicc['spop'].append(spop)
        spop_dicc['n_tot'].append(spop_tot)
        spop_dicc['n_haps'].append(spop_haps)
        spop_dicc['freq'].append(spop_haps / spop_tot)
        # For every population.
        for pop in tgp_dicc[spop]:
            # Determine the total number of chromosomes.
            pop_tot = (tgp_info['POP'] == pop).sum() * 2
            # Determine the number of haplotypes.
            pop_haps = (hap_info[hap_key][hap_info['POP'] == pop]).sum()
            # Update the dictionary.
            pop_dicc['spop'].append(spop)
            pop_dicc['pop'].append(pop)
            pop_dicc['n_tot'].append(pop_tot)
            pop_dicc['n_haps'].append(pop_haps)
            pop_dicc['freq'].append(pop_haps / pop_tot)
    # Convert to dataframes.
    pop_df = pd.DataFrame(pop_dicc)
    spop_df = pd.DataFrame(spop_dicc)
    # Cleanup the column names.
    spop_df.rename(
        columns={
            'spop': 'Super Population',
            'n_tot': 'Total Number of Chromosomes',
            'n_haps': 'Number of Introgressed Haps.',
            'freq': 'Introgressed Hap. Frequency',
        }, inplace=True
    )
    pop_df.rename(
        columns={
            'spop': 'Super Population', 'pop': 'Population',
            'n_tot': 'Total Number of Chromosomes',
            'n_haps': 'Number of Introgressed Haps.',
            'freq': 'Introgressed Hap. Frequency',
        }, inplace=True
    )
    # Export.
    spop_df.to_csv(f'./dataframes/tgp_muc19_introgressed_hap_frequency_per_super_population_{window_size}kb.csv.gz', index=False)
    pop_df.to_csv(f'./dataframes/tgp_muc19_introgressed_hap_frequency_per_population_{window_size}kb.csv.gz', index=False)
    # Determine the total number of chromosomes.
    amr_chroms = spop_df[
        spop_df['Super Population'] == 'AMR'
    ]['Total Number of Chromosomes'].values[0]
    non_amr_chroms = spop_df[
        (spop_df['Super Population'] != 'AMR') & (spop_df['Super Population'] != 'AFR')
    ]['Total Number of Chromosomes'].sum()
    amr_v_non_amr_chroms = np.array([amr_chroms, non_amr_chroms])
    # Determine the total number of haps.
    amr_haps = spop_df[
        spop_df['Super Population'] == 'AMR'
    ]['Number of Introgressed Haps.'].values[0]
    non_amr_haps = spop_df[
        (spop_df['Super Population'] != 'AMR') & (spop_df['Super Population'] != 'AFR')
    ]['Number of Introgressed Haps.'].sum()
    amr_v_non_amr_haps = np.array([amr_haps, non_amr_haps])
    # Perform a z-proprtions test.
    z_stat, z_p = proportions_ztest(
        amr_v_non_amr_haps, amr_v_non_amr_chroms, alternative='larger',
    )
    # Build contigency tables.
    fet_table = [
        [amr_haps, amr_chroms - amr_haps],
        [non_amr_haps, non_amr_chroms - non_amr_haps],
    ]
    # Run a FET.
    odds_ratio, f_p = stats.fisher_exact(fet_table, alternative='greater')
    # Print the Results.
    amr_v_non_amr_freqs = amr_v_non_amr_haps / amr_v_non_amr_chroms
    tract_summary = f"""
    AMR vs Non-AMR Introgressed Haps. Summary {window_size}kb
    ===============================================


    Proportions Z-Test
    ------------------
    Population   Chromosomes   Haps.   Frequency
    ----------   -----------   -----   ---------
    AMR          {amr_chroms}           {amr_haps}      {amr_v_non_amr_freqs[0]}
    Non-AMR      {non_amr_chroms}          {non_amr_haps}      {amr_v_non_amr_freqs[1]}

    Z-statistic: {z_stat}
    P-value:     {z_p}


    Fisher's Exact Test
    -------------------
    Contingency Table:
                        Introgressed   Non-Introgressed
                        Haps.          Haps.
                        ------------   ---------------
                   AMR| {fet_table[0][0]}           {fet_table[0][1]}
               Non-AMR| {fet_table[1][0]}           {fet_table[1][1]}

    Odds Ratio: {odds_ratio}
    P-value:    {f_p}
    """
    print(tract_summary)
    return spop_df, pop_df

# Define a function to summarize the number of pairwise differences between the Denisovan-like haplotypes and the archaics.
def compile_den_like_haps_pw_diffs_summary():
    # Load the haplotype identities.
    _, hap_info = annotate_hap_identities(70, 100, 'tgp')
    # Load the sequence divergence results.
    div_72kb = pd.read_csv('./dataframes/tgp_haplotype_archaic_diplotype_divergence_72kb.csv.gz')
    div_phased = pd.read_csv('./dataframes/tgp_haplotype_late_neanderthal_phased_haplotype_divergence_72kb.csv.gz')
    # Extract the denisovan divergence results.
    den_div_72kb = div_72kb[div_72kb['Archaic'] == 'Denisovan']
    den_hap1_72kb = den_div_72kb['Focal 72kb Region (Pairwise Diffs. Hap. 1)'].values
    den_hap2_72kb = den_div_72kb['Focal 72kb Region (Pairwise Diffs. Hap. 2)'].values
    is_den_hap1_72kb = np.isin(den_div_72kb['Individual'].values, hap_info['IND'][hap_info['IS_HAP1_DEN_72KB']])
    is_den_hap2_72kb = np.isin(den_div_72kb['Individual'].values, hap_info['IND'][hap_info['IS_HAP2_DEN_72KB']])
    den_pwd = np.concatenate((den_hap1_72kb[is_den_hap1_72kb], den_hap2_72kb[is_den_hap2_72kb]))
    # Extract the altai neanderthal divergence results.
    alt_div_72kb = div_72kb[div_72kb['Archaic'] == 'Altai Nean.']
    alt_hap1_72kb = alt_div_72kb['Focal 72kb Region (Pairwise Diffs. Hap. 1)'].values
    alt_hap2_72kb = alt_div_72kb['Focal 72kb Region (Pairwise Diffs. Hap. 2)'].values
    is_alt_hap1_72kb = np.isin(alt_div_72kb['Individual'].values, hap_info['IND'][hap_info['IS_HAP1_DEN_72KB']])
    is_alt_hap2_72kb = np.isin(alt_div_72kb['Individual'].values, hap_info['IND'][hap_info['IS_HAP2_DEN_72KB']])
    alt_pwd = np.concatenate((alt_hap1_72kb[is_alt_hap1_72kb], alt_hap2_72kb[is_alt_hap2_72kb]))
    # Extract the chagyskyra divergence results.
    cha1_div_phased = div_phased[div_phased['Archaic Hap.'] == 'Chagyrskaya Nean. Hap. 1']
    cha1_hap1_72kb = cha1_div_phased['Focal 72kb Region (Pairwise Diffs. Hap. 1)'].values
    cha1_hap2_72kb = cha1_div_phased['Focal 72kb Region (Pairwise Diffs. Hap. 2)'].values
    is_cha1_hap1_72kb = np.isin(cha1_div_phased['Individual'].values, hap_info['IND'][hap_info['IS_HAP1_DEN_72KB']])
    is_cha1_hap2_72kb = np.isin(cha1_div_phased['Individual'].values, hap_info['IND'][hap_info['IS_HAP2_DEN_72KB']])
    cha1_pwd = np.concatenate((cha1_hap1_72kb[is_cha1_hap1_72kb], cha1_hap2_72kb[is_cha1_hap2_72kb]))
    cha2_div_phased = div_phased[div_phased['Archaic Hap.'] == 'Chagyrskaya Nean. Hap. 2']
    cha2_hap1_72kb = cha2_div_phased['Focal 72kb Region (Pairwise Diffs. Hap. 1)'].values
    cha2_hap2_72kb = cha2_div_phased['Focal 72kb Region (Pairwise Diffs. Hap. 2)'].values
    is_cha2_hap1_72kb = np.isin(cha2_div_phased['Individual'].values, hap_info['IND'][hap_info['IS_HAP1_DEN_72KB']])
    is_cha2_hap2_72kb = np.isin(cha2_div_phased['Individual'].values, hap_info['IND'][hap_info['IS_HAP2_DEN_72KB']])
    cha2_pwd = np.concatenate((cha2_hap1_72kb[is_cha2_hap1_72kb], cha2_hap2_72kb[is_cha2_hap2_72kb]))
    # Extract the vindija divergence results.
    vin1_div_phased = div_phased[div_phased['Archaic Hap.'] == 'Vindija Nean. Hap. 1']
    vin1_hap1_72kb = vin1_div_phased['Focal 72kb Region (Pairwise Diffs. Hap. 1)'].values
    vin1_hap2_72kb = vin1_div_phased['Focal 72kb Region (Pairwise Diffs. Hap. 2)'].values
    is_vin1_hap1_72kb = np.isin(vin1_div_phased['Individual'].values, hap_info['IND'][hap_info['IS_HAP1_DEN_72KB']])
    is_vin1_hap2_72kb = np.isin(vin1_div_phased['Individual'].values, hap_info['IND'][hap_info['IS_HAP2_DEN_72KB']])
    vin1_pwd = np.concatenate((vin1_hap1_72kb[is_vin1_hap1_72kb], vin1_hap2_72kb[is_vin1_hap2_72kb]))
    vin2_div_phased = div_phased[div_phased['Archaic Hap.'] == 'Vindija Nean. Hap. 2']
    vin2_hap1_72kb = vin2_div_phased['Focal 72kb Region (Pairwise Diffs. Hap. 1)'].values
    vin2_hap2_72kb = vin2_div_phased['Focal 72kb Region (Pairwise Diffs. Hap. 2)'].values
    is_vin2_hap1_72kb = np.isin(vin2_div_phased['Individual'].values, hap_info['IND'][hap_info['IS_HAP1_DEN_72KB']])
    is_vin2_hap2_72kb = np.isin(vin2_div_phased['Individual'].values, hap_info['IND'][hap_info['IS_HAP2_DEN_72KB']])
    vin2_pwd = np.concatenate((vin2_hap1_72kb[is_vin2_hap1_72kb], vin2_hap2_72kb[is_vin2_hap2_72kb]))
    # Intialize a dictionary to store the results.
    pwd_dicc = {
        'arc': ['Denisovan', 'Altai Nean.', 'Chagyrskaya Nean. Hap. 1', 'Chagyrskaya Nean. Hap. 2', 'Vindija Nean. Hap. 1', 'Vindija Nean. Hap. 2'],
        'mu': [], 'std': [], 'sem': [], 'ci': [],
    }
    # For each archaic.
    for arc_pwd in [den_pwd, alt_pwd, cha1_pwd, cha2_pwd, vin1_pwd, vin2_pwd]:
        # Compute the sem and 95% cis.
        sem, ci = sem_ci_of_mean(arc_pwd)
        # Update the dictionary.
        pwd_dicc['mu'].append(arc_pwd.mean())
        pwd_dicc['std'].append(arc_pwd.std())
        pwd_dicc['sem'].append(sem)
        pwd_dicc['ci'].append(ci)
    # Convert to a dataframe.
    pwd_df = pd.DataFrame(pwd_dicc)
    # Rename the columns to look pretty.
    pwd_df.rename(
        columns={
            'arc': 'Archaic',
            'mu': r'Focal 72kb Region Pairwise Diffs. $\left( \mu\right)$',
            'std': r'Focal 72kb Region Pairwise Diffs. $\left( \sigma \right)$',
            'sem': r'Focal 72kb Region Pairwise Diffs. $\left( SEM \right)$',
            'ci': r'Focal 72kb Region Pairwise Diffs. $\left( \pm CI_{95\%} \right)$',
        }, inplace=True
    )
    # Export the dataframe.
    pwd_df.to_csv('./dataframes/287_denisovan_haps_v_archaics_pairwise_diffs_72kb.csv.gz', index=False)
    return pwd_df

# Define a function to show the results for the recombinant haplotypes at the 72kb region.
def show_rec_haps_seq_div_72kb():
    # Intialize lists.
    rec_hap_ind_list = ['HG01284', 'HG03716', 'HG03920', 'HG04131', 'NA20858', 'NA21141', 'HG03012']
    # Load the sequence divergence results.
    tgp_hap_div_72kb_df = pd.read_csv('./dataframes/tgp_haplotype_archaic_diplotype_divergence_72kb.csv.gz')
    # Subset the focal individuals
    rec_hap_div_72kb_df = tgp_hap_div_72kb_df[np.isin(tgp_hap_div_72kb_df['Individual'].values, rec_hap_ind_list)]
    return rec_hap_div_72kb_df

# Define a function to plot the recombinant haplotypes.
def plot_rec_hap_mat_72kb():
    # Load the phased data information.
    phased_df = pd.read_csv('../meta_data/altai_nean_phased_late_neanderthals_denisovan_genos_72kb.csv.gz')
    # Load the meta information
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the focal indicies.
    mxl_idx = np.where(tgp_df['IND'].values == 'NA19664')[0]
    yri_idx = np.where(tgp_df['IND'].values == 'NA19190')[0]
    # Load the genotype matrix.
    tgp_arcs_no_aa_72kb_gt, _ = load_hap_region('tgp_arcs_masked_no_aa', 12, 40759001, 40831000)
    # Update the dataframe with the focal tgp haplotypes.
    phased_df['MXL_NA19664_HAP'] = calc_ind_alt_freqs(tgp_arcs_no_aa_72kb_gt.take(mxl_idx, axis=1))
    phased_df['YRI_NA19190_HAP'] = calc_ind_alt_freqs(tgp_arcs_no_aa_72kb_gt.take(yri_idx, axis=1))
    # Intialize a list for the recombinant haplotype individuals.
    for pop, ind, hap_idx in [
        ('GIH', 'NA20858', 0), ('GIH', 'NA21141', 0), ('BEB', 'HG04131', 0), ('BEB', 'HG03920', 0),
        ('BEB', 'HG03012', 1), ('ITU', 'HG03716', 0), ('CLM', 'HG01284', 0),
    ]:
        # Determine the individual's index.
        ind_idx = np.where(tgp_df['IND'].values == ind)[0][0]
        # Update the phased information.
        phased_df[f'{pop}_{ind}_HAP'] = tgp_arcs_no_aa_72kb_gt[:, ind_idx, hap_idx]
    # Intialize the haplotype columns.
    hap_cols = phased_df.columns.values[3:]
    # Determine the segregatinge sites among these samples.
    seg_phased_df = phased_df[((phased_df[hap_cols] == 1).any(axis=1)) & ((phased_df[hap_cols] == 0).any(axis=1))]
    # Print the number of segregating sites.
    print(f'Seg. Sites: {seg_phased_df.shape[0]}')
    # Convert the dataframe to a dictionary.
    phased_info = {col: seg_phased_df[col].values for col in seg_phased_df.columns.values}
    # Intialize a haplotype list.
    hap_list = [
        'MXL_NA19664_HAP', 'Chagyrskaya Nean. Hap. 2', 'Vindija Nean. Hap. 2',
        'Denisovan', 'GIH_NA20858_HAP', 'GIH_NA21141_HAP', 'BEB_HG04131_HAP', 'BEB_HG03920_HAP',
        'BEB_HG03012_HAP', 'ITU_HG03716_HAP', 'CLM_HG01284_HAP', 'Altai Nean.',
        'Vindija Nean. Hap. 1', 'Chagyrskaya Nean. Hap. 1', 'YRI_NA19190_HAP',
    ]
    # Intialize a haplotype matrix for plotting.
    hap_mat = np.empty((len(hap_list), seg_phased_df.shape[0]))
    # For every focal haplotype.
    for i, hap in enumerate(hap_list):
        # Update the haplotype matrix.
        hap_mat[i, :] = phased_info[hap]
    # Intialize the haplotype labels.
    hap_labs = np.array([
        'MXL (NA19664)',
        'Chagyrskaya Nean.'+'\n'+'(Hap. 2)     ',
        'Vindija Nean.'+'\n'+'(Hap. 2)     ',
        'Denisovan',
        'GIH (NA20858)', 'GIH (NA21141)', 'BEB (HG04131)', 'BEB (HG03920)', 'BEB (HG03012)', 'ITU (HG03716)', 'CLM (HG01284)',
        'Altai Nean.',
        'Chagyrskaya Nean.'+'\n'+'(Hap. 1)     ',
        'Vindija Nean.'+'\n'+'(Hap. 1)     ',
        'YRI (NA19190)',
    ])
    # Intialize the figure.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.spines.top': True,
        'axes.spines.right': True,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 14,
        'legend.fontsize': 12,
        'legend.title_fontsize': 12,
    })
    # Intialize the figure.
    fig = plt.figure(
        figsize=(12, 6), dpi=300,
        facecolor='white',
    )
    # Intialize the axes.
    ax = fig.add_subplot(111)
    # Intialize color map.
    cmap = ListedColormap(['#0072B2', '#CC79A7'])
    # Plot the haplotype table.
    im = ax.imshow(hap_mat, cmap=cmap, aspect='auto')
    # Label the rows.
    ax.set_yticks(np.arange(hap_labs.size))
    ax.set_yticklabels(hap_labs)
    # Seperate each box and add a grid.
    ax.set_xticks(np.arange(0, hap_mat.shape[1], 1))
    ax.set_yticks(np.arange(0, hap_mat.shape[0], 1))
    ax.set_xticks(np.arange(-0.5, hap_mat.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, hap_mat.shape[0], 1), minor=True)
    ax.grid(which='minor', color='black', linestyle='-', linewidth=0.5)
    # Remove the ticks.
    ax.tick_params(bottom=False, labelbottom=False)
    ax.tick_params(which='minor', left=False, bottom=False, labelbottom=False)
    # Enforce a tight layout.
    plt.tight_layout()
    # Export the plot.
    plt.savefig(
        './supp_figures/png/tgp_recombiant_hap_mat_72kb.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/tgp_recombiant_hap_mat_72kb.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/tgp_recombiant_hap_mat_72kb.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to load the sequence divergence results at papuan introgressed tracts.
def load_pap_v_den_seq_div_at_intro_tracts():
    # Define a list of papuan individuals in the sgdp.
    pap_sgdp = np.array([
        'S_Papuan-1.DG', 'S_Papuan-2.DG', 'S_Papuan-3.DG', 'S_Papuan-4.DG',
        'S_Papuan-5.DG', 'S_Papuan-6.DG', 'S_Papuan-7.DG', 'S_Papuan-8.DG',
        'S_Papuan-9.DG', 'S_Papuan-10.DG', 'S_Papuan-11.DG', 'S_Papuan-12.DG',
        'S_Papuan-13.DG', 'S_Papuan-14.DG', 'B_Papuan-15.DG',
    ])
    # Intialize a list to store the results.
    pap_div = []
    # For every chromosome.
    for chrom in range(1, 23):
        # For every papuan individual.
        for pap in pap_sgdp:
            # For every haplotype.
            for hap in ['hap1', 'hap2']:
                # Load the results.
                pap_div_chrom = np.loadtxt(
                    f'../muc19_results/sgdp_den_masked_no_aa/{pap}_{hap}_den_pw_diffs_seq_div_at_den_intro_tracts_chr{chrom}.txt.gz'
                )
                # Account for the edge case where a papuan individual only has one denisovan tract on the chromosome.
                if pap_div_chrom.size == 2:
                    # Update the list.
                    pap_div.append(np.array([pap_div_chrom]))
                # Account for the edge case where S_Papuan-11.DG_hap1 has no denisovan tracts on chromosome 10 and that S_Papuan-5.DG_hap2 has no denisovan tracts on chromosome 14.
                elif pap_div_chrom.size > 0:
                    # Update the list.
                    pap_div.append(pap_div_chrom)
    # Concatenate the windows.
    pap_div = np.concatenate(pap_div)[:, 1]
    # Compile the summary.
    div_m, div_s = np.nanmean(pap_div), np.nanstd(pap_div)
    div_se, div_ci = sem_ci_of_mean(pap_div)
    # Print a summary.
    print('#***# Papuan v Denisovan Sequence Divergence at Denisovan Introgressed Tracts #***#')
    print(f'Mean: {div_m}')
    print(f'Std Dev: {div_s}')
    print(f'SEM: {div_se}')
    print(f'+/- 95% CIs: {div_ci}')
    return pap_div

# Define a function to plot the distribution of sequence divergence for papuans at denisovan introgressed regions.
def plot_pap_v_den_seq_div_at_intro_tracts():
    # Load the meta data.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    sgdp_df = pd.read_csv(
        '../meta_data/sgdp.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the focal indicies.
    mxl_idx = np.where(tgp_df['IND'].values == 'NA19664')[0][0]
    pap_idx = np.where(sgdp_df['IND'].values == 'B_Papuan-15.DG')[0][0]
    # Load the effective sequence lengths.
    tgp_esl = load_region_esl('tgp_den_masked_no_aa', 72)
    sgdp_esl = load_region_esl('sgdp_den_masked_no_aa', 72)
    # Load the genotype matricies.
    tgp_gt, _ = load_hap_region('tgp_den_masked_no_aa', 12, 40759001, 40831000)
    sgdp_gt, _ = load_hap_region('sgdp_den_masked_no_aa', 12, 40759001, 40831000)
    # Compute the observed sequence divergences.
    mxl_pwd = np.nansum(pwd_per_site(
        calc_ind_alt_freqs(tgp_gt.take([mxl_idx], axis=1)),
        calc_ind_alt_freqs(tgp_gt.take([2347], axis=1)),
    ))
    pap_pwd = np.nansum(pwd_per_site(
        calc_ind_alt_freqs(sgdp_gt.take([pap_idx], axis=1)),
        calc_ind_alt_freqs(sgdp_gt.take([278], axis=1)),
    ))
    mxl_div = mxl_pwd / tgp_esl
    pap_div = pap_pwd / sgdp_esl
    # Load the distribution.
    pap_dist = load_pap_v_den_seq_div_at_intro_tracts()
    # Compute the percentile ranks.
    mxl_pr = stats.percentileofscore(pap_dist, mxl_div, kind='strict')
    pap_pr = stats.percentileofscore(pap_dist, pap_div, kind='strict')
    # Print a summary.
    print(f'NA19664 v Denisovan Pairwise Diffs: {mxl_pwd}')
    print(f'NA19664 v Denisovan Seq Div: {mxl_div}')
    print(f'NA19664 v Denisovan Seq Div Percentile Rank: {mxl_pr}')
    print(f'B_Papuan-15.DG v Denisovan Pairwise Diffs: {pap_pwd}')
    print(f'B_Papuan-15.DG v Denisovan Seq Div: {pap_div}')
    print(f'B_Papuan-15.DG v Denisovan Seq Div Percentile Rank: {pap_pr}')
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the figure and axes.
    fig = plt.figure(
        figsize=(6, 4), dpi=300,
        facecolor='white', constrained_layout=True,
    )
    ax = fig.add_subplot(111)
    # Plot the distribution.
    ax.hist(
        pap_dist, bins=np.linspace(0, pap_dist.max(), 100),
        histtype='stepfilled', color='gray',
    )
    # Add axes labels.
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Papuan v Denisovan Seq. Div. per Introgressed Tract')
    # Annotate the observed value.
    for obs, color in [(mxl_div, '#009E73'), (pap_div, '#F0E442')]:
        ax.annotate(
            '', xy=(obs, 0),
            xytext=(obs, ax.get_ylim()[1] * 0.5),
            arrowprops=dict(color=color, arrowstyle='->', lw=1),
        )
    # Export the plot.
    plt.savefig(
        './supp_figures/png/sequence_divergence_at_denisovan_introgressed_tracts.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/svg/sequence_divergence_at_denisovan_introgressed_tracts.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        './supp_figures/pdf/sequence_divergence_at_denisovan_introgressed_tracts.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to calculate sequence divergence between archaic diplotypes and modern human haplotypes.
def calc_hum_hap_v_arc_dip_diffs(gt, tgp_or_sgdp):
    tgp_or_sgdp_dicc = {
        'tgp': {'meta_path': '../meta_data/tgp_mod.txt', 'arc_idx': 2347},
        'sgdp': {'meta_path': '../meta_data/sgdp.txt', 'arc_idx': 278},
    }
    # Load the meta data.
    meta_data = pd.read_csv(
        tgp_or_sgdp_dicc[tgp_or_sgdp]['meta_path'], sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store the results.
    pwd_dicc = {'hap_1': np.array([]), 'hap_2': np.array([])}
    # Intialize the archaic's allele frequency.
    arc_freq = calc_ind_alt_freqs(gt.take([tgp_or_sgdp_dicc[tgp_or_sgdp]['arc_idx']], axis=1))
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

# Define a function to intialize the phased haplotype divergences.
def init_phased_hap_dist():
    # Intilize a list and a dictioanry.
    hap_list = ['hap_1', 'hap_2']
    dists = {hap: {'pwd': {}, 'div': {}} for hap in hap_list}
    # Load the phased haplotypes.
    phased_df = pd.read_csv('../meta_data/altai_nean_phased_late_neanderthals_denisovan_genos_72kb.csv.gz')
    # Intialize a dictionary to store the archaic haplotypes.
    arc_haps = {
        'CHA_1': phased_df['Chagyrskaya Nean. Hap. 1'].values,
        'CHA_2': phased_df['Chagyrskaya Nean. Hap. 2'].values,
        'VIN_1': phased_df['Vindija Nean. Hap. 1'].values,
        'VIN_2': phased_df['Vindija Nean. Hap. 2'].values,
    }
    # Load the effective sequence lengths.
    arc_esl = load_region_esl('tgp_arcs_masked_no_aa', 72)
    # Adjust the effective sequence lengths to account for the sites that could not be resolved.
    cha_esl = arc_esl[1] - 2
    vin_esl = arc_esl[2] - 3
    # Intialize the effective sequence length dictionary.
    esl = {
        'CHA_1': cha_esl, 'CHA_2': cha_esl,
        'VIN_1': vin_esl, 'VIN_2': vin_esl,
    }
    # Load the genotypes for the focal region.
    gt_72kb, _ = load_hap_region('tgp_arcs_masked_no_aa', 12, 40759001, 40831000)
    # For every archaic.
    for nea, nea_hap in arc_haps.items():
        # Compute the number of pairwise differences.
        pwd = calc_hum_hap_v_cha_vin_hap_diffs(gt_72kb, nea_hap)
        # Compute sequence divergence.
        div = {hap: pwd[hap] / esl[nea] for hap in hap_list}
        # Update the dictionary.
        for hap in hap_list:
            dists[hap]['pwd'][nea] = pwd[hap]
            dists[hap]['div'][nea] = div[hap]
    return dists
    
# Define a function to intialize the haplotype distances.
def init_hap_dist_focal_regions(tgp_or_sgdp):
    # Intilize a lists.
    arc_list = ['DEN', 'ALT', 'CHA', 'VIN']
    hap_list = ['hap_1', 'hap_2']
    # Intialize the effective sequence length dictionary.
    esl_742kb = {
        arc: int(load_region_esl(f'{tgp_or_sgdp}_{arc.lower()}_masked_no_aa', 742)) for arc in arc_list
    }
    esl_72kb = {
        arc: int(load_region_esl(f'{tgp_or_sgdp}_{arc.lower()}_masked_no_aa', 72)) for arc in arc_list
    }
    # Intialize a dictionary.
    dist_742kb = {hap: {'pwd': {}, 'div': {}} for hap in hap_list}
    dist_72kb = {hap: {'pwd': {}, 'div': {}} for hap in hap_list}
    # For every arachaic.
    for arc in arc_list:
        # Load the genotypes for the focal regions.
        gt_742kb, pos_742kb = load_hap_region(f'{tgp_or_sgdp}_{arc.lower()}_masked_no_aa', 12, 40272001, 41014000)
        gt_72kb, pos_72kb = load_hap_region(f'{tgp_or_sgdp}_{arc.lower()}_masked_no_aa', 12, 40759001, 40831000)
        # Compute the number of pairwise differences.
        pwd_742kb = calc_hum_hap_v_arc_dip_diffs(gt_742kb, tgp_or_sgdp)
        pwd_72kb = calc_hum_hap_v_arc_dip_diffs(gt_72kb, tgp_or_sgdp)
        # Compute sequence divergence.
        div_742kb = {hap: pwd_742kb[hap] / esl_742kb[arc] for hap in hap_list}
        div_72kb = {hap: pwd_72kb[hap] / esl_72kb[arc] for hap in hap_list}
        # Update the dictionary.
        for hap in hap_list:
            dist_742kb[hap]['pwd'][arc] = pwd_742kb[hap]
            dist_742kb[hap]['div'][arc] = div_742kb[hap]
            dist_72kb[hap]['pwd'][arc] = pwd_72kb[hap]
            dist_72kb[hap]['div'][arc] = div_72kb[hap]
    return dist_742kb, dist_72kb

# Define a function to plot two distributions.
def plot_tgp_joint_seq_div_72kb_unphased(
    dist_dicc, arc_x, arc_y, is_equal_dims=True,
):
    # Intialize convience dictionaries.
    arc_labs = {
        'DEN': 'Denisovan', 'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.', 'VIN': 'Vindija Nean.',
    }
    color_dicc = {'AMR': '#009E73', 'SAS': '#CC79A7', 'EAS': '#0072B2', 'EUR': '#D55E00', 'AFR': 'grey'}
    # Load the meta information.
    _, info_dicc = annotate_hap_identities(70, 100, 'tgp')
    # Concatenate the distances, population IDs, and haplotype masks.
    x_dist = np.concatenate((
        dist_dicc['hap_1']['div'][arc_x],
        dist_dicc['hap_2']['div'][arc_x],
    ))
    y_dist = np.concatenate((
        dist_dicc['hap_1']['div'][arc_y],
        dist_dicc['hap_2']['div'][arc_y],
    ))
    spop_ids = np.concatenate((
        info_dicc['SUPERPOP'],
        info_dicc['SUPERPOP'],
    ))
    is_den_hap = np.concatenate((
        info_dicc['IS_HAP1_DEN_72KB'],
        info_dicc['IS_HAP2_DEN_72KB'],
    ))
    is_rec_hap = np.concatenate((
        info_dicc['IS_HAP1_REC_72KB'],
        info_dicc['IS_HAP2_REC_72KB'],
    ))
    # Shuffle the indicies.
    np.random.seed(37)
    shuffled_idx = np.random.choice(np.arange(spop_ids.size), size=(spop_ids.size), replace=False)
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the figure.
    fig = plt.figure(
        figsize=(6, 4.5), dpi=300,
        facecolor='white',
    )
    # Define axis positions.
    left = 0.1
    bottom = 0.1
    width = 0.75
    height = 0.75
    gap = 0.025
    margin_size = 0.25
    scatter_ax = fig.add_axes([left, bottom, width, height])  # left, bottom, width, height
    hist_x_ax = fig.add_axes([left, bottom + height + gap, width, margin_size], sharex=scatter_ax)
    hist_y_ax = fig.add_axes([left + width + gap, bottom, margin_size, height], sharey=scatter_ax)
    # For every haplotype.
    for idx in shuffled_idx:
        # If this is a introgressed haplotype.
        if is_den_hap[idx]:
            # Plot the joint distribution.
            scatter_ax.scatter(
                x_dist[idx], y_dist[idx],
                marker='s', facecolor='none',
                edgecolor=color_dicc[spop_ids[idx]],
                alpha=0.5, linewidth=1, s=35,
            )
        # Else-if this is a recomibinant haplotype.
        elif is_rec_hap[idx]:
            # Plot the joint distribution.
            scatter_ax.scatter(
                x_dist[idx], y_dist[idx],
                marker='^', facecolor='none',
                edgecolor=color_dicc[spop_ids[idx]],
                alpha=0.5, linewidth=1, s=35,
            )
        # Else this is a human haplotype.
        else:
            # Plot the joint distribution.
            scatter_ax.scatter(
                x_dist[idx], y_dist[idx],
                marker='o', facecolor='none',
                edgecolor=color_dicc[spop_ids[idx]],
                alpha=0.5, linewidth=1, s=35,
            )
    # If the fucking reviwer really wants the same dimensions.
    if is_equal_dims:
        # Get the current limits of the x and y axes.
        xlim = scatter_ax.get_xlim()
        ylim = scatter_ax.get_ylim()
        # Determine a common limit.
        min_lim = min(xlim[0], ylim[0])
        max_lim = max(xlim[1], ylim[1])
        # Set both axes to use the same limits.
        scatter_ax.set_xlim(min_lim, max_lim)
        scatter_ax.set_ylim(min_lim, max_lim)
    # Add labels.
    scatter_ax.set_xlabel(f'1KG Hap. v {arc_labs[arc_x]} Seq. Div.')
    scatter_ax.set_ylabel(f'1KG Hap. v {arc_labs[arc_y]} Seq. Div.')
    # Plot the marginal distribution for the x-axis.
    counts_x, _, _ = hist_x_ax.hist(
        x_dist, bins='fd',
        histtype='stepfilled', color='grey',
    )
    hist_x_ax.set_ylabel('Num. of Haps.')
    # Plot the marginal distribution for the y-axis.
    counts_y, _, _ = hist_y_ax.hist(
        y_dist, orientation='horizontal', bins='fd',
        histtype='stepfilled', color='grey',
    )
    hist_y_ax.set_xlabel('Num. of Haps.')
    # Determine the maximum count from both histograms.
    max_count = max(counts_x.max(), counts_y.max())
    hist_ticks = np.arange(0, 1001, 250)
    hist_x_ax.set_yticks(hist_ticks)
    hist_x_ax.set_yticklabels(hist_ticks.astype(str))
    hist_y_ax.set_xticks(hist_ticks)
    hist_y_ax.set_xticklabels(hist_ticks.astype(str))
    # Update the count ticks and disable tick labels on the histograms
    hist_x_ax.set_ylim(0, max_count)
    hist_y_ax.set_xlim(0, max_count)
    plt.setp(hist_x_ax.get_xticklabels(), visible=False)
    plt.setp(hist_y_ax.get_yticklabels(), visible=False)
    # Intialize a figure legend.
    legend_handles = [Patch(color=color, label=spop) for spop, color in color_dicc.items()]
    legend_handles.extend([
        Line2D([0], [0], linestyle='none', marker='s', markersize=7.5, markerfacecolor='none', markeredgecolor='black', label=r'$\alpha$ Cluster'),
        Line2D([0], [0], linestyle='none', marker='o', markersize=7.5, markerfacecolor='none', markeredgecolor='black', label=r'$\beta$ Cluster'),
        Line2D([0], [0], linestyle='none', marker='^', markersize=7.5, markerfacecolor='none', markeredgecolor='black', label=r'$\gamma$ Cluster'),
    ])
    # Add the figure legend.
    title_font = fm.FontProperties(weight='bold', size=10)
    fig.legend(
        handles=legend_handles, 
        loc='lower center', bbox_to_anchor=(0.5, -0.175), ncol=4,
        frameon=True, fancybox=True, shadow=True, fontsize=10,
        title='72kb Region', title_fontproperties=title_font,
    )
        # Export the plot.
    plt.savefig(
        f'./supp_figures/png/tgp_joint_seq_div_from_{arc_x.lower()}_and_{arc_y.lower()}_72kb_region.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/tgp_joint_seq_div_from_{arc_x.lower()}_and_{arc_y.lower()}_72kb_region.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/tgp_joint_seq_div_from_{arc_x.lower()}_and_{arc_y.lower()}_72kb_region.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to plot two distributions.
def plot_tgp_joint_seq_div_742kb(
    dist_dicc, arc_x, arc_y, is_equal_dims=True,
):
    # Intialize convience dictionaries.
    arc_labs = {
        'DEN': 'Denisovan', 'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.', 'VIN': 'Vindija Nean.',
    }
    color_dicc = {'AMR': '#009E73', 'SAS': '#CC79A7', 'EAS': '#0072B2', 'EUR': '#D55E00', 'AFR': 'grey'}
    # Load the meta information.
    _, info_dicc = annotate_hap_identities(70, 100, 'tgp')
    # Concatenate the distances, population IDs, and haplotype masks.
    x_dist = np.concatenate((
        dist_dicc['hap_1']['div'][arc_x],
        dist_dicc['hap_2']['div'][arc_x],
    ))
    y_dist = np.concatenate((
        dist_dicc['hap_1']['div'][arc_y],
        dist_dicc['hap_2']['div'][arc_y],
    ))
    spop_ids = np.concatenate((
        info_dicc['SUPERPOP'],
        info_dicc['SUPERPOP'],
    ))
    is_arc_hap = np.concatenate((
        info_dicc['IS_HAP1_ARC_742KB'],
        info_dicc['IS_HAP2_ARC_742KB'],
    ))
    # Shuffle the indicies.
    np.random.seed(37)
    shuffled_idx = np.random.choice(np.arange(spop_ids.size), size=(spop_ids.size), replace=False)
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the figure.
    fig = plt.figure(
        figsize=(6, 4.5), dpi=300,
        facecolor='white',
    )
    # Define axis positions.
    left = 0.1
    bottom = 0.1
    width = 0.75
    height = 0.75
    gap = 0.025
    margin_size = 0.25
    scatter_ax = fig.add_axes([left, bottom, width, height])  # left, bottom, width, height
    hist_x_ax = fig.add_axes([left, bottom + height + gap, width, margin_size], sharex=scatter_ax)
    hist_y_ax = fig.add_axes([left + width + gap, bottom, margin_size, height], sharey=scatter_ax)
    # For every haplotype.
    for idx in shuffled_idx:
        # If this is a introgressed haplotype.
        if is_arc_hap[idx]:
            # Plot the joint distribution.
            scatter_ax.scatter(
                x_dist[idx], y_dist[idx],
                marker='s', facecolor='none',
                edgecolor=color_dicc[spop_ids[idx]],
                alpha=0.5, linewidth=1, s=35,
            )
        # Else this is a human haplotype.
        else:
            # Plot the joint distribution.
            scatter_ax.scatter(
                x_dist[idx], y_dist[idx],
                marker='o', facecolor='none',
                edgecolor=color_dicc[spop_ids[idx]],
                alpha=0.5, linewidth=1, s=35,
            )
    # If the fucking reviwer really wants the same dimensions.
    if is_equal_dims:
        # Get the current limits of the x and y axes.
        xlim = scatter_ax.get_xlim()
        ylim = scatter_ax.get_ylim()
        # Determine a common limit.
        min_lim = min(xlim[0], ylim[0])
        max_lim = max(xlim[1], ylim[1])
        # Set both axes to use the same limits.
        scatter_ax.set_xlim(min_lim, max_lim)
        scatter_ax.set_ylim(min_lim, max_lim)
    # Add labels.
    scatter_ax.set_xlabel(f'1KG Hap. v {arc_labs[arc_x]} Seq. Div.')
    scatter_ax.set_ylabel(f'1KG  Hap. v {arc_labs[arc_y]} Seq. Div.')
    # Plot the marginal distribution for the x-axis.
    counts_x, _, _ = hist_x_ax.hist(
        x_dist, bins='fd',
        histtype='stepfilled', color='grey',
    )
    hist_x_ax.set_ylabel('Num. of Haps.')
    # Plot the marginal distribution for the y-axis.
    counts_y, _, _ = hist_y_ax.hist(
        y_dist, orientation='horizontal', bins='fd',
        histtype='stepfilled', color='grey',
    )
    hist_y_ax.set_xlabel('Num. of Haps.')
    # Determine the maximum count from both histograms.
    max_count = max(counts_x.max(), counts_y.max())
    hist_ticks = np.arange(0, 401, 100)
    hist_x_ax.set_yticks(hist_ticks)
    hist_x_ax.set_yticklabels(hist_ticks.astype(str))
    hist_y_ax.set_xticks(hist_ticks)
    hist_y_ax.set_xticklabels(hist_ticks.astype(str))
    # Update the count ticks and disable tick labels on the histograms
    hist_x_ax.set_ylim(0, max_count)
    hist_y_ax.set_xlim(0, max_count)
    plt.setp(hist_x_ax.get_xticklabels(), visible=False)
    plt.setp(hist_y_ax.get_yticklabels(), visible=False)
    # Intialize a figure legend.
    legend_handles = [Patch(color=color, label=spop) for spop, color in color_dicc.items()]
    legend_handles.extend([
        Line2D([0], [0], linestyle='none', marker='o', markersize=7.5, markerfacecolor='none', markeredgecolor='black', label='Non-Introgressed Hap.'),
        Line2D([0], [0], linestyle='none', marker='s', markersize=7.5, markerfacecolor='none', markeredgecolor='black', label='Introgressed Hap.'),
    ])
    # Add the figure legend.
    title_font = fm.FontProperties(weight='bold', size=10)
    fig.legend(
        handles=legend_handles, 
        loc='lower center', bbox_to_anchor=(0.5, -0.15), ncol=4,
        frameon=True, fancybox=True, shadow=True, fontsize=10,
        title='742kb Region', title_fontproperties=title_font,
    )
        # Export the plot.
    plt.savefig(
        f'./supp_figures/png/tgp_joint_seq_div_from_{arc_x.lower()}_and_{arc_y.lower()}_742kb_region.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/tgp_joint_seq_div_from_{arc_x.lower()}_and_{arc_y.lower()}_742kb_region.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/tgp_joint_seq_div_from_{arc_x.lower()}_and_{arc_y.lower()}_742kb_region.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to plot two distributions.
def plot_tgp_joint_seq_div_72kb_phased(
    unphased_dicc, phased_dicc, arc_x, arc_y, is_equal_dims=True,
):
    # Intialize convience dictionaries.
    arc_labs = {
        'DEN': 'Denisovan', 'ALT': 'Altai Nean.',
        'CHA_1': 'Chagyrskaya Nean. Hap. 1',
        'CHA_2': 'Chagyrskaya Nean. Hap. 2',
        'VIN_1': 'Vindija Nean. Hap. 1',
        'VIN_2': 'Vindija Nean. Hap. 2',
    }
    color_dicc = {'AMR': '#009E73', 'SAS': '#CC79A7', 'EAS': '#0072B2', 'EUR': '#D55E00', 'AFR': 'grey'}
    # Load the meta information.
    _, info_dicc = annotate_hap_identities(70, 100, 'tgp')
    # Concatenate the distances, population IDs, and haplotype masks.
    x_dist = np.concatenate((
        unphased_dicc['hap_1']['div'][arc_x],
        unphased_dicc['hap_2']['div'][arc_x],
    ))
    y_dist = np.concatenate((
        phased_dicc['hap_1']['div'][arc_y],
        phased_dicc['hap_2']['div'][arc_y],
    ))
    spop_ids = np.concatenate((
        info_dicc['SUPERPOP'],
        info_dicc['SUPERPOP'],
    ))
    is_den_hap = np.concatenate((
        info_dicc['IS_HAP1_DEN_72KB'],
        info_dicc['IS_HAP2_DEN_72KB'],
    ))
    is_rec_hap = np.concatenate((
        info_dicc['IS_HAP1_REC_72KB'],
        info_dicc['IS_HAP2_REC_72KB'],
    ))
    # Shuffle the indicies.
    np.random.seed(37)
    shuffled_idx = np.random.choice(np.arange(spop_ids.size), size=(spop_ids.size), replace=False)
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize the figure.
    fig = plt.figure(
        figsize=(6, 4.5), dpi=300,
        facecolor='white',
    )
    # Define axis positions.
    left = 0.1
    bottom = 0.1
    width = 0.75
    height = 0.75
    gap = 0.025
    margin_size = 0.25
    scatter_ax = fig.add_axes([left, bottom, width, height])  # left, bottom, width, height
    hist_x_ax = fig.add_axes([left, bottom + height + gap, width, margin_size], sharex=scatter_ax)
    hist_y_ax = fig.add_axes([left + width + gap, bottom, margin_size, height], sharey=scatter_ax)
    # For every haplotype.
    for idx in shuffled_idx:
        # If this is a introgressed haplotype.
        if is_den_hap[idx]:
            # Plot the joint distribution.
            scatter_ax.scatter(
                x_dist[idx], y_dist[idx],
                marker='s', facecolor='none',
                edgecolor=color_dicc[spop_ids[idx]],
                alpha=0.5, linewidth=1, s=35,
            )
        # Else-if this is a recomibinant haplotype.
        elif is_rec_hap[idx]:
            # Plot the joint distribution.
            scatter_ax.scatter(
                x_dist[idx], y_dist[idx],
                marker='^', facecolor='none',
                edgecolor=color_dicc[spop_ids[idx]],
                alpha=0.5, linewidth=1, s=35,
            )
        # Else this is a human haplotype.
        else:
            # Plot the joint distribution.
            scatter_ax.scatter(
                x_dist[idx], y_dist[idx],
                marker='o', facecolor='none',
                edgecolor=color_dicc[spop_ids[idx]],
                alpha=0.5, linewidth=1, s=35,
            )
    # If the fucking reviwer really wants the same dimensions.
    if is_equal_dims:
        # Get the current limits of the x and y axes.
        xlim = scatter_ax.get_xlim()
        ylim = scatter_ax.get_ylim()
        # Determine a common limit.
        min_lim = min(xlim[0], ylim[0])
        max_lim = max(xlim[1], ylim[1])
        # Set both axes to use the same limits.
        scatter_ax.set_xlim(min_lim, max_lim)
        scatter_ax.set_ylim(min_lim, max_lim)
    # Add labels.
    scatter_ax.set_xlabel(f'1KG Hap. v {arc_labs[arc_x]} Seq. Div.')
    scatter_ax.set_ylabel(f'1KG Hap. v {arc_labs[arc_y]} Seq. Div.')
    # Plot the marginal distribution for the x-axis.
    counts_x, _, _ = hist_x_ax.hist(
        x_dist, bins='fd',
        histtype='stepfilled', color='grey',
    )
    hist_x_ax.set_ylabel('Num. of Haps.')
    # Plot the marginal distribution for the y-axis.
    counts_y, _, _ = hist_y_ax.hist(
        y_dist, orientation='horizontal', bins='fd',
        histtype='stepfilled', color='grey',
    )
    hist_y_ax.set_xlabel('Num. of Haps.')
    # Determine the maximum count from both histograms.
    max_count = max(counts_x.max(), counts_y.max())
    hist_ticks = np.arange(0, 1001, 250)
    hist_x_ax.set_yticks(hist_ticks)
    hist_x_ax.set_yticklabels(hist_ticks.astype(str))
    hist_y_ax.set_xticks(hist_ticks)
    hist_y_ax.set_xticklabels(hist_ticks.astype(str))
    # Update the count ticks and disable tick labels on the histograms
    hist_x_ax.set_ylim(0, max_count)
    hist_y_ax.set_xlim(0, max_count)
    plt.setp(hist_x_ax.get_xticklabels(), visible=False)
    plt.setp(hist_y_ax.get_yticklabels(), visible=False)
    # Intialize a figure legend.
    legend_handles = [Patch(color=color, label=spop) for spop, color in color_dicc.items()]
    legend_handles.extend([
        Line2D([0], [0], linestyle='none', marker='s', markersize=7.5, markerfacecolor='none', markeredgecolor='black', label=r'$\alpha$ Cluster'),
        Line2D([0], [0], linestyle='none', marker='o', markersize=7.5, markerfacecolor='none', markeredgecolor='black', label=r'$\beta$ Cluster'),
        Line2D([0], [0], linestyle='none', marker='^', markersize=7.5, markerfacecolor='none', markeredgecolor='black', label=r'$\gamma$ Cluster'),
    ])
    # Add the figure legend.
    title_font = fm.FontProperties(weight='bold', size=10)
    fig.legend(
        handles=legend_handles, 
        loc='lower center', bbox_to_anchor=(0.5, -0.175), ncol=4,
        frameon=True, fancybox=True, shadow=True, fontsize=10,
        title='72kb Region', title_fontproperties=title_font,
    )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/tgp_joint_seq_div_from_{arc_x.lower()}_and_{arc_y[:3].lower()}_hap{arc_y[-1]}_72kb_region.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/tgp_joint_seq_div_from_{arc_x.lower()}_and_{arc_y[:3].lower()}_hap{arc_y[-1]}_72kb_region.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/tgp_joint_seq_div_from_{arc_x.lower()}_and_{arc_y[:3].lower()}_hap{arc_y[-1]}_72kb_region.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to calculate sequence divergence between archaic diplotypes.
def calc_alt_v_den_diffs(gt, window_size):
    # Intialize a list of archaics and there indicies.
    arc_idx_dicc = {
        'DEN': 3, 'ALT': 0,
    }
    # Load the effective sequence lengths.
    esl = load_region_esl('arcs_masked_no_aa', window_size=window_size)[4]
    # Determine the alternative allele frequencies per archaic.
    arc_1_freq = calc_ind_alt_freqs(gt.take([arc_idx_dicc['DEN']], axis=1))
    arc_2_freq = calc_ind_alt_freqs(gt.take([arc_idx_dicc['ALT']], axis=1))
    # Compute the number of pairwise differences.
    pwd = np.nansum(pwd_per_site(arc_1_freq, arc_2_freq))
    # Compute the sequence divergence.
    div = pwd / esl 
    return {'pwd': pwd, 'div': div}

# Define a function to load the windowed diplotype divergence results.
def load_alt_v_den_windows(window_size):
    # Intialize a dictionary.
    pwd_dicc = {'pwd': [], 'div': []}
    # Load the effective sequence lengths.
    esl_df = load_windows('arcs_masked_no_aa', 'variant', window_size=window_size)
    # For all chromosomes.
    for chrom in range(1, 23):
        # Load the matrix.
        pw_mat = np.loadtxt(f'../muc19_results/arcs_masked_no_aa/den_v_alt_pw_diffs_chr{chrom}_{window_size}kb.txt.gz')
        # Extract the effective sequence lengths.
        esl = esl_df[esl_df['CHR'] == chrom]['DEN-ALT'].values
        # Update the results.
        pwd_dicc['pwd'].append(pw_mat)
        pwd_dicc['div'].append(pw_mat / esl)
    # Concatenate all the windows.
    pwd_dicc['pwd'] = np.concatenate(pwd_dicc['pwd'], axis=0)
    pwd_dicc['div'] = np.concatenate(pwd_dicc['div'], axis=0)
    return pwd_dicc

# Define a function to calculate sequence divergence between africans and the denisovan.
def calc_afr_v_den_diffs(gt, window_size):
    # Load in the meta information as a pandas dataframe.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the super population indicies.
    afr_idx = meta_df[meta_df['SUPERPOP'] == 'AFR'].index.values
    # Load the effective sequence length.
    esl = int(load_region_esl(prefix='tgp_den_masked_no_aa', window_size=window_size))
    # Determine the alternative allele frequencies per archaic.
    den_freq = calc_ind_alt_freqs(gt.take([2347], axis=1))
    afr_freq = calc_alt_freqs(gt.take(afr_idx, axis=1))
    # Compute the number of pairwise differences.
    pwd = np.nansum(pwd_per_site(den_freq, afr_freq))
    # Compute the sequence divergence.
    div = pwd / esl 
    return {'pwd': pwd, 'div': div}

# Define a function to load the windowed diplotype divergence results.
def load_afr_v_den_windows(window_size):
    # Intialize a dictionary.
    pwd_dicc = {'pwd': [], 'div': []}
    # Load the effective sequence lengths.
    esl_df = load_windows('tgp_den_masked_no_aa', 'variant', window_size)
    # For all chromosomes.
    for chrom in range(1, 23):
        # Load the matrix.
        pw_mat = np.loadtxt(f'../muc19_results/tgp_den_masked_no_aa/afr_den_avg_pw_diffs_chr{chrom}_{window_size}kb.txt.gz')
        # Extract the effective sequence lengths.
        esl = esl_df[esl_df['CHR'] == chrom]['DEN'].values
        # Update the results.
        pwd_dicc['pwd'].append(pw_mat)
        pwd_dicc['div'].append(pw_mat / esl)
    # Concatenate all the windows.
    pwd_dicc['pwd'] = np.concatenate(pwd_dicc['pwd'], axis=0)
    pwd_dicc['div'] = np.concatenate(pwd_dicc['div'], axis=0)
    return pwd_dicc

# Define a function to compile the observed values.
def compile_alt_afr_v_den_diffs(arc_gt, tgp_gt, window_size):
    # Intialize a dictionary.
    diff_dicc = {}
    # Load the altai results.
    diff_dicc['ALT'] = calc_alt_v_den_diffs(arc_gt, window_size=window_size)
    diff_dicc['AFR'] = calc_afr_v_den_diffs(tgp_gt, window_size=window_size)
    return diff_dicc

# Define a function to compile the windowed values.
def compile_alt_afr_v_den_windows(window_size):
    # Intialize a dictionary.
    diff_dicc = {}
    # Load the altai results.
    diff_dicc['ALT'] = load_alt_v_den_windows(window_size=window_size)
    diff_dicc['AFR'] = load_afr_v_den_windows(window_size=window_size)
    return diff_dicc

# Define a function to compile the african/altai nean v denisovan results.
def compile_alt_afr_v_den_div_summary(obs_dicc, wind_dicc, window_size):
    # Intialize labels.
    lab_dicc = {
        'ALT': 'Altai Nean.', 'AFR': 'AFR'
    }
    # Intialize dictionaries to store the results.
    df_dicc = {
        'comp': [],
        'obs_pwd': [], 'obs_div': [],
        'wind_m': [], 'wind_s': [],
        'wind_se': [], 'wind_ci': [],
        'wind_p': [],
    }
    # Load the nonoverlapping window indicies of comparable effective sequence length.
    idx_dicc = {
        'ALT': load_esl_qc_windows_idx('arcs_masked_no_aa', window_size=window_size),
        'AFR': load_esl_qc_windows_idx('tgp_den_masked_no_aa', window_size=window_size),
    }
    # For every comparison.
    for comp, lab in lab_dicc.items():
        # Grab the window distribution.
        dist = wind_dicc[comp]['div'][idx_dicc[comp]]
        # Determine the mean, standard deviation, sem, and 95% CIs for non-overlapping windows.
        wind_m = np.nanmean(dist)
        wind_s = np.nanstd(dist)
        wind_se, wind_ci = sem_ci_of_mean(dist)
        # Compute the p-value.
        p_val = (np.count_nonzero(obs_dicc[comp]['div'] <= dist) / np.sum(~np.isnan(dist)))
        # Fill the dictionaries.
        df_dicc['comp'].append(lab)
        df_dicc['obs_pwd'].append(obs_dicc[comp]['pwd'])
        df_dicc['obs_div'].append(obs_dicc[comp]['div'])
        df_dicc['wind_m'].append(wind_m)
        df_dicc['wind_s'].append(wind_s)
        df_dicc['wind_se'].append(wind_se)
        df_dicc['wind_ci'].append(wind_ci)
        df_dicc['wind_p'].append(p_val)
    # Convert the dictionaries to dataframes.
    div_df = pd.DataFrame(df_dicc)
    # Rename all the columns to look pretty.
    div_df.rename(
        columns={
            'comp': 'Comparison',
            'obs_pwd': f'Focal {window_size}kb Region (Pairwise Diffs.)',
            'obs_div': f'Focal {window_size}kb Region (Seq. Div.)',
            'wind_m': fr'{window_size}kb Nonoverlapping Windows $\left( \mu \right)$',
            'wind_s': fr'{window_size}kb Nonoverlapping Windows $\left( \sigma \right)$',
            'wind_se': fr'{window_size}kb Nonoverlapping Windows $\left( SEM \right)$',
            'wind_ci': fr'{window_size}kb Nonoverlapping Windows $\left( \pm CI_{{95\%}} \right)$',
            'wind_p': r'$P-value$',
        }, inplace=True,
    )
    # Export the the dataframe.
    div_df.to_csv(f'./dataframes/altai_nean_afr_v_denisovan_seq_div_{window_size}kb.csv.gz', index=False)
    return div_df

# Define a function to plot the sequence divergence results between the altai neanderthal/africans and the denisovan.
def plot_alt_afr_v_den_seq_div(obs_dicc, wind_dicc, window_size):
    # Intialize labels.
    lab_dicc = {
        'ALT': 'Altai Nean.', 'AFR': 'AFR'
    }
    # Load the nonoverlapping window indicies of comparable effective sequence length.
    idx_dicc = {
        'ALT': load_esl_qc_windows_idx('arcs_masked_no_aa', window_size=window_size),
        'AFR': load_esl_qc_windows_idx('tgp_den_masked_no_aa', window_size=window_size),
    }
    # Intialize the maximum bin value.
    if window_size == 72:
        bin_max = 0.005
    else:
        bin_max = 0.003
    # Determine the significance threshold.
    sig_thresh = 0.05 
    # Intialize a dictionary.
    sig_dicc = {}
    # Intialize the first letter ASCII value.
    c_letter = 65
    # Update the standard matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 2, figsize=(7, 3.5), dpi=300,
        facecolor='white',
        sharex=False, sharey=True,
        constrained_layout=True,
    )
    # For each subplot.
    for i, key in enumerate(lab_dicc):
        # Extract the observed value and windowed distribution.
        obs = obs_dicc[key]['div']
        dist = wind_dicc[key]['div'][idx_dicc[key]]
        # Compute the p-value.
        pval = (np.count_nonzero(obs <= dist) / np.sum(~np.isnan(dist)))
        # If the p-value is significant.
        if pval < sig_thresh:
            # Update the dictionary.
            sig_dicc[i] = (obs, '*', 25, 0.49)
        # Else, it is not significant.
        else:
            # Update the dictionary.
            sig_dicc[i] = (obs, r'$ns$', 10, 0.5225)
        # Plot the distribution.
        axes[i].hist(
            dist, bins=np.linspace(0, bin_max, 100),
            histtype='stepfilled', color='gray',
        )
        # Generate consistent tick labels.
        custom_ticks = np.linspace(0, bin_max, 5)
        axes[i].set_xticks(custom_ticks)
        axes[i].set_xticklabels([f'{x:5f}'.rstrip('0').rstrip('.') for x in custom_ticks])
        # Plot the panel label.
        axes[i].set_title(r'$\bf{'+chr(c_letter+i)+'.}$', loc='left')
        # Add axes labels.
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlabel(f'{lab_dicc[key]} v Denisovan Seq. Div. per {window_size}kb Window')
        # Force labels to appear on all subplots.
        axes[i].xaxis.set_tick_params(which='both', labelbottom=True)
        axes[i].yaxis.set_tick_params(which='both', labelleft=True)
    # For every subplot.
    for i, sig_tup in sig_dicc.items():
        # Unpack.
        obs, label, label_size, constant = sig_tup
        # Annotate the observed value.
        axes[i].annotate(
            '', xy=(obs, 0),
            xytext=(obs, axes[i].get_ylim()[1] * 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
        )
        axes[i].text(
            obs, axes[i].get_ylim()[1] * constant, label,
            ha='center', va='center', fontsize=label_size,
        )
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/altai_nean_afr_v_denisovan_seq_div_{window_size}kb_windows.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/altai_nean_afr_v_denisovan_seq_div_{window_size}kb_windows.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/altai_nean_afr_v_denisovan_seq_div_{window_size}kb_windows.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a PCA function.
def pca(genotype_matrix):
    # Intialize an empty C matrix.
    C_mat = genotype_matrix.T
    # Calculate the column mean vector.
    mu_vec = np.nanmean(C_mat, axis=0)
    # Calculate the column allele frequency vector.
    p_vec = mu_vec / 2
    # Determine the indicies of sites with more than one mutation is observed.
    good_idx = np.where(~((p_vec <= 1 / (C_mat.shape[0] * 2)) | (p_vec == 1)))[0]
    # Calculate the standard deviation vector.
    std_vec = np.sqrt(((C_mat.shape[0] * 2) * p_vec[good_idx] * (1 - p_vec[good_idx])))
    # Convert the C matrix to a zero centered M matrix.
    M_mat = (C_mat[:, good_idx] - mu_vec[good_idx]) / std_vec
    # Set np.nan's to zero.
    M_mat[np.isnan(M_mat)] = 0
    # Compute the covariance matrix.
    X_mat = np.dot(M_mat, M_mat.T) / (M_mat.shape[0] - 1)
    # Compute the eigen -values and -vectors.
    eig_vals, eig_vecs = np.linalg.eigh(X_mat)
    # Sort the -values and -vectors.
    eig_idx = np.argsort(eig_vals)[::-1]
    eig_vals = eig_vals[eig_idx]
    eig_vecs = eig_vecs[:, eig_idx]
    return eig_vals, eig_vecs

# Define a function to plot the TGP and TGP + archaics PCA.
def plot_tgp_arcs_local_pca(tgp_gt, tgp_arcs_gt, region):
    # Intialize the meta data.
    tgp_df, _ = annotate_hap_identities(70, 100, 'tgp')
    spops = tgp_df['SUPERPOP'].values
    tgp_idx = np.arange(tgp_df.shape[0])
    # If this is the 72kb region.
    if region == 72:
        # Create a mask for the recombinant haplotypes.
        is_rec_hap = tgp_df['N_REC_HAPS_72KB'].values == 1
        # Intialize masks for the number of denisovan like haplotypes.
        is_n_haps = {n: (tgp_df['N_DEN_HAPS_72KB'].values == n) & (~is_rec_hap) for n in range(3)}
        is_n_haps[3] = is_rec_hap
        # Intialize dictionaries.
        color_dicc = {
            'AFR': 'grey', 'AMR': '#009E73',
            'SAS': '#CC79A7', 'EAS': '#0072B2', 'EUR': '#D55E00',
            'Denisovan': '#56B4E9', 'Altai Nean.': '#000000',
            'Chagyrskaya Nean.': '#E69F00', 'Vindija Nean.': '#F0E442',
        }
        marker_dicc = {
            0: 'o', 1: '^', 2: 's', 3: 'd',
        }
        arc_dicc = {
            'Altai Nean.': tgp_df.shape[0], 'Chagyrskaya Nean.': tgp_df.shape[0] + 1,
            'Vindija Nean.': tgp_df.shape[0] + 2, 'Denisovan': tgp_df.shape[0] + 3,
        }
        # Intialize the legend.
        legend_handles = [Patch(color=color, label=spop) for spop, color in color_dicc.items()]
        legend_handles.extend([
            Line2D([0], [0], linestyle='none', marker='X', markersize=7.5, markerfacecolor='none', markeredgecolor='#000000', label='Archaic Ind.'),
            Line2D([0], [0], linestyle='none', marker='o', markersize=7.5, markerfacecolor='none', markeredgecolor='#000000', label='0 Introgressed Haps.'),
            Line2D([0], [0], linestyle='none', marker='^', markersize=7.5, markerfacecolor='none', markeredgecolor='#000000', label='1 Introgressed Hap.'),
            Line2D([0], [0], linestyle='none', marker='d', markersize=7.5, markerfacecolor='none', markeredgecolor='#000000', label='1 Recombinant Hap.'),
            Line2D([0], [0], linestyle='none', marker='s', markersize=7.5, markerfacecolor='none', markeredgecolor='#000000', label='2 Introgressed Haps.'),
        ])
    else:
        # Intialize masks for the number of archaic like haplotypes.
        is_n_haps = {n: tgp_df['N_ARC_HAPS_742KB'].values == n for n in range(3)}
        # Intialize dictionaries.
        color_dicc = {
            'AFR': 'grey', 'AMR': '#009E73',
            'SAS': '#CC79A7', 'EAS': '#0072B2', 'EUR': '#D55E00',
            'Denisovan': '#56B4E9', 'Altai Nean.': '#000000',
            'Chagyrskaya Nean.': '#E69F00', 'Vindija Nean.': '#F0E442',
        }
        marker_dicc = {
            0: 'o', 1: '^', 2: 's',
        }
        arc_dicc = {
            'Altai Nean.': tgp_df.shape[0], 'Chagyrskaya Nean.': tgp_df.shape[0] + 1,
            'Vindija Nean.': tgp_df.shape[0] + 2, 'Denisovan': tgp_df.shape[0] + 3,
        }
        # Intialize the legend.
        legend_handles = [Patch(color=color, label=spop) for spop, color in color_dicc.items()]
        legend_handles.extend([
            Line2D([0], [0], linestyle='none', marker='X', markersize=7.5, markerfacecolor='none', markeredgecolor='#000000', label='Archaic Ind.'),
            Line2D([0], [0], linestyle='none', marker='o', markersize=7.5, markerfacecolor='none', markeredgecolor='#000000', label='0 Introgressed Haps.'),
            Line2D([0], [0], linestyle='none', marker='^', markersize=7.5, markerfacecolor='none', markeredgecolor='#000000', label='1 Introgressed Hap.'),
            Line2D([0], [0], linestyle='none', marker='s', markersize=7.5, markerfacecolor='none', markeredgecolor='#000000', label='2 Introgressed Haps.'),
        ])
    # Convert the genotype matrix to an alternative allele count matrix.
    tgp_aac = tgp_gt.to_n_alt(fill=np.nan, dtype=np.float64)
    tgp_arcs_aac = tgp_arcs_gt.to_n_alt(fill=np.nan, dtype=np.float64)
    # Compute the eigen -values and -vectors.
    tgp_e_vals, tgp_e_vecs = pca(tgp_aac)
    tgp_arcs_e_vals, tgp_arcs_e_vecs = pca(tgp_arcs_aac)
    # Compute the variance explained.
    tgp_var_expl = (tgp_e_vals / tgp_e_vals.sum()) * 100
    tgp_arcs_var_expl = (tgp_arcs_e_vals / tgp_arcs_e_vals.sum()) * 100
    # Update the standard matplotlib settings.
    plt.rcParams.update({
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'ytick.labelsize': 8,
        'xtick.labelsize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 2, figsize=(7, 3.5), dpi=300,
        facecolor='white',
        sharex=False, sharey=False,
        constrained_layout=True,
    )
    # For every haplotype group.
    for hap_key, hap_mask in is_n_haps.items():
        # For every individual.
        for ind_idx, color_key in zip(tgp_idx[hap_mask], spops[hap_mask]):
            # Plot the PC1 v PC2 results.
            axes[0].scatter(
                tgp_e_vecs[ind_idx, 0], tgp_e_vecs[ind_idx, 1],
                marker=marker_dicc[hap_key], facecolor='none', edgecolor=color_dicc[color_key],
                alpha=0.25, linewidth=1, s=25,
            )
            axes[1].scatter(
                tgp_arcs_e_vecs[ind_idx, 0], tgp_arcs_e_vecs[ind_idx, 1],
                marker=marker_dicc[hap_key], facecolor='none', edgecolor=color_dicc[color_key],
                alpha=0.25, linewidth=1, s=25,
            )
    # For every archaic.
    for color_key, ind_idx in arc_dicc.items():
        # Plot the PC1 v PC2 results.
        axes[1].scatter(
            tgp_arcs_e_vecs[ind_idx, 0], tgp_arcs_e_vecs[ind_idx, 1],
            marker='X', facecolor='none', edgecolor=color_dicc[color_key],
            alpha=0.75, linewidth=1, s=25,
        )   
    # Label the axes.
    axes[0].set_xlabel(f'PC1 ({tgp_var_expl[0]:.2f}%)')
    axes[0].set_ylabel(f'PC2 ({tgp_var_expl[1]:.2f}%)')
    axes[1].set_xlabel(f'PC1 ({tgp_arcs_var_expl[0]:.2f}%)')
    axes[1].set_ylabel(f'PC2 ({tgp_arcs_var_expl[1]:.2f}%)')
    # Add labels to the subplots.
    axes[0].set_title(r'$\bf{A.}$', loc='left')
    axes[1].set_title(r'$\bf{B.}$', loc='left')
    # Add the legend.
    title_font = fm.FontProperties(weight='bold', size=10)
    fig.legend(
        handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5375, -0.005),
        ncol=5, frameon=True, fancybox=True, shadow=True, prop={'size': 8},
        title=f'{region}kb Region', title_fontproperties=title_font,
    )
    # If this is the 72kb region.
    if region == 72:
        # Intialize the clusters.
        tgp_ellipses = [
            ((-0.007, 0.0065), 0.0325, 0.09, r'$\delta$'),
            ((0.035, 0.00375), 0.01, 0.045, r'$\epsilon$'),
            ((0.051, -0.0025), 0.02, 0.055, r'$\zeta$'),
            ((0.11, -0.0015), 0.0175, 0.0175, r'$\eta$'),
        ]
        # For every cluster.
        for center, width, height, numeral in tgp_ellipses:
            # Outline the cluster.
            axes[0].add_patch(Ellipse(
                center, width=width, height=height, linestyle='dashed',
                edgecolor='lightgrey', facecolor='none', linewidth=1, 
            ))
            # Annotate the cluster.
            label_position = (center[0], center[1] + height/2 + 0.001)
            axes[0].annotate(numeral, label_position, color='black', fontsize=10, ha='center', va='bottom')
        # Intialize the clusters.
        tgp_arc_ellipses = [
            ((-0.0065, 0.007), 0.035, 0.085, r'$\delta$'),
            ((0.0355, 0.0031), 0.01, 0.045, r'$\epsilon$'),
            ((0.051, -0.0025), 0.0185, 0.055, r'$\zeta$'),
            ((0.1075, -0.0035), 0.0175, 0.0175, r'$\eta$'),
        ]
        # For every cluster.
        for center, width, height, numeral in tgp_arc_ellipses:
            # Outline the cluster.
            axes[1].add_patch(Ellipse(
                center, width=width, height=height, linestyle='dashed',
                edgecolor='lightgrey', facecolor='none', linewidth=1,
            ))
            # Annotate the cluster.
            label_position = (center[0], center[1] + height/2 + 0.001)
            axes[1].annotate(numeral, label_position, color='black', fontsize=10, ha='center', va='bottom')
    # Export the plot.
    plt.savefig(
        f'./supp_figures/png/tgp_arcs_pca_{region}kb_region.png', format='png',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/svg/tgp_arcs_pca_{region}kb_region.svg', format='svg',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    plt.savefig(
        f'./supp_figures/pdf/tgp_arcs_pca_{region}kb_region.pdf', format='pdf',
        facecolor='white', bbox_inches='tight', dpi=500,
    )
    # Show the plot.
    plt.show()
    return

# Define a function to find fixed differences.
def find_tgp_denisovan_fixed_diffs(tgp_den_gt):
    # Compute allele frequencies.
    tgp_aaf = calc_alt_freqs(tgp_den_gt.take(np.arange(2347), axis=1))
    arc_aaf = calc_ind_alt_freqs(tgp_den_gt.take([2347], axis=1))
    # Determine where the tgp is fixed.
    is_tgp_fxd_ref = tgp_aaf == 0
    is_tgp_fxd_alt = tgp_aaf == 1
    # Determine the archaic genotypes.
    is_arc_hom_ref = arc_aaf == 0
    is_arc_hom_alt = arc_aaf == 1
    is_arc_het = arc_aaf == 0.5
    # Find the fixed and half differences.
    is_arc_fxd_tgp_ref = is_tgp_fxd_ref & is_arc_hom_alt
    is_arc_hlf_tgp_ref = is_tgp_fxd_ref & is_arc_het
    is_arc_fxd_tgp_alt = is_tgp_fxd_alt & is_arc_hom_ref
    is_arc_hlf_tgp_alt = is_tgp_fxd_alt & is_arc_het
    is_arc_fxd = is_arc_fxd_tgp_ref | is_arc_fxd_tgp_alt
    is_arc_hlf = is_arc_hlf_tgp_ref | is_arc_hlf_tgp_alt
    # Print a summary.
    print(f'Fixed Differences: {is_arc_fxd.sum()}')
    print(f'Half Differences: {is_arc_hlf.sum()}')
    return

# Define a function to show the san sequence divergence results from the Denisovan.
def show_san_denisovan_seq_div_72kb():
    # Load the meta information
    sgdp_df = pd.read_csv(
        '../meta_data/sgdp.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the indicies.
    san_idx = np.where(sgdp_df['IND'].values == 'S_Khomani_San-2.DG')[0][0]
    # Load the genotype matrix.
    sgdp_den_gt_72kb, _ = load_hap_region('sgdp_den_masked_no_aa', 12, 40759001, 40831000)
    # Load the effective sequence length.
    esl = int(load_region_esl('sgdp_den_masked_no_aa', 72))
    # Extract the san haplotypes.
    san_hap1 = sgdp_den_gt_72kb[:, san_idx, 0] # The first haplotype carries the missense mutations.
    # Calculate the denisovan's allele frequencies.
    den_aaf = calc_ind_alt_freqs(sgdp_den_gt_72kb.take([278], axis=1))
    # Compute the numbper of pairwise differences.
    san_hap1_pwd = np.nansum(pwd_per_site(san_hap1, den_aaf))
    # Print a summary.
    print('#***# S_Khomani_San-2.DG Hap. 1 v Denisovan #***#')
    print(f'Effective Sequence Length: {esl}')
    print(f'Pairwise Differences: {san_hap1_pwd}')
    print(f'Sequence Divergence: {round(san_hap1_pwd / esl, 6)}')
    return

# Define a function to annotate haplotype identities at the 72kb region.
def annotate_tgp_hap_identities_72kb(den_thresh, rec_thresh):
    # Intialize the meta information.
    meta_data = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the inidcies.
    all_inds = meta_data['IND'].values
    # Load the sequence divergence results.
    div_72kb = pd.read_csv('./dataframes/tgp_haplotype_archaic_diplotype_divergence_72kb.csv.gz')
    # Extract the denisovan divergence results.
    den_div_72kb = div_72kb[div_72kb['Archaic'] == 'Denisovan']
    den_hap1_72kb = den_div_72kb['Focal 72kb Region (Pairwise Diffs. Hap. 1)'].values
    den_hap2_72kb = den_div_72kb['Focal 72kb Region (Pairwise Diffs. Hap. 2)'].values
    # Intialize masks for the unphased haplotype identities at the 72kb region.
    # Note: this threshold was chosen from the S-curves.
    is_den_hap1_72kb = np.isin(all_inds, den_div_72kb[den_hap1_72kb < den_thresh]['Individual'].values)
    is_den_hap2_72kb = np.isin(all_inds, den_div_72kb[den_hap2_72kb < den_thresh]['Individual'].values)
    is_rec_hap1_72kb = np.isin(
        all_inds, den_div_72kb[(den_hap1_72kb > den_thresh) & (den_hap1_72kb < rec_thresh)]['Individual'].values,
    )
    is_rec_hap2_72kb = np.isin(
        all_inds, den_div_72kb[(den_hap2_72kb > den_thresh) & (den_hap2_72kb < rec_thresh)]['Individual'].values,
    )
    is_hum_hap1_72kb = ~(is_den_hap1_72kb | is_rec_hap1_72kb)
    is_hum_hap2_72kb = ~(is_den_hap2_72kb | is_rec_hap2_72kb)
    # Update the dataframe.
    meta_data['IS_HAP1_DEN_72KB'] = is_den_hap1_72kb
    meta_data['IS_HAP2_DEN_72KB'] = is_den_hap2_72kb
    meta_data['IS_HAP1_HUM_72KB'] = is_hum_hap1_72kb
    meta_data['IS_HAP2_HUM_72KB'] = is_hum_hap2_72kb
    meta_data['IS_HAP1_REC_72KB'] = is_rec_hap1_72kb
    meta_data['IS_HAP2_REC_72KB'] = is_rec_hap2_72kb
    meta_data['N_DEN_HAPS_72KB'] = is_den_hap1_72kb.astype(int) + is_den_hap2_72kb.astype(int)
    meta_data['N_HUM_HAPS_72KB'] = is_hum_hap1_72kb.astype(int) + is_hum_hap2_72kb.astype(int)
    meta_data['N_REC_HAPS_72KB'] = is_rec_hap1_72kb.astype(int) + is_rec_hap2_72kb.astype(int)
    return meta_data, {key: np.array(value) for key, value in meta_data.to_dict(orient='list').items()}