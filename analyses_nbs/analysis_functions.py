# Import packages.
import allel
from functools import reduce
import gzip
import itertools
import matplotlib
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import numcodecs
import pandas as pd
from scipy import stats as stats
from scipy.stats import norm
import sys
import zarr

# Intialize the matplolib styling.
plt.rcParams.update({
    'figure.constrained_layout.use': True,
    'figure.dpi': 300,
    'figure.facecolor': 'white',
    'figure.titlesize': 10,
    'font.family': 'serif',
    'font.size': 8,
    'axes.titlesize': 8,
    'axes.labelsize': 8,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'legend.frameon': False,
    'legend.fontsize': 8,
    'lines.markersize': 3,
    'lines.linewidth': 0.8,
    'savefig.dpi': 500,
    'savefig.format': 'svg',
    'ytick.labelsize': 6,
    'xtick.labelsize': 6,
})

# Intialize the numpy warning preferences.
np.seterr(invalid='ignore')

# Intialize the pandas preferences.
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)



##########################
### GENERAL FUNCTIONS ###
##########################

# Define a function to load genotyope and positions arrays.
def load_gt_pos(prefix, chrom):
    # Intialize the file path.
    path = '../vcf_data/zarr_data/{0}.zarr'.format(prefix)
    # Load the zarr array.
    zarr_array = zarr.open_group(path, mode='r')
    # Extract the genotype callset.
    callset = zarr_array['{0}/calldata/GT'.format(chrom)]
    # Convert the genotype callset to an array.
    gt = allel.GenotypeArray(callset)
    # Load the positions.
    pos = allel.SortedIndex(zarr_array['{0}/variants/POS'.format(chrom)])
    return gt, pos

# Define a function to load genotyope and positions arrays.
def load_callset_pos(prefix, chrom):
    # Intialize the file path.
    path = '../vcf_data/zarr_data/{0}_chr{1}.zarr'.format(prefix, chrom)
    # Load the zarr array.
    zarr_array = zarr.open_group(path, mode='r')
    # Extract the genotype callset.
    callset = zarr_array['{0}/calldata/GT'.format(chrom)]
    # Load the positions.
    pos = allel.SortedIndex(zarr_array['{0}/variants/POS'.format(chrom)])
    return callset, pos

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

# Define a function to calculate alternative allele frequencies for the archaics
def calc_arc_alt_freqs(arc_gt):
    # Compute the frequency for each site.
    raw_freqs = np.nansum(arc_gt, axis=2).flatten() / 2
    # Set missing data to Nan
    alt_freqs = np.where(raw_freqs == -1, np.nan, raw_freqs)
    return alt_freqs

# Define a function to load the QC'ed windows.
def load_windows(tgp_or_arc, var_or_invar, window_size=72):
    # Load the data frame.
    wind_df = pd.read_csv(f'../windowing/{tgp_or_arc}/{window_size}kb_nonoverlapping_{var_or_invar}_windows.csv.gz')
    return wind_df

# Define a function to load the window indicies that passed effective sequence length QC.
def load_esl_qc_windows_idx(tgp_or_arc, window_size=72):
    # Load the QC'ed windows indicies.
    qc_winds_idx = np.loadtxt(
        f'../windowing/{tgp_or_arc}/{window_size}kb_esl_qced_nonoverlapping_variant_windows.csv.gz',
        delimiter=',', dtype=int,
    )
    return qc_winds_idx

# Define a function to load the haplotype effective sequence lengths.
def load_hap_esl_mat(prefix='arc_and_tgp', window_size=72):
    # Load the effective sequence length matrix.
    hap_esl_mat = np.loadtxt(
            f'../meta_data/{prefix}_{window_size}kb_eff_seq_len.csv',
            delimiter=',', dtype=int,
        )
    return hap_esl_mat

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

# Define a function to polarize a genotype matrix.
def polarize_gt(gt):
    # Extract and determine the ancestral sequence.
    anc_seq = np.mean(gt[:, -1, :], axis=1)
    # Intialize an empty genotype matrix.
    p_gt = np.empty_like(gt[:, :-1, :])
    # For every sample...
    for samp in range(gt[:, :-1, :].shape[1]):
        # Extract the haplotypes.
        hap_1 = gt[:, samp, 0]
        hap_2 = gt[:, samp, 1]
        # Polarize.
        p_hap_1 = np.where(anc_seq == 1, np.abs(hap_1 - 1), hap_1)
        p_hap_2 = np.where(anc_seq == 1, np.abs(hap_2 - 1), hap_2)
        # Build the polarized genotype matrix.
        p_gt[:, samp, 0] = np.where(p_hap_1 == 2, -1, p_hap_1)
        p_gt[:, samp, 1] = np.where(p_hap_2 == 2, -1, p_hap_2)
    # Intialize the polarized genotype matrix into sci-kit allele.
    polarized_gt = allel.GenotypeArray(p_gt)
    return polarized_gt


##################
### WINDOWING ###
#################

# Define a function to visualize distributions with respect to each archaic,
# including an observed value.
def plot_arc_dist(arc_dicc, winds_df, title, obs_label, x_label):
    # Intialize a list of populations to plot.
    pop_list = ['DEN', 'ALT', 'CHA', 'VIN']
    # Intialize a row and columns list.
    rows = [0, 0, 1, 1]
    cols = [0, 1, 0, 1]
    # Intialize the figure and axes.
    fig, axes = plt.subplots(
        2, 2, dpi=300,
        sharex=True, sharey=True,
    )
    # For every population...
    for idx in range(len(pop_list)):
        # Extract the population.
        pop = pop_list[idx]
        # Plot the distribution.
        axes[rows[idx], cols[idx]].hist(
            winds_df[pop].values,
            bins='fd',
            histtype='stepfilled', color='tab:blue',
        )
        # Add a subplot title.
        axes[rows[idx], cols[idx]].set_title(pop)
        # If this is the last plot...
        if (idx == 3):
            # Plot the observed value with a label.
            axes[rows[idx], cols[idx]].axvline(
                arc_dicc[pop], 0, 1,
                color='black', linestyle='dashed',
                label=obs_label,
            )
        # Else...
        else:
            # Plot the observed value.
            axes[rows[idx], cols[idx]].axvline(
                arc_dicc[pop], 0, 1,
                color='black', linestyle='dashed',
            )
    # Add a figure legend.
    fig.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    # Title the plot.
    fig.suptitle(title)
    # Label the super-axes.
    fig.supylabel('Frequency', size=10)
    fig.supxlabel(x_label, size=10)
    # Show the plot!
    plt.show()
    return


############
### PBS ###
###########

# Define a function to plot the distribution of Indigenous American ancestry in MXL.
def plot_mxl_nat_anc_dist(anc_array, cut_off):
    # Intialize the figure.
    fig = plt.figure(
        figsize=(6, 3), dpi=300.0,
    )
    # Intialize the axes.
    ax = fig.add_subplot(111)
    # Plot the distribution.
    ax.hist(
        anc_array, bins=np.arange(0, 1.1, 0.1),
        histtype='stepfilled', color='tab:blue',
    )
    # Plot the observed value.
    ax.axvline(
        cut_off, 0, 1,
        color='black', linestyle='dashed',
    )
    # Add a subplot title.
    ax.set_title('MXL')
    # Label the axes.
    ax.set_ylabel('Frequency')
    ax.set_xlabel(f'Percent Indigenous American Ancetry')
    # Show the plot!
    plt.show()
    return

# Define a function to calculate PBS.
def calc_pbs(gt, pop_a, pop_b, pop_c):
    # Determine allele counts.
    a_ac = gt.take(pop_a, axis=1).count_alleles()
    b_ac = gt.take(pop_b, axis=1).count_alleles()
    c_ac = gt.take(pop_c, axis=1).count_alleles()
    # Calculate the numerator and denominator for Hudson's Fst estimator.
    a_b_num, a_b_den = allel.hudson_fst(a_ac, b_ac)
    a_c_num, a_c_den = allel.hudson_fst(a_ac, c_ac)
    c_b_num, c_b_den = allel.hudson_fst(c_ac, b_ac)
    # Calculate Fst.
    a_b_raw_fst = a_b_num / a_b_den
    a_c_raw_fst = a_c_num / a_c_den
    c_b_raw_fst = c_b_num / c_b_den
    # Correct for Fst values of 1 that will lead to inf.
    a_b_fst_raw = np.where(a_b_raw_fst == 1, 0.99999, a_b_raw_fst)
    a_c_fst_raw = np.where(a_c_raw_fst == 1, 0.99999, a_c_raw_fst)
    c_b_fst_raw = np.where(c_b_raw_fst == 1, 0.99999, c_b_raw_fst)
    # Correct for negative Fst values.
    a_b_fst = np.where(a_b_fst_raw < 0, 0, a_b_fst_raw)
    a_c_fst = np.where(a_c_fst_raw < 0, 0, a_c_fst_raw)
    c_b_fst = np.where(c_b_fst_raw < 0, 0, c_b_fst_raw)
    # Calculate PBS.
    pbs = (
        ((np.log(1.0 - a_b_fst) * -1.0) +\
         (np.log(1.0 - a_c_fst) * -1.0) -\
         (np.log(1.0 - c_b_fst) * -1.0)) / float(2)
    )
    return pbs, np.nanmean(pbs)

# Define a function to load the per-SNP PBS values.
def load_pbs_chromosomes(anc):
    # Intialize an array to store the results.
    all_pbs = np.array([])
    # For all chromosomes...
    for chrom in range(1, 23):
        # Load the pbs data.
        chrom_pbs = np.loadtxt(
            f'../pbs/chroms/mxl_{anc}_chb_ceu_pbs_chr{chrom}.csv.gz', delimiter=',',
        )
        # Append the results.
        all_pbs = np.append(all_pbs, chrom_pbs)
    return all_pbs

# Define a function to load the Denisovan-specific SNP PBS values.
def load_den_sites_pbs_chromosomes():
    # Intialize an array to store the results.
    den_pbs = np.array([])
    # For all chromosomes...
    for chrom in range(1, 23):
        # Load the pbs data.
        chrom_pbs = np.loadtxt(
            f'../pbs/chroms/mxl_chb_ceu_den_sites_pbs_all_chr{chrom}.csv.gz', delimiter=',',
        )
        # Append the results.
        den_pbs = np.append(den_pbs, chrom_pbs)
    return den_pbs

# Define a function to load the windowed average PBS values.
def load_pbs_windows(window_size=748):
    # Intialize an array to store the results.
    pbs_winds = np.array([])
    # For all chromosomes...
    for chrom in range(1, 23):
        # Load the pbs data.
        chrom_pbs = np.loadtxt(
            f'../pbs/windows/mxl_pop_chb_ceu_pbs_chr{chrom}_{window_size}kb.csv.gz', delimiter=',',
        )
        # Append the results.
        pbs_winds = np.append(pbs_winds, chrom_pbs)
    return pbs_winds

# Define a function to exctract a archaic specific snps.
def extract_arc_der_sites(gt, pos, pop_dicc):
    # Calculate alternative allele frequencies.
    afr_alt_freq = calc_alt_freqs(gt.take(pop_dicc['AFR'], axis=1))
    ooa_alt_freq = calc_alt_freqs(gt.take(pop_dicc['OOA'], axis=1))
    alt_alt_freq = calc_alt_freqs(gt.take(pop_dicc['ALT'], axis=1))
    cha_alt_freq = calc_alt_freqs(gt.take(pop_dicc['CHA'], axis=1))
    vin_alt_freq = calc_alt_freqs(gt.take(pop_dicc['VIN'], axis=1))
    den_alt_freq = calc_alt_freqs(gt.take(pop_dicc['DEN'], axis=1))
    anc_freq = calc_alt_freqs(gt.take([-1], axis=1))
    # Polarize the samples.
    afr_der_freq = np.where(anc_freq == 1, np.abs(afr_alt_freq - 1), afr_alt_freq)
    ooa_der_freq = np.where(anc_freq == 1, np.abs(ooa_alt_freq - 1), ooa_alt_freq)
    alt_der_freq = np.where(anc_freq == 1, np.abs(alt_alt_freq - 1), alt_alt_freq)
    cha_der_freq = np.where(anc_freq == 1, np.abs(cha_alt_freq - 1), cha_alt_freq)
    vin_der_freq = np.where(anc_freq == 1, np.abs(vin_alt_freq - 1), vin_alt_freq)
    den_der_freq = np.where(anc_freq == 1, np.abs(den_alt_freq - 1), den_alt_freq)
    # Intialize conditions.
    c_hum = (afr_der_freq < 0.01) & (ooa_der_freq > 0)
    c_den = (den_der_freq == 1)
    c_not_den = (den_der_freq != 1)
    c_alt = (alt_der_freq == 1)
    c_not_alt = (alt_der_freq != 1)
    c_cha = (cha_der_freq == 1)
    c_not_cha = (cha_der_freq != 1)
    c_vin = (vin_der_freq == 1)
    c_not_vin = (vin_der_freq != 1)
    c_nea = c_alt | c_cha | c_vin
    c_not_nea = c_not_alt & c_not_cha & c_not_vin
    # Determine the archaic specific sites.
    den_only_idx = np.where(c_hum & c_not_nea & c_den)[0]
    nea_only_idx = np.where(c_hum & c_not_den & c_nea)[0]
    # Construct the masks.
    den_mask = np.in1d(np.arange(pos.size), den_only_idx)
    nea_mask = np.in1d(np.arange(pos.size), nea_only_idx)
    # Compile a dictionary with the results.
    arc_sites_dicc = {
        'DEN': {'IDX': den_only_idx, 'MASK': den_mask},
        'NEA': {'IDX': nea_only_idx, 'MASK': nea_mask},
    }
    return arc_sites_dicc

# Define a function to exctract a archaic specific snps.
def extract_arc_anc_sites(gt, pos, pop_dicc):
    # Calculate alternative allele frequencies.
    afr_alt_freq = calc_alt_freqs(gt.take(pop_dicc['AFR'], axis=1))
    ooa_alt_freq = calc_alt_freqs(gt.take(pop_dicc['OOA'], axis=1))
    alt_alt_freq = calc_alt_freqs(gt.take(pop_dicc['ALT'], axis=1))
    cha_alt_freq = calc_alt_freqs(gt.take(pop_dicc['CHA'], axis=1))
    vin_alt_freq = calc_alt_freqs(gt.take(pop_dicc['VIN'], axis=1))
    den_alt_freq = calc_alt_freqs(gt.take(pop_dicc['DEN'], axis=1))
    anc_freq = calc_alt_freqs(gt.take([-1], axis=1))
    # Polarize the samples.
    afr_der_freq = np.where(anc_freq == 1, np.abs(afr_alt_freq - 1), afr_alt_freq)
    ooa_der_freq = np.where(anc_freq == 1, np.abs(ooa_alt_freq - 1), ooa_alt_freq)
    alt_der_freq = np.where(anc_freq == 1, np.abs(alt_alt_freq - 1), alt_alt_freq)
    cha_der_freq = np.where(anc_freq == 1, np.abs(cha_alt_freq - 1), cha_alt_freq)
    vin_der_freq = np.where(anc_freq == 1, np.abs(vin_alt_freq - 1), vin_alt_freq)
    den_der_freq = np.where(anc_freq == 1, np.abs(den_alt_freq - 1), den_alt_freq)
    # Intialize conditions.
    c_hum = ((1 - afr_der_freq) < 0.01) & ((1 - ooa_der_freq) > 0)
    c_den = (den_der_freq == 0)
    c_not_den = (den_der_freq != 0)
    c_alt = (alt_der_freq == 0)
    c_not_alt = (alt_der_freq != 0)
    c_cha = (cha_der_freq == 0)
    c_not_cha = (cha_der_freq != 0)
    c_vin = (vin_der_freq == 0)
    c_not_vin = (vin_der_freq != 0)
    c_nea = c_alt | c_cha | c_vin
    c_not_nea = c_not_alt & c_not_cha & c_not_vin
    # Determine the archaic specific sites.
    den_only_idx = np.where(c_hum & c_not_nea & c_den)[0]
    nea_only_idx = np.where(c_hum & c_not_den & c_nea)[0]
    # Construct the masks.
    den_mask = np.in1d(np.arange(pos.size), den_only_idx)
    nea_mask = np.in1d(np.arange(pos.size), nea_only_idx)
    # Compile a dictionary with the results.
    arc_sites_dicc = {
        'DEN': {'IDX': den_only_idx, 'MASK': den_mask},
        'NEA': {'IDX': nea_only_idx, 'MASK': nea_mask},
    }
    return arc_sites_dicc

# Define a function to plot the windowed PBS results.
def plot_avg_pbs_winds(gt, pop_a, pop_b, pop_c, obs_label, window_size=748):
    # Load the observed, nonoverlapping windows, and window indicies of,
    # comparable effective sequence length.
    _, obs_pbs = calc_pbs(gt=gt, pop_a=pop_a, pop_b=pop_b, pop_c=pop_c)
    wind_pbs = load_pbs_windows(window_size=window_size)
    wind_idx = load_esl_qc_windows_idx('tgp', window_size=window_size)
    # Determine the significance threshold.
    sig_threshold = 0.05
    # Compute the p-value.
    wind_p = (np.count_nonzero(obs_pbs <= wind_pbs[wind_idx])
              / np.sum(~np.isnan(wind_pbs[wind_idx])))
    # Intialize the figure.
    fig = plt.figure(
        figsize=(6, 3), dpi=300.0,
    )
    # Intialize the axes.
    ax = fig.add_subplot(111)
    # Plot the distribution.
    ax.hist(
        wind_pbs[wind_idx], bins='fd',
        histtype='stepfilled', color='tab:blue',
    )
    # Plot the observed value.
    ax.axvline(
        obs_pbs, 0, 1,
        color='black', linestyle='dashed',
    )
    # If the p-value is significant...
    if wind_p < sig_threshold:
        # If the p-value is 0.
        if wind_p == 0:
            # Construct the title.
            title = r'$P-value=$'+'<3.137e-04'+r'$^{*}$'
        # Else...
        else:
            # Construct the title.
            title = r'$P-value=$'+'{:.3e}'.format(wind_p)+r'$^{*}$'
    # Else...
    else:
        # If the p-value is 1.
        if round(wind_p, 3) == 1:
            # Construct the title.
            title = r'$P-value=$'+'>0.99969'+r'$^{ns}$'
        # Else...
        else:
            # Construct the title.
            title = r'$P-value=$'+'{:.3e}'.format(wind_p)+r'$^{ns}$'
    # Add a subplot title.
    ax.set_title(title)
    # Construct the legend
    legend_elements = [
        Line2D([0], [0], color='black', linestyle='dashed', label=obs_label),
    ]
    # Add a figure lgend.
    ax.legend(
        handles=legend_elements, loc='center left',
        bbox_to_anchor=(1.0, 0.5), frameon=False,
    )
    # Label the axes.
    ax.set_ylabel('Frequency')
    ax.set_xlabel(f'PBS(MXL:CHB:CEU) in {window_size}kb Windows')
    # Show the plot!
    plt.show()
    return


##################
### Haplotype ###
#################

# Define a function to count the number of archaic specific derived snps.
def count_arc_der_sites(gt, pop_dicc):
    # Calculate alternative allele frequencies.
    afr_alt_freq = calc_alt_freqs(gt.take(pop_dicc['AFR'], axis=1))
    ooa_alt_freq = calc_alt_freqs(gt.take(pop_dicc['OOA'], axis=1))
    alt_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['ALT'], axis=1))
    cha_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['CHA'], axis=1))
    vin_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['VIN'], axis=1))
    den_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['DEN'], axis=1))
    anc_freq = calc_alt_freqs(gt.take([-1], axis=1))
    # Polarize the samples.
    afr_der_freq = np.where(anc_freq == 1, np.abs(afr_alt_freq - 1), afr_alt_freq)
    ooa_der_freq = np.where(anc_freq == 1, np.abs(ooa_alt_freq - 1), ooa_alt_freq)
    alt_der_freq = np.where(anc_freq == 1, np.abs(alt_alt_freq - 1), alt_alt_freq)
    cha_der_freq = np.where(anc_freq == 1, np.abs(cha_alt_freq - 1), cha_alt_freq)
    vin_der_freq = np.where(anc_freq == 1, np.abs(vin_alt_freq - 1), vin_alt_freq)
    den_der_freq = np.where(anc_freq == 1, np.abs(den_alt_freq - 1), den_alt_freq)
    # Intialize conditions.
    c_hum = (afr_der_freq < 0.01) & (ooa_der_freq > 0)
    c_den = (den_der_freq == 1)
    c_not_den = (den_der_freq != 1)
    c_alt = (alt_der_freq == 1)
    c_not_alt = (alt_der_freq != 1)
    c_cha = (cha_der_freq == 1)
    c_not_cha = (cha_der_freq != 1)
    c_vin = (vin_der_freq == 1)
    c_not_vin = (vin_der_freq != 1)
    c_nea = c_alt | c_cha | c_vin
    c_not_nea = c_not_alt & c_not_cha & c_not_vin
    # Determine the archaic specific sites.
    den_only_idx = np.where(c_hum & c_not_nea & c_den)[0]
    nea_only_idx = np.where(c_hum & c_not_den & c_nea)[0]
    return [den_only_idx.size, nea_only_idx.size]

# Define a function to count the number of archaic specific ancestral snps.
def count_arc_anc_sites(gt, pop_dicc):
    # Calculate alternative allele frequencies.
    afr_alt_freq = calc_alt_freqs(gt.take(pop_dicc['AFR'], axis=1))
    ooa_alt_freq = calc_alt_freqs(gt.take(pop_dicc['OOA'], axis=1))
    alt_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['ALT'], axis=1))
    cha_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['CHA'], axis=1))
    vin_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['VIN'], axis=1))
    den_alt_freq = calc_arc_alt_freqs(gt.take(pop_dicc['DEN'], axis=1))
    anc_freq = calc_alt_freqs(gt.take([-1], axis=1))
    # Polarize the samples.
    afr_der_freq = np.where(anc_freq == 1, np.abs(afr_alt_freq - 1), afr_alt_freq)
    ooa_der_freq = np.where(anc_freq == 1, np.abs(ooa_alt_freq - 1), ooa_alt_freq)
    alt_der_freq = np.where(anc_freq == 1, np.abs(alt_alt_freq - 1), alt_alt_freq)
    cha_der_freq = np.where(anc_freq == 1, np.abs(cha_alt_freq - 1), cha_alt_freq)
    vin_der_freq = np.where(anc_freq == 1, np.abs(vin_alt_freq - 1), vin_alt_freq)
    den_der_freq = np.where(anc_freq == 1, np.abs(den_alt_freq - 1), den_alt_freq)
    # Intialize conditions.
    c_hum = ((1 - afr_der_freq) < 0.01) & ((1 - ooa_der_freq) > 0)
    c_den = (den_der_freq == 0)
    c_not_den = (den_der_freq != 0)
    c_alt = (alt_der_freq == 0)
    c_not_alt = (alt_der_freq != 0)
    c_cha = (cha_der_freq == 0)
    c_not_cha = (cha_der_freq != 0)
    c_vin = (vin_der_freq == 0)
    c_not_vin = (vin_der_freq != 0)
    c_nea = c_alt | c_cha | c_vin
    c_not_nea = c_not_alt & c_not_cha & c_not_vin
    # Determine the archaic specific sites.
    den_only_idx = np.where(c_hum & c_not_nea & c_den)[0]
    nea_only_idx = np.where(c_hum & c_not_den & c_nea)[0]
    return [den_only_idx.size, nea_only_idx.size]

# Define a function to compute archaic snp density for muc19.
def tgp_arc_snp_denisty(gt):
    # Load in the meta information as a pandas dataframe.
    tgp_meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary to store all sample indicies.
    samp_idx_dicc = {
        'ALT': np.array([2347]), 'CHA': np.array([2348]),
        'VIN': np.array([2349]), 'DEN': np.array([2350]),
    }
    # Append the dictionary with sample AFR and OOA indicies.
    samp_idx_dicc['AFR'] = tgp_meta_df[tgp_meta_df['SUPERPOP'] == 'AFR'].index.values
    samp_idx_dicc['OOA'] = tgp_meta_df[tgp_meta_df['SUPERPOP'] != 'AFR'].index.values
    # Compute the number of derived and ancestral archaic snps.
    c_der = count_arc_der_sites(gt=gt, pop_dicc=samp_idx_dicc)
    c_anc = count_arc_anc_sites(gt=gt, pop_dicc=samp_idx_dicc)
    # Fill the dictionary.
    snp_dicc = {
        'DEN': {'DER': c_der[0], 'ANC': c_anc[0]},
        'NEA': {'DER': c_der[1], 'ANC': c_anc[1]},
    }
    return snp_dicc

# Define a function to load the archaic snp density windowed results.
def load_tgp_arc_snp_density_windows(window_size=72):
    # Intialize a dictionary to store the results.
    snp_dicc = {
        'DEN': {'DER': np.array([]), 'ANC': np.array([])},
        'NEA': {'DER': np.array([]), 'ANC': np.array([])},
    }
    # For all chromosomes...
    for chrom in range(1, 23):
        # Load the archaic snp density.
        arc_snps = np.loadtxt(
            f'../arc_snp_density/tgp/windows/chr{chrom}_{window_size}kb.csv.gz',
            delimiter=',',
        )
        # Fill the dictionary.
        snp_dicc['DEN']['DER'] = np.append(snp_dicc['DEN']['DER'], arc_snps[:, 0])
        snp_dicc['NEA']['DER'] = np.append(snp_dicc['NEA']['DER'], arc_snps[:, 1])
        snp_dicc['DEN']['ANC'] = np.append(snp_dicc['DEN']['ANC'], arc_snps[:, 2])
        snp_dicc['NEA']['ANC'] = np.append(snp_dicc['NEA']['ANC'], arc_snps[:, 3])
    return snp_dicc

# Define a function to compile the archaic snp denisty results.
def compile_tgp_arc_snp_denisty_summary(gt, window_size=72, export=True):
    # Intialize a dictionary to store the results.
    df_dicc = {
        'arc': [], 'type': [],
        'muc19_a': [], 'wind_m': [],
        'wind_s': [], 'wind_p': [],
    }
    # Intialize label dictionary.
    label_dicc = {
        'DEN': 'Denisovan',
        'NEA': 'Neanderthal',
        'DER': 'Derived',
        'ANC': 'Ancestral',
    }
    # Load the observed, nonoverlapping windows, and window indicies of,
    # comparable effective sequence length.
    muc19_dicc = tgp_arc_snp_denisty(gt=gt)
    wind_dicc = load_tgp_arc_snp_density_windows(window_size=window_size)
    wind_idx = load_esl_qc_windows_idx('tgp', window_size=window_size)
    # For every archaich...
    for arc in ['DEN', 'NEA']:
        # For every snp type.
        for snp in ['DER', 'ANC']:
            # Determine the mean and standard deviation for non-overlapping windows.
            wind_m = np.nanmean(wind_dicc[arc][snp][wind_idx])
            wind_s = np.nanstd(wind_dicc[arc][snp][wind_idx])
            # Compute the p-values.
            wind_p = (np.count_nonzero(muc19_dicc[arc][snp] <= wind_dicc[arc][snp][wind_idx])
                      / np.sum(~np.isnan(wind_dicc[arc][snp][wind_idx])))
            # Fill the dictionary.
            df_dicc['arc'].append(label_dicc[arc])
            df_dicc['type'].append(label_dicc[snp])
            df_dicc['muc19_a'].append(muc19_dicc[arc][snp])
            df_dicc['wind_m'].append(wind_m)
            df_dicc['wind_s'].append(wind_s)
            df_dicc['wind_p'].append(wind_p)
    # Convert the dictionary to a dataframe.
    snp_df = pd.DataFrame(data=df_dicc)
    # Adjust p-values of 0 by the number of permutations.
    snp_df['wind_p'] = np.where(snp_df['wind_p'] == 0, '<3.242e-05', snp_df['wind_p'])
    snp_df['wind_p'] = np.where(snp_df['wind_p'] == '1.0', '>0.999968', snp_df['wind_p'])
    # Rename all the columns to look pretty.
    snp_df.rename(
        columns={
            'arc': 'Archaic',
            'type': 'SNP Type',
            'muc19_a': r'$MUC19$ (Archaic SNP Denisty)',
            'wind_m': r'Nonoverlapping Windows $\left( \mu \right)$',
            'wind_s': r'Nonoverlapping Windows $\left( \sigma \right)$',
            'wind_p': r'$P-value$',
        }, inplace=True,
    )
    # If export is true...
    if export:
        # Export the dataframe as a csv.
        snp_df.to_csv(f'./dataframes/tgp_archaic_snp_denisty_{window_size}kb.csv', index=False)
    return snp_df

# Define a function to plot archaic snp density.
def plot_tgp_arc_snp_denisty_summary(gt, arc, obs_label, window_size=72, export=True):
    # Load the observed, nonoverlapping windows, and window indicies of,
    # comparable effective sequence length.
    muc19_dicc = tgp_arc_snp_denisty(gt=gt)
    wind_dicc = load_tgp_arc_snp_density_windows(window_size=window_size)
    wind_idx = load_esl_qc_windows_idx('tgp', window_size=window_size)
    # Determine the significance threshold.
    sig_threshold = 0.05 / 2
    # Define a snp types.
    snp_type = ['DER', 'ANC']
    # Intialize label dictionary.
    label_dicc = {
        'DEN': 'Denisovan',
        'NEA': 'Neanderthal',
    }
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 2, figsize=(8, 4), dpi=300,
        sharex=False, sharey=False,
    )
    # For every condition...
    for idx, snp in enumerate(snp_type):
        # Compute the p-values.
        wind_p = (np.count_nonzero(muc19_dicc[arc][snp] <= wind_dicc[arc][snp][wind_idx])
                  / np.sum(~np.isnan(wind_dicc[arc][snp][wind_idx])))
        # Plot the distribution.
        axes[idx].hist(
            wind_dicc[arc][snp][wind_idx],
            bins=np.arange(min(wind_dicc[arc][snp][wind_idx]), max(wind_dicc[arc][snp][wind_idx]) + 1, 1),
            histtype='stepfilled', color='tab:blue',
        )
        # Plot the observed value.
        axes[idx].axvline(
            muc19_dicc[arc][snp], 0, 1,
            color='black', linestyle='dashed',
        )
        # If the p-value is significant after correcting for multiple comparisons...
        if wind_p < sig_threshold:
            # If the p-value is 0.
            if wind_p == 0:
                # Construct the title.
                title = r'$P-value=$'+'<3.242e-05'+r'$^{*}$'
            # Else...
            else:
                # Construct the title.
                title = r'$P-value=$'+'{:.3e}'.format(wind_p)+r'$^{*}$'
        # Else...
        else:
            # If the p-value is 1.
            if round(wind_p, 3) == 1:
                # Construct the title.
                title = r'$P-value=$'+'>0.99968'+r'$^{ns}$'
            # Else...
            else:
                # Construct the title.
                title = r'$P-value=$'+'{:.3e}'.format(wind_p)+r'$^{ns}$'
        # Add a subplot title.
        axes[idx].set_title(title)
    # Construct the legend
    legend_elements = [
        Line2D([0], [0], color='black', linestyle='dashed', label=obs_label),
    ]
    # Add a figure lgend.
    fig.legend(
        handles=legend_elements, loc='center left',
        bbox_to_anchor=(1.0, 0.5), frameon=False,
    )    
    # Label the axes.
    axes[0].set_ylabel('Frequency')
    axes[0].set_xlabel(f'Derived {label_dicc[arc]}-Specifc Sites in {window_size}kb Windows')
    axes[1].set_xlabel(f'Ancestral {label_dicc[arc]}-Specifc Sites in {window_size}kb Windows')
    # If export is set to true...
    if export:
        # Export the plot.
        plt.savefig(
            f'./supp_figures/tgp_{arc.lower()}_snp_denisty_{window_size}kb.png', format='png',
            facecolor='white', bbox_inches='tight', dpi=500,
        )
    # Show the plot!
    plt.show()
    return

# Define a function to compute sequence divergence.
def sequence_divergence(px, py):
    # Calculate sequence divergence.
    seq_div = np.nansum(((px * (1 - py)) + (py * (1 - px))))
    return seq_div

# Define a function to calculate sequence divergence between archaic diplotypes and tgp haplotypes.
def tgp_haplotype_distances(p_gt, window_size=72, export=True):
    # Load the meta data.
    meta_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Load the effective sequence length matrix.
    esl_mat = load_hap_esl_mat(prefix='arc_and_tgp', window_size=window_size)
    # Intialize archaic information dictionaries.
    freq_dicc = {
        'DEN': p_gt.take([2350], axis=1).count_alleles().to_frequencies()[:, 1],
        'ALT': p_gt.take([2347], axis=1).count_alleles().to_frequencies()[:, 1],
        'CHA': p_gt.take([2348], axis=1).count_alleles().to_frequencies()[:, 1],
        'VIN': p_gt.take([2349], axis=1).count_alleles().to_frequencies()[:, 1],
    }
    dist_dicc = {
        'DEN': {'hap_1': np.array([]), 'hap_2': np.array([])},
        'ALT': {'hap_1': np.array([]), 'hap_2': np.array([])},
        'CHA': {'hap_1': np.array([]), 'hap_2': np.array([])},
        'VIN': {'hap_1': np.array([]), 'hap_2': np.array([])},
    }
    esl_dicc = {
        'DEN': esl_mat[1, 3], 'ALT': esl_mat[1, 0],
        'CHA': esl_mat[1, 1], 'VIN': esl_mat[1, 2],
    }
    # Intialize data frame dictionaries.
    abs_dicc = {
        'Sample': meta_df['IND'].values,
        'Population': meta_df['POP'].values,
        'Super Population': meta_df['SUPERPOP'].values,
        'Chr. 1 (DEN)': [], 'Chr. 2 (DEN)': [],
        'Chr. 1 (ALT)': [], 'Chr. 2 (ALT)': [],
        'Chr. 1 (CHA)': [], 'Chr. 2 (CHA)': [],
        'Chr. 1 (VIN)': [], 'Chr. 2 (VIN)': [],
    }
    norm_dicc = {
        'Sample': meta_df['IND'].values,
        'Population': meta_df['POP'].values,
        'Super Population': meta_df['SUPERPOP'].values,
        'Chr. 1 (DEN)': [], 'Chr. 2 (DEN)': [],
        'Chr. 1 (ALT)': [], 'Chr. 2 (ALT)': [],
        'Chr. 1 (CHA)': [], 'Chr. 2 (CHA)': [],
        'Chr. 1 (VIN)': [], 'Chr. 2 (VIN)': [],
    }
    # Extract haplotype arrays.
    tgp_hap_1 = p_gt[:, :-4, 0]
    tgp_hap_2 = p_gt[:, :-4, 1]
    # For every archaic individual...
    for arc in ['DEN', 'ALT', 'CHA', 'VIN']:
        # Determine the number of called sites and update the dictionary.
        n_sites = esl_dicc[arc]
        dist_dicc[arc]['n_sites'] = n_sites
        # For every tgp sample...
        for samp in range(p_gt[:, :-4, :].shape[1]):
            # Extract the sample's haplotypes.
            hap_1 = tgp_hap_1[:, samp]
            hap_2 = tgp_hap_2[:, samp]
            # Compute the distance from the archaic diplotype.
            dist_1 = sequence_divergence(hap_1, freq_dicc[arc])
            dist_2 = sequence_divergence(hap_2, freq_dicc[arc])
            # Update dictionaries.
            dist_dicc[arc]['hap_1'] = np.append(dist_dicc[arc]['hap_1'], dist_1)
            dist_dicc[arc]['hap_2'] = np.append(dist_dicc[arc]['hap_2'], dist_2)
            abs_dicc['Chr. 1 ({0})'.format(arc)].append(dist_1)
            abs_dicc['Chr. 2 ({0})'.format(arc)].append(dist_2)
            norm_dicc['Chr. 1 ({0})'.format(arc)].append(dist_1/n_sites)
            norm_dicc['Chr. 2 ({0})'.format(arc)].append(dist_2/n_sites)
    # Intailize the dataframes.
    abs_df = pd.DataFrame(data=abs_dicc)
    norm_df = pd.DataFrame(data=norm_dicc)
    # If export is true...
    if export:
        # Export the dataframe as a csv.
        norm_df.to_csv(f'./dataframes/tgp_arc_hap_seq_div_summary_{window_size}kb.csv', index=False)
    return abs_df, norm_df, dist_dicc

# Define a function to calculate sequence divergence between archaic diplotypes and papuan haplotypes.
def pap_haplotype_distances(p_gt, window_size=72):
    # Load the meta data.
    meta_df = pd.read_csv(
        '../meta_data/sgdp_pap.txt', sep='\t',
        names=['IDX', 'IND', 'POP', 'SUPERPOP'],
    )
    # Load the effective sequence length matrix.
    esl_mat = load_hap_esl_mat(prefix='sgdp', window_size=window_size)
    # Intialize archaic information dictionaries.
    freq_dicc = {
        'DEN': p_gt.take([281], axis=1).count_alleles().to_frequencies()[:, 1],
        'ALT': p_gt.take([278], axis=1).count_alleles().to_frequencies()[:, 1],
        'CHA': p_gt.take([279], axis=1).count_alleles().to_frequencies()[:, 1],
        'VIN': p_gt.take([280], axis=1).count_alleles().to_frequencies()[:, 1],
    }
    dist_dicc = {
        'DEN': {'hap_1': np.array([]), 'hap_2': np.array([])},
        'ALT': {'hap_1': np.array([]), 'hap_2': np.array([])},
        'CHA': {'hap_1': np.array([]), 'hap_2': np.array([])},
        'VIN': {'hap_1': np.array([]), 'hap_2': np.array([])},
    }
    esl_dicc = {
        'DEN': esl_mat[3], 'ALT': esl_mat[0],
        'CHA': esl_mat[1], 'VIN': esl_mat[2],
    }
    # Intialize data frame dictionaries.
    abs_dicc = {
        'Sample': meta_df['IND'].values,
        'Population': meta_df['POP'].values,
        'Super Population': meta_df['SUPERPOP'].values,
        'Chr. 1 (DEN)': [], 'Chr. 2 (DEN)': [],
        'Chr. 1 (ALT)': [], 'Chr. 2 (ALT)': [],
        'Chr. 1 (CHA)': [], 'Chr. 2 (CHA)': [],
        'Chr. 1 (VIN)': [], 'Chr. 2 (VIN)': [],
    }
    norm_dicc = {
        'Sample': meta_df['IND'].values,
        'Population': meta_df['POP'].values,
        'Super Population': meta_df['SUPERPOP'].values,
        'Chr. 1 (DEN)': [], 'Chr. 2 (DEN)': [],
        'Chr. 1 (ALT)': [], 'Chr. 2 (ALT)': [],
        'Chr. 1 (CHA)': [], 'Chr. 2 (CHA)': [],
        'Chr. 1 (VIN)': [], 'Chr. 2 (VIN)': [],
    }
    # Extract haplotype arrays.
    pap_hap_1 = p_gt[:, :-4, 0]
    pap_hap_2 = p_gt[:, :-4, 1]
    # For every archaic individual...
    for arc in ['DEN', 'ALT', 'CHA', 'VIN']:
        # Determine the number of called sites and update the dictionary.
        n_sites = esl_dicc[arc]
        dist_dicc[arc]['n_sites'] = n_sites
        # For every pap sample...
        for samp in meta_df['IDX'].values:
            # Extract the sample's haplotypes.
            hap_1 = pap_hap_1[:, samp]
            hap_2 = pap_hap_2[:, samp]
            # Compute the distance from the archaic diplotype.
            dist_1 = sequence_divergence(hap_1, freq_dicc[arc])
            dist_2 = sequence_divergence(hap_2, freq_dicc[arc])
            # Update dictionaries.
            dist_dicc[arc]['hap_1'] = np.append(dist_dicc[arc]['hap_1'], dist_1)
            dist_dicc[arc]['hap_2'] = np.append(dist_dicc[arc]['hap_2'], dist_2)
            abs_dicc['Chr. 1 ({0})'.format(arc)].append(dist_1)
            abs_dicc['Chr. 2 ({0})'.format(arc)].append(dist_2)
            norm_dicc['Chr. 1 ({0})'.format(arc)].append(dist_1/n_sites)
            norm_dicc['Chr. 2 ({0})'.format(arc)].append(dist_2/n_sites)
    # Intailize the dataframes.
    abs_df = pd.DataFrame(data=abs_dicc)
    norm_df = pd.DataFrame(data=norm_dicc)
    return abs_df, norm_df, dist_dicc

# Define a function to plot haplostrips s curves for the tgp.
def tgp_plot_s_curves(hap_dist_dicc, window_size=72):
    # Intialize a dictionary to store haplotype distances.
    dist_dicc = {}
    # For every archaic...
    for arc in ['DEN', 'ALT', 'CHA', 'VIN']:
        # Fill the dictionary.
        dist_dicc[arc] = np.concatenate((hap_dist_dicc[arc]['hap_1'], hap_dist_dicc[arc]['hap_2']))
    # Concatenate all distances.
    all_dists = np.concatenate((dist_dicc['DEN'], dist_dicc['ALT'], dist_dicc['CHA'], dist_dicc['VIN']))
    # Intialize the figure and axes.
    fig, axes = plt.subplots(
        2, 2, figsize=(6, 4), dpi=300,
        sharex=True, sharey=True,
    )
    # Define lists needed for plotting.
    arcs = ['DEN', 'ALT', 'CHA', 'VIN']
    titles = [
        'Denisovan', 'Altai Nean.',
        'Chagyrskaya Nean.', 'Vindija Nean.',
    ]
    rows = [0, 0, 1, 1]
    cols = [0, 1, 0, 1]
    # For all arachics...
    for i in range(len(arcs)):
        # Plot the results.
        axes[rows[i], cols[i]].plot(
            np.sort(dist_dicc[arcs[i]]), np.arange(dist_dicc[arcs[i]].size),
            color='tab:blue', linestyle='', marker='x',
            markersize=5, markeredgewidth=0.5,
        )
        # Label the subplot.
        axes[rows[i], cols[i]].set_title(titles[i]+'\n'+'({0} sites)'.format(hap_dist_dicc[arcs[i]]['n_sites']))
        # Set the axes ticks.
        axes[rows[i], cols[i]].set_xticks(np.arange(all_dists.min(), all_dists.max() + 1, 25))
        axes[rows[i], cols[i]].set_xticklabels(
            np.arange(all_dists.min(), all_dists.max() + 1, 25).astype(str),
            size=6, rotation=45, ha='right', rotation_mode='anchor',
        )
        axes[rows[i], cols[i]].set_yticks(np.arange(0, dist_dicc[arcs[i]].size, 500))
        axes[rows[i], cols[i]].set_yticklabels(
            np.arange(0, dist_dicc[arc].size, 500).astype(str), size=6,
        )
    # Set the axes labels.
    axes[0, 0].set_ylabel('Haplotypes')
    axes[1, 0].set_ylabel('Haplotypes')
    axes[1, 0].set_xlabel('Number of Differences')
    axes[1, 1].set_xlabel('Number of Differences')
    # Show the plot.
    plt.show()
    return

# Define a function to plot haplostrips s curves for the Papuans.
def pap_plot_s_curves(hap_dist_dicc, window_size=72):
    # Intialize a dictionary to store haplotype distances.
    dist_dicc = {}
    # For every archaic...
    for arc in ['DEN', 'ALT', 'CHA', 'VIN']:
        # Fill the dictionary.
        dist_dicc[arc] = np.concatenate((hap_dist_dicc[arc]['hap_1'], hap_dist_dicc[arc]['hap_2']))
    # Concatenate all distances.
    all_dists = np.concatenate((dist_dicc['DEN'], dist_dicc['ALT'], dist_dicc['CHA'], dist_dicc['VIN']))
    # Intialize the figure and axes.
    fig, axes = plt.subplots(
        2, 2, figsize=(6, 4), dpi=300,
        sharex=True, sharey=True,
    )
    # Define lists needed for plotting.
    arcs = ['DEN', 'ALT', 'CHA', 'VIN']
    titles = [
        'Denisovan', 'Altai Nean.',
        'Chagyrskaya Nean.', 'Vindija Nean.',
    ]
    rows = [0, 0, 1, 1]
    cols = [0, 1, 0, 1]
    # For all arachics...
    for i in range(len(arcs)):
        # Plot the results.
        axes[rows[i], cols[i]].plot(
            np.sort(dist_dicc[arcs[i]]), np.arange(dist_dicc[arcs[i]].size),
            color='tab:blue', linestyle='', marker='x',
            markersize=5, markeredgewidth=0.5,
        )
        # Label the subplot.
        axes[rows[i], cols[i]].set_title(titles[i]+'\n'+'({0} sites)'.format(hap_dist_dicc[arcs[i]]['n_sites']))
        # Set the axes ticks.
        axes[rows[i], cols[i]].set_xticks(np.arange(all_dists.min(), all_dists.max() + 1, 25))
        axes[rows[i], cols[i]].set_xticklabels(
            np.arange(all_dists.min(), all_dists.max() + 1, 25).astype(str),
            size=6, rotation=45, ha='right', rotation_mode='anchor',
        )
        axes[rows[i], cols[i]].set_yticks(np.arange(0, dist_dicc[arcs[i]].size+1, 5))
        axes[rows[i], cols[i]].set_yticklabels(
            np.arange(0, dist_dicc[arc].size+1, 5).astype(str), size=6,
        )
    # Set the axes labels.
    axes[0, 0].set_ylabel('Haplotypes')
    axes[1, 0].set_ylabel('Haplotypes')
    axes[1, 0].set_xlabel('Number of Differences')
    axes[1, 1].set_xlabel('Number of Differences')
    # Show the plot.
    plt.show()
    return

# Define a function to plot haplotype distances focusing on the denisovan.
def plot_den_hap_dist_summary(dist_df, window_size=72, export=True):
    # Intialize figures and axes.
    fig, axes = plt.subplots(
         4, 5, figsize=(12, 9),
        sharex=True, sharey='row', dpi=300,
    )
    ### TGP Results ###
    # Intialize a super population list.
    s_pop_list = ['AFR', 'SAS', 'EAS', 'EUR', 'AMR']
    # Intialize a columns list.
    cols = [
        0, 1, 2, 3, 4,
    ]
    # Determine the max and min distance from the Denisovan.
    max_dist = dist_df[['Chr. 1 (DEN)', 'Chr. 2 (DEN)']].max().max()
    min_dist = dist_df[['Chr. 1 (DEN)', 'Chr. 2 (DEN)']].min().min()
    # For every super population...
    for idx in range(len(s_pop_list)):
        # Extract the super population.
        s_pop = s_pop_list[idx]
        # Extract the haplotype distances for the Denisovan and Neanderthals.
        den_1 = dist_df[dist_df['Super Population'].isin([s_pop])]['Chr. 1 (DEN)'].values
        den_2 = dist_df[dist_df['Super Population'].isin([s_pop])]['Chr. 2 (DEN)'].values
        den_dists = np.concatenate((den_1, den_2))
        alt_1 = dist_df[dist_df['Super Population'].isin([s_pop])]['Chr. 1 (ALT)'].values
        alt_2 = dist_df[dist_df['Super Population'].isin([s_pop])]['Chr. 2 (ALT)'].values
        alt_dists = np.concatenate((alt_1, alt_2))
        cha_1 = dist_df[dist_df['Super Population'].isin([s_pop])]['Chr. 1 (CHA)'].values
        cha_2 = dist_df[dist_df['Super Population'].isin([s_pop])]['Chr. 2 (CHA)'].values
        cha_dists = np.concatenate((cha_1, cha_2))
        vin_1 = dist_df[dist_df['Super Population'].isin([s_pop])]['Chr. 1 (VIN)'].values
        vin_2 = dist_df[dist_df['Super Population'].isin([s_pop])]['Chr. 2 (VIN)'].values
        vin_dists = np.concatenate((vin_1, vin_2))
        # Plot the haplotype distribution.
        axes[0, cols[idx]].hist(
            den_dists, bins=np.round(np.linspace(min_dist, max_dist, 25), 4),
            color='tab:blue', histtype='bar', edgecolor='black', linewidth=0.5,
        )
        # Plot the distance from the Denisovan vs distance from Neanderthals.
        axes[1, cols[idx]].scatter(
            den_dists, alt_dists, color='tab:blue',
            marker='x', linewidths=0.5,
        )
        axes[2, cols[idx]].scatter(
            den_dists, cha_dists, color='tab:blue',
            marker='x', linewidths=0.5,
        )
        axes[3, cols[idx]].scatter(
            den_dists, vin_dists, color='tab:blue',
            marker='x', linewidths=0.5,
        )
        # Add a title.
        axes[0, cols[idx]].set_title(s_pop)
        # Add x-axes labels.
        axes[3, cols[idx]].set_xlabel('Sequence Divergence'+'\n'+'from the Denisovan')
        # Despine the axes.
        axes[0, cols[idx]].spines['right'].set_visible(False)
        axes[0, cols[idx]].spines['top'].set_visible(False)
        axes[1, cols[idx]].spines['right'].set_visible(False)
        axes[1, cols[idx]].spines['top'].set_visible(False)
        axes[2, cols[idx]].spines['right'].set_visible(False)
        axes[2, cols[idx]].spines['top'].set_visible(False)
        axes[3, cols[idx]].spines['right'].set_visible(False)
        axes[3, cols[idx]].spines['top'].set_visible(False)
    # Add the y-axes labels.
    axes[0, 0].set_ylabel('Number of Haplotypes')
    axes[1, 0].set_ylabel('Sequence Divergence'+'\n'+'from the Altai Nean.')
    axes[2, 0].set_ylabel('Sequence Divergence'+'\n'+'from the Chagyrskaya Nean.')
    axes[3, 0].set_ylabel('Sequence Divergence'+'\n'+'from the Vindija Nean.')
    # If export is set to true...
    if export:
        # Export the plot.
        plt.savefig(
            f'./supp_figures/tgp_arc_hap_dist_summary_{window_size}kb.png', format='png',
            facecolor='white', bbox_inches='tight', dpi=500,
        )
    # Show the plot.
    plt.show()
    return

# Define a function to summarize tgp haplotype frequencies.
def tgp_hap_freq_summary(dist_dicc, arc, threshold, window_size=72):
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod_arc.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize an ordered dictionary.
    tgp_dicc = {
        'AFR': ['LWK', 'GWD', 'MSL', 'ESN', 'YRI'],
        'SAS': ['BEB', 'STU', 'ITU', 'PJL', 'GIH'],
        'EAS': ['CHB', 'KHV', 'CHS', 'JPT', 'CDX'],
        'EUR': ['TSI', 'CEU', 'IBS', 'GBR', 'FIN'],
        'AMR': ['PEL', 'MXL', 'CLM', 'PUR'],
    }
    # Intialize a samples dictionary.
    idx_dicc = {}
    # For every superpopulation.
    for s_pop in tgp_dicc.keys():
        # Fill the dictionary.
        idx_dicc[s_pop] = tgp_df[tgp_df['SUPERPOP'] == s_pop].index.values
        # For every population.
        for pop in tgp_dicc[s_pop]:
            # Fill the dictionary.
            idx_dicc[pop] = tgp_df[tgp_df['POP'] == pop].index.values
    # Intialize dictionaries.
    hap_dicc = {}
    hap_pop_dicc = {
        'Super Population': [],
        'Population': [],
        'Total Haplotypes': [],
        'Haplotype Frequency': [],
        'Homozygous Individuals': [],
        'Heterozygous Individuals': [],
    }
    hap_spop_dicc = {
        'Super Population': [],
        'Total Haplotypes': [],
        'Haplotype Frequency': [],
        'Homozygous Individuals': [],
        'Heterozygous Individuals': [],
    }
    # Extract the sample indicies for the Denisovan like haplotype.
    arc_hap_1_idx = np.where(dist_dicc[arc]['hap_1'] < threshold)[0]
    arc_hap_2_idx = np.where(dist_dicc[arc]['hap_2'] < threshold)[0]
    # For every super population...
    for s_pop in tgp_dicc.keys():
        # Extract the population indicies.
        s_all_idx = idx_dicc[s_pop]
        # Intialize arrays to store results.
        s_u_hap_idx = np.array([])
        s_hom_hap_idx = np.array([])
        s_het_hap_idx = np.array([])
        # Intialize a totla haplotype counter.
        s_all_hap = 0
        # For every population...
        for pop in tgp_dicc[s_pop]:
            # Extract the population indicies.
            all_idx = idx_dicc[pop]
            # Determine if there are any samples who carry the archaic like haplotype.
            hap_1_idx = np.intersect1d(all_idx, arc_hap_1_idx)
            hap_2_idx = np.intersect1d(all_idx, arc_hap_2_idx)
            all_hap_idx = np.concatenate((hap_1_idx, hap_2_idx))
            # Update the super population count.
            s_all_hap += all_hap_idx.size
            # Determine the unique haplotypes.
            u_hap_idx, hap_counts = np.unique(all_hap_idx, return_counts=True)
            # Create a mask for homozygous samples.
            hom_mask = hap_counts > 1
            # Determine the homozygous and hetyerozygous samples.
            hom_hap_idx = u_hap_idx[hom_mask]
            het_hap_idx = u_hap_idx[~hom_mask]
            # Append to the super population.
            s_u_hap_idx = np.append(s_u_hap_idx, u_hap_idx)
            s_hom_hap_idx = np.append(s_hom_hap_idx, hom_hap_idx)
            s_het_hap_idx = np.append(s_het_hap_idx, het_hap_idx)
            # Intialize and build the poulation dictionary.
            hap_dicc[pop] = {}
            hap_dicc[pop]['int_idx'] = u_hap_idx
            hap_dicc[pop]['hom_idx'] = hom_hap_idx
            hap_dicc[pop]['het_idx'] = het_hap_idx
            # Build the dataframe dictionary.
            hap_pop_dicc['Super Population'].append(s_pop)
            hap_pop_dicc['Population'].append(pop)
            hap_pop_dicc['Total Haplotypes'].append(all_hap_idx.size)
            hap_pop_dicc['Haplotype Frequency'].append(all_hap_idx.size/(all_idx.size*2))
            hap_pop_dicc['Homozygous Individuals'].append(hom_hap_idx.size)
            hap_pop_dicc['Heterozygous Individuals'].append(het_hap_idx.size)
        # Intialize and build the superpopulation dictionary.
        hap_dicc[s_pop] = {}
        hap_dicc[s_pop]['int_idx'] = s_u_hap_idx
        hap_dicc[s_pop]['hom_idx'] = s_hom_hap_idx
        hap_dicc[s_pop]['het_idx'] = s_het_hap_idx
        # Build the dataframe dictionary.
        hap_spop_dicc['Super Population'].append(s_pop)
        hap_spop_dicc['Total Haplotypes'].append(s_all_hap)
        hap_spop_dicc['Haplotype Frequency'].append(s_all_hap/(s_all_idx.size*2))
        hap_spop_dicc['Homozygous Individuals'].append(s_hom_hap_idx.size)
        hap_spop_dicc['Heterozygous Individuals'].append(s_het_hap_idx.size)
    # Convert the dictionaries to data frames.
    hap_pop_summary_df = pd.DataFrame(data=hap_pop_dicc)
    hap_spop_summary_df = pd.DataFrame(data=hap_spop_dicc)
    return hap_pop_summary_df, hap_spop_summary_df, hap_dicc

# Define a function to summarize pap haplotype frequencies.
def pap_hap_freq_summary(dist_dicc, arc, threshold, window_size=72):
    # Intialize a dictionary.
    hap_pop_dicc = {
        'Super Population': [],
        'Population': [],
        'Total Haplotypes': [],
        'Haplotype Frequency': [],
        'Homozygous Individuals': [],
        'Heterozygous Individuals': [],
    }
    # Extract the sample indicies for the Denisovan like haplotype.
    arc_hap_1_idx = np.where(dist_dicc[arc]['hap_1'] < threshold)[0]
    arc_hap_2_idx = np.where(dist_dicc[arc]['hap_2'] < threshold)[0]
    # Extract the population indicies.
    all_idx = np.arange(dist_dicc[arc]['hap_1'].size)
    # Determine if there are any samples who carry the archaic like haplotype.
    hap_1_idx = np.intersect1d(all_idx, arc_hap_1_idx)
    hap_2_idx = np.intersect1d(all_idx, arc_hap_2_idx)
    all_hap_idx = np.concatenate((hap_1_idx, hap_2_idx))
    # Determine the unique haplotypes.
    u_hap_idx, hap_counts = np.unique(all_hap_idx, return_counts=True)
    # Create a mask for homozygous samples.
    hom_mask = hap_counts > 1
    # Determine the homozygous and hetyerozygous samples.
    hom_hap_idx = u_hap_idx[hom_mask]
    het_hap_idx = u_hap_idx[~hom_mask]
    # Build the dataframe dictionary.
    hap_pop_dicc['Super Population'].append('MSEA')
    hap_pop_dicc['Population'].append('PNG')
    hap_pop_dicc['Total Haplotypes'].append(all_hap_idx.size)
    hap_pop_dicc['Haplotype Frequency'].append(all_hap_idx.size/(all_idx.size*2))
    hap_pop_dicc['Homozygous Individuals'].append(hom_hap_idx.size)
    hap_pop_dicc['Heterozygous Individuals'].append(het_hap_idx.size)
    # Convert the dictionaries to data frames.
    hap_pop_summary_df = pd.DataFrame(data=hap_pop_dicc)
    return hap_pop_summary_df

# Define a function to calculate sequence divergence between archaic diplotypes and tgp haplotypes.
def tgp_haplotype_divergence(p_gt, window_size=748):
    # Load the effective sequence length matrix.
    esl_mat = load_hap_esl_mat(prefix='arc_and_tgp', window_size=window_size)
    # Intialize archaic information dictionaries.
    freq_dicc = {
        'DEN': calc_arc_alt_freqs(p_gt.take([2350], axis=1)),
        'ALT': calc_arc_alt_freqs(p_gt.take([2347], axis=1)),
        'CHA': calc_arc_alt_freqs(p_gt.take([2348], axis=1)),
        'VIN': calc_arc_alt_freqs(p_gt.take([2349], axis=1)),
    }
    esl_dicc = {
        'DEN': esl_mat[1, 3], 'ALT': esl_mat[1, 0],
        'CHA': esl_mat[1, 1], 'VIN': esl_mat[1, 2],
    }
    pwd_dicc = {
        'ABS': {
            'DEN': {'hap_1': np.array([]), 'hap_2': np.array([])}, 
            'ALT': {'hap_1': np.array([]), 'hap_2': np.array([])},
            'CHA': {'hap_1': np.array([]), 'hap_2': np.array([])},
            'VIN': {'hap_1': np.array([]), 'hap_2': np.array([])},
        },
        'NORM': {
            'DEN': {'hap_1': np.array([]), 'hap_2': np.array([])}, 
            'ALT': {'hap_1': np.array([]), 'hap_2': np.array([])},
            'CHA': {'hap_1': np.array([]), 'hap_2': np.array([])},
            'VIN': {'hap_1': np.array([]), 'hap_2': np.array([])},
        },
    }
    # Extract haplotype arrays.
    tgp_hap_1 = p_gt[:, :-4, 0]
    tgp_hap_2 = p_gt[:, :-4, 1]
    # For every archaic individual...
    for arc in ['DEN', 'ALT', 'CHA', 'VIN']:
        # For every tgp sample...
        for samp in range(p_gt[:, :-4, :].shape[1]):
            # Extract the sample's haplotypes.
            hap_1 = tgp_hap_1[:, samp]
            hap_2 = tgp_hap_2[:, samp]
            # Compute the distance from the archaic diplotype.
            pwd_1 = sequence_divergence(hap_1, freq_dicc[arc])
            pwd_2 = sequence_divergence(hap_2, freq_dicc[arc])
            # Update dictionaries.
            pwd_dicc['ABS'][arc]['hap_1'] = np.append(pwd_dicc['ABS'][arc]['hap_1'], pwd_1)
            pwd_dicc['ABS'][arc]['hap_2'] = np.append(pwd_dicc['ABS'][arc]['hap_2'], pwd_2)
            pwd_dicc['NORM'][arc]['hap_1'] = np.append(pwd_dicc['NORM'][arc]['hap_1'], pwd_1 / esl_dicc[arc])
            pwd_dicc['NORM'][arc]['hap_2'] = np.append(pwd_dicc['NORM'][arc]['hap_2'], pwd_2 / esl_dicc[arc])
    return pwd_dicc

# Define a function to load the windowed haplotype divergence results.
def load_tgp_hap_div_windows(window_size=748):
    # Intialize dictionaries to store the results.
    pwd_dicc = {
        'ABS': {
            'DEN': {}, 'ALT': {}, 'CHA': {}, 'VIN': {},
        },
        'NORM': {
            'DEN': {}, 'ALT': {}, 'CHA': {}, 'VIN': {},
        },
    }
    # Load the effective sequence lengths.
    esl_df = load_windows('tgp', 'variant', window_size=window_size)
    # For every archaic...
    for arc in ['DEN', 'ALT', 'CHA', 'VIN']:
        # For every haplotype.
        for hap in ['hap_1', 'hap_2']:
            # For all chromosomes...
            for chrom in range(1, 23):
                # If this is the first chromosome...
                if chrom == 1:
                    # Load the matricies.
                    pw_mat = np.loadtxt(
                        f'../sequence_divergence/tgp/windows/{arc.lower()}_{hap}_pw_diffs_chr{chrom}_{window_size}kb.csv.gz',
                        delimiter=',',
                    )
                    # Extract the effective sequence lengths.
                    esl = esl_df[esl_df['CHR'] == chrom][arc].values
                    # Intialize the results.
                    pwd_dicc['ABS'][arc][hap] = pw_mat
                    pwd_dicc['NORM'][arc][hap] = pw_mat / esl[:, np.newaxis]
                # Else...
                else:
                    # Load the matricies.
                    pw_mat = np.loadtxt(
                        f'../sequence_divergence/tgp/windows/{arc.lower()}_{hap}_pw_diffs_chr{chrom}_{window_size}kb.csv.gz',
                        delimiter=',',
                    )
                    # Extract the effective sequence lengths and heterozygous sites.
                    esl = esl_df[esl_df['CHR'] == chrom][arc].values
                    # Append the results.
                    pwd_dicc['ABS'][arc][hap] = np.concatenate((pwd_dicc['ABS'][arc][hap], pw_mat), axis=0)
                    pwd_dicc['NORM'][arc][hap] = np.concatenate((pwd_dicc['NORM'][arc][hap], pw_mat / esl[:, np.newaxis]), axis=0)
    return pwd_dicc

# Define a function to compile the haplotype divergence summary.
def compile_hap_div_summary(p_gt, focal_ind, window_size=748, export=True):
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Extract the sample and population arrays.
    inds = tgp_df['IND'].values
    pops = tgp_df['POP'].values
    # Intialize archaic labels.
    arc_dicc = {
        'DEN': 'Denisovan', 'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.', 'VIN': 'Vindija Nean.',
    }
    # Intialize dictionaries to store the results.
    pwd_dicc = {
        'ind': [],
        'pop': [], 'arc': [],
        'muc19_a1': [], 'muc19_a2': [],
        'muc19_n1': [], 'muc19_n2': [],
        'wind_m1': [], 'wind_m2': [], 
        'wind_s1': [], 'wind_s2': [],
        'wind_p1': [], 'wind_p2': [],
    }
    
    # Load the observed, nonoverlapping windows, and window indicies of,
    # comparable effective sequence length.
    muc19_pwd = tgp_haplotype_divergence(p_gt=p_gt, window_size=window_size)
    wind_idx = load_esl_qc_windows_idx('tgp', window_size=window_size)
    wind_pwd = load_tgp_hap_div_windows(window_size=window_size)
    # For every archaic...
    for arc in ['DEN', 'ALT', 'CHA', 'VIN']:
        # For every individual...
        for idx, ind in enumerate(inds):
            # Determine the mean and standard deviation for non-overlapping windows.
            wind_pwd_m1 = np.nanmean(wind_pwd['NORM'][arc]['hap_1'][:, idx][wind_idx])
            wind_pwd_m2 = np.nanmean(wind_pwd['NORM'][arc]['hap_2'][:, idx][wind_idx])
            wind_pwd_s1 = np.nanstd(wind_pwd['NORM'][arc]['hap_1'][:, idx][wind_idx])
            wind_pwd_s2 = np.nanstd(wind_pwd['NORM'][arc]['hap_2'][:, idx][wind_idx])
            # Compute the p-values.
            pwd_p1 = (np.count_nonzero(muc19_pwd['NORM'][arc]['hap_1'][idx] >= wind_pwd['NORM'][arc]['hap_1'][:, idx][wind_idx])
                        / np.sum(~np.isnan(wind_pwd['NORM'][arc]['hap_1'][:, idx][wind_idx])))
            pwd_p2 = (np.count_nonzero(muc19_pwd['NORM'][arc]['hap_2'][idx] >= wind_pwd['NORM'][arc]['hap_2'][:, idx][wind_idx])
                        / np.sum(~np.isnan(wind_pwd['NORM'][arc]['hap_2'][:, idx][wind_idx])))
            # Fill the dictionaries.
            pwd_dicc['ind'].append(ind)
            pwd_dicc['pop'].append(pops[idx])
            pwd_dicc['arc'].append(arc_dicc[arc])
            pwd_dicc['muc19_a1'].append(muc19_pwd['ABS'][arc]['hap_1'][idx])
            pwd_dicc['muc19_a2'].append(muc19_pwd['ABS'][arc]['hap_2'][idx])
            pwd_dicc['muc19_n1'].append(muc19_pwd['NORM'][arc]['hap_1'][idx])
            pwd_dicc['muc19_n2'].append(muc19_pwd['NORM'][arc]['hap_2'][idx])
            pwd_dicc['wind_m1'].append(wind_pwd_m1)
            pwd_dicc['wind_m2'].append(wind_pwd_m2)
            pwd_dicc['wind_s1'].append(wind_pwd_s1)
            pwd_dicc['wind_s2'].append(wind_pwd_s2)
            pwd_dicc['wind_p1'].append(pwd_p1)
            pwd_dicc['wind_p2'].append(pwd_p2)
    # Convert the dictionaries to dataframes.
    pwd_df = pd.DataFrame(data=pwd_dicc)
    # Rename all the columns to look pretty.
    pwd_df.rename(
        columns={
            'ind': 'Individual',
            'pop': 'Population', 'arc': 'Archaic',
            'muc19_a1': f'{window_size}kb (Pairwise Differences Chr. 1)',
            'muc19_a2': f'{window_size}kb (Pairwise Differences Chr. 2)',
            'muc19_n1': f'{window_size}kb (Sequence Divergence Chr. 1)',
            'muc19_n2': f'{window_size}kb (Sequence Divergence Chr. 2)',
            'wind_m1': r'Nonoverlapping Windows ($\mu$ Chr. 1)',
            'wind_m2': r'Nonoverlapping Windows ($\mu$ Chr. 2)',
            'wind_s1': r'Nonoverlapping Windows ($\sigma$ Chr. 1)',
            'wind_s2': r'Nonoverlapping Windows ($\sigma$ Chr. 2)',
            'wind_p1': r'$P-value$ (Chr. 1)',
            'wind_p2': r'$P-value$ (Chr. 2)',
        }, inplace=True,
    )
    # Show the results for the focal individual of interest.
    pwd_df = pwd_df[pwd_df['Individual'] == focal_ind].reset_index(drop=True)
    # If export is true...
    if export:
        # Export the dataframe as a csv.
        pwd_df.to_csv(f'./dataframes/{focal_ind.lower()}_longest_tract_haplotype_divergence_{window_size}kb.csv', index=False)
    return pwd_df


######################
### HETEROZYGOSITY ###
######################

# Define a function to compute the number of heterozygous sites for the archaics.
def arc_het(gt, window_size=72):
    # Intialize a list of archaics and there indicies.
    arc_idx_dicc = {
        'DEN': 3, 'ALT': 0,
        'CHA': 1, 'VIN': 2,
    }
    # Intialize a dictionary for heterozygosity results.
    het_dicc = {}
    # Count the number of heterozygous sites per archaic.
    het_counts = gt.count_het(axis=0)[:-1]
    # For every archaic...
    for arc in arc_idx_dicc.keys():
        # Fill the dictionary.
        het_dicc[arc] = het_counts[arc_idx_dicc[arc]]
    return het_dicc

# Define a function to load the heterozygosity results from the window data.
def load_arc_het_windows(window_size=72):
    # Intialize a list of archaics and there indicies.
    arc_idx_dicc = {
        'DEN': 3, 'ALT': 0,
        'CHA': 1, 'VIN': 2,
    }
    # Intialize a dictionary for heterozygosity results.
    het_dicc = {
        'DEN': np.array([]), 'ALT': np.array([]), 'CHA': np.array([]), 'VIN': np.array([]),
    }
    # For all chromosomes...
    for chrom in range(1, 23):
        # Load the window data.
        het_counts = np.loadtxt(
            f'../heterozygosity/arc/windows/arc_het_counts_chr{chrom}_{window_size}kb.csv',
            delimiter=',',
        )
        # For every archaic...
        for arc in arc_idx_dicc.keys():
            # Fill the dictionary.
            het_dicc[arc] = np.append(het_dicc[arc], het_counts[:, arc_idx_dicc[arc]])
    return het_dicc

# Define a function to compile and summarize the archaic heterozygosity results.
def compile_arc_het_summary(gt, window_size=72, export=True):
    # Intialize a list of archaics.
    arc_dicc = {
        'DEN': 'Denisovan', 'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.', 'VIN': 'Vindija Nean.',
    }
    # Intialize a dictionary to store the results.
    df_dicc = {
        'arc': [], 'muc19_a': [],
        'wind_m': [], 'wind_s': [],
        'wind_p': [],
    }
    # Load the observed, nonoverlapping windows, and window indicies of,
    # comparable effective sequence length.
    muc19_dicc = arc_het(gt, window_size=window_size)
    wind_dicc = load_arc_het_windows(window_size=window_size)
    wind_idx = load_esl_qc_windows_idx('arc', window_size=window_size)
    # For every archaic...
    for arc in arc_dicc.keys():
        # Determine the mean and standard deviation for non-overlapping windows.
        wind_m = np.nanmean(wind_dicc[arc][wind_idx])
        wind_s = np.nanstd(wind_dicc[arc][wind_idx])
        # Compute the p-value.
        wind_p = (np.count_nonzero(muc19_dicc[arc] <= wind_dicc[arc][wind_idx])
                  / np.sum(~np.isnan(wind_dicc[arc][wind_idx])))
        # Append the results.
        df_dicc['arc'].append(arc_dicc[arc])
        df_dicc['muc19_a'].append(muc19_dicc[arc])
        df_dicc['wind_m'].append(wind_m)
        df_dicc['wind_s'].append(wind_s)
        df_dicc['wind_p'].append(wind_p)
    # Convert the dictionary to a dataframe.
    het_df = pd.DataFrame(data=df_dicc)
    # Adjust p-values of 0 by the number of permutations.
    het_df['wind_p'] = np.where(het_df['wind_p'] == 0, '<3.242e-05', het_df['wind_p'])
    # Rename all the columns to look pretty.
    het_df.rename(
        columns={
            'arc': 'Archaic',
            'muc19_a': r'$MUC19$ (Heterozygous Sites)',
            'wind_m': r'Nonoverlapping Windows $\left( \mu \right)$',
            'wind_s': r'Nonoverlapping Windows $\left( \sigma \right)$',
            'wind_p': r'$P-value$',
        }, inplace=True,
    )
    # If export is true...
    if export:
        # Export the dataframe as a csv.
        het_df.to_csv(f'./dataframes/archaic_heterozygous_sites_{window_size}kb.csv', index=False)
    return het_df

# Define a function to plot the distribution of heterozygous sites amongst the archaics.
def plot_arc_het(gt, obs_label, window_size=72, export=True):
    # Intialize a list of archaics.
    arc_dicc = {
        'DEN': 'Denisovan', 'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.', 'VIN': 'Vindija Nean.',
    }
    arc_list = ['DEN', 'ALT', 'CHA', 'VIN']
    # Intialize a row and columns list.
    rows = [0, 0, 1, 1]
    cols = [0, 1, 0, 1]
    # Load the observed, nonoverlapping windows, and window indicies of,
    # comparable effective sequence length.
    muc19_dicc = arc_het(gt, window_size=window_size)
    wind_dicc = load_arc_het_windows(window_size=window_size)
    wind_idx = load_esl_qc_windows_idx('arc', window_size=window_size)
    # Determine the significance threshold.
    sig_threshold = 0.05 / 4
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        2, 2, figsize=(6, 4), dpi=300,
        sharex=True, sharey=True,
    )
    # For every population...
    for idx in range(len(arc_list)):
        # Extract the archaics 
        arc = arc_list[idx]
        # Compute the p-value.
        wind_p = (np.count_nonzero(muc19_dicc[arc] <= wind_dicc[arc][wind_idx])
                  / np.sum(~np.isnan(wind_dicc[arc][wind_idx])))
        # Plot the distribution of heterozygous sites.
        axes[rows[idx], cols[idx]].hist(
            wind_dicc[arc][wind_idx],
            bins=np.arange(0, 200, 1),
            histtype='stepfilled', color='tab:blue',
        )
        # Plot the observed value for muc19.
        axes[rows[idx], cols[idx]].axvline(
            muc19_dicc[arc], 0, 1,
            color='black', linestyle='dashed',
        )
        # If the p-value is significant after correcting for multiple comparisons...
        if wind_p < sig_threshold:
            # If the rounded p-value is 0...
            if wind_p == 0:
                # Construct the title.
                title = arc_dicc[arc]+'\n'+r'$P-value=$'+'<3.242e-05'+r'$^{*}$'
            else:
                # Construct the title.
                title = arc_dicc[arc]+'\n'+r'$P-value=$'+'{:.3e}'.format(wind_p)+r'$^{*}$'
        # Else...
        else:
            # Construct the title.
            title = arc_dicc[arc]+'\n'+r'$P-value=$'+'{:.3e}'.format(wind_p)+r'$^{ns}$'
        # Add the subplot title.
        axes[rows[idx], cols[idx]].set_title(title)
    # Configure the legend.
    legend_elements = [
        Line2D([0], [0], color='black', linestyle='dashed', label=obs_label),
    ]
    # Add a figure legend.
    fig.legend(
        handles=legend_elements, loc='center left',
        bbox_to_anchor=(1.0, 0.5), frameon=False,
    )
    # Label the super-axes.
    fig.supylabel('Frequency')
    fig.supxlabel(f'Heterozygous Sites in {window_size}kb Windows')
    # If export is set to true...
    if export:
        # Export the plot.
        plt.savefig(
            f'./supp_figures/archaic_heterozygous_sites_{window_size}kb.png', format='png',
            facecolor='white', bbox_inches='tight', dpi=500,
        )
    # Show the plot!
    plt.show()
    return

# Define a function to compute heterozygosity for the combined dataset.
def tgp_het(gt, window_size=72):
    # Intialize a dictionary for introgressed individuals.
    int_dicc = {
        'HET': np.loadtxt(f'../meta_data/{window_size}kb_all_het_int_idx.csv', delimiter=',', dtype=int),
        'HOM': np.loadtxt(f'../meta_data/{window_size}kb_all_hom_int_idx.csv', delimiter=',', dtype=int),
    }
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod_arc.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize an ordered population list.
    pop_list = [
        'LWK', 'GWD', 'MSL', 'ESN', 'YRI', # AFR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
        'PEL', 'MXL', 'CLM', 'PUR', # AMR.
        'DEN', 'ALT', 'CHA', 'VIN', #ARC.
    ]
    # Intialize a dictionary to store the results.
    het_dicc = {}
    # Determine the number of heterozygous sites in MUC19 per ind.
    het_array = gt.count_het(axis=0)[:-1]
    # For every population...
    for pop in pop_list:
        # Extract the sample indicies and ids.
        pop_idx = tgp_df[tgp_df['POP'] == pop].index.values
        pop_ind = tgp_df[tgp_df['POP'] == pop]['IND'].values
        # Intialize the subdictionary.
        het_dicc[pop] = {}
        # Determine the het sites for this population.
        hets = het_array[pop_idx]
        # Determine the average number of het sites.
        het_dicc[pop]['AVG'] = np.nanmean(hets)
        # For every individual in the target population...
        for idx in range(pop_idx.size):
            # Extract the individual's index and id.
            samp_idx = pop_idx[idx]
            samp_id = pop_ind[idx]
            # Intialize the subdictionaries.
            het_dicc[pop][samp_id] = {}
            # Determine the individual's number of heterozygous sites.
            ind_het = hets[idx]
            # Fill the dictionaries.
            het_dicc[pop][samp_id]['HET'] = ind_het
            # If the individual carries two Denisovan-like haplotypes...
            if samp_idx in int_dicc['HOM']:
                # Fill the dictionaries.
                het_dicc[pop][samp_id]['INT'] = 2
            # Else-if the individual carries one Denisovan-like haplotype...
            elif samp_idx in int_dicc['HET']:
                # Fill the dictionaries.
                het_dicc[pop][samp_id]['INT'] = 1
            # Else...
            else:
            # Fill the dictionaries.
                het_dicc[pop][samp_id]['INT'] = 0
    return het_dicc

# Define a function to load the heterozygosity results from the window data.
def load_tgp_het_windows(window_size=72):
    # Intialize a dictionary for introgressed individuals.
    int_dicc = {
        'HET': np.loadtxt(f'../meta_data/{window_size}kb_all_het_int_idx.csv', delimiter=',', dtype=int),
        'HOM': np.loadtxt(f'../meta_data/{window_size}kb_all_hom_int_idx.csv', delimiter=',', dtype=int),
    }
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod_arc.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize an ordered population list.
    pop_list = [
        'LWK', 'GWD', 'MSL', 'ESN', 'YRI', # AFR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
        'PEL', 'MXL', 'CLM', 'PUR', # AMR.
        'DEN', 'ALT', 'CHA', 'VIN', #ARC.
    ]
    # Intialize a dictionary to store the results.
    het_dicc = {}
    # For every population...
    for pop in pop_list:
        # Extract the sample indicies and ids.
        pop_idx = tgp_df[tgp_df['POP'] == pop].index.values
        pop_ind = tgp_df[tgp_df['POP'] == pop]['IND'].values
        # Intialize the subdictionaries.
        het_dicc[pop] = {}
        het_dicc[pop]['AVG'] = np.array([])
        # For every individual in the target population...
        for idx in range(pop_idx.size):
            # Extract the individual's index and id.
            samp_idx = pop_idx[idx]
            samp_id = pop_ind[idx]
            # Intialize the subdictionaries.
            het_dicc[pop][samp_id] = {}
            het_dicc[pop][samp_id]['HET'] = np.array([])
            # If the individual carries two Denisovan-like haplotypes...
            if samp_idx in int_dicc['HOM']:
                # Fill the dictionaries.
                het_dicc[pop][samp_id]['INT'] = 2
            # Else-if the individual carries one Denisovan-like haplotype...
            elif samp_idx in int_dicc['HET']:
                # Fill the dictionaries.
                het_dicc[pop][samp_id]['INT'] = 1
            # Else...
            else:
                # Fill the dictionaries.
                het_dicc[pop][samp_id]['INT'] = 0
    # For all chromosomes...
    for chrom in range(1, 23):
        # Load the heterozygosity data.
        hets = np.loadtxt(
            f'../heterozygosity/tgp/windows/tgp_arc_het_counts_chr{chrom}_{window_size}kb.csv',
            delimiter=',',
        )
        # For every population...
        for pop in pop_list:
            # Extract the sample indicies and ids.
            pop_idx = tgp_df[tgp_df['POP'] == pop].index.values
            pop_ind = tgp_df[tgp_df['POP'] == pop]['IND'].values
            # Extract the population results.
            pop_het = hets[:, pop_idx]
            # Append the population results.
            het_dicc[pop]['AVG'] = np.append(het_dicc[pop]['AVG'], np.nanmean(pop_het, axis=1))
            # For every individual in the target population...
            for idx in range(pop_idx.size):
                # Extract the individual's id.
                samp_id = pop_ind[idx]
                # Fill the dictionaries.
                het_dicc[pop][samp_id]['HET'] = np.append(het_dicc[pop][samp_id]['HET'], pop_het[:, idx])
    return het_dicc

# Define a function to compile and summarize the heterozygosity results.
def compile_afr_hap_het_summary(gt, window_size=72, export=True):
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod_arc.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary for introgressed individuals.
    int_dicc = {
        'HET': np.loadtxt(f'../meta_data/{window_size}kb_all_het_int_idx.csv', delimiter=',', dtype=int),
        'HOM': np.loadtxt(f'../meta_data/{window_size}kb_all_hom_int_idx.csv', delimiter=',', dtype=int),
    }
    # Intialize a dictionary of group labels.
    group_labels = {
        'AFR': 'African Inds.',
        'HET': 'Heterozygous Inds.',
        'HOM': 'Homozygous Inds.',
    }
    # Intialize ordered population lists.
    ooa_pops = [
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
        'PEL', 'MXL', 'CLM', 'PUR', # AMR.
    ]
    afr_pops = ['LWK', 'GWD', 'MSL', 'ESN', 'YRI']
    # Load the effective sequence lengths.
    esl_df = load_windows('tgp', 'variant', window_size=window_size)
    # Intialize dictionaries to store the results.
    mat_dicc = {
        'AFR': {
            'M': np.array([]),
            'W': np.empty((esl_df.shape[0], tgp_df[tgp_df['SUPERPOP'] == 'AFR'].shape[0])),
        },
        'HET': {
            'M': np.array([]),
            'W': np.empty((esl_df.shape[0], int_dicc['HET'].size)),
        },
        'HOM': {
            'M': np.array([]),
            'W': np.empty((esl_df.shape[0], int_dicc['HOM'].size)),
        },
    }
    df_dicc = {
        'group': [],
        'muc19_a': [],
        'wind_m': [], 'wind_s': [], 'wind_p': [],
    }
    # Load the observed, nonoverlapping windows, and window indicies of,
    # comparable effective sequence length.
    muc19_dicc = tgp_het(gt, window_size=window_size)
    wind_dicc = load_tgp_het_windows(window_size=window_size)
    wind_idx = load_esl_qc_windows_idx('tgp', window_size=window_size)
    # Intialize a counter for the AFR individuals.
    afr_c = 0
    # For every AFR population...
    for pop in afr_pops:
        # For every individual...
        for key in muc19_dicc[pop].keys():
            # If the key isn't the population average key...
            if key != 'AVG':
                # Fill the dictionary.
                mat_dicc['AFR']['M'] = np.append(mat_dicc['AFR']['M'], muc19_dicc[pop][key]['HET'])
                mat_dicc['AFR']['W'][:, afr_c] = wind_dicc[pop][key]['HET']
                # Move the Africa counter forward.
                afr_c += 1
    # Intialize a counter for the het and hom individuals.
    het_c = 0
    hom_c = 0
    # For every AFR population...
    for pop in ooa_pops:
        # For every individual...
        for key in muc19_dicc[pop].keys():
            # If the key isn't the population average key...
            if key != 'AVG':
                # If this is a het individual...
                if muc19_dicc[pop][key]['INT'] == 1:
                    # Fill the dictionary.
                    mat_dicc['HET']['M'] = np.append(mat_dicc['HET']['M'], muc19_dicc[pop][key]['HET'])
                    mat_dicc['HET']['W'][:, het_c] = wind_dicc[pop][key]['HET']
                    # Move the het counter forward.
                    het_c += 1
                # Else-if this is a hom individual.
                elif muc19_dicc[pop][key]['INT'] == 2:
                    # Fill the dictionary.
                    mat_dicc['HOM']['M'] = np.append(mat_dicc['HOM']['M'], muc19_dicc[pop][key]['HET'])
                    mat_dicc['HOM']['W'][:, hom_c] = wind_dicc[pop][key]['HET']
                    # Move the hom counter forward.
                    hom_c += 1
    # For every focal group...
    for key in mat_dicc.keys():
        # Compute the means.
        muc19_mean = np.nanmean(mat_dicc[key]['M'])
        wind_mean = np.nanmean(mat_dicc[key]['W'], axis=1)
        # Determine the mean and standard deviation for non-overlapping windows.
        wind_m = np.nanmean(wind_mean[wind_idx])
        wind_s = np.nanstd(wind_mean[wind_idx])
        # If this is the homozygous individuals...
        if key == 'HOM':
            # Compute the p-value.
            wind_p = (np.count_nonzero(muc19_mean >= wind_mean[wind_idx])
                      / np.sum(~np.isnan(wind_mean[wind_idx])))
        # Else...
        else:
            # Compute the p-value.
            wind_p = (np.count_nonzero(muc19_mean <= wind_mean[wind_idx])
                      / np.sum(~np.isnan(wind_mean[wind_idx])))
        # Append the population results.
        df_dicc['group'].append(group_labels[key])
        df_dicc['muc19_a'].append(muc19_mean)
        df_dicc['wind_m'].append(wind_m)
        df_dicc['wind_s'].append(wind_s)
        df_dicc['wind_p'].append(wind_p)
    # Convert the dictionary to a dataframe.
    group_df = pd.DataFrame(data=df_dicc)
    # Adjust p-values of 0 by the number of permutations.
    group_df['wind_p'] = np.where(group_df['wind_p'] == 0, '<3.242e-05', group_df['wind_p'])
    # Rename all the columns to look pretty.
    group_df.rename(
        columns={
            'group': 'Group',
            'muc19_a': r'$MUC19$ (Average Number of Heterozygous Sites)',
            'wind_m': r'Nonoverlapping Windows $\left( \mu \right)$',
            'wind_s': r'Nonoverlapping Windows $\left( \sigma \right)$',
            'wind_p': r'$P-value$',
        }, inplace=True,
    )
    # If export is true...
    if export:
        # Export the dataframe as a csv.
        group_df.to_csv(f'./dataframes/african_and_introgressed_inds_heterozygous_sites_{window_size}kb.csv', index=False)
    return group_df

# Define a function to plot the distribution of heterozygous sites amongst the focal groups.
def plot_afr_hap_het(gt, obs_label, window_size=72, export=True):
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod_arc.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize a dictionary for introgressed individuals.
    int_dicc = {
        'HET': np.loadtxt(f'../meta_data/{window_size}kb_all_het_int_idx.csv', delimiter=',', dtype=int),
        'HOM': np.loadtxt(f'../meta_data/{window_size}kb_all_hom_int_idx.csv', delimiter=',', dtype=int),
    }
    # Intialize a dictionary of group labels.
    group_labels = {
        'AFR': 'African Inds.',
        'HET': 'Heterozygous Inds.',
        'HOM': 'Homozygous Inds.',
    }
    # Intialize ordered population lists.
    ooa_pops = [
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
        'PEL', 'MXL', 'CLM', 'PUR', # AMR.
    ]
    afr_pops = ['LWK', 'GWD', 'MSL', 'ESN', 'YRI']
    # Load the effective sequence lengths.
    esl_df = load_windows('tgp', 'variant', window_size=window_size)
    # Intialize dictionaries to store the results.
    mat_dicc = {
        'AFR': {
            'M': np.array([]),
            'W': np.empty((esl_df.shape[0], tgp_df[tgp_df['SUPERPOP'] == 'AFR'].shape[0])),
        },
        'HET': {
            'M': np.array([]),
            'W': np.empty((esl_df.shape[0], int_dicc['HET'].size)),
        },
        'HOM': {
            'M': np.array([]),
            'W': np.empty((esl_df.shape[0], int_dicc['HOM'].size)),
        },
    }
    # Load the observed, nonoverlapping windows, and window indicies of,
    # comparable effective sequence length.
    muc19_dicc = tgp_het(gt, window_size=window_size)
    wind_dicc = load_tgp_het_windows(window_size=window_size)
    wind_idx = load_esl_qc_windows_idx('tgp', window_size=window_size)
    # Intialize a counter for the AFR individuals.
    afr_c = 0
    # For every AFR population...
    for pop in afr_pops:
        # For every individual...
        for key in muc19_dicc[pop].keys():
            # If the key isn't the population average key...
            if key != 'AVG':
                # Fill the dictionary.
                mat_dicc['AFR']['M'] = np.append(mat_dicc['AFR']['M'], muc19_dicc[pop][key]['HET'])
                mat_dicc['AFR']['W'][:, afr_c] = wind_dicc[pop][key]['HET']
                # Move the Africa counter forward.
                afr_c += 1
    # Intialize a counter for the het and hom individuals.
    het_c = 0
    hom_c = 0
    # For every AFR population...
    for pop in ooa_pops:
        # For every individual...
        for key in muc19_dicc[pop].keys():
            # If the key isn't the population average key...
            if key != 'AVG':
                # If this is a het individual...
                if muc19_dicc[pop][key]['INT'] == 1:
                    # Fill the dictionary.
                    mat_dicc['HET']['M'] = np.append(mat_dicc['HET']['M'], muc19_dicc[pop][key]['HET'])
                    mat_dicc['HET']['W'][:, het_c] = wind_dicc[pop][key]['HET']
                    # Move the het counter forward.
                    het_c += 1
                # Else-if this is a hom individual.
                elif muc19_dicc[pop][key]['INT'] == 2:
                    # Fill the dictionary.
                    mat_dicc['HOM']['M'] = np.append(mat_dicc['HOM']['M'], muc19_dicc[pop][key]['HET'])
                    mat_dicc['HOM']['W'][:, hom_c] = wind_dicc[pop][key]['HET']
                    # Move the hom counter forward.
                    hom_c += 1
    # Determine the significance threshold.
    sig_threshold = 0.05 / 3
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 3, figsize=(6, 3), dpi=300,
        sharex=True, sharey=True,
    )
    # For every population...
    for idx in range(len(list(group_labels.keys()))):
        # Extract the the focal group.
        group = list(group_labels.keys())[idx]
        # Compute the means.
        muc19_mean = np.nanmean(mat_dicc[group]['M'])
        wind_mean = np.nanmean(mat_dicc[group]['W'], axis=1)
        # If this is the homozygous individuals...
        if group == 'HOM':
            # Compute the p-value.
            wind_p = (np.count_nonzero(muc19_mean >= wind_mean[wind_idx])
                      / np.sum(~np.isnan(wind_mean[wind_idx])))
        # Else...
        else:
            # Compute the p-value.
            wind_p = (np.count_nonzero(muc19_mean <= wind_mean[wind_idx])
                      / np.sum(~np.isnan(wind_mean[wind_idx])))
        # Plot the distribution of heterozygous sites.
        axes[idx].hist(
            wind_mean,
            bins=np.arange(0, 220, 1),
            histtype='stepfilled', color='tab:blue',
        )
        # Plot the observed value for muc19.
        axes[idx].axvline(
            muc19_mean, 0, 1,
            color='black', linestyle='dashed',
        )
        # If the p-value is significant after correcting for multiple comparisons...
        if wind_p < sig_threshold:
            # If the p-value is 0.
            if wind_p == 0:
                # Construct the title.
                title = group_labels[group]+'\n'+r'$P-value=$'+'<3.242e-05'+r'$^{*}$'
            # Else...
            else:
                # Construct the title.
                title = group_labels[group]+'\n'+r'$P-value=$'+'{:.3e}'.format(wind_p)+r'$^{*}$'
        # Else...
        else:
            # Construct the title.
            title = group_labels[group]+'\n'+r'$P-value=$'+'{:.3e}'.format(wind_p)+r'$^{ns}$'
        # Add the subplot title.
        axes[idx].set_title(title)
    # Configure the legend.
    legend_elements = [
        Line2D([0], [0], color='black', linestyle='dashed', label=obs_label),
    ]
    # Add a figure legend.
    fig.legend(
        handles=legend_elements, loc='center left',
        bbox_to_anchor=(1.0, 0.5), frameon=False,
    )
    # Label the super-axes.
    fig.supylabel('Frequency')
    fig.supxlabel(f'Average Number of Heterozygous Sites in {window_size}kb Windows')
    # If export is set to true...
    if export:
        # Export the plot.
        plt.savefig(
            f'./supp_figures/african_and_introgressed_inds_heterozygous_sites_{window_size}kb.png', format='png',
            facecolor='white', bbox_inches='tight', dpi=500,
        )
    # Show the plot!
    plt.show()
    return


############
### PAP ###
###########

# Define a function to calculate pap scores.
def psuedo_ancestry_painting(het_p_gt, target, source_1, source_2):
    # Determine the derived allele counts.
    t_dac = het_p_gt[:, target[0], :].sum(axis=1)
    s1_dac = het_p_gt[:, source_1[0], :].sum(axis=1)
    s2_dac = het_p_gt[:, source_2[0], :].sum(axis=1)
    # Determine the sites that can be painted.
    s1_anc_s2_der = np.where((s1_dac == 0) & (s2_dac == 2) & (t_dac == 1))[0]
    s1_der_s2_anc = np.where((s1_dac == 2) & (s2_dac == 0) & (t_dac == 1))[0]
    # Determine pap sites and its complement.
    pap_idx = np.union1d(s1_anc_s2_der, s1_der_s2_anc)
    com_idx = np.setdiff1d(np.arange(0, het_p_gt.shape[0]), pap_idx)
    return pap_idx, com_idx


#####################
### SITE PATTERNS ###
#####################

# Define a function to calculate site pattern counts.
def site_patterns(p1, p2, p3):
    # Calculate site pattern counts.
    abba = np.sum((1 - p1) * (p2) * (p3))
    baba = np.sum((p1) * (1 - p2) * (p3))
    bbaa = np.sum((p1) * (p2) * (1 - p3))
    baaa = np.sum((p1) * (1 - p2) * (1 - p3))
    abaa = np.sum((1 - p1) * (p2) * (1 - p3))
    aaba = np.sum((1 - p1) * (1 - p2) * (p3))
    return abba, baba, bbaa, baaa, abaa, aaba

# Define a function to calculate site patterns between the TGP and archaics.
def tgp_arc_site_patterns(
    gt,
    p1_idx, p2_idx, p3_idx,
):
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
        var_mask = gt.take(samp_list, axis=1).compress(called_mask, axis=0).count_alleles().is_variant()
        # If there are no variable sites...
        if (var_mask.sum() == 0):
            # Set the results to 0 since we are iterating over QC'ed regions.
            results = np.zeros(6)
        # Else...
        else:
            # Calculate the alternative allele frequencies.
            p1_alt_freqs = calc_alt_freqs(gt.take(p1_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p2_alt_freqs = calc_alt_freqs(gt.take(p2_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p3_alt_freqs = calc_alt_freqs(gt.take(p3_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            anc_freqs = calc_alt_freqs(gt.take([-1], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            # Polarize the samples.
            p1_der_freqs = np.where(anc_freqs == 1, np.abs(p1_alt_freqs - 1), p1_alt_freqs)
            p2_der_freqs = np.where(anc_freqs == 1, np.abs(p2_alt_freqs - 1), p2_alt_freqs)
            p3_der_freqs = np.where(anc_freqs == 1, np.abs(p3_alt_freqs - 1), p3_alt_freqs)
            # Calculate the site pattern counts.
            abba, baba, bbaa, baaa, abaa, aaba = site_patterns(
                p1_der_freqs, p2_der_freqs, p3_der_freqs,
            )
            # Intialize a results array.
            results = np.array([abba, baba, bbaa, baaa, abaa, aaba])
    return results

# Define a function to filter and calculate site pattern counts for the archaics.
def arc_anc_site_patterns(gt, p1, p2, p3):
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
        var_mask = gt.take(samp_list, axis=1).compress(called_mask, axis=0).count_alleles().is_variant()
        # If there are no variable sites...
        if (var_mask.sum() == 0):
            # Set the results to 0 since we are iterating over QC'ed regions.
            results = np.zeros(6)
        # Else...
        else:
            # Calculate the alternative allele frequencies.
            p1_alt_freqs = calc_alt_freqs(gt.take([idx_dicc[p1]], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p2_alt_freqs = calc_alt_freqs(gt.take([idx_dicc[p2]], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p3_alt_freqs = calc_alt_freqs(gt.take([idx_dicc[p3]], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            anc_freqs = calc_alt_freqs(gt.take([-1], axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            # Polarize the samples.
            p1_der_freqs = np.where(anc_freqs == 1, np.abs(p1_alt_freqs - 1), p1_alt_freqs)
            p2_der_freqs = np.where(anc_freqs == 1, np.abs(p2_alt_freqs - 1), p2_alt_freqs)
            p3_der_freqs = np.where(anc_freqs == 1, np.abs(p3_alt_freqs - 1), p3_alt_freqs)
            # Calculate the site pattern counts.
            abba, baba, bbaa, baaa, abaa, aaba = site_patterns(
                p1_der_freqs, p2_der_freqs, p3_der_freqs,
            )
            # Intialize a results array.
            results = np.array([abba, baba, bbaa, baaa, abaa, aaba])
    return results

# Define a function to load site pattern windowed results.
def load_site_patterns_windows(prefix, p1, p2, p3, window_size=72):
    # Intialize a dictionary to store results.
    results_dicc = {
        'ABBA': np.array([]), 'BABA': np.array([]), 'BBAA': np.array([]),
        'BAAA': np.array([]), 'ABAA': np.array([]), 'AABA': np.array([]),
        'D+': np.array([]),
    }
    # Intialize the file path prefix.
    path_prefix = '../site_patterns/{0}/windows/{1}_{2}_{3}_'.format(
        prefix, p1.lower(), p2.lower(), p3.lower(),
    )
    # For all chromosomes...
    for chromosome in range(1, 23):
        # Load the site pattern results.
        site_patterns = np.loadtxt(
            path_prefix+'chr{0}_{1}kb.csv.gz'.format(chromosome, window_size), delimiter=',',
        )
        # Extract the site pattern arrays.
        abba = site_patterns[:, 0].flatten()
        baba = site_patterns[:, 1].flatten()
        bbaa = site_patterns[:, 2].flatten()
        baaa = site_patterns[:, 3].flatten()
        abaa = site_patterns[:, 4].flatten()
        aaba = site_patterns[:, 5].flatten()
        # Calculate D+.
        dplus = (((abba - baba) + (baaa - abaa)) / ((abba + baba) + (baaa + abaa)))
        # Append the dictionary.
        results_dicc['ABBA'] = np.append(results_dicc['ABBA'], abba)
        results_dicc['BABA'] = np.append(results_dicc['BABA'], baba)
        results_dicc['BBAA'] = np.append(results_dicc['BBAA'], bbaa)
        results_dicc['BAAA'] = np.append(results_dicc['BAAA'], baaa)
        results_dicc['ABAA'] = np.append(results_dicc['ABAA'], abaa)
        results_dicc['AABA'] = np.append(results_dicc['AABA'], aaba)
        results_dicc['D+'] = np.append(results_dicc['D+'], dplus)
    return results_dicc

# Define a funtion to load the tgp window reuslts.
def build_tgp_arc_scenario_window_dicc(
    prefix, p1_list,
    p2_list, p3_list, window_size=72,
):
    # Intialize a dictionary to store all of the results.
    wind_dicc = {}
    # For every p1...
    for p1 in p1_list:
        # For every p2...
        for p2 in p2_list:
            # For every p3...
            for p3 in p3_list:
                # Intialize the template.
                template = '(({0}, {1}), {2})'.format(p1, p2, p3)
                # Load the window results.
                wind_dicc[template] = load_site_patterns_windows(
                    prefix, p1, p2, p3, window_size,
                )
    return wind_dicc

# Define a funtion to load the archaic window reuslts.
def build_arc_anc_scenario_window_dicc(
    prefix, config_list, window_size=72,
):
    # Intialize a dictionary to store all of the results.
    wind_dicc = {}
    # For every configuration...
    for config in config_list:
        # Unpack the configuration.
        p1, p2, p3 = config
        # Intialize the template.
        template = '(({0}, {1}), {2})'.format(p1, p2, p3)
        # Load the window results.
        wind_dicc[template] = load_site_patterns_windows(
                prefix, p1, p2, p3, window_size,
            )
    return wind_dicc

# Define a function to load the tgp scenarios.
def calculate_tgp_arc_scenarios(
    gt, wind_dicc,
    p1_list, p2_list, p3_list, window_size=72,
):
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Intialize an ordered population list.
    tgp_pop_list = [
        'LWK', 'GWD', 'MSL', 'ESN', 'YRI', # AFR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
        'PEL', 'MXL', 'CLM', 'PUR', # AMR.
    ]
    # Intialize a samples dictionary.
    pop_idx_dicc = {
        'DEN': np.array([2350]), 'ALT': np.array([2347]),
        'CHA': np.array([2348]), 'VIN': np.array([2349]),
    }
    # For every population...
    for pop in tgp_pop_list:
        # Fill the sample dictionary.
        pop_idx_dicc[pop] = tgp_df[tgp_df['POP'] == pop].index.values
    # Load the window indicies that passed QC.
    wind_idx = load_esl_qc_windows_idx('tgp', window_size=window_size)
    # Intialize a dictionary to store all of the results.
    muc19_dicc = {}
    # For every p1...
    for p1 in p1_list:
        # For every p2...
        for p2 in p2_list:
            # For every p3...
            for p3 in p3_list:
                # Intialize the template.
                template = '(({0}, {1}), {2})'.format(p1, p2, p3)
                # Calculate the site pattern counts.
                site_patterns = tgp_arc_site_patterns(
                    gt=gt,
                    p1_idx=pop_idx_dicc[p1],
                    p2_idx=pop_idx_dicc[p2],
                    p3_idx=pop_idx_dicc[p3],
                )
                # Calculate D+.
                dplus_num = ((site_patterns[0] - site_patterns[1]) + (site_patterns[3] - site_patterns[4]))
                dplus_den = ((site_patterns[0] + site_patterns[1]) + (site_patterns[3] + site_patterns[4]))
                dplus = (dplus_num / dplus_den)
                # Calculate the p-value from the bootstrapped distribution.
                pval = norm.sf(x=abs(dplus), loc=0, scale=np.nanstd(wind_dicc[template]['D+'][wind_idx]))
                # Append the D+ results to the site pattern results.
                site_patterns = np.append(site_patterns, dplus)
                # Append the p-value to the site pattern results.
                site_patterns = np.append(site_patterns, pval)
                # Fill the dictionary.
                muc19_dicc[template] = site_patterns
    return muc19_dicc

# Define a function to load the archiac scenarios.
def calculate_arc_anc_scenarios(
    gt, config_list, wind_dicc, window_size=72,
):
    # Load the window indicies that passed QC.
    wind_idx = load_esl_qc_windows_idx('arc', window_size=window_size)
    # Intialize a dictionary to store all of the results.
    muc19_dicc = {}
    # For every configuration...
    for config in config_list:
        # Unpack the configuration.
        p1, p2, p3 = config
        # Intialize the template.
        template = '(({0}, {1}), {2})'.format(p1, p2, p3)
        # Calculate the site pattern counts.
        site_patterns = arc_anc_site_patterns(
            gt=gt, p1=p1, p2=p2, p3=p3,
        )
        # Calculate D+.
        dplus_num = ((site_patterns[0] - site_patterns[1]) + (site_patterns[3] - site_patterns[4]))
        dplus_den = ((site_patterns[0] + site_patterns[1]) + (site_patterns[3] + site_patterns[4]))
        dplus = (dplus_num / dplus_den)
        # Calculate the p-value from the bootstrapped distribution.
        pval = norm.sf(x=abs(dplus), loc=0, scale=np.nanstd(wind_dicc[template]['D+'][wind_idx]))
        # Append the D+ results to the site pattern results.
        site_patterns = np.append(site_patterns, dplus)
        # Append the p-value to the site pattern results.
        site_patterns = np.append(site_patterns, pval)
        # Fill the dictionary.
        muc19_dicc[template] = site_patterns
    return muc19_dicc

# Define a function to generate tgp data frames.
def build_tgp_arc_scenario_tables(
    muc19_dicc,
    wind_dicc,
    p1_list,
    p2_list,
    p3_list,
    window_size=72,
    export=True,
):
    # Load the window indicies that passed QC.
    wind_idx = load_esl_qc_windows_idx('tgp', window_size=window_size)
    # Intialize a dataframe dictionary
    results_dicc = {
        'P1': [], 'P2': [], 'P3': [],
        'ABBA': [], 'BABA': [],
        'BAAA' : [], 'ABAA': [],
        'BBAA': [],  'AABA': [],
        r'$D+$': [],
        r'Nonoverlapping Windows $\left( \mu \right)$': [],
        r'Nonoverlapping Windows $\left( \sigma \right)$': [],
        r'$P-value$': [],
        
    }
    # For every p1...
    for p1 in p1_list:
        # For every p2...
        for p2 in p2_list:
            # For every p3...
            for p3 in p3_list:
                # Intialize the dictionaries key.
                key = '(({0}, {1}), {2})'.format(p1, p2, p3)
                # Fill the dataframe dictionary.
                results_dicc['P1'].append(p1)
                results_dicc['P2'].append(p2)
                results_dicc['P3'].append(p3)
                results_dicc['ABBA'].append(muc19_dicc[key][0])
                results_dicc['BABA'].append(muc19_dicc[key][1])
                results_dicc['BBAA'].append(muc19_dicc[key][2])
                results_dicc['BAAA'].append(muc19_dicc[key][3])
                results_dicc['ABAA'].append(muc19_dicc[key][4])
                results_dicc['AABA'].append(muc19_dicc[key][5])
                results_dicc[r'$D+$'].append(muc19_dicc[key][-2])
                results_dicc[r'Nonoverlapping Windows $\left( \mu \right)$'].append(np.nanmean(wind_dicc[key]['D+'][wind_idx]))
                results_dicc[r'Nonoverlapping Windows $\left( \sigma \right)$'].append(np.nanstd(wind_dicc[key]['D+'][wind_idx]))
                results_dicc[r'$P-value$'].append(muc19_dicc[key][-1])
    # Compile the dataframes
    results_df = pd.DataFrame(data=results_dicc)
    # If export is true...
    if export:
        # Export the dataframe as a csv.
        results_df.to_csv(f'./dataframes/p1_yri_p2_ooa_p3_{p3_list[0].lower()}_{window_size}kb.csv', index=False)
    return results_df

# Define a function to generate tgp data frames.
def build_arc_anc_scenario_tables(
    muc19_dicc,
    wind_dicc,
    config_list,
    window_size=72,
    export=True,
):
    # Load the window indicies that passed QC.
    wind_idx = load_esl_qc_windows_idx('arc', window_size=window_size)
    # Intialize a dataframe dictionary
    results_dicc = {
        'P1': [], 'P2': [], 'P3': [],
        'ABBA': [], 'BABA': [],
        'BAAA' : [], 'ABAA': [],
        'BBAA': [],  'AABA': [],
        r'$D+$': [],
        r'Nonoverlapping Windows $\left( \mu \right)$': [],
        r'Nonoverlapping Windows $\left( \sigma \right)$': [],
        r'$P-value$': [],
        
    }
    # For every configuration...
    for config in config_list:
        # Unpack the configuration.
        p1, p2, p3 = config
        # Intialize the dictionaries key.
        key = '(({0}, {1}), {2})'.format(p1, p2, p3)
        # Fill the dataframe dictionary.
        results_dicc['P1'].append(p1)
        results_dicc['P2'].append(p2)
        results_dicc['P3'].append(p3)
        results_dicc['ABBA'].append(muc19_dicc[key][0])
        results_dicc['BABA'].append(muc19_dicc[key][1])
        results_dicc['BBAA'].append(muc19_dicc[key][2])
        results_dicc['BAAA'].append(muc19_dicc[key][3])
        results_dicc['ABAA'].append(muc19_dicc[key][4])
        results_dicc['AABA'].append(muc19_dicc[key][5])
        results_dicc[r'$D+$'].append(muc19_dicc[key][-2])
        results_dicc[r'Nonoverlapping Windows $\left( \mu \right)$'].append(np.nanmean(wind_dicc[key]['D+'][wind_idx]))
        results_dicc[r'Nonoverlapping Windows $\left( \sigma \right)$'].append(np.nanstd(wind_dicc[key]['D+'][wind_idx]))
        results_dicc[r'$P-value$'].append(muc19_dicc[key][-1])
    # Compile the dataframes
    results_df = pd.DataFrame(data=results_dicc)
    # If export is true...
    if export & (p1 == 'ALT'):
        # Export the dataframe as a csv.
        results_df.to_csv(f'./dataframes/p1_alt_p2_nea_p3_den_{window_size}kb.csv', index=False)
    elif export & (p1 == 'CHA'):
        # Export the dataframe as a csv.
        results_df.to_csv(f'./dataframes/p1_cha_p2_vin_p3_arc_{window_size}kb.csv', index=False)
    return results_df

def plot_tgp_arc_scenarios(p1, p3, wind_dicc, muc19_dicc, obs_label, window_size=72, export=True):
    # Intialize a list of p3 populations.
    p2_list = [
        'PEL', 'MXL', 'CLM', 'PUR', 'NA', # AMR.
        'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
    ]
    # Intialize the template.
    template = '(({0},'.format(p1)
    # Intialize a list of archaics and there indicies.
    arc_dicc = {
        'DEN': 'Denisovan', 'ALT': 'Altai Nean.',
        'CHA': 'Chagyrskaya Nean.', 'VIN': 'Vindija Nean.',
    }
    # Intialize a row and columns list.
    rows = [
        0, 0, 0, 0, 0,
        1, 1, 1, 1, 1,
        2, 2, 2, 2, 2,
        3, 3, 3, 3, 3,
    ]
    cols = [
        0, 1, 2, 3, 4,
        0, 1, 2, 3, 4,
        0, 1, 2, 3, 4,
        0, 1, 2, 3, 4,
    ]
    # Load the window indicies that passed QC.
    wind_idx = load_esl_qc_windows_idx('tgp', window_size=window_size)
    # Determine the significance threshold.
    sig_threshold = 0.05 / 19
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        4, 5, figsize=(8, 4), dpi=300,
        sharex=True, sharey=True, 
    )
    # For every p2 population.
    for idx in range(len(p2_list)):
        # If we are in the last column in the first row...
        if idx == 4:
            # Turn off the axes.
            axes[0, 4].axis('off')
        # Else...
        else:
            # Identify the p2 population and fill in the template.
            key = template+' {0}), {1})'.format(p2_list[idx], p3)
            # Plot the distribution of D+ values.
            axes[rows[idx], cols[idx]].hist(
                wind_dicc[key]['D+'][wind_idx], np.arange(-1, 1, 0.01),
                histtype='stepfilled', color='tab:blue',
            )
            # Plot the observed value for muc19.
            axes[rows[idx], cols[idx]].axvline(
                muc19_dicc[key][-2], 0, 1,
                color='black', linestyle='dashed', linewidth=1,
            )
            # Grab and round the p-value.
            wind_p = muc19_dicc[key][-1]
            # If the p-value is significant after correcting for multiple comparisons...
            if wind_p < sig_threshold:
                # Construct the title.
                title = 'P2 = {0}'.format(p2_list[idx])+'\n'+r'$P-value=$'+'{:.3e}'.format(wind_p)+r'$^{*}$'
            # Else...
            else:
                # Construct the title.
                title = 'P2 = {0}'.format(p2_list[idx])+'\n'+r'$P-value=$'+'{:.3e}'.format(wind_p)+r'$^{ns}$'
            # Add the subplot title.
            axes[rows[idx], cols[idx]].set_title(title)
    # Configure the legend.
    legend_elements = [
        Line2D([0], [0], color='black', linestyle='dashed', label=obs_label),
    ]
    # Add a figure legend.
    fig.legend(
        handles=legend_elements, loc='center left',
        bbox_to_anchor=(1.0, 0.5), frameon=False,
    )
    # Label the super-axes.
    fig.supylabel('Frequency')
    fig.supxlabel(r'$D+$'+f'(({p1}, P2), {arc_dicc[p3]}) in {window_size}kb Windows')
    # If export is set to true...
    if export:
        # Export the plot.
        plt.savefig(
            f'./supp_figures/p1_yri_p2_ooa_p3_{p3.lower()}_{window_size}kb.png', format='png',
            facecolor='white', bbox_inches='tight', dpi=500,
        )
    plt.show()
    return

def plot_arc_anc_scenarios(p3_list, wind_dicc, muc19_dicc, obs_label, window_size=72, export=True):
    # If the Altain Nean. is in the p3 list...
    if 'ALT' in p3_list:
        # Intialize a list of configurations.
        config_list = [
            ('CHA', 'VIN', 'DEN'),
            ('CHA', 'VIN', 'ALT'),
        ]
    # Else....
    else:
        # Intialize a list of configurations.
        config_list = [
            ('ALT', 'CHA', 'DEN'),
            ('ALT', 'VIN', 'DEN'),
        ]
    # Intialize a row and columns list.
    cols = [0, 1]
    # Load the window indicies that passed QC.
    wind_idx = load_esl_qc_windows_idx('arc', window_size=window_size)
    # Determine the significance threshold.
    sig_threshold = 0.05 / 2
    # Intialize figures and axes.
    fig, axes = plt.subplots(
        1, 2, figsize=(6, 3), dpi=300,
        sharex=True, sharey=True, 
    )
    # For every configuration
    for idx in range(len(config_list)):
        # Unpack the configuration.
        p1, p2, p3 = config_list[idx]
        # Intialize the dictionaries key.
        key = '(({0}, {1}), {2})'.format(p1, p2, p3)
        # Plot the distribution of D+ values.
        axes[cols[idx]].hist(
            wind_dicc[key]['D+'][wind_idx], bins=np.arange(-1, 1, 0.01),
            histtype='stepfilled', color='tab:blue',
        )
        # Plot the observed value for muc19.
        axes[cols[idx]].axvline(
            muc19_dicc[key][-2], 0, 1,
            color='black', linestyle='dashed',
        )
        # Grab and round the p-value.
        wind_p = muc19_dicc[key][-1]
        # If the p-value is significant after correcting for multiple comparisons...
        if wind_p < sig_threshold:
            # Construct the title.
            title = key+'\n'+r'$P-value=$'+'{:.3e}'.format(wind_p)+r'$^{*}$'
        # Else...
        else:
            # Construct the title.
            title = key+'\n'+r'$P-value=$'+'{:.3e}'.format(wind_p)+r'$^{ns}$'
        # Add the subplot title.
        axes[cols[idx]].set_title(title)
    # Configure the legend.
    legend_elements = [
        Line2D([0], [0], color='black', linestyle='dashed', label=obs_label),
    ]
    # Add a figure legend.
    fig.legend(
        handles=legend_elements, loc='center left',
        bbox_to_anchor=(1.0, 0.5), frameon=False,
    )
    # Label the super-axes.
    fig.supylabel('Frequency')
    fig.supxlabel(r'$D+$'+f'((P1, P2), P3) in {window_size}kb Windows')
    # If export is set to true...
    if export & (p1 == 'ALT'):
        # Export the plot.
        plt.savefig(
            f'./supp_figures/p1_alt_p2_nea_p3_den_{window_size}kb.png', format='png',
            facecolor='white', bbox_inches='tight', dpi=500,
        )
    elif export & (p1 == 'CHA'):
        # Export the plot.
        plt.savefig(
            f'./supp_figures/p1_cha_p2_vin_p3_arc_{window_size}kb.png', format='png',
            facecolor='white', bbox_inches='tight', dpi=500,
        )
    # Show the plot!
    plt.show()
    return