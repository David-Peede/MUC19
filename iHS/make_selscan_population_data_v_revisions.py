# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

### sys.argv[1] = chromosome ###
### sys.argv[2] = population ###


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

# Define a function to calculate alternative allele frequencies for a single individual.
def calc_ind_alt_freqs(gt):
    # Compute the frequency for each site.
    raw_freqs = np.nansum(gt, axis=2).flatten() / 2
    # Set missing data to Nan
    alt_freqs = np.where(raw_freqs == -1, np.nan, raw_freqs)
    return alt_freqs

# Define a function to generate data for running selscan.
def generate_selscan_pop_data(chrom, pop):
    # Load in the meta information as a pandas dataframe.
    tgp_df = pd.read_csv(
        '../meta_data/tgp_mod.txt', sep='\t',
        names=['IND', 'POP', 'SUPERPOP'],
    )
    # Grab the MXL indicies and individuals
    pop_idx = tgp_df[tgp_df['POP'] == pop.upper()].index.values
    pop_inds = tgp_df[tgp_df['POP'] == pop.upper()]['IND'].values
    # Load the recombination maps.
    hap_map_df = pd.read_csv(
        f'./plink_maps/plink.chr{chrom}.GRCh37.map', sep='\t',
        names=['CHROM', 'ID', 'G_POS', 'P_POS'],
    )
    # Extract the physical and genetic positions.
    hap_p_pos = hap_map_df['P_POS'].values
    hap_g_pos = hap_map_df['G_POS'].values
    # Load the genotype matrix and positions.
    chrom_gt, chrom_pos = load_gt_pos('tgp_mod_aa', chrom)
    # Extract the mxl genotype matrix.
    pop_gt = chrom_gt.take(pop_idx, axis=1)
    # Generate a mask for the segregating sites.
    pop_seg_mask = pop_gt.count_alleles().is_segregating()
    # Subset the positions array.
    pop_seg_pos = chrom_pos[pop_seg_mask]
    # Subset the mxl genotype matrix.
    pop_seg_gt = pop_gt.compress(pop_seg_mask, axis=0)
    # Extract the ancestral sequence.
    anc_seq = calc_ind_alt_freqs(chrom_gt.take([-1], axis=1).compress(pop_seg_mask, axis=0))
    # Determine what positions need to be interpolated.
    inter_pos_hap = np.setdiff1d(pop_seg_pos, hap_p_pos)
    # Create masks for the physical positions.
    hap_inter_mask = np.isin(pop_seg_pos, inter_pos_hap)
    # Create a mask for the positions in the recombination maps.
    pos_in_hap = np.isin(hap_p_pos, pop_seg_pos)
    # Interperlote the genetic positions.
    hap_inter_g_pos = np.interp(inter_pos_hap, hap_p_pos, hap_g_pos)
    # Intialize genetic position arrays.
    pop_hap_g_pos = np.zeros(pop_seg_pos.size)
    # Fill the genetic positions that were already known.
    pop_hap_g_pos[~hap_inter_mask] = hap_g_pos[pos_in_hap]
    # Fill the genetic positions that were interpolated.
    pop_hap_g_pos[hap_inter_mask] = hap_inter_g_pos
    # Generate locus id's.
    loc_ids = [f'{chrom}_{pos}' for pos in pop_seg_pos]
    # Generate the dataframe.
    hap_map = pd.DataFrame({
        'CHROM': np.full(pop_seg_pos.size, chrom),
        'ID': loc_ids,
        'G_POS': pop_hap_g_pos,
        'P_POS': pop_seg_pos,
    })
    # Export the map files.
    hap_map.to_csv(f'./maps/{pop}_selscan_chr{chrom}.map', sep='\t', header=False, index=False)
    # Intialize the first nine header columns.
    header_cols = '\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'])
    # Intialize the mxl sample columns.
    pop_samp_cols = '\t'.join(pop_inds)
    # Complete the selscan vcf header.
    selscan_vcf_header = header_cols+'\t'+pop_samp_cols+'\n'
    # Intialize a buffer list for writting the vcf.
    selscan_vcf_lines = [selscan_vcf_header]
    # For every segregating site amongst mxl.
    for i, pos in enumerate(pop_seg_pos):
        # Intialize the current vcf line.
        c_line = [
            f'{chrom}', f'{pos}', f'{chrom}_{pos}', # CHROM, POS, LOCUS_ID.
            '0', '1', '.', '.', '.', '.',           # ANC, DER, QUAL, FILTER, INFO, FORMAT.
        ]
        # Extract the phased genotypes for this region.
        pop_genos = pop_seg_gt[i, :, :]
        # If the alternative allele is the ancestral allele
        if anc_seq[i] == 1:
            # Polarize the data.
            pop_genos = np.abs(pop_genos - 1)
        # Add the polarized phased haplotypes to the current line.
        c_line.extend([f'{geno[0]}|{geno[1]}' for geno in pop_genos])
        # Update the vcf buffer.
        selscan_vcf_lines.append('\t'.join(c_line)+'\n')
        # If 50_000 lines have been accumulated.
        if len(selscan_vcf_lines) == 50_000:
            # Write the lines to the selscan vcf file in chunks to improve performance.
            sys.stdout.writelines(selscan_vcf_lines)
            # Clear the written vcf lines.
            selscan_vcf_lines.clear()
    # If there are remaining lines to be written.
    if selscan_vcf_lines:
        # Write the remaining qc lines.
        sys.stdout.writelines(selscan_vcf_lines)
    return

# Generate the selscan data.
generate_selscan_pop_data(chrom=int(sys.argv[1]), pop=str(sys.argv[2]))
