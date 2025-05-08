# Dryad Data

**Citation:** Please cite [Villanea and Peede et al. 2024](https://doi.org/10.1101/2023.09.25.559202) (doi: https://doi.org/10.1101/2023.09.25.559202) when using this data.

To use the data, you will first need to extract the individual data archives (`.tar.gz`). For example, run `tar -xf {data_archive}.tar.gz`. After extraction, you may need to edit the code on [GitHub](https://github.com/David-Peede/MUC19) so that the paths are appropriate for the computer on which the code is run. The paths provided in this README reflect my original project directory structure.

Below is a quick reference of the naming conventions I used throughout the project.

**Naming Conventions**
- Archaic Individuals
	- `den`: Altai Denisovan
	- `alt`: Altai Neanderthal
	- `cha`: Chagyrskaya Neanderthal
		- This individual is considered a late Neanderthal.
	- `vin`: Vindija Neanderthal
		- This individual is considered a late Neanderthal.
- Modern Human Datasets
	- `tgp`: 1000 Genomes Project
	- `sgdp`: Simon's Genome Diversity Project
- 1000 Genomes Project Superpopulations & Populations
	- `afr`: Africans
		- Populations: `yri`, `lwk`, `gwd`, `msl`, `esn`
	- `amr`: Admixed Americans
		- Populations: `mxl`, `pel`, `clm`, `pur`
	- `sas`: South Asians
		- Populations: `beb`, `stu`, `itu`, `pjl`, `gih`
	- `eas`: East Asians
		- Populations: `chb`, `khv`, `chs`, `jpt`, `cdx`
	- `eur`: Europeans
		- Populations: `tsi`, `ceu`, `ibs`, `gbr`, `fin`

## Data Description

**Note on Paths:** The file paths listed below correspond to the directory structure used in the original study. After extracting the archives, you may need to organize the files accordingly or adjust paths in the analysis scripts.

### `amr_lai`

#### `72kb_amr_beds.tar.gz` & `short_read_repeat_amr_beds.tar.gz`
Both `72kb_amr_beds.tar.gz` and `short_read_repeat_amr_beds.tar.gz` contain BED files, which are the output from the `bedtools` code found [here](https://github.com/David-Peede/MUC19/tree/main/amr_lai).

**Paths**
- `72kb_amr_beds.tar.gz`
	- `./muc19/amr_lai/region_beds/72kb/{amr_ind}_{A,B}.bed`
- `short_read_repeat_amr_beds.tar.gz`
	- `./muc19/amr_lai/region_beds/short_read_repeat/{amr_ind}_{A,B}.bed`

##### `{amr_ind}_{A,B}.bed`
Each admixed American individual (`{amr_ind}`) has two corresponding BED files, one per haplotype (`{A,B}`). Each BED file has the following columns:
- Columns 1-6: Original columns from `RFmix`.
- Columns 7-9: Coordinates for the intersecting region of interest.
	- The BED files for the intersecting region of interest can be found [here](https://github.com/David-Peede/MUC19/tree/main/amr_lai/region_beds).
- Column 10: Number of base pairs of overlap with respect to the intersecting region of interest.

### `arc_snp_density`

#### `classify_tgp_snps_chromosome.tar.gz`
`classify_tgp_snps_chromosome.tar.gz` is the output from `classify_tgp_snps_chromosome_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/arc_snp_density). The original path is `./muc19/muc19_results/tgp_arcs_masked_no_aa/{non_afr_pop,tgp}_{snp_type}_chr{1..22}.txt.gz`.

##### ‌`{non_afr_pop,tgp}_{snp_type}_chr{1..22}.txt.gz`
For each non-African population (`{non_afr_pop}`) and the entire 1000 Genomes Project (`{tgp}`), each file has one row listing the positions for the specified SNP classification (`{snp_type}`), stratified by autosome (`chr{1..22}`).

#### `tgp_archaic_snp_denisty_windows.tar.gz`
`tgp_archaic_snp_denisty_windows.tar.gz` is the output from `tgp_archaic_snp_denisty_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/arc_snp_density). The original path is `./muc19/muc19_results/tgp_arcs_masked_no_aa/{non_afr_pop}_arc_snp_denisty_chr{1..22}_{742,72}kb.txt.gz`.

##### `{non_afr_pop}_arc_snp_denisty_chr{1..22}_{742,72}kb.txt.gz`
For each non-African population (`{non_afr_pop}`), each row contains archaic SNP partition counts for a non-overlapping 742kb or 72kb window (`{742,72}kb`), stratified by autosome (`chr{1..22}`), with the following columns:
- Column 1: Denisovan-specific SNPs.
- Column 2: Neanderthal-specific SNPs.
- Column 3: Shared archaic SNPs.
- Column 4: Archaic SNPs.

### `heterozygosity`

#### `{den,alt,cha,vin}_heterozygosity_chromosome.tar.gz`
`den_heterozygosity_chromosome.tar.gz`, `alt_heterozygosity_chromosome.tar.gz`, `cha_heterozygosity_chromosome.tar.gz`, and `vin_heterozygosity_chromosome.tar.gz` are the output from `archaic_heterozygosity_chromosome_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/heterozygosity).

**Paths**
- `den_heterozygosity_chromosome.tar.gz`
	- `./muc19/muc19_results/den_masked_no_aa/archaic_het_sites_heterozygosity_eff_seq_len_chr{1..22}.txt.gz`
- `alt_heterozygosity_chromosome.tar.gz`
	- `./muc19/muc19_results/alt_masked_no_aa/archaic_het_sites_heterozygosity_eff_seq_len_chr{1..22}.txt.gz`
- `cha_heterozygosity_chromosome.tar.gz`
	- `./muc19/muc19_results/cha_masked_no_aa/archaic_het_sites_heterozygosity_eff_seq_len_chr{1..22}.txt.gz`
- `vin_heterozygosity_chromosome.tar.gz`
	- `./muc19/muc19_results/vin_masked_no_aa/archaic_het_sites_heterozygosity_eff_seq_len_chr{1..22}.txt.gz`

##### `archaic_het_sites_heterozygosity_eff_seq_len_chr{1..22}.txt.gz`
For each archaic individual (`{den,alt,cha,vin}`), there is one corresponding file per autosome (`chr{1..22}`) with the following columns:
- Column 1: Number of heterozygous sites.
- Column 2: Gene diversity.
- Column 3: Effective sequence length.

#### `{den,alt,cha,vin}_heterozygosity_windows.tar.gz`
`den_heterozygosity_windows.tar.gz`, `alt_heterozygosity_windows.tar.gz`, `cha_heterozygosity_windows.tar.gz`, and `vin_heterozygosity_windows.tar.gz` are the output from `archaic_heterozygosity_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/heterozygosity).

**Paths**
- `den_heterozygosity_windows.tar.gz`
	- `./muc19/muc19_results/den_masked_no_aa/archaic_het_sites_heterozygosity_chr{1..22}_72kb.txt.gz`
- `alt_heterozygosity_windows.tar.gz`
	- `./muc19/muc19_results/alt_masked_no_aa/archaic_het_sites_heterozygosity_chr{1..22}_72kb.txt.gz`
- `cha_heterozygosity_windows.tar.gz`
	- `./muc19/muc19_results/cha_masked_no_aa/archaic_het_sites_heterozygosity_chr{1..22}_72kb.txt.gz`
- `vin_heterozygosity_windows.tar.gz`
	- `./muc19/muc19_results/vin_masked_no_aa/archaic_het_sites_heterozygosity_chr{1..22}_72kb.txt.gz`

##### `archaic_het_sites_heterozygosity_chr{1..22}_72kb.txt.gz`
For each archaic individual (`{den,alt,cha,vin}`), there is one corresponding file per autosome (`chr{1..22}`). Each row corresponds to a non-overlapping 72kb window and has the following columns:
- Column 1: Number of heterozygous sites.
- Column 2: Gene diversity.

#### `tgp_heterozygosity_chromosome.tar.gz`
`tgp_heterozygosity_chromosome.tar.gz` is the output from `tgp_heterozygosity_chromosome_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/heterozygosity).

**Paths**
- `tgp_heterozygosity_chromosome.tar.gz`
	- `./muc19/muc19_results/tgp_mod_no_aa/tgp_het_sites_chr{1..22}.txt.gz`
	- `./muc19/muc19_results/tgp_mod_no_aa/afr_heterozygosity_eff_seq_len_chr{1..22}.txt.gz`

##### `tgp_het_sites_chr{1..22}.txt.gz`
For each autosome (`chr{1..22}`), the columns correspond to a 1000 Genomes Project individual, in the same order as individuals appear in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data). Each entry is the number of heterozygous sites for that individual in that chromosome.

##### `afr_heterozygosity_eff_seq_len_chr{1..22}.txt.gz`
For the African superpopulation (`afr`) in the 1000 Genomes Project, there is one corresponding file per autosome (`chr{1..22}`) with the following columns:
- Column 1: Gene diversity.
- Column 2: Effective sequence length.

#### `tgp_heterozygosity_windows.tar.gz`
`tgp_heterozygosity_windows.tar.gz` is the output from `tgp_heterozygosity_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/heterozygosity). The original path is `./muc19/muc19_results/tgp_mod_no_aa/{afr,het,hom}_het_sites_chr{1..22}_72kb.txt.gz`.

##### `{afr,het,hom}_het_sites_chr{1..22}_72kb.txt.gz`
For each focal group (`{afr,het,hom}`), there is one corresponding file per autosome (`chr{1..22}`). Each row corresponds to a non-overlapping 72kb window, and the columns correspond to the number of heterozygous sites per 1000 Genomes Project individual. The order of individuals is as specified in the following metadata files:
- `afr`: African individuals.
	- [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)
- `het`: Individuals with one introgressed haplotype at the 72kb region.
	- [`./muc19/meta_data/tgp_het_den_like_ind_idx_72kb.txt.gz`](https://github.com/David-Peede/MUC19/tree/main/meta_data)
- `hom`: Individuals with two introgressed haplotypes at the 72kb region.
	- [`./muc19/meta_data/tgp_hom_den_like_ind_idx_72kb.txt.gz`](https://github.com/David-Peede/MUC19/tree/main/meta_data)

### `iHS`

#### `tgp_selscan_maps.tar.gz` & `tgp_selscan_vcfs.tar.gz`
`tgp_selscan_maps.tar.gz` and `tgp_selscan_vcfs.tar.gz` are the outputs from `make_selscan_population_data_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/iHS).

**Paths**
- `tgp_selscan_maps.tar.gz`
	- `./muc19/iHS/maps/{pop}_selscan_chr{1..22}.map`
- `tgp_selscan_vcfs.tar.gz`
	- `./muc19/iHS/selscan_vcfs/{pop}_selscan_chr{1..22}.vcf`

##### ‌‌`{pop}_selscan_chr{1..22}.map`
For every population in the 1000 Genomes Project (`{pop}`), there is one corresponding map file per autosome (`chr{1..22}`) with the following columns:
- Column 1: Chromosome.
- Column 2: Locus ID.
- Column 3: Genetic position.
- Column 4: Physical position.

##### ‌‌`{pop}_selscan_chr{1..22}.vcf`
For every population in the 1000 Genomes Project (`{pop}`), there is one corresponding `selscan` formatted VCF file per autosome (`chr{1..22}`) with the following columns:
- Column 1: Chromosome.
- Column 2: Physical position.
- Column 3: Locus ID.
- Column 4: Ancestral allele.
- Column 5: Derived allele.
- Column 6: Quality field.
- Column 7: Filter field.
- Column 8: Info field.
- Column 9: Format field.
- Subsequent columns contain genotype information for individuals in the same order as they appear in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data).

#### `tgp_selscan_ihs_output.tar.gz`
`tgp_selscan_ihs_output.tar.gz` is the output from the `selscan` code found [here](https://github.com/David-Peede/MUC19/tree/main/iHS).

**Paths**
- `tgp_selscan_ihs_output.tar.gz`
	- `./muc19/muc19_results/tgp_mod_aa/{pop}_chr{1..22}.ihs.out`
	- `./muc19/muc19_results/tgp_mod_aa/{pop}_chr{1..22}.ihs.out.100bins.norm`

##### `{pop}_chr{1..22}.ihs.out` & `{pop}_chr{1..22}.ihs.out.100bins.norm`
For every population in the 1000 Genomes Project (`{pop}`), there is a corresponding unnormalized (`.ihs.out`) and normalized (`.ihs.out.100bins.norm`) *iHS* file per autosome (`chr{1..22}`). The columns are described in the `selscan` documentation, which you can find [here](https://github.com/szpiech/selscan).

#### `tgp_ihs_tables.tar.gz`
`tgp_ihs_tables.tar.gz` is the output from `compile_ihs_tables_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/iHS).

**Paths**
- `tgp_ihs_tables.tar.gz`:
	- `./muc19/muc19_results/tgp_mod_aa/{pop}_ihs_genome_wide.csv.gz`
	- `./muc19/muc19_results/tgp_mod_aa/{pop}_ihs_windows_{742,72}kb.csv.gz`

##### `{pop}_ihs_genome_wide.csv.gz`
For every population in the 1000 Genomes Project (`{pop}`), this file contains genome-wide *iHS* data with the following columns:
- `ID`: Locus ID.
- `CHR`: Chromosome.
- `POS`: Physical position.
- `DAF`: Derived allele frequency.
- `IHH1`: Integrated haplotype homozygosity for the derived allele.
- `IHH0`: Integrated haplotype homozygosity for the ancestral allele.
- `U_IHS`: Unnormalized integrated haplotype score.
- `N_IHS`: Normalized integrated haplotype score.
- `CRIT`: Is the absolute value of the normalized *iHS* greater than two? (numeric, e.g., 1 if true, 0 if false).
- `DEN`: Is this a Denisovan-specific SNP? (boolean).
- `NEA`: Is this a Neanderthal-specific SNP? (boolean).
- `SHR`: Is this a shared archaic SNP? (boolean).
- `ARC`: Is this an archaic SNP? (boolean).
- `HOM`: Is this a shared Hominin SNP? (boolean).
- `HUM`: Is this a human-specific SNP? (boolean).
- `ALL_CRIT`: Is the absolute value of the normalized *iHS* greater than two? (boolean).
- `DEN_CRIT`: Is this a Denisovan-specific SNP AND is the absolute value of the normalized *iHS* greater than two? (boolean).
- `NEA_CRIT`: Is this a Neanderthal-specific SNP AND is the absolute value of the normalized *iHS* greater than two? (boolean).
- `SHR_CRIT`: Is this a shared archaic SNP AND is the absolute value of the normalized *iHS* greater than two? (boolean).
- `ARC_CRIT`: Is this an archaic SNP AND is the absolute value of the normalized *iHS* greater than two? (boolean).
- `HOM_CRIT`: Is this a shared Hominin SNP AND is the absolute value of the normalized *iHS* greater than two? (boolean).
- `HUM_CRIT`: Is this a human-specific SNP AND is the absolute value of the normalized *iHS* greater than two? (boolean).

##### `{pop}_ihs_windows_{742,72}kb.csv.gz`
For every population in the 1000 Genomes Project (`{pop}`), this file contains windowed *iHS* data. Each row corresponds to a non-overlapping window (`{742,72}kb`) with the following columns:
- `CHR`: Chromosome.
- `START`: Start position (inclusive).
- `STOP`: Stop position (exclusive).
- `N_ALL_SNPS`: Number of SNPs.
- `N_DEN_SNPS`: Number of Denisovan-specific SNPs.
- `N_NEA_SNPS`: Number of Neanderthal-specific SNPs.
- `N_SHR_SNPS`: Number of shared archaic SNPs.
- `N_ARC_SNPS`: Number of archaic SNPs.
- `N_HOM_SNPS`: Number of shared Hominin SNPs.
- `N_HUM_SNPS`: Number of human-specific SNPs.
- `N_ALL_CRIT`: Number of SNPs with an absolute value of the normalized *iHS* greater than two.
- `PROP_ALL_CRIT`: Proportion of SNPs with an absolute value of the normalized *iHS* greater than two.
- `N_DEN_CRIT`: Number of Denisovan-specific SNPs with an absolute value of the normalized *iHS* greater than two.
- `PROP_DEN_CRIT`: Proportion of Denisovan-specific SNPs with an absolute value of the normalized *iHS* greater than two.
- `N_NEA_CRIT`: Number of Neanderthal-specific SNPs with an absolute value of the normalized *iHS* greater than two.
- `PROP_NEA_CRIT`: Proportion of Neanderthal-specific SNPs with an absolute value of the normalized *iHS* greater than two.
- `N_SHR_CRIT`: Number of shared archaic SNPs with an absolute value of the normalized *iHS* greater than two.
- `PROP_SHR_CRIT`: Proportion of shared archaic SNPs with an absolute value of the normalized *iHS* greater than two.
- `N_ARC_CRIT`: Number of archaic SNPs with an absolute value of the normalized *iHS* greater than two.
- `PROP_ARC_CRIT`: Proportion of archaic SNPs with an absolute value of the normalized *iHS* greater than two.
- `N_HOM_CRIT`: Number of shared Hominin SNPs with an absolute value of the normalized *iHS* greater than two.
- `PROP_HOM_CRIT`: Proportion of shared Hominin SNPs with an absolute value of the normalized *iHS* greater than two.
- `N_HUM_CRIT`: Number of human-specific SNPs with an absolute value of the normalized *iHS* greater than two.
- `PROP_HUM_CRIT`: Proportion of human-specific SNPs with an absolute value of the normalized *iHS* greater than two.

### `introgression`

#### `archaic_site_patterns_windows.tar.gz`
`archaic_site_patterns_windows.tar.gz` is the output from `archaic_site_patterns_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/introgression).

**Paths**
- `archaic_site_patterns_windows.tar.gz`
	- `./muc19/muc19_results/arcs_masked_aa/alt_{cha,vin}_den_chr{1..22}_72kb.txt.gz`

##### `alt_{cha,vin}_den_chr{1..22}_72kb.txt.gz`
For each late Neanderthal (`{cha,vin}`), there is one corresponding file per autosome (`chr{1..22}`). Each row corresponds to a non-overlapping 72kb window with the following columns:
- Column 1: Number of ABBA sites.
- Column 2: Number of BABA sites.
- Column 3: Number of BBAA sites.
- Column 4: Number of BAAA sites.
- Column 5: Number of ABAA sites.
- Column 6: Number of AABA sites.

#### `yri_na19664_archaic_site_patterns_windows.tar.gz`
`yri_na19664_archaic_site_patterns_windows.tar.gz` is the output from `mxl_archaic_site_patterns_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/introgression).

**Paths**
- `yri_na19664_archaic_site_patterns_windows.tar.gz`:
	- `./muc19/muc19_results/tgp_den_masked_aa/yri_na19664_den_chr{1..22}_{742,72}kb.txt.gz`
	- `./muc19/muc19_results/tgp_alt_masked_aa/yri_na19664_alt_chr{1..22}_{742,72}kb.txt.gz`
	- `./muc19/muc19_results/tgp_cha_masked_aa/yri_na19664_cha_chr{1..22}_{742,72}kb.txt.gz`
	- `./muc19/muc19_results/tgp_vin_masked_aa/yri_na19664_vin_chr{1..22}_{742,72}kb.txt.gz`

##### `yri_na19664_{den,alt,cha,vin}_chr{1..22}_{742,72}kb.txt.gz`
For each archaic individual (`{den,alt,cha,vin}`), there is one corresponding file per autosome (`chr{1..22}`). Each row corresponds to a non-overlapping window (`{742,72}kb`) with the following columns:
- Column 1: Number of ABBA sites.
- Column 2: Number of BABA sites.
- Column 3: Number of BBAA sites.
- Column 4: Number of BAAA sites.
- Column 5: Number of ABAA sites.
- Column 6: Number of AABA sites.

#### `tgp_q95_u30_afr_b_den_windows.tar.gz`
`tgp_q95_u30_afr_b_den_windows.tar.gz` is the output from `tgp_q95_u30_afr_b_den_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/introgression).

**Paths**
- `tgp_q95_u30_afr_b_den_windows.tar.gz`
	- `./muc19/muc19_results/tgp_den_masked_no_aa/q95_u20_u30_afr_{non_afr_pop}_den_chr{1..22}_{742,72}kb.txt.gz`

##### `q95_u20_u30_afr_{non_afr_pop}_den_chr{1..22}_{742,72}kb.txt.gz`
For each non-African population (`{non_afr_pop}`), there is one corresponding file per autosome (`chr{1..22}`). Each row corresponds to a non-overlapping window (`{742,72}kb`) with the following columns:
- Column 1: The *Q95* value.
- Column 2: The *U20* value.
- Column 3: The *U30* value.

#### `tgp_u30_afr_b_den_windows.tar.gz`
`tgp_u30_afr_b_den_windows.tar.gz` is the output from `tgp_u30_afr_b_den_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/introgression).

**Paths**
- `tgp_u30_afr_b_den_windows.tar.gz`
	- `./muc19/muc19_results/tgp_den_masked_no_aa/u30_afr_b_den_chr{1..22}_{742,72}kb.txt.gz`

##### `u30_afr_b_den_chr{1..22}_{742,72}kb.txt.gz`
For each autosome (`chr{1..22}`), there is one corresponding file. Each row corresponds to a non-overlapping window (`{742,72}kb`), and the columns correspond to the non-African population used to compute *U30* in the following order: MXL, PEL, CLM, PUR, BEB, STU, ITU, PJL, GIH, CHB, KHV, CHS, JPT, CDX, TSI, CEU, IBS, GBR, FIN.

### `pbs`

#### `amr_asn_eur_pbs_windows.tar.gz`
`amr_asn_eur_pbs_windows.tar.gz` is the output from `amr_asn_eur_pbs_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/pbs).

**Paths**
- `amr_asn_eur_pbs_windows.tar.gz`
	- `./muc19/muc19_results/tgp_mod_no_aa/{amr_pop}_{asn_pop}_{eur_pop}_chr{1..22}_{742,72}kb.txt.gz`

##### `{amr_pop}_{asn_pop}_{eur_pop}_chr{1..22}_{742,72}kb.txt.gz`
*PBS(A, B, C)* values for all unique combinations of *A* = AMR population (`{amr_pop}`), *B*  = EAS/SAS population (`{asn_pop}`), and *C* = EUR population (`{eur_pop}`). There is one corresponding file per autosome (`chr{1..22}`). Each file contains a single row, where each column in that row represents the *PBS(A,B,C)* value for a non-overlapping window (`{742,72}kb`).

#### `mxl_chb_ceu_pbs_chromsome.tar.gz`
`mxl_chb_ceu_pbs_chromsome.tar.gz` is the output from `mxl_chb_ceu_pbs_chromsome_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/pbs).

**Paths**
- `mxl_chb_ceu_pbs_chromsome.tar.gz`
	- `./muc19/muc19_results/tgp_mod_no_aa/mxl_chb_ceu_pbs_partitions_chr{1..22}.txt.gz`

##### `mxl_chb_ceu_pbs_partitions_chr{1..22}.txt.gz`
For each autosome (`chr{1..22}`), there is one corresponding file. Each row corresponds to a genomic position and has the following columns:
- Column 1: *PBS* using all MXL individuals.
- Column 2: *PBS* using only MXL individuals with more than 50% Indigenous American ancestry.
- Column 3: *PBS* using only MXL individuals with less than 50% Indigenous American ancestry.

#### `mxl_chb_ceu_pbs_windows.tar.gz`
`mxl_chb_ceu_pbs_windows.tar.gz` is the output from `mxl_chb_ceu_pbs_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/pbs).

**Paths**
- `mxl_chb_ceu_pbs_windows.tar.gz`
	- `./muc19/muc19_results/tgp_mod_no_aa/mxl_chb_ceu_pbs_partitions_chr{1..22}_{742,72}kb.txt.gz`

##### `mxl_chb_ceu_pbs_partitions_chr{1..22}_{742,72}kb.txt.gz`
For each autosome (`chr{1..22}`), there is one corresponding file. Each row corresponds to a non-overlapping window (`{742,72}kb`) with the following columns:
- Column 1: *PBS* using all MXL individuals.
- Column 2: *PBS* using only MXL individuals with more than 50% Indigenous American ancestry.
- Column 3: *PBS* using only MXL individuals with less than 50% Indigenous American ancestry.

#### `sprime_sites_mxl_chb_ceu_pbs_windows.tar.gz`
`sprime_sites_mxl_chb_ceu_pbs_windows.tar.gz` is the output from `sprime_sites_mxl_chb_ceu_pbs_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/pbs).

**Paths**
- `sprime_sites_mxl_chb_ceu_pbs_windows.tar.gz`
	- `./muc19/muc19_results/tgp_mod_no_aa/mxl_chb_ceu_pbs_sprime_all_arc_chr{1..22}_{742,72}kb.txt.gz`

##### `mxl_chb_ceu_pbs_sprime_all_arc_chr{1..22}_{742,72}kb.txt.gz`
For each autosome (`chr{1..22}`), there is one corresponding file. Each row corresponds to a non-overlapping window (`{742,72}kb`) with the following columns:
- Column 1: *PBS* for `SPrime` sites in MXL that match with either the Altai Denisovan or Altai Neanderthal.
- Column 2: Number of `SPrime` sites in MXL that match with either the Altai Denisovan or Altai Neanderthal.

### `pbs_sims`

#### `{neutral,negative,positive_neutral_s1,positive_neutral_s01,positive_neutral_s0015}_recomb_map_slimulated_vcfs.tar.gz`
These archives (listed below) contain the VCF files from the SLiMulations. The `SLiM` scripts can be found [here](https://github.com/David-Peede/MUC19/tree/main/pbs_sims).
- `neutral_recomb_map_slimulated_vcfs.tar.gz` is output from `sim_sneutral_rmap.slim`.
- `negative_recomb_map_slimulated_vcfs.tar.gz` is output from `sim_snegative_rmap.slim`.
- `positive_neutral_s1_recomb_map_slimulated_vcfs.tar.gz` is output from `sim_spositive_neutral_1_rmap.slim`.
- `positive_neutral_s01_recomb_map_slimulated_vcfs.tar.gz` is output from `sim_spositive_neutral_01_rmap.slim`.
- `positive_neutral_s0015_recomb_map_slimulated_vcfs.tar.gz` is output from `sim_spositive_neutral_0015_rmap.slim`.

**Paths**
- `neutral_recomb_map_slimulated_vcfs.tar.gz`
	- `./muc19/pbs_sims/smodel/neutral/recomb_map/VCF_sneutral_{sim_id}_{rep_id}.txt.gz`
- `negative_recomb_map_slimulated_vcfs.tar.gz`
	- `./muc19/pbs_sims/smodel/negative/recomb_map/VCF_snegative_{sim_id}_{rep_id}.txt.gz`
- `positive_neutral_s1_recomb_map_slimulated_vcfs.tar.gz`
	- `./muc19/pbs_sims/smodel/positive_neutral_1/VCF_spositive_neutral_{sim_id}_{rep_id}.txt`
- `positive_neutral_s01_recomb_map_slimulated_vcfs.tar.gz`
	- `./muc19/pbs_sims/smodel/positive_neutral_01/VCF_spositive_neutral_{sim_id}_{rep_id}.txt`
- `positive_neutral_s0015_recomb_map_slimulated_vcfs.tar.gz`
	- `./muc19/pbs_sims/smodel/positive_neutral_0015/VCF_spositive_neutral_{sim_id}_{rep_id}.txt`

##### `VCF_{sneutral,snegative,spositive_neutral}_{sim_id}_{rep_id}.{txt,txt.gz}`
There are 10,000 simulated replicates (`{sim_id}_{rep_id}`), each as a bgzipped VCF file (`.txt.gz`), for the neutral (`{sneutral}`) and heterosis (`{snegative}`) SLiMulations. There are 1,000 simulated replicates, each as a VCF file (`.txt`), for the positive selection (`{spositive_neutral}`) SLiMulations. In each VCF file, the first 99 samples are simulated European individuals, the next 103 samples are simulated East Asian individuals, and the last 64 samples are simulated Mexican individuals. The VCF format specification can be found [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

#### `mxl_slimulations.tar.gz`
`mxl_slimulations.tar.gz` is the output from `extract_neutral_negative_recomb_map_slimulated_results.py` and `extract_positive_slimulations_info.py` found [here](https://github.com/David-Peede/MUC19/tree/main/pbs_sims).

**Paths**
- `mxl_slimulations.tar.gz`
	- `./muc19/muc19_results/mxl_slimulations/{neutral,negative}_per_10k_replicates.csv.gz`
	- `./muc19/muc19_results/mxl_slimulations/{neutral,negative}_{all,arc}_snp_freqs_{742,72}kb_per_snp.txt.gz`
	- `./muc19/muc19_results/mxl_slimulations/{neutral,negative}_pbs_{all,arc}_snps_{742,72}kb_per_snp.txt.gz`
	- `./muc19/muc19_results/mxl_slimulations/positive_s{s1,s01,s0015}_per_1k_replicates.csv.gz`

##### `{neutral,negative}_per_10k_replicates.csv.gz`
For each selection scenario (`{neutral,negative}`), there is one corresponding file. Each row corresponds to one simulated replicate, with the following columns:
- `geq_all_snps_742kb`: Number of SNPs in MXL segregating at a frequency greater than or equal to 0.3 in the 742kb region.
- `pbs_all_snps_742kb`: *PBS(A=MXL,B=EAS,C=EUR)* value for the 742kb region using all SNPs.
- `prop_out_pbs_all_snps_742kb`: Proportion of per-SNP *PBS(A=MXL,B=EAS,C=EUR)* values greater than the empirical 99.95th percentile for the 742kb region.
- `geq_arc_snps_742kb`: Number of archaic SNPs in MXL segregating at a frequency greater than or equal to 0.3 in the 742kb region.
- `pbs_arc_snps_742kb`: *PBS(A=MXL,B=EAS,C=EUR)* value for the 742kb region using archaic SNPs.
- `prop_out_pbs_arc_snps_742kb`: Proportion of per-archaic SNP *PBS(A=MXL,B=EAS,C=EUR)* values greater than the empirical 99.95th percentile for the 742kb region.
- `geq_all_snps_72kb`: Number of SNPs in MXL segregating at a frequency greater than or equal to 0.3 in the 72kb region.
- `pbs_all_snps_72kb`: *PBS(A=MXL,B=EAS,C=EUR)* value for the 72kb region using all SNPs.
- `prop_out_pbs_all_snps_72kb`: Proportion of per-SNP *PBS(A=MXL,B=EAS,C=EUR)* values greater than the empirical 99.95th percentile for the 72kb region.
- `geq_arc_snps_72kb`: Number of archaic SNPs in MXL segregating at a frequency greater than or equal to 0.3 in the 72kb region.
- `pbs_arc_snps_72kb`: *PBS(A=MXL,B=EAS,C=EUR)* value for the 72kb region using archaic SNPs.
- `prop_out_pbs_arc_snps_72kb`: Proportion of per-archaic SNP *PBS(A=MXL,B=EAS,C=EUR)* values greater than the empirical 99.95th percentile for the 72kb region.

##### `{neutral,negative}_{all,arc}_snp_freqs_{742,72}kb_per_snp.txt.gz`
For each selection scenario (`{neutral,negative}`), there is one corresponding file for all (`{all}`) and archaic (`{arc}`) SNPs per focal region (`{742,72}kb`). Each file has one row, where the values represent the allele frequency for sites segregating in MXL across all 10,000 simulated replicates.

##### `{neutral,negative}_pbs_{all,arc}_snps_{742,72}kb_per_snp.txt.gz`
For each selection scenario (`{neutral,negative}`), there is one corresponding file for all (`{all}`) and archaic (`{arc}`) SNPs per focal region (`{742,72}kb`). Each file has one row, where the values represent  the per-SNP *PBS(A=MXL,B=EAS,C=EUR)* values across all 10,000 simulated replicates.

##### `positive_s{s1,s01,s0015}_per_1k_replicates.csv.gz`
For each selection coefficient (`{s1,s01,s0015}`), there is one corresponding file. Each row corresponds to one simulated replicate, with the following columns:
- `vcf_file`: Corresponding VCF file.
- `rng_seed`: Corresponding SLiM rng seed.
- `org_rep_id`: Arbitrary replicate ID.
- `n_den_742kb`: Number of Denisovan SNPs in the 742kb region.
- `n_nea_742kb`: Number of Neanderthal SNPs in the 742kb region.
- `n_arc_742kb`: Number of archaic SNPs in the 742kb region.
- `n_den_72kb`: Number of Denisovan SNPs in the 72kb region.
- `n_nea_72kb`: Number of Neanderthal SNPs in the 72kb region.
- `n_arc_72kb`: Number of archaic SNPs in the 72kb region.
- `sel_pos`: Position of the selected mutation.
- `sel_origin`: Population origin of the selected mutation.
- `sel_mxl_daf`: Frequency of the selected allele in the sampled MXL individuals.
- `sel_eas_daf`: Frequency of the selected allele in the sampled EAS individuals.
- `sel_eur_daf`: Frequency of the selected allele in the sampled EUR individuals.
- `MXB_before`: Frequency of the mutation in the MXB population before positive selection.
- `EAS_before`: Frequency of the mutation in the EAS population before positive selection.
- `EUR_before`: Frequency of the mutation in the EUR population before positive selection.
- `MXB_after`: Frequency of the selected allele in the MXB population at the end of the SLiMulation.
- `MXL_after`: Frequency of the selected allele in the MXL population at the end of the SLiMulation.
- `EAS_after`: Frequency of the selected allele in the EAS population at the end of the SLiMulation.
- `EUR_after`: Frequency of the selected allele in the EUR population at the end of the SLiMulation.
- `n_seg_all_snps_742kb`: Number of segregating sites in the 742kb region.
- `pbs_all_snps_742kb`: *PBS(A=MXL,B=EAS,C=EUR)* value for the 742kb region using all SNPs.
- `fst_mxl_eas_all_snps_742kb`: *Fst(MXL,EAS)* value for the 742kb region using all SNPs.
- `fst_mxl_eur_all_snps_742kb`: *Fst(MXL,EUR)* value for the 742kb region using all SNPs.
- `fst_eas_eur_all_snps_742kb`: *Fst(EAS,EUR)* value for the 742kb region using all SNPs.
- `n_seg_arc_snps_742kb`: Number of archaic SNPs segregating in the 742kb region.
- `pbs_arc_snps_742kb`: *PBS(A=MXL,B=EAS,C=EUR)* value for the 742kb region using archaic SNPs.
- `fst_mxl_eas_arc_snps_742kb`: *Fst(MXL,EAS)* value for the 742kb region using archaic SNPs.
- `fst_mxl_eur_arc_snps_742kb`: *Fst(MXL,EUR)* value for the 742kb region using archaic SNPs.
- `fst_eas_eur_arc_snps_742kb`: *Fst(EAS,EUR)* value for the 742kb region using archaic SNPs.
- `n_seg_all_snps_72kb`: Number of segregating sites in the 722kb region.
- `pbs_all_snps_72kb`: *PBS(A=MXL,B=EAS,C=EUR)* value for the 72kb region using all SNPs.
- `fst_mxl_eas_all_snps_72kb`: *Fst(MXL,EAS)* value for the 72kb region using all SNPs.
- `fst_mxl_eur_all_snps_72kb`: *Fst(MXL,EUR)* value for the 72kb region using all SNPs.
- `fst_eas_eur_all_snps_72kb`: *Fst(EAS,EUR)* value for the 72kb region using all SNPs.
- `n_seg_arc_snps_72kb`: Number of archaic SNPs segregating in the 72kb region.
- `pbs_arc_snps_72kb`: *PBS(A=MXL,B=EAS,C=EUR)* value for the 72kb region using archaic SNPs.
- `fst_mxl_eas_arc_snps_72kb`: *Fst(MXL,EAS)* value for the 72kb region using archaic SNPs.
- `fst_mxl_eur_arc_snps_72kb`: *Fst(MXL,EUR)* value for the 72kb region using archaic SNPs.
- `fst_eas_eur_arc_snps_72kb`: *Fst(EAS,EUR)* value for the 72kb region using archaic SNPs.

### `psuedo_ancestry_painting`

#### `archaic_psuedo_ancestry_painting_windows.tar.gz`
`archaic_psuedo_ancestry_painting_windows.tar.gz` is the output from `archaic_psuedo_ancestry_painting_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/psuedo_ancestry_painting).

**Paths**
- `archaic_psuedo_ancestry_painting_windows.tar.gz`
	- `./muc19/muc19_results/arcs_masked_no_aa/{cha,vin}_den_alt_pap_counts_chr{1..22}_72kb.txt.gz`

##### `{cha,vin}_den_alt_pap_counts_chr{1..22}_72kb.txt.gz`
For each late Neanderthal (`{cha,vin}`), there is one corresponding file per autosome (`chr{1..22}`). Each row corresponds to a non-overlapping 72kb window with the following columns:
- Column 1: Number of *PAP* sites.
- Column 2: Number of heterozygous sites.

#### `tgp_archaic_psuedo_ancestry_painting_windows.tar.gz`
`tgp_archaic_psuedo_ancestry_painting_windows.tar.gz` is the output from `tgp_archaic_psuedo_ancestry_painting_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/psuedo_ancestry_painting).

**Paths**
- `tgp_archaic_psuedo_ancestry_painting_windows.tar.gz`
	- `./muc19/muc19_results/tgp_arcs_masked_no_aa/{den,alt,cha,vin}_mxl_yri_pap_counts_chr{1..22}_72kb.txt.gz`

##### `{cha,vin}_den_alt_pap_counts_chr{1..22}_72kb.txt.gz`
For each archaic individual (`{den,alt,cha,vin}`), there is one corresponding file per autosome (`chr{1..22}`). Each row corresponds to a non-overlapping 72kb window with the following columns:
- Column 1: Number of *PAP* sites.
- Column 2: Number of heterozygous sites.

### `sequence_divergence`

#### `den_v_alt_divergence_windows.tar.gz`
`den_v_alt_divergence_windows.tar.gz` is the output from `arc_dip_v_arc_dip_divergence_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/sequence_divergence).

**Paths**
- `den_v_alt_divergence_windows.tar.gz`
	- `./muc19/muc19_results/arcs_masked_no_aa/den_v_alt_pw_diffs_chr{1..22}_72kb.txt.gz`

##### `den_v_alt_pw_diffs_chr{1..22}_72kb.txt.gz`
For each autosome (`chr{1..22}`), there is one corresponding file. This file contains a single row, where each column in that row corresponds to the average number of pairwise differences between the Altai Denisovan's and Altai Neanderthal's two chromosomes per non-overlapping 72kb window.

#### `sgdp_denisovan_sequence_divergence_at_denisovan_intro_tracts_in_papuans.tar.gz`
`sgdp_denisovan_sequence_divergence_at_denisovan_intro_tracts_in_papuans.tar.gz` is the output from `sgdp_denisovan_sequence_divergence_at_denisovan_intro_tracts_in_papuans_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/sequence_divergence).

**Paths**
- `sgdp_denisovan_sequence_divergence_at_denisovan_intro_tracts_in_papuans.tar.gz`
	- `./muc19/muc19_results/sgdp_den_masked_no_aa/{pap_ind}_hap{1,2}_den_pw_diffs_seq_div_at_den_intro_tracts_chr{1..22}.txt.gz`

##### `{pap_ind}_hap{1,2}_den_pw_diffs_seq_div_at_den_intro_tracts_chr{1..22}.txt.gz`
For each Papuan individual (`{pap_ind}`), there is one corresponding file per autosome (`chr{1..22}`) for each haploid genome (`hap{1,2}`). Each row corresponds to a Denisovan introgressed tract for the corresponding haploid genome (`hap{1,2}`) in the respective Papuan individual (`{pap_ind}`) with the following columns:
- Column 1: The average number of pairwise differences between the Papuan individual's haplotype and the Altai Denisovan's two chromosomes.
- Column 2: The corresponding sequence divergence.

#### `tgp_hap_v_arc_dip_divergence_windows.tar.gz`
`tgp_hap_v_arc_dip_divergence_windows.tar.gz` is the output from `tgp_hap_v_arc_dip_divergence_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/sequence_divergence).

**Paths**
- `tgp_hap_v_arc_dip_divergence_windows.tar.gz`
	- `./muc19/muc19_results/tgp_den_masked_no_aa/den_hap_{1,2}_pw_diffs_chr{1..22}_{742,72}kb.txt.gz`
	- `./muc19/muc19_results/tgp_alt_masked_no_aa/alt_hap_{1,2}_pw_diffs_chr{1..22}_{742,72}kb.txt.gz`
	- `./muc19/muc19_results/tgp_cha_masked_no_aa/cha_hap_{1,2}_pw_diffs_chr{1..22}_{742,72}kb.txt.gz`
	- `./muc19/muc19_results/tgp_vin_masked_no_aa/vin_hap_{1,2}_pw_diffs_chr{1..22}_{742,72}kb.txt.gz`

##### `{den,alt,cha,vin}_hap_{1,2}_pw_diffs_chr{1..22}_{742,72}kb.txt.gz`
For each archaic individual (`{den,alt,cha,vin}`), there is one corresponding file per autosome (`chr{1..22}`) for each haploid genome (`hap_{1,2}`). Each row corresponds to a non-overlapping window (`{742,72}kb`). The columns correspond to the average number of pairwise differences between a 1000 Genomes Project individual's haplotype (`hap_{1,2}`) and the archaic individual's (`{den,alt,cha,vin}`) two chromosomes, in the same order as individuals appear in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data).

#### `tgp_hap_v_arc_dip_divergence_windows.tar.gz`
`tgp_hap_v_arc_dip_divergence_windows.tar.gz` is the output from `tgp_hap_v_arc_psuedo_hap_divergence_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/sequence_divergence).

**Paths**
- `tgp_hap_v_arc_psuedo_hap_divergence_windows.tar.gz`
	- `./muc19/muc19_results/tgp_cha_masked_no_aa/cha_hap_{1,2}_psuedo_hap_diffs_chr{1..22}_72kb.txt.gz`
	- `./muc19/muc19_results/tgp_vin_masked_no_aa/vin_hap_{1,2}_psuedo_hap_diffs_chr{1..22}_72kb.txt.gz`

##### `{cha,vin}_hap_{1,2}_psuedo_hap_diffs_chr{1..22}_72kb.txt.gz`
For each late Neanderthal (`{cha,vin}`), there is one corresponding file per autosome (`chr{1..22}`) for each haploid genome (`hap_{1,2}`). Each row corresponds to a non-overlapping 72kb window. The columns correspond to the average number of pairwise differences between a 1000 Genomes Project individual's haplotype (`hap_{1,2}`) and the late Neanderthal's (`{cha,vin}`) pseudo-haplotype, in the same order as individuals appear in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data).

#### `afr_v_den_divergence_windows.tar.gz`
`afr_v_den_divergence_windows.tar.gz` is the output from `tgp_spop_v_arc_dip_divergence_windows_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/sequence_divergence).

**Paths**
- `afr_v_den_divergence_windows.tar.gz`
	- `./muc19/muc19_results/tgp_den_masked_no_aa/afr_den_avg_pw_diffs_chr{1..22}_72kb.txt.gz`

##### `afr_den_avg_pw_diffs_chr{1..22}_72kb.txt.gz`
For each autosome (`chr{1..22}`), there is one corresponding file. This file contains one row, where each column corresponds to the average number of pairwise differences between all African chromosomes in the 1000 Genomes Project and the Altai Denisovan's two chromosomes per non-overlapping 72kb window.

### `vcf_data`
Since all of the VCF files and bookkeeping information are well over 3TB, we provide the converted Zarr arrays, which are the output of `vcf_to_zarr.py` found [here](https://github.com/David-Peede/MUC19/tree/main/vcf_data). Additionally, in the following `windowing` subsubsection, we provide the corresponding bookkeeping files used in our analyses. The `{prefix}.tar.gz` files are named by the dataset (`{prefix}`) and, unless otherwise noted, contain the corresponding Zarr arrays for all autosomes (`chr{1..22}`).

**Datasets & Paths**
- `arcs_masked_no_aa.tar.gz`
	- All archaic individuals without ancestral allele calls. Individuals are in the following order: Altai Neanderthal, Chagyrskaya Neanderthal, Vindija Neanderthal, Altai Denisovan.
	- `./muc19/zarr_data/arcs_masked_no_aa_chr{1..22}.zarr`
- `arcs_masked_aa.tar.gz`
	- All archaic individuals with ancestral allele calls. Individuals are in the following order: Altai Neanderthal, Chagyrskaya Neanderthal, Vindija Neanderthal, Altai Denisovan, and a placeholder individual "Ancestor" who is always homozygous for the ancestral allele.
	- `./muc19/zarr_data/arcs_masked_no_aa_chr{1..22}.zarr`
- `{den,alt,cha,vin}_masked_no_aa.tar.gz`
	- `{den}`: Altai Denisovan without ancestral allele calls.
		- `./muc19/zarr_data/den_masked_no_aa_chr{1..22}.zarr`
	- `{alt}`: Altai Denisovan without ancestral allele calls.
		- `./muc19/zarr_data/alt_masked_no_aa_chr{1..22}.zarr`
	- `{cha}`: Chagyrskaya Neanderthal without ancestral allele calls.
		- `./muc19/zarr_data/cha_masked_no_aa_chr{1..22}.zarr`
	- `{vin}`: Vindija Neanderthal without ancestral allele calls.
		- `./muc19/zarr_data/vin_masked_no_aa_chr{1..22}.zarr`
- `tgp_arcs_masked_no_aa.tar.gz`
	- 1000 Genomes Project and all archaic individuals without ancestral allele calls. Individuals are in the following order: 1000 Genomes Project individuals (order as in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)), Altai Neanderthal, Chagyrskaya Neanderthal, Vindija Neanderthal, Altai Denisovan.
	- `./muc19/zarr_data/tgp_arcs_masked_no_aa_chr{1..22}.zarr`
- `tgp_{den,alt,cha,vin}_masked_no_aa.tar.gz`
	- `{den}`: 1000 Genomes Project (order as in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)) and the Altai Denisovan, without ancestral allele calls.
		- `./muc19/zarr_data/tgp_den_masked_no_aa_chr{1..22}.zarr`
	- `{alt}`: 1000 Genomes Project (order as in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)) and the Altai Neanderthal, without ancestral allele calls.
		- `./muc19/zarr_data/tgp_alt_masked_no_aa_chr{1..22}.zarr`
	- `{cha}`: 1000 Genomes Project (order as in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)) and the Chagyrskaya Neanderthal, without ancestral allele calls.
		- `./muc19/zarr_data/tgp_cha_masked_no_aa_chr{1..22}.zarr`
	- `{vin}`: 1000 Genomes Project (order as in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)) and the Vindija Neanderthal, without ancestral allele calls.
		- `./muc19/zarr_data/tgp_vin_masked_no_aa_chr{1..22}.zarr`
- `tgp_{den,alt,cha,vin}_masked_aa.tar.gz`
	- `{den}`: 1000 Genomes Project (order as in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)), the Altai Denisovan, and a placeholder individual "Ancestor" who is always homozygous for the ancestral allele.
		- `./muc19/zarr_data/tgp_den_masked_aa_chr{1..22}.zarr`
	- `{alt}`: 1000 Genomes Project (order as in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)), the Altai Neanderthal, and a placeholder individual "Ancestor" who is always homozygous for the ancestral allele.
		- `./muc19/zarr_data/tgp_alt_masked_aa_chr{1..22}.zarr`
	- `{cha}`: 1000 Genomes Project (order as in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)), the Chagyrskaya Neanderthal, and a placeholder individual "Ancestor" who is always homozygous for the ancestral allele.
		- `./muc19/zarr_data/tgp_cha_masked_aa_chr{1..22}.zarr`
	- `{vin}`: 1000 Genomes Project (order as in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)), the Vindija Neanderthal, and a placeholder individual "Ancestor" who is always homozygous for the ancestral allele.
		- `./muc19/zarr_data/tgp_vin_masked_aa_chr{1..22}.zarr`
- `sgdp_arcs_masked_no_aa.tar.gz`
	- Simon's Genome Diversity Project and all archaic individuals without ancestral allele calls. Individuals are in the following order: Simon's Genome Diversity Project individuals (order as in [`./muc19/meta_data/sgdp.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)), Altai Neanderthal, Chagyrskaya Neanderthal, Vindija Neanderthal, Altai Denisovan.
	- `./muc19/zarr_data/sgdp_arcs_masked_no_aa_chr{1..22}.zarr`
- `sgdp_{den,alt,cha,vin}_masked_no_aa.tar.gz`
	- `{den}`: Simon's Genome Diversity Project (order as in [`./muc19/meta_data/sgdp.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)) and the Altai Denisovan, without ancestral allele calls.
		- `./muc19/zarr_data/sgdp_den_masked_no_aa_chr{1..22}.zarr`
	- `{alt}`: Simon's Genome Diversity Project (order as in [`./muc19/meta_data/sgdp.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)) and the Altai Neanderthal, without ancestral allele calls.
		- `./muc19/zarr_data/sgdp_alt_masked_no_aa_chr{1..22}.zarr`
	- `{cha}`: Simon's Genome Diversity Project (order as in [`./muc19/meta_data/sgdp.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)) and the Chagyrskaya Neanderthal, without ancestral allele calls.
		- `./muc19/zarr_data/sgdp_cha_masked_no_aa_chr{1..22}.zarr`
	- `{vin}`: Simon's Genome Diversity Project (order as in [`./muc19/meta_data/sgdp.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data)) and the Vindija Neanderthal, without ancestral allele calls.
		- `./muc19/zarr_data/sgdp_vin_masked_no_aa_chr{1..22}.zarr`
- `tgp_mod_no_aa.tar.gz`
	- 1000 Genomes Project without ancestral allele calls. Individuals are in the same order as they appear in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data).
	- `./muc19/zarr_data/tgp_mod_no_aa_chr{1..22}.zarr`
- `tgp_mod_aa.tar.gz`
	- 1000 Genomes Project with ancestral allele calls. Individuals are in the same order as they appear in [`./muc19/meta_data/tgp_mod.txt`](https://github.com/David-Peede/MUC19/tree/main/meta_data), and a placeholder individual "Ancestor" who is always homozygous for the ancestral allele.
	- `./muc19/zarr_data/tgp_mod_aa_chr{1..22}.zarr`
- `cha_phased_ref_panel_all_inds.tar.gz`
	- Corresponding output from `BEAGLE` before read-based phasing for the focal 72kb region, without ancestral allele calls. Individuals are in the following order: Altai Neanderthal, Chagyrskaya Neanderthal, Altai Denisovan.
	- `./muc19/zarr_data/cha_phased_ref_panel_all_inds.zarr`
		- Note that this is a single Zarr array for the focal 72kb region.
- `vin_phased_ref_panel_all_inds.tar.gz`
	- Corresponding output from `BEAGLE` before read-based phasing for the focal 72kb region, without ancestral allele calls. Individuals are in the following order: Altai Neanderthal, Vindija Neanderthal, Altai Denisovan.
	- `./muc19/zarr_data/vin_phased_ref_panel_all_inds.zarr`
		- Note that this is a single Zarr array for the focal 72kb region.

### `windowing`

#### `arcs_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
`arcs_masked_no_aa_window_info_and_eff_seq_len.tar.gz` is the output from `consolidate_all_archaics_tgp_sgdp_windows_v_revisions.py` and `compute_region_effective_sequence_lengths_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/windowing).

**Paths**
- `arcs_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/arcs_masked_no_aa/72kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/arcs_masked_no_aa/72kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
	- `./muc19/windowing/arcs_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`

##### `72kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
One file for invariant (`_invariant_`) and one for variant (`_variant_`) 72kb non-overlapping windows, with the following columns:
- `IDX`: Window index with respect to the `CHR` column.
- `CHR`: Chromosome.
- `START`: Start position (inclusive).
- `STOP`: Stop position (inclusive).
- `ALT`: Altai Neanderthal effective sequence length.
- `CHA`: Chagyrskaya Neanderthal effective sequence length.
- `VIN`: Vindija Neanderthal effective sequence length.
- `DEN`: Altai Denisovan effective sequence length.
- `DEN-ALT`: Effective sequence length for Altai Denisovan-Altai Neanderthal comparisons.
- `DEN-CHA`: Effective sequence length for Altai Denisovan-Chagyrskaya Neanderthal comparisons.
- `DEN-VIN`: Effective sequence length for Altai Denisovan-Vindija Neanderthal comparisons.
- `ALT-CHA`: Effective sequence length for Altai Neanderthal-Chagyrskaya Neanderthal comparisons.
- `ALT-VIN`: Effective sequence length for Altai Neanderthal-Vindija Neanderthal comparisons.
- `CHA-VIN`: Effective sequence length for Chagyrskaya Neanderthal-Vindija Neanderthal comparisons.
- `S`: Number of segregating sites.
- `QC`: QC condition (numeric, e.g., 1 if true, 0 if false).

**Note**: If only column headers are present in the invariant file (`_invariant_`), no invariant windows passed initial QC.

##### `72kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
Indices of 72kb variant windows of comparable effective sequence length, used for subsetting `72kb_nonoverlapping_variant_windows.csv.gz`.

##### `{742,72}kb_eff_seq_len.txt.gz`
Effective sequence lengths for the focal regions (`{742,72}kb`) with the following columns:
- Column 1: Altai Neanderthal effective sequence length.
- Column 2: Chagyrskaya Neanderthal effective sequence length.
- Column 3: Vindija Neanderthal effective sequence length.
- Column 4: Altai Denisovan effective sequence length.
- Column 5: Effective sequence length for Altai Denisovan-Altai Neanderthal comparisons.
- Column 6: Effective sequence length for Altai Denisovan-Chagyrskaya Neanderthal comparisons.
- Column 7: Effective sequence length for Altai Denisovan-Vindija Neanderthal comparisons.
- Column 8: Effective sequence length for Altai Neanderthal-Chagyrskaya Neanderthal comparisons.
- Column 9: Effective sequence length for Altai Neanderthal-Vindija Neanderthal comparisons.
- Column 10: Effective sequence length for Chagyrskaya Neanderthal-Vindija Neanderthal comparisons.

#### `arcs_masked_aa_window_info_and_eff_seq_len.tar.gz`
`arcs_masked_aa_window_info_and_eff_seq_len.tar.gz` is the output from `consolidate_all_archaics_tgp_sgdp_windows_v_revisions.py` and `compute_region_effective_sequence_lengths_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/windowing).

**Paths**
- `arcs_masked_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/arcs_masked_aa/72kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/arcs_masked_aa/72kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
	- `./muc19/windowing/arcs_masked_aa/{742,72}kb_eff_seq_len.txt.gz`

##### `72kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
One file for invariant (`_invariant_`) and one for variant (`_variant_`) 72kb non-overlapping windows, with the following columns:
- `IDX`: Window index with respect to the `CHR` column.
- `CHR`: Chromosome.
- `START`: Start position (inclusive).
- `STOP`: Stop position (inclusive).
- `ALT`: Altai Neanderthal effective sequence length.
- `CHA`: Chagyrskaya Neanderthal effective sequence length.
- `VIN`: Vindija Neanderthal effective sequence length.
- `DEN`: Altai Denisovan effective sequence length.
- `DEN-ALT`: Effective sequence length for Altai Denisovan-Altai Neanderthal comparisons.
- `DEN-CHA`: Effective sequence length for Altai Denisovan-Chagyrskaya Neanderthal comparisons.
- `DEN-VIN`: Effective sequence length for Altai Denisovan-Vindija Neanderthal comparisons.
- `ALT-CHA`: Effective sequence length for Altai Neanderthal-Chagyrskaya Neanderthal comparisons.
- `ALT-VIN`: Effective sequence length for Altai Neanderthal-Vindija Neanderthal comparisons.
- `CHA-VIN`: Effective sequence length for Chagyrskaya Neanderthal-Vindija Neanderthal comparisons.
- `S`: Number of segregating sites.
- `QC`: QC condition (numeric, e.g., 1 if true, 0 if false).

**Note**: If only column headers are present in the invariant file (`_invariant_`), no invariant windows passed initial QC.

##### `72kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
Indices of 72kb variant windows of comparable effective sequence length, used for subsetting `72kb_nonoverlapping_variant_windows.csv.gz`.

##### `{742,72}kb_eff_seq_len.txt.gz`
Effective sequence lengths for the focal regions (`{742,72}kb`) with the following columns:
- Column 1: Altai Neanderthal effective sequence length.
- Column 2: Chagyrskaya Neanderthal effective sequence length.
- Column 3: Vindija Neanderthal effective sequence length.
- Column 4: Altai Denisovan effective sequence length.
- Column 5: Effective sequence length for Altai Denisovan-Altai Neanderthal comparisons.
- Column 6: Effective sequence length for Altai Denisovan-Chagyrskaya Neanderthal comparisons.
- Column 7: Effective sequence length for Altai Denisovan-Vindija Neanderthal comparisons.
- Column 8: Effective sequence length for Altai Neanderthal-Chagyrskaya Neanderthal comparisons.
- Column 9: Effective sequence length for Altai Neanderthal-Vindija Neanderthal comparisons.
- Column 10: Effective sequence length for Chagyrskaya Neanderthal-Vindija Neanderthal comparisons.

#### `{den,alt,cha,vin}_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
`arcs_masked_no_aa_window_info_and_eff_seq_len.tar.gz` is the output from `consolidate_all_archaics_tgp_sgdp_windows_v_revisions.py` and `compute_region_effective_sequence_lengths_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/windowing).

**Paths**
- `den_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/den_masked_no_aa/72kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/den_masked_no_aa/72kb_esl_qced_nonoverlapping_{invariant,variant}_windows.txt.gz`
	- `./muc19/windowing/den_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`
- `alt_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/alt_masked_no_aa/72kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/alt_masked_no_aa/72kb_esl_qced_nonoverlapping_{invariant,variant}_windows.txt.gz`
	- `./muc19/windowing/alt_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`
- `cha_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/cha_masked_no_aa/72kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/cha_masked_no_aa/72kb_esl_qced_nonoverlapping_{invariant,variant}_windows.txt.gz`
	- `./muc19/windowing/cha_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`
- `vin_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/vin_masked_no_aa/72kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/vin_masked_no_aa/72kb_esl_qced_nonoverlapping_{invariant,variant}_windows.txt.gz`
	- `./muc19/windowing/vin_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`

##### `72kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
One file for invariant (`_invariant_`) and one for variant (`_variant_`) 72kb non-overlapping windows, with the following columns:
- `IDX`: Window index with respect to the `CHR` column.
- `CHR`: Chromosome.
- `START`: Start position (inclusive).
- `STOP`: Stop position (inclusive).
- `{ARC}`: Effective sequence length with respect to the given archaic individual (e.g., `DEN` for Altai Denisovan).
- `S`: Number of segregating sites.
- `QC`: QC condition (numeric, e.g., 1 if true, 0 if false).

**Note**: If only column headers are present in the invariant file (`_invariant_`), no invariant windows passed initial QC.

##### `72kb_esl_qced_nonoverlapping_{invariant,variant}_windows.txt.gz`
Indices of non-overlapping 72kb invariant (`_invariant_`) and variant (`_variant_`) windows of comparable effective sequence length, used for subsetting the corresponding `72kb_nonoverlapping_{invariant,variant}_windows.csv.gz` file.

##### `{742,72}kb_eff_seq_len.txt.gz`
Effective sequence lengths for the focal regions (`{742,72}kb`) with respect to the given archaic individual—singe value per file.

#### `tgp_arcs_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
`tgp_arcs_masked_no_aa_window_info_and_eff_seq_len.tar.gz` is the output from `consolidate_all_archaics_tgp_sgdp_windows_v_revisions.py` and `compute_region_effective_sequence_lengths_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/windowing).

**Paths**
- `tgp_arcs_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/tgp_arcs_masked_no_aa/{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/tgp_arcs_masked_no_aa/{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
	- `./muc19/windowing/tgp_arcs_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`

##### `{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
One file for invariant (`_invariant_`) and one for variant (`_variant_`) per window size (`{742,72}kb`), with the following columns:
- `IDX`: Window index with respect to the `CHR` column.
- `CHR`: Chromosome.
- `START`: Start position (inclusive).
- `STOP`: Stop position (inclusive).
- `ALT`:  Effective sequence length for Altai Neanderthal comparisons.
- `CHA`: Effective sequence length for Chagyrskaya Neanderthal comparisons.
- `VIN`: Effective sequence length for Vindija Neanderthal comparisons.
- `DEN`: Effective sequence length for Altai Denisovan comparisons.
- `DEN-ALT`: Effective sequence length for Altai Denisovan-Altai Neanderthal comparisons.
- `DEN-CHA`: Effective sequence length for Altai Denisovan-Chagyrskaya Neanderthal comparisons.
- `DEN-VIN`: Effective sequence length for Altai Denisovan-Vindija Neanderthal comparisons.
- `ALT-CHA`: Effective sequence length for Altai Neanderthal-Chagyrskaya Neanderthal comparisons.
- `ALT-VIN`: Effective sequence length for Altai Neanderthal-Vindija Neanderthal comparisons.
- `CHA-VIN`: Effective sequence length for Chagyrskaya Neanderthal-Vindija Neanderthal comparisons.
- `S`: Number of segregating sites.
- `QC`: QC condition (numeric, e.g., 1 if true, 0 if false).

**Note**: If only column headers are present in the invariant file (`_invariant_`), no invariant windows passed initial QC.

##### `{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
Indices of variant windows per window size (`{742,72}kb`) of comparable effective sequence length, used for subsetting `{742,72}kb_nonoverlapping_variant_windows.csv.gz`.

##### `{742,72}kb_eff_seq_len.txt.gz`
Effective sequence lengths for the focal regions (`{742,72}kb`) with the following columns:
- Column 1: Effective sequence length for Altai Neanderthal comparisons.
- Column 2: Effective sequence length for Chagyrskaya Neanderthal comparisons.
- Column 3: Effective sequence length for Vindija Neanderthal comparisons.
- Column 4: Effective sequence length for Altai Denisovan comparisons.
- Column 5: Effective sequence length for Altai Denisovan-Altai Neanderthal comparisons.
- Column 6: Effective sequence length for Altai Denisovan-Chagyrskaya Neanderthal comparisons.
- Column 7: Effective sequence length for Altai Denisovan-Vindija Neanderthal comparisons.
- Column 8: Effective sequence length for Altai Neanderthal-Chagyrskaya Neanderthal comparisons.
- Column 9: Effective sequence length for Altai Neanderthal-Vindija Neanderthal comparisons.
- Column 10: Effective sequence length for Chagyrskaya Neanderthal-Vindija Neanderthal comparisons.

#### `tgp_{den,alt,cha,vin}_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
`tgp_{den,alt,cha,vin}_masked_no_aa_window_info_and_eff_seq_len.tar.gz` is the output from `consolidate_all_archaics_tgp_sgdp_windows_v_revisions.py` and `compute_region_effective_sequence_lengths_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/windowing).

**Paths**
- `tgp_den_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/tgp_den_masked_no_aa/{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/tgp_den_masked_no_aa/{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
	- `./muc19/windowing/tgp_den_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`
- `tgp_alt_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/tgp_alt_masked_no_aa/{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/tgp_alt_masked_no_aa/{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
	- `./muc19/windowing/tgp_alt_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`
- `tgp_cha_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/tgp_cha_masked_no_aa/{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/tgp_cha_masked_no_aa/{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
	- `./muc19/windowing/tgp_cha_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`
- `tgp_vin_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/tgp_vin_masked_no_aa/{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/tgp_vin_masked_no_aa/{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
	- `./muc19/windowing/tgp_vin_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`

##### `{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
One file for invariant (`_invariant_`) and one for variant (`_variant_`) per window size (`{742,72}kb`), with the following columns:
- `IDX`: Window index with respect to the `CHR` column.
- `CHR`: Chromosome.
- `START`: Start position (inclusive).
- `STOP`: Stop position (inclusive).
- `{ARC}`: Effective sequence length with respect to the given archaic individual (e.g., `DEN` for Altai Denisovan).
- `S`: Number of segregating sites.
- `QC`: QC condition (numeric, e.g., 1 if true, 0 if false).

**Note**: If only column headers are present in the invariant file (`_invariant_`), no invariant windows passed initial QC.

##### `{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
Indices of variant windows per window size (`{742,72}kb`) of comparable effective sequence length, used for subsetting `{742,72}kb_nonoverlapping_variant_windows.csv.gz`.

##### `{742,72}kb_eff_seq_len.txt.gz`
Effective sequence lengths for the focal regions (`{742,72}kb`) with respect to the given archaic individual—singe value per file.

#### `tgp_{den,alt,cha,vin}_masked_aa_window_info_and_eff_seq_len.tar.gz`
`tgp_{den,alt,cha,vin}_masked_aa_window_info_and_eff_seq_len.tar.gz` is the output from `consolidate_all_archaics_tgp_sgdp_windows_v_revisions.py` and `compute_region_effective_sequence_lengths_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/windowing).

**Paths**
- `tgp_den_masked_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/tgp_den_masked_aa/{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/tgp_den_masked_aa/{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
	- `./muc19/windowing/tgp_den_masked_aa/{742,72}kb_eff_seq_len.txt.gz`
- `tgp_alt_masked_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/tgp_alt_masked_aa/{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/tgp_alt_masked_aa/{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
	- `./muc19/windowing/tgp_alt_masked_aa/{742,72}kb_eff_seq_len.txt.gz`
- `tgp_cha_masked_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/tgp_cha_masked_aa/{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/tgp_cha_masked_aa/{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
	- `./muc19/windowing/tgp_cha_masked_aa/{742,72}kb_eff_seq_len.txt.gz`
- `tgp_vin_masked_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/tgp_vin_masked_aa/{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/tgp_vin_masked_aa/{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
	- `./muc19/windowing/tgp_vin_masked_aa/{742,72}kb_eff_seq_len.txt.gz`

##### `{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
One file for invariant (`_invariant_`) and one for variant (`_variant_`) per window size (`{742,72}kb`), with the following columns:
- `IDX`: Window index with respect to the `CHR` column.
- `CHR`: Chromosome.
- `START`: Start position (inclusive).
- `STOP`: Stop position (inclusive).
- `{ARC}`: Effective sequence length with respect to the given archaic individual (e.g., `DEN` for Altai Denisovan).
- `S`: Number of segregating sites.
- `QC`: QC condition (numeric, e.g., 1 if true, 0 if false).

**Note**: If only column headers are present in the invariant file (`_invariant_`), no invariant windows passed initial QC.

##### `{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
Indices of variant windows per window size (`{742,72}kb`) of comparable effective sequence length, used for subsetting `{742,72}kb_nonoverlapping_variant_windows.csv.gz`.

##### `{742,72}kb_eff_seq_len.txt.gz`
Effective sequence lengths for the focal regions (`{742,72}kb`) with respect to the given archaic individual—singe value per file.

#### `sgdp_{den,alt,cha,vin}_masked_no_aa_eff_seq_len.tar.gz`
`sgdp_{den,alt,cha,vin}_masked_no_aa_window_info_and_eff_seq_len.tar.gz` is the output from `compute_region_effective_sequence_lengths_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/windowing).

**Paths**
- `sgdp_den_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/sgdp_den_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`
- `sgdp_alt_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/sgdp_alt_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`
- `sgdp_cha_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/sgdp_cha_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`
- `sgdp_vin_masked_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/sgdp_vin_masked_no_aa/{742,72}kb_eff_seq_len.txt.gz`

##### `{742,72}kb_eff_seq_len.txt.gz`
Effective sequence lengths for the focal regions (`{742,72}kb`) with respect to the given archaic individual—singe value per file.

#### `tgp_mod_no_aa_window_info_and_eff_seq_len.tar.gz`
`tgp_mod_no_aa_window_info_and_eff_seq_len.tar.gz` is the output from `consolidate_all_archaics_tgp_sgdp_windows_v_revisions.py` and `compute_region_effective_sequence_lengths_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/windowing).

**Paths**
- `tgp_mod_no_aa_window_info_and_eff_seq_len.tar.gz`
	- `./muc19/windowing/tgp_mod_no_aa/{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
	- `./muc19/windowing/tgp_mod_no_aa/{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
	- `./muc19/windowing/tgp_mod_no_aa/{742,72}kb_eff_seq_len.txt.gz`

##### `{742,72}kb_nonoverlapping_{invariant,variant}_windows.csv.gz`
One file for invariant (`_invariant_`) and one for variant (`_variant_`) per window size (`{742,72}kb`), with the following columns:
- `IDX`: Window index with respect to the `CHR` column.
- `CHR`: Chromosome.
- `START`: Start position (inclusive).
- `STOP`: Stop position (inclusive).
- `HUM`: Effective sequence length.
- `S`: Number of segregating sites.
- `QC`: QC condition (numeric, e.g., 1 if true, 0 if false).

**Note**: If only column headers are present in the invariant file (`_invariant_`), no invariant windows passed initial QC.

##### `{742,72}kb_esl_qced_nonoverlapping_variant_windows.txt.gz`
Indices of variant windows per window size (`{742,72}kb`) of comparable effective sequence length, used for subsetting `{742,72}kb_nonoverlapping_variant_windows.csv.gz`.

##### `{742,72}kb_eff_seq_len.txt.gz`
Effective sequence lengths for the focal regions (`{742,72}kb`) with respect to the given archaic individual—singe value per file.

#### `sgdp_denisovan_intro_tracts_in_papuans_qc_info.tar.gz`
`sgdp_denisovan_intro_tracts_in_papuans_qc_info.tar.gz` is the output from `consolidate_sgdp_denisovan_intro_tracts_in_papuans_v_revisions.py` found [here](https://github.com/David-Peede/MUC19/tree/main/windowing).

**Paths**
- `sgdp_denisovan_intro_tracts_in_papuans_qc_info.tar.gz`
	- `./muc19/windowing/sgdp_den_masked_no_aa/{pap_ind}_hap{1,2}_den_intro_{invariant,variant}_tracts.csv.gz`

##### `sgdp_denisovan_intro_tracts_in_papuans_qc_info.tar.gz`
For each Papuan individual (`{pap_ind}`), there is one file for invariant (`_invariant_`) and one for variant (`_variant_`) for each haploid genome (`hap{1,2}`), with the following columns:
- `IDX`: Denisovan introgressed tract index with respect to the `CHR` column.
- `CHR`: Chromosome.
- `START`: Start position (inclusive).
- `STOP`: Stop position (inclusive).
- `DEN`: Effective sequence length for Altai Denisovan comparisons.
- `S`: Number of segregating sites.
- `QC`: QC condition (numeric, e.g., 1 if true, 0 if false). 



























