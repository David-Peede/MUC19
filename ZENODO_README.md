# Zenodo Data

**Citation:** Please cite [Villanea and Peede et al. 2025](https://doi.org/10.1126/science.adl0882) (https://doi.org/10.1126/science.adl0882) when using this data.

Note that for the archived datasets you will first need to extract the individual data archives (`.tar.gz`). For example, run `tar -xf {data_archive}.tar.gz`.

## Dataset Descriptions

### `papuan_phased_hmmmix_haploid_tracts_per_ind.tar.gz`
`papuan_phased_hmmmix_haploid_tracts_per_ind.tar.gz` contains the unfiltered, decoded, and annotated `hmmix` autosomal tracts for all Papuans in the Simons Genome Diversity Project. 

#### `decoded.{pap_ind}.hap{1,2}.txt`
For each Papuan individual (`{pap_ind}`), there is one corresponding file for each haploid genome (`hap{1,2}`). The columns are described in the `hmmix` documentation, which you can find [here](https://github.com/LauritsSkov/Introgression-detection).

### `tgp_hmmmix_haploid_tracts_per_ind.tar.gz`
`tgp_hmmmix_haploid_tracts_per_ind.tar.gzz` contains the unfiltered, decoded, and annotated `hmmix` autosomal tracts for all non-African individuals in the 1000 Genomes Project.

#### `decoded.{non_afr_ind}.autosome.hap{1,2}.txt`
For each non-African individual (`{pap_ind}`), there is one corresponding file for each haploid genome (`hap{1,2}`). The columns are described in the `hmmix` documentation, which you can find [here](https://github.com/LauritsSkov/Introgression-detection).

### `den_tracts_in_pap.tar.gz`
`den_tracts_in_pap.tar.gz` contains the coordinates for the putative Denisovan tracts for the Papuan individuals.

#### `den_intro_tracts_in_{pap_ind}_hap{1,2}.csv.gz`
For each Papuan individual (`{pap_ind}`), there is one corresponding file for each haploid genome (`hap{1,2}`) with the following columns:
- `chrom`: Chromosome.
- `start`: Start position (inclusive).
- `stop`: Stop position (inclusive).

### `tgp_hmmix_haplotype_tracts_chr12.csv.gz`
This file contains the initially filtered `hmmix` chromosome 12 archaic tracts for all non-African individuals in the 1000 Genomes Project with the following columns:
- `IND`: The non-African individual and corresponding haploid genome (`{non_afr_ind}-hap{1,2}`).
- `CHR`: Chromosome.
- `START`: Start position (inclusive).
- `END`: Stop position (exclusive).
- `LENGTH`: The length of the tract (`END` - `START`).
- `STATE`: Inferred hidden ancestry state from `hmmix`.
- `N_SNPS`: Number of derived SNPs found in `IND` as reported by `hmmix`.
- `N_ARC_SNPS`: Number of derived SNPs found in `IND` that are shared with at least one of the four high-coverage individuals as reported by `hmmix`.
- `POP`: Population ID.
- `SUPERPOP`: Superpopulation ID.
- `PROB`: Posterior probability of the hidden ancestry state from `hmmix`.
- `N_ALT_SNPS`: Number of derived SNPs found in `IND` that are shared with the Altai Neanderthal as reported by `hmmix`.
- `N_DEN_SNPS`: Number of derived SNPs found in `IND` that are shared with the Altai Denisovan as reported by `hmmix`.
- `N_VIN_SNPS`: Number of derived SNPs found in `IND` that are shared with the Vindija Neanderthal as reported by `hmmix`.
- `N_CHA_SNPS`: Number of derived SNPs found in `IND` that are shared with the Chagyrskaya Neanderthal as reported by `hmmix`.

### `tgp_hmmix_haplotype_tracts_muc19.csv.gz`
This file contains the initially filtered `hmmix` archaic tracts that overlap *MUC19* by at least one base pair for all non-African individuals in the 1000 Genomes Project, with the following columns:
- `IND`: The non-African individual and corresponding haploid genome (`{non_afr_ind}-hap{1,2}`).
- `CHR`: Chromosome.
- `START`: Start position (inclusive).
- `END`: Stop position (exclusive).
- `LENGTH`: The length of the tract (`END` - `START`).
- `STATE`: Inferred hidden ancestry state from `hmmix`.
- `N_SNPS`: Number of derived SNPs found in `IND` as reported by `hmmix`.
- `N_ARC_SNPS`: Number of derived SNPs found in `IND` that are shared with at least one of the four high-coverage individuals as reported by `hmmix`.
- `POP`: Population ID.
- `SUPERPOP`: Superpopulation ID.
- `PROB`: Posterior probability of the hidden ancestry state from `hmmix`.
- `N_ALT_SNPS`: Number of derived SNPs found in `IND` that are shared with the Altai Neanderthal as reported by `hmmix`.
- `N_DEN_SNPS`: Number of derived SNPs found in `IND` that are shared with the Altai Denisovan as reported by `hmmix`.
- `N_VIN_SNPS`: Number of derived SNPs found in `IND` that are shared with the Vindija Neanderthal as reported by `hmmix`.
- `N_CHA_SNPS`: Number of derived SNPs found in `IND` that are shared with the Chagyrskaya Neanderthal as reported by `hmmix`.

### `papuan_hmmix_haplotype_tracts.csv.gz`
This file contains the initially filtered `hmmix` autosomal archaic tracts for all Papuan individuals in the Simons Genome Diversity Project, with the following columns:
- `IND`: The Papuan individual and corresponding haploid genome (`{non_afr_ind}-hap{1,2}`).
- `CHR`: Chromosome.
- `START`: Start position (inclusive).
- `END`: Stop position (exclusive).
- `LENGTH`: The length of the tract (`END` - `START`).
- `STATE`: Inferred hidden ancestry state from `hmmix`.
- `N_SNPS`: Number of derived SNPs found in `IND` as reported by `hmmix`.
- `N_ARC_SNPS`: Number of derived SNPs found in `IND` that are shared with at least one of the four high-coverage individuals as reported by `hmmix`.
- `POP`: Population ID.
- `SUPERPOP`: Superpopulation ID.
- `PROB`: Posterior probability of the hidden ancestry state from `hmmix`.
- `N_ALT_SNPS`: Number of derived SNPs found in `IND` that are shared with the Altai Neanderthal as reported by `hmmix`.
- `N_DEN_SNPS`: Number of derived SNPs found in `IND` that are shared with the Altai Denisovan as reported by `hmmix`.
- `N_VIN_SNPS`: Number of derived SNPs found in `IND` that are shared with the Vindija Neanderthal as reported by `hmmix`.
- `N_CHA_SNPS`: Number of derived SNPs found in `IND` that are shared with the Chagyrskaya Neanderthal as reported by `hmmix`.

### `tgp_muc19_short_read_vntr_info.csv.gz`
This file contains the number of *MUC19* VNTR copies estimated from high-coverage short-read data and the number of `hmmix` introgressed tracts overlapping the *MUC19* VNTR region (hg19, Chr12:40876395-40885001) for individuals in the 1000 Genomes Project, with the following columns:
- `Individual`: Individual ID.
- `Population`: Population ID.
- `Super Population`: Super population ID.
- `Repeat Copies`: Number of 30bp repeat units for the corresponding `Individual`.
- `Number of Tracts Overlapping the Repeat Region`: Number of `hmmix` introgressed tracts overlapping the VNTR region (hg19, Chr12:40876395-40885001).
	- Note that we inferred introgressed tracts per haploid genome, so an individual can have at most two introgressed tracts overlapping the VNTR region (i.e., one per haploid genome since humans are diploid).
- `Repeat Copies > 95th Percentile`: Is an individual's `Repeat Copies` value greater than 487 (boolean, e.g., True or False)?

**Note**: The number of *MUC19* VNTR copies estimated from short-read data for the high-coverage archaic individuals is not included in this CSV file, but is reported here for the sake of completeness:
- Altai Denisovan: 296 copies.
- Altai Neanderthal: 379 copies.
- Vindija Neanderthal: 268 copies.
- Chagyrskaya Neanderthal: 293 copies.

### `hprc_hgsv_muc19_long_read_vntr_info.csv.gz`
This file contains the repeat length of the *MUC19* VNTR estimated from high-coverage long-read data for individuals in the Human Pangenome Reference Consortium and Human Genome Structural Variant Consortium, with the following columns:
- `Super Population`: Super population ID.
- `Population`: Population ID.
- `Individual`: Individual ID.
- `Haplotype`: Haploid assembly ID for the corresponding `Individual`.
- `Repeat Length`: The total number of bases across all repeat units for the corresponding `Haplotype`.

**Note**: Individuals from the Human Pangenome Reference Consortium have the `mat` and `pat` midfix, and individuals in the Human Genome Structural Variant Consortium have either the `clr` or `css` midfix in their corresponding haploid assembly ID, respectively.

### `altai_nean_phased_late_neanderthals_denisovan_genos_72kb.csv.gz`
This file contains the phased haplotypes for the late Neanderthals and the unphased genotypes for the Altai Denisovan and Altai Neanderthal for the 1669 variable sites with respect to the 1000 Genomes Project and the four high-coverage archaic individuals at the focal 72kb *MUC19* region, with the following columns:
- `POS`: Genomic coordinate with respect to the hg19 reference assembly.
- `REF`: Reference allele with respect to the hg19 reference assembly.
- `ALT`: Alternative allele with respect to the hg19 reference assembly.
- `Altai Nean.`: Alternative allele frequency with respect to the Altai Neanderthal—i.e., homozygous reference = 0, heterozygous = 0.5, homozygous alternative = 1, failed QC = missing.
- `Chagyrskaya Nean. Hap. 1`: The allele present on the Chagyrskaya Neanderthal's *Neanderthal-like* haplotype—i.e., reference = 0, alternative = 1, failed QC = missing.
- `Chagyrskaya Nean. Hap. 2`: The allele present on the Chagyrskaya Neanderthal's *Denisovan-like* haplotype—i.e., reference = 0, alternative = 1, failed QC = missing.
- `Vindija Nean. Hap. 1`: The allele present on the Vindija Neanderthal's *Neanderthal-like* haplotype—i.e., reference = 0, alternative = 1, failed QC = missing.
- `Vindija Nean. Hap. 2`: The allele present on the Vindija Neanderthal's *Denisovan-like* haplotype—i.e., reference = 0, alternative = 1, failed QC = missing.
- `Denisovan`: Alternative allele frequency with respect to the Altai Denisovan—i.e., homozygous reference = 0, heterozygous = 0.5, homozygous alternative = 1, failed QC = missing.

### `dataset_1.csv.gz`
This file contains the sequence divergence results between the haplotypes in the 1000 Genomes Project and the four high-coverage archaic individuals at the focal 72kb *MUC19* region, with the following columns:
- `Individual`: Individual ID.
- `Super Population`: Super population ID.
- `Population`: Population ID.
- `Archaic`: The archaic individual used to compute sequence divergence.
- `Focal 72kb Region (Pairwise Diffs. Hap. 1)`: The average number of pairwise differences between the `Individual`'s first haplotype and the `Archaic`'s two unphased chromosomes.
- `Focal 72kb Region (Pairwise Diffs. Hap. 2)`: The average number of pairwise differences between the `Individual`'s second haplotype and the `Archaic`'s two unphased chromosomes.
- `Focal 72kb Region (Seq. Div. Hap. 1)`: `Focal 72kb Region (Pairwise Diffs. Hap. 1)` normalized by the effective sequence length for the corresponding `Archaic`.
- `Focal 72kb Region (Seq. Div. Hap. 2)`: `Focal 72kb Region (Pairwise Diffs. Hap. 2)` normalized by the effective sequence length for the corresponding `Archaic`.
- `72kb Non-overlapping Windows $\left( \mu \right)$`: The mean of the `Individual`'s sequence divergence genomic background distribution for the corresponding `Archaic` used to compute the *P-value*.
- `72kb Non-overlapping Windows $\left( \sigma \right)$`: The standard deviation of the `Individual`'s sequence divergence genomic background distribution for the corresponding `Archaic` used to compute the *P-value*.
- `72kb Non-overlapping Windows $\left( SEM \right)$`: The standard error of the mean of the `Individual`'s sequence divergence genomic background distribution for the corresponding `Archaic` used to compute the *P-value*.
- `72kb Non-overlapping Windows $\left( \pm CI_{95\%} \right)$`: The 95% confidence intervals of the `Individual`'s sequence divergence genomic background distribution for the corresponding `Archaic` used to compute the *P-value*.
- `$P-value$ (Hap. 1)`: The proportion of non-overlapping 72kb windows of comparable effective sequence length where the sequence divergence is less than or equal to `Focal 72kb Region (Seq. Div. Hap. 1)`.
	- After correcting for two multiple comparisons—i.e., one per haplotype—a *P-value* less than 0.025 is considered statistically significant.
- `$P-value$ (Hap. 2)`: The proportion of non-overlapping 72kb windows of comparable effective sequence length where the sequence divergence is less than or equal to `Focal 72kb Region (Seq. Div. Hap. 2)`.
	- After correcting for two multiple comparisons—i.e., one per haplotype—a *P-value* less than 0.025 is considered statistically significant.

**Note**: At the focal 72kb region, the effective sequence length with respect to the 1000 Genomes Project is 48435bp, 48654bp, 48121bp, and 48447bp for the Denisovan, Altai Neanderthal, Chagyrskaya Neanderthal, and Vindija Neanderthal, respectively.

### `dataset_2.csv.gz`
This file contains the sequence divergence results between the haplotypes in the 1000 Genomes Project and the four high-coverage archaic individuals at the focal 742kb *MUC19* region, with the following columns:
- `Individual`: Individual ID.
- `Super Population`: Super population ID.
- `Population`: Population ID.
- `Archaic`: The archaic individual used to compute sequence divergence.
- `Focal 742kb Region (Pairwise Diffs. Hap. 1)`: The average number of pairwise differences between the `Individual`'s first haplotype and the `Archaic`'s two unphased chromosomes.
- `Focal 742kb Region (Pairwise Diffs. Hap. 2)`: The average number of pairwise differences between the `Individual`'s second haplotype and the `Archaic`'s two unphased chromosomes.
- `Focal 742kb Region (Seq. Div. Hap. 1)`: `Focal 742kb Region (Pairwise Diffs. Hap. 1)` normalized by the effective sequence length for the corresponding `Archaic`.
- `Focal 742kb Region (Seq. Div. Hap. 2)`: `Focal 742kb Region (Pairwise Diffs. Hap. 2)` normalized by the effective sequence length for the corresponding `Archaic`.
- `742kb Non-overlapping Windows $\left( \mu \right)$`: The mean of the `Individual`'s sequence divergence genomic background distribution for the corresponding `Archaic` used to compute the *P-value*.
- `742kb Non-overlapping Windows $\left( \sigma \right)$`: The standard deviation of the `Individual`'s sequence divergence genomic background distribution for the corresponding `Archaic` used to compute the *P-value*.
- `742kb Non-overlapping Windows $\left( SEM \right)$`: The standard error of the mean of the `Individual`'s sequence divergence genomic background distribution for the corresponding `Archaic` used to compute the *P-value*.
- `742kb Non-overlapping Windows $\left( \pm CI_{95\%} \right)$`: The 95% confidence intervals of the `Individual`'s sequence divergence genomic background distribution for the corresponding `Archaic` used to compute the *P-value*.
- `$P-value$ (Hap. 1)`: The proportion of non-overlapping 742kb windows of comparable effective sequence length where the sequence divergence is less than or equal to `Focal 742kb Region (Seq. Div. Hap. 1)`.
	- After correcting for two multiple comparisons—i.e., one per haplotype—a *P-value* less than 0.025 is considered statistically significant.
- `$P-value$ (Hap. 2)`: The proportion of non-overlapping 742kb windows of comparable effective sequence length where the sequence divergence is less than or equal to `Focal 742kb Region (Seq. Div. Hap. 2)`.
	- After correcting for two multiple comparisons—i.e., one per haplotype—a *P-value* less than 0.025 is considered statistically significant.

**Note**: At the focal 742kb region, the effective sequence length with respect to the 1000 Genomes Project is 48435bp, 48654bp, 48121bp, and 48447bp for the Denisovan, Altai Neanderthal, Chagyrskaya Neanderthal, and Vindija Neanderthal, respectively.

### `dataset_3.csv.gz`
This file contains the sequence divergence results between the haplotypes in the 1000 Genomes Project and the phased late Neanderthal haplotypes at the focal 72kb *MUC19* region, with the following columns:
- `Individual`: Individual ID.
- `Super Population`: Super population ID.
- `Population`: Population ID.
- `Archaic Hap.`: The late Neanderthal haplotype used to compute sequence divergence.
	- The *Denisovan-like* in each of the late Neanderthal's corresponds to their second haplotype (i.e., Chagyrskaya Nean. Hap. 2 and Vindija Nean. Hap. 2).
- `Focal 72kb Region (Pairwise Diffs. Hap. 1)`: The number of pairwise differences between the `Individual`'s first haplotype and the `Archaic Hap.`.
- `Focal 72kb Region (Pairwise Diffs. Hap. 2)`: The number of pairwise differences between the `Individual`'s second haplotype and the `Archaic Hap.`.
- `Focal 72kb Region (Seq. Div. Hap. 1)`: `Focal 72kb Region (Pairwise Diffs. Hap. 1)` normalized by the effective sequence length for the corresponding `Archaic Hap.`.
- `Focal 72kb Region (Seq. Div. Hap. 2)`: `Focal 72kb Region (Pairwise Diffs. Hap. 2)` normalized by the effective sequence length for the corresponding `Archaic Hap.` 
- `72kb Non-overlapping Windows $\left( \mu \right)$`: The mean of the `Individual`'s pseudo-haplotype divergence genomic background distribution for the corresponding late Neanderthal used to compute the *P-value*.
- `72kb Non-overlapping Windows $\left( \sigma \right)$`: The standard deviation of the `Individual`'s pseudo-haplotype divergence genomic background distribution for the corresponding late Neanderthal used to compute the *P-value*.
- `72kb Non-overlapping Windows $\left( SEM \right)$`: The standard error of the `Individual`'s pseudo-haplotype divergence genomic background distribution for the corresponding late Neanderthal used to compute the *P-value*.
- `72kb Non-overlapping Windows $\left( \pm CI_{95\%} \right)$`: The 95% confidence intervals of the `Individual`'s pseudo-haplotype divergence genomic background distribution for the corresponding late Neanderthal used to compute the *P-value*.
- `$P-value$ (Hap. 1)`: The proportion of non-overlapping 72kb windows of comparable effective sequence length where the sequence divergence is less than or equal to `Focal 72kb Region (Seq. Div. Hap. 1)`.
	- After correcting for two multiple comparisons—i.e., two per each modern human haplotype—a *P-value* less than 0.0125 is considered statistically significant.
- `$P-value$ (Hap. 2)`: The proportion of non-overlapping 72kb windows of comparable effective sequence length where the sequence divergence is less than or equal to `Focal 72kb Region (Seq. Div. Hap. 2)`.
	- After correcting for two multiple comparisons—i.e., two per each modern human haplotype—a *P-value* less than 0.0125 is considered statistically significant.

**Note**: At the phased focal 72kb region, the effective sequence length with respect to the 1000 Genomes Project is 48119bp and 48444bp for the  Chagyrskaya and Vindija Neanderthal, respectively.

### `dataset_4.csv.gz`
This file contains the genotype information for the high-coverage long-read sequenced individuals aligned to hg38 used to validate the presence of the 135 Denisvoan-specific alleles identified from short-read data aligned to hg19, with the following columns:
- `Chr12 Position (Hg19)`: Genomic coordinate with respect to the hg19 reference assembly.
- `Ref. Allele (Hg19)`: Reference allele with respect to the hg19 reference assembly.
	- Note that these alleles are identical to those in `Ref. Allele (Hg38)`.
- `Denisovan Allele (Hg19)`: Denisovan allele with respect to the hg19 reference assembly.
	- Note that these alleles are always alternative with respect to the hg19 reference assembly and identical to those in `Denisovan Allele (Hg38)`.
- `Chr12 Position (Hg38)`: Genomic coordinate with respect to the hg38 reference assembly.
- `Ref. Allele (Hg38)`: Reference allele with respect to the hg38 reference assembly.
	- Note that these alleles are identical to those in `Ref. Allele (Hg19)`.
- `Denisovan Allele (Hg38)`: Denisovan allele with respect to the hg38 reference assembly.
	- Note that these alleles are always alternative with respect to the hg38 reference assembly and identical to those in `Denisovan Allele (Hg19)`.
- `HG00864 (PAV)`: Phased genotypes called using the Phased Assembly Variant caller (PAV) for the CDX individual (HG00864) from the Human Genome Structural Variation Consortium.
	- Note that this individual has one *Denisovan-like* haplotype at the 72kb region.
- `HG03009 (PAV)`: Phased genotypes called using the Phased Assembly Variant caller (PAV) for the BEB individual (HG03009) from the Human Genome Structural Variation Consortium.
	- Note that this individual has one *Denisovan-like* haplotype at the 72kb region.
- `NA20847 (PAV)`: Phased genotypes called using the Phased Assembly Variant caller (PAV) for the GIH individual (NA20847) from the Human Genome Structural Variation Consortium.
	- Note that this individual has one *Denisovan-like* haplotype at the 72kb region.
- `HG01122 (Clair3)`: Unphased genotypes called using Clair3 for the CLM individual (HG01122) from the 1000 Genomes Project Oxford Nanopore Technologies Sequencing Consortium.
	- Note that this individual has one *Denisovan-like* haplotype at the 72kb region.
- `HG01122 (PMDV)`: Unphased genotypes called using PEPPER-Margin-DeepVariant (PMDV) for the CLM individual (HG01122) from the 1000 Genomes Project Oxford Nanopore Technologies Sequencing Consortium.
- `HG02252 (Clair3)`: Unphased genotypes called using Clair3 for the PEL individual (HG02252) from the 1000 Genomes Project Oxford Nanopore Technologies Sequencing Consortium.
	- Note that this individual has one *Denisovan-like* haplotype at the 72kb region.
- `HG02252 (PMDV)`: Unphased genotypes called using PEPPER-Margin-DeepVariant (PMDV) for the PEL individual (HG02252) from the 1000 Genomes Project Oxford Nanopore Technologies Sequencing Consortium.
	- Note that this individual has one *Denisovan-like* haplotype at the 72kb region.
- `HG02262 (Clair3)`: Unphased genotypes called using Clair3 for the PEL individual (HG02262) from the 1000 Genomes Project Oxford Nanopore Technologies Sequencing Consortium.
	- Note that this individual has two *Denisovan-like* haplotypes at the 72kb region.
- `HG02262 (PMDV)`: Unphased genotypes called using PEPPER-Margin-DeepVariant (PMDV) for the PEL individual (HG02262) from the 1000 Genomes Project Oxford Nanopore Technologies Sequencing Consortium.
	- Note that this individual has two *Denisovan-like* haplotypes at the 72kb region.












