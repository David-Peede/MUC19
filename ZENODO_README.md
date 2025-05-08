# Zenodo Data

**Citation:** Please cite [Villanea and Peede et al. 2024](https://doi.org/10.1101/2023.09.25.559202) (doi: https://doi.org/10.1101/2023.09.25.559202) when using this data.

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

### `dataset_1.csv.gz`
This file contains the sequence divergence results between the haplotypes in the 1000 Genomes Project and the four high-coverage archaic individuals at the focal 72kb *MUC19* region, with the following columns:
- `Individual`: Individual ID.
- `Super Population`: Superpopulation ID.
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
- `Super Population`: Superpopulation ID.
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
- `Super Population`: Superpopulation ID.
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