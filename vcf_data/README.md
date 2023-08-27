# `vcf_data`

This directory contains the code to replicate the VCF data used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`BCFtoools v1.13`](https://samtools.github.io/bcftools/bcftools.html)
- [`Tabix v0.2.6`](http://www.htslib.org/doc/tabix.html)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Data

All data is publicly available and can be downloaded from the following locations:

- [EPO Ancestral Allele Calls](http://ftp.ensembl.org/pub/release-74/fasta/ancestral_alleles/)
- [Denisovan Genome](http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Denisova/)
- [Altai Neanderthal Genome](http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Altai/)
- [Chagyrskaya Neandertal Genome](http://ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF/)
- [Vindija Neanderthal Genome](http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/)
- [Simons Genome Diversity Project Phased Dataset](https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data2021/)
- [Phase 3 Release of the 1000 Genomes Project](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)

## Code

### Archaic Dataset

__Merge the unfiltered data.__

```bash
# Intialize file paths.
ALT=./Altai
CHA=./Chagyrskaya
VIN=./Vindija33.19
DEN=./Denisova

# Merge all the VCF files.
for CHR in {1..22}; do
bcftools merge -Oz -o all_archaics_merged_raw_chr${CHR}.vcf.gz ${ALT}/chr${CHR}_mq25_mapab100.vcf.gz ${CHA}/chr${CHR}.noRB.vcf.gz ${VIN}/chr${CHR}_mq25_mapab100.vcf.gz ${DEN}/chr${CHR}_mq25_mapab100.vcf.gz
done
```

__Filter the merged data.__

```bash
# Filter the merged VCF files.
for CHR in {1..22}; do
python filter_all_archaics_merged_raw_vcfs.py all_archaics_merged_raw_chr${CHR}.vcf.gz ${CHR} | bgzip > all_archaics_merged_filtered_all_sites_no_aa_chr${CHR}.vcf.gz
done
```

__Remove monomorphic sites.__

```bash
# Filter out monomorphic sites.
for CHR in {1..22}; do
bcftools view -m2 -M2 -Oz -o all_archaics_merged_filtered_variable_sites_no_aa_chr${CHR}.vcf.gz all_archaics_merged_filtered_all_sites_no_aa_chr${CHR}.vcf.gz
done
```

__Append the ancestral allele calls to the filtered variable sites data.__

```bash
# Add ancestral allele calls.
for CHR in {1..22}; do
python arc_append_ancestral_allele_calls.py all_archaics_merged_filtered_variable_sites_no_aa_chr${CHR}.vcf.gz ${CHR} | bgzip > all_archaics_merged_filtered_variable_sites_aa_calls_chr${CHR}.vcf.gz
done
```

__Compile an all sites report for the merged data.__

```bash
# Create an all sites VCF file and summary report.
for CHR in {1..22}; do
python all_archaics_merged_all_sites.py all_archaics_merged_filtered_all_sites_no_aa_chr${CHR}.vcf.gz ${CHR} | bgzip > all_archaics_merged_filtered_all_sites_aa_calls_chr${CHR}.vcf.gz
done
```

__Index the filtered variable sites data.__

```bash
# Index all variable VCF files.
for CHR in {1..22}; do
tabix -p vcf all_archaics_merged_filtered_variable_sites_aa_calls_chr${CHR}.vcf.gz
done
```

__Encode the filtered variable sites VCF files as `zarr` arrays to minimize file storage.__

```bash
# Convert the VCF files to zarr arrays.
for CHR in {1..22}; do
python vcf_to_zarr.py all_archaics_merged_filtered_variable_sites_aa_calls_chr${CHR} arc_anc_chr${CHR} ${CHR}
done
```

### Simons Genome Diversity Project + Archaic Dataset

__Merge the unfiltered data.__

```bash
# Intialize file paths.
ALT=./Altai
CHA=./Chagyrskaya
VIN=./Vindija33.19
DEN=./Denisova
SGDP=./sgdp

# Merge all the VCF files.
for CHR in {1..22}; do
bcftools merge -Oz -o sgdp_all_archaics_merged_raw_chr${CHR}.vcf.gz ${SGDP}/sgdp.phased.unfiltered.chr${CHR}.vcf.gz ${ALT}/chr${CHR}_mq25_mapab100.vcf.gz ${CHA}/chr${CHR}.noRB.vcf.gz ${VIN}/chr${CHR}_mq25_mapab100.vcf.gz ${DEN}/chr${CHR}_mq25_mapab100.vcf.gz
done
```

__Filter the merged data.__

```bash
# Filter the merged VCF files.
for CHR in {1..22}; do
python filter_sgdp_all_archaics_merged_raw_vcfs.py sgdp_all_archaics_merged_raw_chr${CHR}.vcf.gz ${CHR} | bgzip > sgdp_all_archaics_merged_filtered_all_sites_no_aa_chr${CHR}.vcf.gz
done
```

__Remove monomorphic sites.__

```bash
# Filter out monomorphic sites.
for CHR in {1..22}; do
bcftools view -m2 -M2 -Oz -o sgdp_all_archaics_merged_filtered_variable_sites_no_aa_chr${CHR}.vcf.gz sgdp_all_archaics_merged_filtered_all_sites_no_aa_chr${CHR}.vcf.gz
done
```

__Append the ancestral allele calls to the filtered variable sites data.__

```bash
# Add ancestral allele calls.
for CHR in {1..22}; do
python sgdp_append_ancestral_allele_calls.py sgdp_all_archaics_merged_filtered_variable_sites_no_aa_chr${CHR}.vcf.gz ${CHR} | bgzip > sgdp_all_archaics_merged_filtered_variable_sites_aa_calls_chr${CHR}.vcf.gz
done
```

__Compile an all sites report for the merged data.__

```bash
# Create an all sites VCF file and summary report.
for CHR in {1..22}; do
python sgdp_all_archaics_merged_all_sites.py sgdp_all_archaics_merged_filtered_all_sites_no_aa_chr${CHR}.vcf.gz ${CHR} | bgzip > sgdp_all_archaics_merged_filtered_all_sites_aa_calls_chr${CHR}.vcf.gz
done
```

__Index the filtered variable sites data.__

```bash
# Index all variable VCF files.
for CHR in {1..22}; do
tabix -p vcf sgdp_all_archaics_merged_filtered_variable_sites_aa_calls_chr${CHR}.vcf.gz
done
```

__Encode the filtered variable sites VCF files as `zarr` arrays to minimize file storage.__

```bash
# Convert the VCF files to zarr arrays.
for CHR in {1..22}; do
python vcf_to_zarr.py sgdp_all_archaics_merged_filtered_variable_sites_aa_calls_chr${CHR} sgdp_arc_anc_chr${CHR} ${CHR}
done
```

### 1000 Genomes Project + Archaic Dataset

__Merge the unfiltered data.__

```bash
# Intialize file paths.
ALT=./Altai
CHA=./Chagyrskaya
VIN=./Vindija33.19
DEN=./Denisova
TGP=./tgp

# Merge all the VCF files.
for CHR in {1..22}; do
bcftools merge -Oz -o tgp_all_archaics_merged_raw_chr${CHR}.vcf.gz ${TGP}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${ALT}/chr${CHR}_mq25_mapab100.vcf.gz ${CHA}/chr${CHR}.noRB.vcf.gz ${VIN}/chr${CHR}_mq25_mapab100.vcf.gz ${DEN}/chr${CHR}_mq25_mapab100.vcf.gz
done
```

__Filter the merged data.__

```bash
# Filter the merged VCF files.
for CHR in {1..22}; do
python filter_tgp_all_archaics_merged_raw_vcfs.py tgp_all_archaics_merged_raw_chr${CHR}.vcf.gz ${CHR} | bgzip > tgp_all_archaics_merged_filtered_all_sites_no_aa_chr${CHR}.vcf.gz
done
```

__Remove monomorphic sites.__

```bash
# Filter out monomorphic sites.
for CHR in {1..22}; do
bcftools view -m2 -M2 -Oz -o tgp_all_archaics_merged_filtered_variable_sites_chr${CHR}.vcf.gz tgp_all_archaics_merged_filtered_all_sites_chr${CHR}.vcf.gz
done
```

__Append the ancestral allele calls to the filtered variable sites data.__

```bash
# Add ancestral allele calls.
for CHR in {1..22}; do
python tgp_append_ancestral_allele_calls.py tgp_all_archaics_merged_filtered_variable_sites_no_aa_chr${CHR}.vcf.gz ${CHR} | bgzip > tgp_all_archaics_merged_filtered_variable_sites_aa_calls_chr${CHR}.vcf.gz
done
```

__Compile an all sites report for the merged data.__

```bash
# Create an all sites VCF file and summary report.
for CHR in {1..22}; do
python tgp_all_archaics_merged_all_sites.py tgp_all_archaics_merged_filtered_all_sites_no_aa_chr${CHR}.vcf.gz ${CHR} | bgzip > tgp_all_archaics_merged_filtered_all_sites_aa_calls_chr${CHR}.vcf.gz
done
```

__Index the filtered variable sites data.__

```bash
# Index all variable VCF files.
for CHR in {1..22}; do
tabix -p vcf tgp_all_archaics_merged_filtered_variable_sites_aa_calls_chr${CHR}.vcf.gz
done
```

__Remove individuals from the ASW and ACB populations.__

```bash
# Subset all the TGP data to not include ASW or ACB.
for CHR in {1..22}; do
bcftools view -S ../meta_data/tgp_mod_arc_anc_samps.txt -Oz -o tgp_mod_all_archaics_merged_filtered_variable_sites_aa_calls_chr${CHR}.vcf.gz tgp_all_archaics_merged_filtered_variable_sites_aa_calls_chr${CHR}.vcf.gz
done
```

__Remove duplicate positions.__

Due to the way the 1000 Genomes Projects encodes different variant types there are duplicated positions that need to be removed. Note that to replicate you will need to first decompress the `../meta_data/tgp_snp_dups` directory as GitHub only allows 100 Mb of storage.

```bash 
# Filter out duplicate records.
for CHR in {1..22}; do
bcftools view -T ^../meta_data/tgp_snp_dups/tgp_dup_snps_chr${CHR}.txt -Oz -o tgp_mod_all_archaics_merged_filtered_variable_sites_dedupped_aa_calls_chr${CHR}.vcf.gz tgp_mod_all_archaics_merged_filtered_variable_sites_aa_calls_chr${CHR}.vcf.gz
done
```

__Exclude duplicated positions from the all sites report.__

```bash
# Remove duplicated records from the all sites report.
for CHR in {1..22}; do
python dedup_all_sites_report.py ${CHR}
done
```

__Index the modified filtered variable sites data.__

```bash
# Index all the modified variable VCF files.
for CHR in {1..22}; do
tabix -p vcf tgp_mod_all_archaics_merged_filtered_variable_sites_dedupped_aa_calls_chr${CHR}.vcf.gz
done
```

__Encode the filtered variable sites VCF files as `zarr` arrays to minimize file storage.__

```bash
# Convert the VCF files to zarr arrays.
for CHR in {1..22}; do
python vcf_to_zarr.py tgp_mod_all_archaics_merged_filtered_variable_sites_dedupped_aa_calls_chr${CHR} tgp_mod_arc_anc_chr${CHR} ${CHR}
done
```

