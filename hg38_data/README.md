# `hg38_data`

This directory contains the code to replicate the long-read validation data used in this study.

## Data

The data is publicly available and can be downloaded from the following location:

- [Human Genome Structural Variation Consortium Data](https://www.internationalgenome.org/data-portal/data-collection/hgsvc2)
- [1000 Genomes Project Oxford Nanopore Technologies Sequencing Consortium Data](https://s3.amazonaws.com/1000g-ont/index.html?prefix=ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100)

## Code

__Download the HGSV2 data.__
```bash
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_snv_snv_alt.vcf.gz
```

__Download the ONT data using Clair3 to call variants.__
```bash
wget https://s3.amazonaws.com/1000g-ont/ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100/IN-HOUSE_MINIMAP2/HG38/HG02262-ONT-hg38-R9-LSK110-guppy-sup-5mC/HG02262-ONT-hg38-R9-LSK110-guppy-sup-5mC.clair3.notPhased.vcf.gz
wget https://s3.amazonaws.com/1000g-ont/ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100/IN-HOUSE_MINIMAP2/HG38/HG02252-ONT-hg38-R9-LSK110-guppy-sup-5mC/HG02252-ONT-hg38-R9-LSK110-guppy-sup-5mC.clair3.notPhased.vcf.gz
wget https://s3.amazonaws.com/1000g-ont/ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100/IN-HOUSE_MINIMAP2/HG38/HG01122-ONT-hg38-R9-LSK110-guppy-sup-5mC/HG01122-ONT-hg38-R9-LSK110-guppy-sup-5mC.clair3.notPhased.vcf.gz
```

__Download the ONT data using PMDV to call variants.__
```bash
wget https://s3.amazonaws.com/1000g-ont/ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100/NAPU_PIPELINE/HG38/HG02262-ONT-hg38-R9-LSK110-guppy-sup-5mC/HG02262-ONT-hg38-R9-LSK110-guppy-sup-5mC.PMDV_FINAL.vcf.gz
wget https://s3.amazonaws.com/1000g-ont/ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100/NAPU_PIPELINE/HG38/HG02252-ONT-hg38-R9-LSK110-guppy-sup-5mC/HG02252-ONT-hg38-R9-LSK110-guppy-sup-5mC.PMDV_FINAL.vcf.gz
wget https://s3.amazonaws.com/1000g-ont/ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100/NAPU_PIPELINE/HG38/HG01122-ONT-hg38-R9-LSK110-guppy-sup-5mC/HG01122-ONT-hg38-R9-LSK110-guppy-sup-5mC.PMDV_FINAL.vcf.gz
```