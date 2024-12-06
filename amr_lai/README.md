# `amr_lai`

This directory contains the code to replicate the AMR ancestry proportions data used in this study.

## Packages

The package is publicly available and its documentation can be viewed at:

- [`bedtools v2.31.0`](https://bedtools.readthedocs.io/en/latest/)

## Data

The data is publicly available and can be downloaded from the following location:

- [RFMIX Local Ancestry Calls](https://personal.broadinstitute.org/armartin/tgp_admixture)

## Code

__Download the genome-wide ancestry proportions.__
```bash
wget https://personal.broadinstitute.org/armartin/tgp_admixture/AMR_lai.txt
```


__Subset the AMR bed files for the regions of interest.__
```bash
cat tgp_mod_amr_inds.txt | while read IND; do for REGION in short_read_repeat 72kb; do
bedtools intersect -a ./rfmix_beds/${IND}_A.bed -b ./region_beds/${REGION}_region.bed -wo > ./region_beds/${REGION}/${IND}_A.bed
bedtools intersect -a ./rfmix_beds/${IND}_B.bed -b ./region_beds/${REGION}_region.bed -wo > ./region_beds/${REGION}/${IND}_B.bed
done; done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `amr_lai/region_beds/short_read_repeat_amr_beds.tar.gz` and `amr_lai/region_beds/72kb_amr_beds.tar.gz`.