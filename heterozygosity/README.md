# `heterozygosity`

This directory contains the code to replicate the heterozygosity results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__Compute the number of heterozygous sites per TGP individuals per site.__
```bash
for CHR in {1..22}; do
python tgp_heterozygosity_chromosome_v_revisions.py ${CHR}
done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/tgp_mod_no_aa/tgp_heterozygosity_chromosome.tar.gz`.


__Compute the number of heterozygous sites amongst TGP individuals in non-overlapping windows.__
```bash
for CHR in {1..22}; do
python tgp_heterozygosity_windows_v_revisions.py ${CHR} 72
done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/tgp_mod_no_aa/tgp_heterozygosity_windows.tar.gz`.


__Compute the number of heterozygous sites per Archaic individual per site.__
```bash
for CHR in {1..22}; do for ARC in DEN ALT CHA VIN; do
python archaic_heterozygosity_chromosome_v_revisions.py ${CHR} ${ARC}
done; done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/den_masked_no_aa/archaic_heterozygosity_chromosome.tar.gz`, `muc19_results/alt_masked_no_aa/archaic_heterozygosity_chromosome.tar.gz` `muc19_results/cha_masked_no_aa/archaic_heterozygosity_chromosome.tar.gz`, and `muc19_results/vin_masked_no_aa/archaic_heterozygosity_chromosome.tar.gz`.


__Compute the number of heterozygous sites amongst archaic individuals in non-overlapping windows.__
```bash
for CHR in {1..22}; do for ARC in DEN ALT CHA VIN; do
python archaic_heterozygosity_windows_v_revisions.py ${CHR} 72 ${ARC}
done; done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/den_masked_no_aa/archaic_heterozygosity_windows.tar.gz`, `muc19_results/alt_masked_no_aa/archaic_heterozygosity_windows.tar.gz` `muc19_results/cha_masked_no_aa/archaic_heterozygosity_windows.tar.gz`, and `muc19_results/vin_masked_no_aa/archaic_heterozygosity_windows.tar.gz`.