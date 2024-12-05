# `introgression`

This directory contains the code to replicate the introgression results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__Compute _U{AFR,B,Denisovan}(1%, 30%, 100%)_ per non-AFR population (i.e., _B_) per NCBI RefSeq gene.__
```bash
for CHR in {1..22}; do
python tgp_u30_afr_b_den_genes_v_revisions.py ${CHR}
done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/tgp_u30_afr_b_den_genes.tar.gz`.


__Compute _U{AFR,B,Denisovan}(1%, 30%, 100%)_ per non-AFR population (i.e., _B_) in non-overlapping windows.__
```bash
for CHR in {1..22}; do for WIND in 72 742; do
python tgp_u30_afr_b_den_windows_v_revisions.py ${CHR} ${WIND}
done; done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/tgp_den_masked_no_aa/tgp_u30_afr_b_den_windows.tar.gz`.


__Jointly compute _Q95{AFR,B,Denisovan}(1%, 30%, 100%)_ _U{AFR,B,Denisovan}(1%, 30%, 100%)_ per non-AFR population (i.e., _B_) in non-overlapping windows.__
```bash
for CHR in {1..22}; do for WIND in 72 742; do
python tgp_q95_u30_afr_b_den_windows_v_revisions.py ${CHR} ${WIND} DEN
done; done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/tgp_den_masked_no_aa/tgp_q95_u30_afr_b_den_windows.tar.gz`.


__Compute _D+((YRI, NA19664), Archaic))_ in non-overlapping windows.__
```bash
for CHR in {1..22}; do for WIND in 72 742; do for ARC in DEN ALT CHA VIN; do
python mxl_archaic_site_patterns_windows_v_revisions.py ${CHR} ${WIND} ${ARC}
done; done; done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/tgp_den_masked_aa/mxl_archaic_site_patterns_windows.tar.gz`, `muc19_results/tgp_alt_masked_aa/mxl_archaic_site_patterns_windows.tar.gz` `muc19_results/tgp_cha_masked_aa/mxl_archaic_site_patterns_windows.tar.gz`, and `muc19_results/tgp_vin_masked_aa/mxl_archaic_site_patterns_windows.tar.gz`.


__Compute _D+((Altai Neanderthal, Late Neanderthals), Denisovan))_ in non-overlapping windows.__
```bash
for CHR in {1..22}; do
python archaic_site_patterns_windows_v_revisions.py ${CHR} 72
done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/arcs_masked_aa`.