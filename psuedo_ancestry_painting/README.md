# `psuedo_ancestry_painting`

This directory contains the code to replicate the pseudo ancestry painting results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__Compute the number of PAP sites for TGP + archaic individuals in non-overlapping windows.__
```bash
for CHR in {1..22}; do
python tgp_archaic_psuedo_ancestry_painting_windows_v_revisions.py ${CHR} 72
done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/tgp_arcs_masked_no_aa/tgp_archaic_psuedo_ancestry_painting_windows.tar.gz`.


__Compute the number of PAP sites for archaic individuals in non-overlapping windows.__
```bash
for CHR in {1..22}; do
python archaic_psuedo_ancestry_painting_windows_v_revisions.py ${CHR} 72
done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/arcs_masked_no_aa/archaic_psuedo_ancestry_painting_windows.tar.gz`.