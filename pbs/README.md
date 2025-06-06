# `pbs`

This directory contains the code to replicate the empirical PBS results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Data

The data is publicly available and can be downloaded from the following location:

- [SPrime Introgression Tracts](https://data.mendeley.com/datasets/y7hyt83vxr/1)

## Code

__Compute PBS for MXL:CHB:CEU partitioned by ancestry per non-overlapping window.__
```bash
for CHR in {1..22}; do for WIND in 72 742; do
python mxl_chb_ceu_pbs_windows_v_revisions.py ${CHR} ${WIND}
done; done
```


__Compute PBS for MXL:CHB:CEU partitioned by ancestry per SNP.__
```bash
for CHR in {1..22}; do
python mxl_chb_ceu_pbs_chromsome_v_revisions.py ${CHR}
done
```


__Compute PBS for MXL:CHB:CEU for SPrime sites per non-overlapping window.__
```bash
for CHR in {1..22}; do for WIND in 72 742; do
python sprime_sites_mxl_chb_ceu_pbs_windows_v_revisions.py ${CHR} ${WIND}
done; done
```


__Compute PBS for AMR:ASN:EUR per non-overlapping window.__
```bash
for CHR in {1..22}; do for WIND in 72 742; do
python amr_asn_eur_pbs_windows_v_revisions.py ${CHR} ${WIND}
done; done
```