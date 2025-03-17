# `arc_snp_density`

This directory contains the code to replicate the archaic SNP density results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__Compute archaic SNP density in non-overlapping windows across the autosomes.__
```bash
for CHR in {1..22}; do for WIND in 72 742; do
python tgp_archaic_snp_denisty_windows_v_revisions.py ${CHR} ${WIND}
done; done
```


__Classify all snp types for all chromsomes in all non-AFR population and globally in the TGP.__
```bash
for CHR in {1..22}; do
python classify_tgp_snps_chromosome_v_revisions.py ${CHR}
done
```