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


__Compute the number of heterozygous sites amongst TGP individuals in non-overlapping windows.__
```bash
for CHR in {1..22}; do
python tgp_heterozygosity_windows_v_revisions.py ${CHR} 72
done
```


__Compute the number of heterozygous sites per Archaic individual per site.__
```bash
for CHR in {1..22}; do for ARC in DEN ALT CHA VIN; do
python archaic_heterozygosity_chromosome_v_revisions.py ${CHR} ${ARC}
done; done
```

__Compute the number of heterozygous sites amongst archaic individuals in non-overlapping windows.__
```bash
for CHR in {1..22}; do for ARC in DEN ALT CHA VIN; do
python archaic_heterozygosity_windows_v_revisions.py ${CHR} 72 ${ARC}
done; done
```