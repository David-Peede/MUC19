# `site_patterns`

This directory contains the code to replicate the site pattern results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__Compute _D+_((YRI, OOA Population), Archaic)).__

```bash
# Compute site patterns for P1=YRI, P2=OOA, P3=DEN/ALT/CHA/VIN in non-overlapping windows.
for CHR in {1..22}; do for P2 in BEB STU ITU PJL GIH CHB KHV CHS JPT CDX TSI CEU IBS GBR FIN PEL MXL CLM PUR; do for P3 in DEN ALT CHA VIN; do
python tgp_site_patterns_windows.py ${CHR} 72 YRI ${P2} ${P3}
done; done; done
```

__Compute _D+_((Altai Neanderthal, Late Neanderthals), Denisovan)).__

```bash
# Compute site patterns for P1=ALT, P2=CHA/VIN, P3=DEN in non-overlapping windows.
for CHR in {1..22}; do
python arc_site_patterns_windows.py ${CHR} 72 ALT CHA DEN
python arc_site_patterns_windows.py ${CHR} 72 ALT VIN DEN
done
```

__Compute _D+_((Chagyrskaya Neanderthal, Vindija Neanderthal), Altai Neanderthal/Denisovan)).__

```bash
# Compute site patterns for P1=CHA, P2=VIN, P3=ALT/DEN in non-overlapping windows.
for CHR in {1..22}; do
python arc_site_patterns_windows.py ${CHR} 72 CHA VIN ALT
python arc_site_patterns_windows.py ${CHR} 72 CHA VIN DEN
done
```

