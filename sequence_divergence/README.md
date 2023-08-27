# `sequence_divergence`

This directory contains the code to replicate the sequence divergence results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__Compute sequence divergence between the Denisovan and Neanderthal diplotypes.__

```bash
# Compute sequence divergence between archaic diplotypes in non-overlapping windows.
for CHR in {1..22}; do
python den_nea_sequence_divergence_windows.py ${CHR} 72
done
```

__Compute sequence divergence between individuals in the 1000 Genomes Project and the archaic diplotypes.__

```bash
# Compute sequence divergence from the archaics in non-overlapping windows.
for CHR in {1..22}; do for ARC in DEN ALT CHA VIN; do for TGP in LWK GWD MSL ESN YRI BEB STU ITU PJL GIH CHB KHV CHS JPT CDX TSI CEU IBS GBR FIN PEL MXL CLM PUR; do
python tgp_arc_sequence_divergence_windows.py ${CHR} 72 ${TGP} ${ARC}
done; done; done
```

__Compute sequence divergence between haplotype groups in the 1000 Genomes Project and the archaic diplotypes.__

```bash
# Compute the haplotype group sequence divergence from the archaics in non-overlapping windows.
for CHR in {1..22}; do
python tgp_haplotype_group_arc_diplotype_sequence_divergence_windows.py ${CHR} 72
done
```

__Compute sequence divergence between haplotypes in the 1000 Genomes Project and the archaic diplotypes.__

```bash
# Compute sequence divergence between modern haplotypes and archaic diplotypes in non-overlapping windows.
for CHR in {1..22}; do for WIND in 72 748; do
python tgp_haplotype_arc_diplotype_sequence_divergence_windows.py ${CHR} ${WIND}
done; done
```

__Compute the number of ((African, Altai Neanderthal), Denisovan) shared divergent sites.__

```bash
# Compute the shared sequence divergence between AFR-ALT and DEN in non-overlapping windows.
for CHR in {1..22}; do
python afr_alt_den_sequence_divergence_windows.py ${CHR} 72 2350
done
```

