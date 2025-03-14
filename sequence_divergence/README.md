# `sequence_divergence`

This directory contains the code to replicate the sequence divergence results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__Compute the sequence divergence between haplotypes in the 1000 Genomes Project and archaic diplotypes in non-overlapping windows.__
```bash
for CHR in {1..22}; do for WIND in 72 742; do for ARC in DEN ALT CHA VIN; do
python tgp_hap_v_arc_dip_divergence_windows_v_revisions.py ${CHR} ${WIND} ${ARC}
done; done; done
```


__Compile the sequence divergence summary between all pairwise combinations 1000 Genomes Project haplotypes and archaic diplotypes.__
```bash
for WIND in 72 742; do
python tgp_hap_v_arc_dip_divergence_summary_v_revisions.py ${WIND}
done
```


__Compute the sequence divergence between haplotypes in the 1000 Genomes Project and late Neanderthal pseudo-haplotypes in non-overlapping windows.__
```bash
for CHR in {1..22}; do for ARC in CHA VIN; do
python tgp_hap_v_arc_psuedo_hap_divergence_windows_v_revisions.py ${CHR} 72 ${ARC}
done; done
```


__Compile the sequence divergence summary between all pairwise combinations 1000 Genomes Project haplotypes and phased late Neanderthal haplotypes.__
```bash
python tgp_hap_v_cha_vin_phased_hap_divergence_summary_v_revisions.py
```


__Compute the sequence divergence between Papuan haplotypes in the Simons Genome Diversity Project and the Denisovan per Denisovan introgressed tract.__
```bash
for CHR in {1..22}; do for PAP in S_Papuan-1.DG S_Papuan-2.DG S_Papuan-3.DG S_Papuan-4.DG S_Papuan-5.DG S_Papuan-6.DG S_Papuan-7.DG S_Papuan-8.DG S_Papuan-9.DG S_Papuan-10.DG S_Papuan-11.DG S_Papuan-12.DG S_Papuan-13.DG S_Papuan-14.DG B_Papuan-15.DG; do
python sgdp_denisovan_sequence_divergence_at_denisovan_intro_tracts_in_papuans_v_revisions.py ${CHR} ${PAP}
done; done
```


__Compute the sequence divergence between the Denisovan and Altai Neanderthal diplotypes in non-overlapping windows.__
```bash
for CHR in {1..22}; do
python arc_dip_v_arc_dip_divergence_windows_v_revisions.py ${CHR} 72
done
```


__Compute the average sequence divergence between all  African individuals in the 1000 Genomes Project and Denisovan diplotype in non-overlapping windows.__
```bash
for CHR in {1..22}; do
python tgp_spop_v_arc_dip_divergence_windows_v_revisions.py ${CHR} 72 AFR DEN
done
```


