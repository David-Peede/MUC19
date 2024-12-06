# `iHS`

This directory contains the code to replicate the _iHS_ results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`selscan v2.0.0`](https://github.com/szpiech/selscan)
- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__Generate the input population files for selscan.__
```bash
for CHR in {1..22}; do for POP in mxl pel clm pur chb ceu yri lwk gwd msl esn beb stu itu pjl gih khv chs jpt cdx tsi ibs gbr fin; do
python make_selscan_population_data_v_revisions.py ${CHR} ${POP} > ./selscan_vcfs/${POP}_selscan_chr${CHR}.vcf
done; done
```
Note that the recombination maps can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `iHS/plink_maps.tar.gz` and the corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `iHS/maps/tgp_selscan_maps.tar.gz` and `iHS/selscan_vcfs/tgp_selscan_vcfs.tar.gz`.


__Compute the unstanderized _iHS_ values.__
```bash
for CHR in {1..22}; do for POP in mxl pel clm pur chb ceu yri lwk gwd msl esn beb stu itu pjl gih khv chs jpt cdx tsi ibs gbr fin; do
selscan-2.0.0/bin/linux/./selscan --ihs --vcf ./selscan_vcfs/${POP}_selscan_chr${CHR}.vcf --map ./maps/${POP}_selscan_chr${CHR}.map --out ../muc19_results/tgp_mod_aa/${POP}_chr${CHR}
done; done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/tgp_mod_aa/tgp_unstanderized_ihs.tar.gz`.


__Compute the normalized _iHS_ values by frequency bins.__
```bash
for CHR in {1..22}; do for POP in mxl pel clm pur chb ceu yri lwk gwd msl esn beb stu itu pjl gih khv chs jpt cdx tsi ibs gbr fin; do
selscan-2.0.0/bin/linux/./norm --ihs --files ../muc19_results/tgp_mod_aa/${POP}_chr1.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr2.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr3.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr4.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr5.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr6.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr7.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr8.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr9.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr10.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr11.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr12.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr13.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr14.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr15.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr16.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr17.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr18.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr19.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr20.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr21.ihs.out ../muc19_results/tgp_mod_aa/${POP}_chr22.ihs.out
done; done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/tgp_mod_aa/tgp_normalized_ihs.tar.gz`.


__Compile the _iHS_ results.__
```bash
for POP in MXL PEL CLM PUR CHB CEU YRI LWK GWD MSL ESN BEB STU ITU PJL GIH KHV CHS JPT CDX TSI IBS GBR FIN; do
python compile_ihs_tables_v_revisions.py ${POP}
done
```
The corresponding outputs can be found on [Google Drive](https://drive.google.com/drive/folders/1w1uz1a0-l9LwR6x3CKWPgPtT02F1uKzv?usp=sharing) at `muc19_results/tgp_mod_aa/compile_ihs_tables.tar.gz`.