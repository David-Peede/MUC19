# `vcf_data`

This directory contains the code to replicate the VCF data used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`BCFtoools v1.13`](https://samtools.github.io/bcftools/bcftools.html)
- [`Tabix v0.2.6`](http://www.htslib.org/doc/tabix.html)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)
- [`SnpEff v5.1`](https://pcingola.github.io/SnpEff/)
- [`Beagle v5.4`](https://faculty.washington.edu/browning/beagle/b5_4.html)

## Data

All data is publicly available and can be downloaded from the following locations:

- [EPO Ancestral Allele Calls](http://ftp.ensembl.org/pub/release-74/fasta/ancestral_alleles/)
- [Altai Denisovan Data](http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Denisova/)
- [Altai Neanderthal Data](http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Altai/)
- [Chagyrskaya Neandertal Data](http://ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF/)
- [Vindija Neanderthal Data](http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/)
- [Simons Genome Diversity Project Phased Dataset](https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data2021/)
- [Phase 3 Release of the 1000 Genomes Project](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)

## Code

Note that because the intermediate VCF and bookkeeping files take up roughly 3TB of storage.

### Single Archaic Datasets

__Mask out low-quality regions.__
```bash
ALT=./Altai
CHA=./Chagyrskaya
VIN=./Vindija33.19
DEN=./Denisova

for CHR in {1..22}; do
bcftools view -R ./masks/alt/chr${CHR}_mask.bed.gz -Oz -o ./arcs/alt_masked_chr${CHR}.vcf.gz ${ALT}/chr${CHR}_mq25_mapab100.vcf.gz
bcftools view -R ./masks/cha/chr${CHR}_mask.bed.gz -Oz -o ./arcs/cha_masked_chr${CHR}.vcf.gz ${CHA}/chr${CHR}.noRB.vcf.gz
bcftools view -R ./masks/vin/chr${CHR}_mask.bed.gz -Oz -o ./arcs/vin_masked_chr${CHR}.vcf.gz ${VIN}/chr${CHR}_mq25_mapab100.vcf.gz
bcftools view -R ./masks/den/chr${CHR}_mask.bed.gz -Oz -o ./arcs/den_masked_chr${CHR}.vcf.gz ${DEN}/chr${CHR}_mq25_mapab100.vcf.gz
done
```


__Index the masked VCF files.__
```bash
for CHR in {1..22}; do for ARC in alt cha vin den; do
tabix -p vcf ./arcs/${ARC}_masked_chr${CHR}.vcf.gz
done; done
```


__Filter each archaic excluding the ancestral allele.__
```bash
for CHR in {1..22}; do for ARC in alt cha vin den; do
python filter_archaic_no_aa_call.py ./arcs/${ARC}_masked_chr${CHR}.vcf.gz ${CHR} ${ARC}_masked | bgzip > ./filtered_merge/${ARC}_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz
done; done
```


__Convert the filtered VCF files to `Zarr` arrays.__
```bash
for CHR in {1..22}; do for ARC in alt cha vin den; do
python vcf_to_zarr.py ./filtered_merge/${ARC}_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz ${ARC}_masked_no_aa_chr${CHR} ${CHR}
done; done
```


### Combined Archaic Datasets

__Merge the unfiltered data.__
```bash
for CHR in {1..22}; do
bcftools merge -m all -Oz -o ./init_merge/all_archaics_masked_init_merge_chr${CHR}.vcf.gz ./arcs/alt_masked_chr${CHR}.vcf.gz ./arcs/cha_masked_chr${CHR}.vcf.gz ./arcs/vin_masked_chr${CHR}.vcf.gz ./arcs/den_masked_chr${CHR}.vcf.gz
done
```


__Identify duplicate records.__
```bash
for CHR in {1..22}; do
python identify_duplicate_records.py ./init_merge/all_archaics_masked_init_merge_chr${CHR}.vcf.gz ${CHR} all_archaics_masked | bgzip > ./dups/all_archaics_masked_dup_records_chr${CHR}.vcf.gz
done
```


__Filter the dataset excluding the ancestral allele.__
```bash
for CHR in {1..22}; do
python all_archaics_init_merge_filter_no_aa_call.py ./init_merge/all_archaics_masked_init_merge_chr${CHR}.vcf.gz ${CHR} all_archaics_masked | bgzip > ./filtered_merge/all_archaics_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz
done
```


__Filter the dataset including the ancestral allele.__
```bash
for CHR in {1..22}; do
python all_archaics_init_merge_filter_aa_call.py ./init_merge/all_archaics_masked_init_merge_chr${CHR}.vcf.gz ${CHR} all_archaics_masked | bgzip > ./filtered_merge/all_archaics_masked_var_sites_aa_calls_chr${CHR}.vcf.gz
done
```


__Convert the filtered VCF files to `Zarr` arrays.__
```bash
for CHR in {1..22}; do
python vcf_to_zarr.py ./filtered_merge/all_archaics_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz arcs_masked_no_aa_chr${CHR} ${CHR}
python vcf_to_zarr.py ./filtered_merge/all_archaics_masked_var_sites_aa_calls_chr${CHR}.vcf.gz arcs_masked_aa_chr${CHR} ${CHR}
done
```

### Modern Human Datasets

__Identify duplicate records.__
```bash
TGP=./tgp
SGDP=./sgdp

for CHR in {1..22}; do
python identify_duplicate_records.py ${TGP}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${CHR} tgp | bgzip > ./dups/tgp_dup_records_chr${CHR}.vcf.gz
python identify_duplicate_records.py ${SGDP}/sgdp.phased.unfiltered.chr${CHR}.vcf.gz ${CHR} sgdp | bgzip > ./dups/sgdp_dup_records_chr${CHR}.vcf.gz
done
```


__Filter each dataset excluding the ancestral allele.__
```bash
for CHR in {1..22}; do
python tgp_sgdp_init_merge_filter_no_aa_call.py ${TGP}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${CHR} tgp | bgzip > ./filtered_merge/tgp_var_sites_no_aa_calls_chr${CHR}.vcf.gz
python tgp_sgdp_init_merge_filter_no_aa_call.py ${SGDP}/sgdp.phased.unfiltered.chr${CHR}.vcf.gz ${CHR} sgdp | bgzip > ./filtered_merge/sgdp_var_sites_no_aa_calls_chr${CHR}.vcf.gz
done
```


__Filter the TGP dataset including the ancestral allele.__
```bash
for CHR in {1..22}; do
python tgp_sgdp_init_merge_filter_aa_call.py ${TGP}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${CHR} tgp | bgzip > ./filtered_merge/tgp_var_sites_aa_calls_chr${CHR}.vcf.gz
done
```


__Index the filtered TGP data.__
```bash
for CHR in {1..22}; do for PREFIX in tgp_var_sites_aa_calls tgp_var_sites_no_aa_calls; do
tabix -p vcf ./filtered_merge/${PREFIX}_chr${CHR}.vcf.gz
done; done
```


__Subset the filtered TGP data to not include ASW or ACB.__
```bash
for CHR in {1..22}; do
bcftools view -S ../meta_data/tgp_mod_samps.txt -Oz -o ./filtered_merge/tgp_mod_var_sites_no_aa_calls_chr${CHR}.vcf.gz ./filtered_merge/tgp_var_sites_no_aa_calls_chr${CHR}.vcf.gz
bcftools view -S ../meta_data/tgp_mod_anc_samps.txt -Oz -o ./filtered_merge/tgp_mod_var_sites_aa_calls_chr${CHR}.vcf.gz ./filtered_merge/tgp_var_sites_aa_calls_chr${CHR}.vcf.gz
done
```


__Convert the filtered VCF files to `Zarr` arrays.__
```bash
for CHR in {1..22}; do
python vcf_to_zarr.py ./filtered_merge/tgp_mod_var_sites_no_aa_calls_chr${CHR}.vcf.gz tgp_mod_no_aa_chr${CHR} ${CHR}
python vcf_to_zarr.py ./filtered_merge/tgp_mod_var_sites_aa_calls_chr${CHR}.vcf.gz tgp_mod_aa_chr${CHR} ${CHR}
python vcf_to_zarr.py ./filtered_merge/sgdp_var_sites_no_aa_calls_chr${CHR}.vcf.gz sgdp_no_aa_chr${CHR} ${CHR}
done
```


### Single Archaic + Modern Human Datasets

**‌Merge the unfiltered data.**
```bash
TGP=./tgp
SGDP=./sgdp

for CHR in {1..22}; do for ARC in alt cha vin den; do
bcftools merge -m all -Oz -o ./init_merge/tgp_${ARC}_masked_init_merge_chr${CHR}.vcf.gz ${TGP}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ./arcs/${ARC}_masked_chr${CHR}.vcf.gz
bcftools merge -m all -Oz -o ./init_merge/sgdp_${ARC}_masked_init_merge_chr${CHR}.vcf.gz ${SGDP}/sgdp.phased.unfiltered.chr${CHR}.vcf.gz ./arcs/${ARC}_masked_chr${CHR}.vcf.gz
done; done
```


**‌Identify duplicate records.**
```bash
for CHR in {1..22}; do for ARC in alt cha vin den; do
python identify_duplicate_records.py ./init_merge/tgp_${ARC}_masked_init_merge_chr${CHR}.vcf.gz ${CHR} tgp_${ARC}_masked | bgzip > ./dups/tgp_${ARC}_masked_dup_records_chr${CHR}.vcf.gz
python identify_duplicate_records.py ./init_merge/sgdp_${ARC}_masked_init_merge_chr${CHR}.vcf.gz ${CHR} sgdp_${ARC}_masked | bgzip > ./dups/sgdp_${ARC}_masked_dup_records_chr${CHR}.vcf.gz
done; done
```


**‌Filter each dataset excluding the ancestral allele.**
```bash
for CHR in {1..22}; do for ARC in alt cha vin den; do
python tgp_sgdp_single_archaic_init_merge_filter_no_aa_call.py ./init_merge/tgp_${ARC}_masked_init_merge_chr${CHR}.vcf.gz ${CHR} tgp_${ARC}_masked | bgzip > ./filtered_merge/tgp_${ARC}_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz
python tgp_sgdp_single_archaic_init_merge_filter_no_aa_call.py ./init_merge/sgdp_${ARC}_masked_init_merge_chr${CHR}.vcf.gz ${CHR} sgdp_${ARC}_masked | bgzip > ./filtered_merge/sgdp_${ARC}_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz
done; done
```


__Filter each combined TGP dataset including the ancestral allele.__
```bash
for CHR in {1..22}; do for ARC in alt cha vin den; do
python tgp_sgdp_single_archaic_init_merge_filter_aa_call.py ./init_merge/tgp_${ARC}_masked_init_merge_chr${CHR}.vcf.gz ${CHR} tgp_${ARC}_masked | bgzip > ./filtered_merge/tgp_${ARC}_masked_var_sites_aa_calls_chr${CHR}.vcf.gz
done; done
```


__Index the filtered TGP data.__
```bash
for CHR in {1..22}; do for PREFIX in tgp_alt_masked_var_sites_aa_calls tgp_alt_masked_var_sites_no_aa_calls tgp_cha_masked_var_sites_aa_calls tgp_cha_masked_var_sites_no_aa_calls tgp_vin_masked_var_sites_aa_calls tgp_vin_masked_var_sites_no_aa_calls tgp_den_masked_var_sites_aa_calls tgp_den_masked_var_sites_no_aa_calls; do
tabix -p vcf ./filtered_merge/${PREFIX}_chr${CHR}.vcf.gz
done; done
```


__Subset the filtered TGP data to not include ASW or ACB.__
```bash
for CHR in {1..22}; do for ARC in alt cha vin den; do
bcftools view -S ../meta_data/tgp_mod_${ARC}_samps.txt -Oz -o ./filtered_merge/tgp_mod_${ARC}_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz ./filtered_merge/tgp_${ARC}_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz
bcftools view -S ../meta_data/tgp_mod_${ARC}_anc_samps.txt -Oz -o ./filtered_merge/tgp_mod_${ARC}_masked_var_sites_aa_calls_chr${CHR}.vcf.gz ./filtered_merge/tgp_${ARC}_masked_var_sites_aa_calls_chr${CHR}.vcf.gz
done; done
```


__Convert the filtered VCF files to `Zarr` arrays.__
```bash
for CHR in {1..22}; do for ARC in alt cha vin den; do
python vcf_to_zarr.py ./filtered_merge/tgp_mod_${ARC}_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz tgp_${ARC}_masked_no_aa_chr${CHR} ${CHR}
python vcf_to_zarr.py ./filtered_merge/tgp_mod_${ARC}_masked_var_sites_aa_calls_chr${CHR}.vcf.gz tgp_${ARC}_masked_aa_chr${CHR} ${CHR}
python vcf_to_zarr.py ./filtered_merge/sgdp_${ARC}_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz sgdp_${ARC}_masked_no_aa_chr${CHR} ${CHR}
done; done
```


### All Archaics + Modern Human Combined Datasets

**‌Merge the unfiltered data.**
```bash
TGP=./tgp

for CHR in {1..22}; do 
bcftools merge -m all -Oz -o ./init_merge/tgp_all_archaics_masked_init_merge_chr${CHR}.vcf.gz $ {TGP}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ./arcs/alt_masked_chr${CHR}.vcf.gz ./arcs/cha_masked_chr${CHR}.vcf.gz ./arcs/vin_masked_chr${CHR}.vcf.gz ./arcs/den_masked_chr${CHR}.vcf.gz
done
```


**‌Identify duplicate records.**
```bash
for CHR in {1..22}; do
python identify_duplicate_records.py ./init_merge/tgp_all_archaics_masked_init_merge_chr${CHR}.vcf.gz ${CHR} tgp_all_archaics_masked | bgzip > ./dups/tgp_all_archaics_masked_dup_records_chr${CHR}.vcf.gz
done
```


**‌Filter each dataset excluding the ancestral allele.**
```bash
for CHR in {1..22}; do
python tgp_sgdp_all_archaics_init_merge_filter_no_aa_call.py ./init_merge/tgp_all_archaics_masked_init_merge_chr${CHR}.vcf.gz ${CHR} tgp_all_archaics_masked | bgzip > ./filtered_merge/tgp_all_archaics_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz
done
```


__Index the filtered TGP data.__
```bash
for CHR in {1..22}; do for PREFIX in tgp_all_archaics_masked_var_sites_no_aa_calls; do
tabix -p vcf ./filtered_merge/${PREFIX}_chr${CHR}.vcf.gz
done; done
```


__Subset the filtered TGP data to not include ASW or ACB.__
```bash
for CHR in {1..22}; do
bcftools view -S ../meta_data/tgp_mod_arc_samps.txt -Oz -o ./filtered_merge/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz ./filtered_merge/tgp_all_archaics_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz
done
```


__Subset the focal 72kb region for phasing.__
```bash
bcftools view -r 12:40759001-40831000 -Oz -o ./muc19/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb.vcf.gz ./filtered_merge/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_chr12.vcf.gz
```


__Index the focal 72kb region for phasing.__
```bash
tabix -p vcf ./muc19/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb.vcf.gz
```


__Convert the filtered VCF files to `Zarr` arrays.__
```bash
for CHR in {1..22}; do
python vcf_to_zarr.py ./filtered_merge/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz tgp_arcs_masked_no_aa_chr${CHR} ${CHR}
done
```


**Annotate the VCF files.**
```bash
for CHR in {1..22}; do
java -Xmx8g -jar ./snpEff/snpEff.jar ann -csvStats ./ann_summary/tgp_arcs_masked_no_aa_chr${CHR} -onlyTr ../annotations/hg19_genes/ncbi_refseq_transcripts.txt GRCh37.p13 ./filtered_merge/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_chr${CHR}.vcf.gz | bgzip > ./filtered_merge_ann/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_all_annotations_chr${CHR}.vcf.gz
done
```


**Extract the coding information.**
```bash
for CHR in {1..22}; do for SUFFIX in syn non; do
python subset_filtered_annotated_vcf_by_annotation.py ${CHR} tgp_mod_all_archaics_masked_var_sites_no_aa_calls ${SUFFIX} | bgzip > ./filtered_merge_ann/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_${SUFFIX}_muts_chr${CHR}.vcf.gz
done; done
```


### Phasing

**Subset the focal samples and reference panels.**
```bash
bcftools view -S ^../meta_data/tgp_mod_samps.txt -Oz -o ./phasing/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_focal_samps.vcf.gz ./muc19/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb.vcf.gz
bcftools view -S ../meta_data/tgp_mod_samps.txt -Oz -o ./phasing/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_ref_panel_all_inds.vcf.gz ./muc19/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb.vcf.gz
```


**Index the focal samples and reference panels.**
```bash
for PREFIX in tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_focal_samps tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_ref_panel_all_inds; do
tabix -p vcf ./phasing/${PREFIX}.vcf.gz
done
```


**‌Subset the reference panels for the different phasing datasets.**
```bash
for DATA in syn cha vin; do
bcftools view -R ./phasing/${DATA}_to_phase_sites.txt -Oz -o ./phasing/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_ref_panel_all_inds_${DATA}_to_phase_sites.vcf.gz ./phasing/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_ref_panel_all_inds.vcf.gz
done
```


**‌Subset the focal samples for the different phasing datasets.**
```bash
bcftools view -R ./phasing/syn_to_phase_sites.txt -Oz -o ./phasing/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_focal_samps_syn_to_phase_sites.vcf.gz ./phasing/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_focal_samps.vcf.gz
bcftools view -R ./phasing/cha_to_phase_sites.txt -s AltaiNeandertal,Chagyrskaya-Phalanx,Denisova -Oz -o ./phasing/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_focal_samps_cha_to_phase_sites.vcf.gz ./phasing/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_focal_samps.vcf.gz
bcftools view -R ./phasing/vin_to_phase_sites.txt -s AltaiNeandertal,Vindija33.19,Denisova -Oz -o ./phasing/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_focal_samps_vin_to_phase_sites.vcf.gz ./phasing/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_focal_samps.vcf.gz
```


**Remove the header from the focal samples VCF file for generating the synthetic Neanderthal.**
```bash
zcat tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_focal_samps_syn_to_phase_sites.vcf.gz | grep -v '^##' | gzip > tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_focal_samps_syn_to_phase_sites.txt.gz
```


**Extract the header from the focal samples VCF file for generating the synthetic Neanderthal.**
```bash
zcat tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_focal_samps_syn_to_phase_sites.vcf.gz | grep '^##' > tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_focal_samps_syn_to_phase_sites_header.txt
```


**Phase the synthetic Neanderthal.**
```bash
java -Xmx2g -jar beagle.22Jul22.46e.jar gt=./phasing/alt_syn_nean_den_unphased.vcf ref=./phasing/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_ref_panel_all_inds_syn_to_phase_sites.vcf.gz out=./phasing/alt_syn_nean_den_phased_ref_panel_all_inds map=./phasing/plink_maps/plink.chr12.GRCh37.map chrom=12:40759001-40831000 impute=False seed=42
```


**Phase the late Neanderthals.**
```bash
for DATA in cha vin; do
java -Xmx2g -jar beagle.22Jul22.46e.jar gt=./phasing/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_focal_samps_${DATA}_to_phase_sites.vcf.gz  ref=./phasing/tgp_mod_all_archaics_masked_var_sites_no_aa_calls_muc19_72kb_ref_panel_all_inds_${DATA}_to_phase_sites.vcf.gz out=./phasing/${DATA}_phased_ref_panel_all_inds map=./phasing/plink_maps/plink.chr12.GRCh37.map chrom=12:40759001-40831000 impute=False seed=42
done
```


__Convert the phased late Neanderthal VCF files to `Zarr` arrays.__
```bash
for DATA in cha vin; do
python vcf_to_zarr.py ./phasing/${DATA}_phased_ref_panel_all_inds.vcf.gz ${DATA}_phased_ref_panel_all_inds 12
done
```
