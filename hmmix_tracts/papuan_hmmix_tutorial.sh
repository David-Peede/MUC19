#!/bin/sh

# ----------------------------------------------------
# install hmmix and dependencies
# ----------------------------------------------------
# pip install hmmix 
# conda install -c bioconda vcftools bcftools

# ----------------------------------------------------
# initialize folders
# ----------------------------------------------------
mkdir external_files
mkdir external_files/anc  
mkdir external_files/arc  
mkdir external_files/bcf_outgroup
mkdir external_files/bcf_ingroup
mkdir external_files/ref

mkdir obs_files
mkdir trained_files
mkdir decoded_files

# ----------------------------------------------------
# Download external data 
# ----------------------------------------------------
# Ingroup and outgroup sequences (bcf_ingroup/ and bcf_outgroup/), 
# callability mask, and ancestral information (anc/).
# Variants from high coverage archaic sequences (arc/) and Reference genome 
#  information (ref/) are included in this run because ingroup and outgroup 
#  samples are from different datasets.  Please see 
#  'Finding snps which are derived in the outgroup' section of hmmix README
#  https://github.com/LauritsSkov/Introgression-detection/blob/master/README.md#finding-snps-which-are-derived-in-the-outgroup 
# 
# Download links are provided in section 'Getting data' of hmmix README
# https://github.com/LauritsSkov/Introgression-detection/blob/master/README.md 

# Ancestral information (anc/)
wget -P external_files/anc ftp.ensembl.org/pub/release-74/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2
tar -xf external_files/anc/homo_sapiens_ancestor_GRCh37_e71.tar.bz2 --strip-components=1 --include="homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_[123456789]*.fa" -C external_files/anc
# rm external_files/anc/homo_sapiens_ancestor_GRCh37_e71.tar.bz2

# Archaic variants (Altai, Vindija, Chagyrskaya and Denisova in hg19; arc/)
wget -P external_files/arc https://zenodo.org/api/records/7246376/files-archive
tar -xf external_files/arc/files-archive -C external_files/arc
# rm external_files/arc/files-archive

# Variants for the ingroup and outgroups (bcf_ingroup/ and bcf_outgroup/)
for i in {1..22}
do
	# 1kG autosome for outgroup
	wget -P external_files/bcf_outgroup ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/bcf_files/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf
	# SGDP autosome for ingroup
	wget -P external_files/bcf_ingroup https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data2021/chr.sgdp.pub.$i.bcf
done

# 1kG callability mask (remember to remove chr in the beginning of each line to make it compatible with hg19 e.g. chr1 > 1)
wget ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20141020.strict_mask.whole_genome.bed
sed 's/^chr\|%$//g' 20141020.strict_mask.whole_genome.bed | awk '{print $1"\t"$2"\t"$3}' > external_files/strickmask.bed
rm 20141020.strict_mask.whole_genome.bed

# Reference genome (ref/)
wget -P external_files/ref ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar -xf external_files/ref/chromFa.tar.gz --include="chr[123456789].fa" -C external_files/ref
tar -xf external_files/ref/chromFa.tar.gz --include="chr[123456789][123456789].fa" -C external_files/ref
# rm external_files/ref/chromFa.tar.gz


# ----------------------------------------------------
# Call archaic segments in SGDP Papuan samples, using 1kG outgroup 
# ----------------------------------------------------

# .json file listing ingroup and outgroup samples

# text file listing sample ids for ingroup samples

# ----------------------------------------------------
# Run hmmix set up 
# ----------------------------------------------------

# hmmix Step 1: create outgroup
hmmix create_outgroup -ind=papuan_individuals.json -vcf='external_files/bcf/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,external_files/bcf/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf' -ancestral='external_files/anc/homo_sapiens_ancestor_10.fa,external_files/anc/homo_sapiens_ancestor_11.fa,external_files/anc/homo_sapiens_ancestor_12.fa,external_files/anc/homo_sapiens_ancestor_13.fa,external_files/anc/homo_sapiens_ancestor_14.fa,external_files/anc/homo_sapiens_ancestor_15.fa,external_files/anc/homo_sapiens_ancestor_16.fa,external_files/anc/homo_sapiens_ancestor_17.fa,external_files/anc/homo_sapiens_ancestor_18.fa,external_files/anc/homo_sapiens_ancestor_19.fa,external_files/anc/homo_sapiens_ancestor_1.fa,external_files/anc/homo_sapiens_ancestor_20.fa,external_files/anc/homo_sapiens_ancestor_21.fa,external_files/anc/homo_sapiens_ancestor_22.fa,external_files/anc/homo_sapiens_ancestor_2.fa,external_files/anc/homo_sapiens_ancestor_3.fa,external_files/anc/homo_sapiens_ancestor_4.fa,external_files/anc/homo_sapiens_ancestor_5.fa,external_files/anc/homo_sapiens_ancestor_6.fa,external_files/anc/homo_sapiens_ancestor_7.fa,external_files/anc/homo_sapiens_ancestor_8.fa,external_files/anc/homo_sapiens_ancestor_9.fa' -refgenome='external_files/ref/chr10.fa,external_files/ref/chr11.fa,external_files/ref/chr12.fa,external_files/ref/chr13.fa,external_files/ref/chr14.fa,external_files/ref/chr15.fa,external_files/ref/chr16.fa,external_files/ref/chr17.fa,external_files/ref/chr18.fa,external_files/ref/chr19.fa,external_files/ref/chr1.fa,external_files/ref/chr20.fa,external_files/ref/chr21.fa,external_files/ref/chr22.fa,external_files/ref/chr2.fa,external_files/ref/chr3.fa,external_files/ref/chr4.fa,external_files/ref/chr5.fa,external_files/ref/chr6.fa,external_files/ref/chr7.fa,external_files/ref/chr8.fa,external_files/ref/chr9.fa' -weights=external_files/strickmask.bed -out=outgroup.txt

# hmmix Step 2: calculate mutation rate
hmmix mutation_rate -outgroup=outgroup.txt -weights=external_files/strickmask.bed -window_size=100000 -out=mutationrate.bed

# ----------------------------------------------------
# Run hmmix on each sample (loop example)
# ----------------------------------------------------
# We perform hmmix Steps 3, 4, and 5 on each ingroup sample.
# Note that Step 5 includes explicitly includes variants from the high-coverage 
# archaic samples, as recommended in 
# 'Finding snps which are derived in the outgroup' section of hmmix README when 
# ingroup and outgroup samples are not fromthe same dataset.
# Loop through haplotype samples in papuan_ingroup_haps.txt.

# -----start loop-------
while read ingroup_hap; do
	echo "Ingroup_hap: "$ingroup_hap

	# hmmix Step 3: create ingroup
	echo "create_ingroup is starting"
	ingroup_file="obs.${ingroup_hap}.txt"
	hmmix create_ingroup -ind=$ingroup_hap -vcf='external_files/bcf_ingroup/chr.sgdp.pub.10.bcf,external_files/bcf_ingroup/chr.sgdp.pub.11.bcf,external_files/bcf_ingroup/chr.sgdp.pub.12.bcf,external_files/bcf_ingroup/chr.sgdp.pub.13.bcf,external_files/bcf_ingroup/chr.sgdp.pub.14.bcf,external_files/bcf_ingroup/chr.sgdp.pub.15.bcf,external_files/bcf_ingroup/chr.sgdp.pub.16.bcf,external_files/bcf_ingroup/chr.sgdp.pub.17.bcf,external_files/bcf_ingroup/chr.sgdp.pub.18.bcf,external_files/bcf_ingroup/chr.sgdp.pub.19.bcf,external_files/bcf_ingroup/chr.sgdp.pub.1.bcf,external_files/bcf_ingroup/chr.sgdp.pub.20.bcf,external_files/bcf_ingroup/chr.sgdp.pub.21.bcf,external_files/bcf_ingroup/chr.sgdp.pub.22.bcf,external_files/bcf_ingroup/chr.sgdp.pub.2.bcf,external_files/bcf_ingroup/chr.sgdp.pub.3.bcf,external_files/bcf_ingroup/chr.sgdp.pub.4.bcf,external_files/bcf_ingroup/chr.sgdp.pub.5.bcf,external_files/bcf_ingroup/chr.sgdp.pub.6.bcf,external_files/bcf_ingroup/chr.sgdp.pub.7.bcf,external_files/bcf_ingroup/chr.sgdp.pub.8.bcf,external_files/bcf_ingroup/chr.sgdp.pub.9.bcf' -ancestral=external_files/anc/homo_sapiens_ancestor_10.fa,external_files/anc/homo_sapiens_ancestor_11.fa,external_files/anc/homo_sapiens_ancestor_12.fa,external_files/anc/homo_sapiens_ancestor_13.fa,external_files/anc/homo_sapiens_ancestor_14.fa,external_files/anc/homo_sapiens_ancestor_15.fa,external_files/anc/homo_sapiens_ancestor_16.fa,external_files/anc/homo_sapiens_ancestor_17.fa,external_files/anc/homo_sapiens_ancestor_18.fa,external_files/anc/homo_sapiens_ancestor_19.fa,external_files/anc/homo_sapiens_ancestor_1.fa,external_files/anc/homo_sapiens_ancestor_20.fa,external_files/anc/homo_sapiens_ancestor_21.fa,external_files/anc/homo_sapiens_ancestor_22.fa,external_files/anc/homo_sapiens_ancestor_2.fa,external_files/anc/homo_sapiens_ancestor_3.fa,external_files/anc/homo_sapiens_ancestor_4.fa,external_files/anc/homo_sapiens_ancestor_5.fa,external_files/anc/homo_sapiens_ancestor_6.fa,external_files/anc/homo_sapiens_ancestor_7.fa,external_files/anc/homo_sapiens_ancestor_8.fa,external_files/anc/homo_sapiens_ancestor_9.fa -weights=external_files/strickmask.bed -outgroup=outgroup.txt
	
	# hmmix Step 4: train HMM
	echo "train is starting"
	trained_file="trained_files/trained.${ingroup_hap}.haploid.json"
	hmmix train -obs=$ingroup_file -weights=external_files/strickmask.bed -mutrates=mutationrate.bed -out=$trained_file -haploid

	# hmmix Step 5: decode archaic states
	echo "decode is starting"
	decoded_file="decoded_files/decoded.${ingroup_hap}.autosome"
	hmmix decode -obs=$ingroup_file -weights=external_files/strickmask.bed -mutrates=mutationrate.bed -param=$trained_file -out=$decoded_file -haploid -admixpop='external_files/arc/highcov_ind_10.bcf,external_files/arc/highcov_ind_11.bcf,external_files/arc/highcov_ind_12.bcf,external_files/arc/highcov_ind_13.bcf,external_files/arc/highcov_ind_14.bcf,external_files/arc/highcov_ind_15.bcf,external_files/arc/highcov_ind_16.bcf,external_files/arc/highcov_ind_17.bcf,external_files/arc/highcov_ind_18.bcf,external_files/arc/highcov_ind_19.bcf,external_files/arc/highcov_ind_1.bcf,external_files/arc/highcov_ind_20.bcf,external_files/arc/highcov_ind_21.bcf,external_files/arc/highcov_ind_22.bcf,external_files/arc/highcov_ind_2.bcf,external_files/arc/highcov_ind_3.bcf,external_files/arc/highcov_ind_4.bcf,external_files/arc/highcov_ind_5.bcf,external_files/arc/highcov_ind_6.bcf,external_files/arc/highcov_ind_7.bcf,external_files/arc/highcov_ind_8.bcf,external_files/arc/highcov_ind_9.bcf'

	# tidy up
	mv $ingroup_file obs_files/
	
done <papuan_ingroup_haps.txt
# -----end loop-------

# ----------------------------------------------------
# * Run hmmix on each sample (SLURM example)
# ----------------------------------------------------
# We perform hmmix Steps 3, 4, and 5 on each ingroup sample.
# Note that Step 5 includes explicitly includes variants from the high-coverage 
# archaic samples, as recommended in 
# 'Finding snps which are derived in the outgroup' section of hmmix README when 
# ingroup and outgroup samples are not fromthe same dataset.
# Submit SLURM array job for each haplotype sample in papuan_ingroup_haps.txt.

# !!!: to do this in practice:
# 		[ ] move below code to a separate file
#       [ ] uncomment one level
# 		[ ] set email address
# 		[ ] confirm upper limit of `--array` matches length of `papuan_ingroup_haps.txt`
# 		[ ] adjust the `source` and `module` commands to fit your installation 

# ---------vvvv SLURM batch submission file vvvvv------

# #SBATCH --mail-type=FAIL,END
# #SBATCH --mail-user=<your@email.here>

# #SBATCH -e arrayjob-%a.err
# #SBATCH -o arrayjob-%a.out

# #SBATCH -t 4:00:00
# #SBATCH --mem 4Gb

# #SBATCH --array=1-14

# source ~/hmmix_env/bin/activate
# module load python bcftools vcftools

# ingroup_hap=`cat papuan_ingroup_haps.txt | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

# echo "Ingroup_hap: "$ingroup_hap

# # hmmix Step 3: create ingroup
# echo "create_ingroup is starting"
# ingroup_file="obs.${ingroup_hap}.txt"
# hmmix create_ingroup -ind=$ingroup_hap -vcf='external_files/bcf_ingroup/chr.sgdp.pub.10.bcf,external_files/bcf_ingroup/chr.sgdp.pub.11.bcf,external_files/bcf_ingroup/chr.sgdp.pub.12.bcf,external_files/bcf_ingroup/chr.sgdp.pub.13.bcf,external_files/bcf_ingroup/chr.sgdp.pub.14.bcf,external_files/bcf_ingroup/chr.sgdp.pub.15.bcf,external_files/bcf_ingroup/chr.sgdp.pub.16.bcf,external_files/bcf_ingroup/chr.sgdp.pub.17.bcf,external_files/bcf_ingroup/chr.sgdp.pub.18.bcf,external_files/bcf_ingroup/chr.sgdp.pub.19.bcf,external_files/bcf_ingroup/chr.sgdp.pub.1.bcf,external_files/bcf_ingroup/chr.sgdp.pub.20.bcf,external_files/bcf_ingroup/chr.sgdp.pub.21.bcf,external_files/bcf_ingroup/chr.sgdp.pub.22.bcf,external_files/bcf_ingroup/chr.sgdp.pub.2.bcf,external_files/bcf_ingroup/chr.sgdp.pub.3.bcf,external_files/bcf_ingroup/chr.sgdp.pub.4.bcf,external_files/bcf_ingroup/chr.sgdp.pub.5.bcf,external_files/bcf_ingroup/chr.sgdp.pub.6.bcf,external_files/bcf_ingroup/chr.sgdp.pub.7.bcf,external_files/bcf_ingroup/chr.sgdp.pub.8.bcf,external_files/bcf_ingroup/chr.sgdp.pub.9.bcf' -ancestral=external_files/anc/homo_sapiens_ancestor_10.fa,external_files/anc/homo_sapiens_ancestor_11.fa,external_files/anc/homo_sapiens_ancestor_12.fa,external_files/anc/homo_sapiens_ancestor_13.fa,external_files/anc/homo_sapiens_ancestor_14.fa,external_files/anc/homo_sapiens_ancestor_15.fa,external_files/anc/homo_sapiens_ancestor_16.fa,external_files/anc/homo_sapiens_ancestor_17.fa,external_files/anc/homo_sapiens_ancestor_18.fa,external_files/anc/homo_sapiens_ancestor_19.fa,external_files/anc/homo_sapiens_ancestor_1.fa,external_files/anc/homo_sapiens_ancestor_20.fa,external_files/anc/homo_sapiens_ancestor_21.fa,external_files/anc/homo_sapiens_ancestor_22.fa,external_files/anc/homo_sapiens_ancestor_2.fa,external_files/anc/homo_sapiens_ancestor_3.fa,external_files/anc/homo_sapiens_ancestor_4.fa,external_files/anc/homo_sapiens_ancestor_5.fa,external_files/anc/homo_sapiens_ancestor_6.fa,external_files/anc/homo_sapiens_ancestor_7.fa,external_files/anc/homo_sapiens_ancestor_8.fa,external_files/anc/homo_sapiens_ancestor_9.fa -weights=external_files/strickmask.bed -outgroup=outgroup.txt

# # hmmix Step 4: train HMM
# echo "train is starting"
# trained_file="trained_files/trained.${ingroup_hap}.haploid.json"
# hmmix train -obs=$ingroup_file -weights=external_files/strickmask.bed -mutrates=mutationrate.bed -out=$trained_file -haploid

# # hmmix Step 5: decode archaic states
# echo "decode is starting"
# decoded_file="decoded_files/decoded.${ingroup_hap}.autosome"
# hmmix decode -obs=$ingroup_file -weights=external_files/strickmask.bed -mutrates=mutationrate.bed -param=$trained_file -out=$decoded_file -haploid

# # tidy up
# mv $ingroup_file obs_files/
