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
mkdir obs_files
mkdir trained_files
mkdir decoded_files

# ----------------------------------------------------
# download external data (.bcf, .fa, and callability mask)
# ----------------------------------------------------

# download links from section Getting data of hmmix README
# https://github.com/LauritsSkov/Introgression-detection/blob/master/README.md 

# .bcf file for chromosome 12
wget -P external_files/ ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/bcf_files/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf

# callability mask (remember to remove chr in the beginning of each line to make it compatible with hg19 e.g. chr1 > 1)
wget ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20141020.strict_mask.whole_genome.bed
sed 's/^chr\|%$//g' 20141020.strict_mask.whole_genome.bed | awk '{print $1"\t"$2"\t"$3}' > external_files/strickmask.bed
# rm 20141020.strict_mask.whole_genome.bed

# Ancestral information
# !!!: may take a while to download and extract
wget -P external_files/ ftp.ensembl.org/pub/release-74/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2
tar -xvjf external_files/homo_sapiens_ancestor_GRCh37_e71.tar.bz2 --strip-components=1 -C external_files homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_12.fa
# rm external_files/homo_sapiens_ancestor_GRCh37_e71.tar.bz2

# ----------------------------------------------------
# choose samples for which to call archaic segments
# ----------------------------------------------------

# .json file listing ingroup and outgroup samples
wget https://raw.githubusercontent.com/ramachandran-lab/archaic-chrX-sexbias/main/hmmix/individuals.json

# text file listing sample ids for ingroup samples
wget https://raw.githubusercontent.com/ramachandran-lab/archaic-chrX-sexbias/main/hmmix/ingroup_haps.txt

# ----------------------------------------------------
# Run hmmix set up 
# ----------------------------------------------------

# hmmix Step 1: create outgroup
hmmix create_outgroup -ind=individuals.json -vcf=external_files/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf -ancestral=external_files/homo_sapiens_ancestor_12.fa -weights=external_files/strickmask.bed -out=outgroup.txt

# hmmix Step 2: calculate mutation rate
hmmix mutation_rate -outgroup=outgroup.txt -weights=external_files/strickmask.bed -window_size=100000 -out=mutationrate.bed

# ----------------------------------------------------
# Run hmmix on each sample (loop example)
# ----------------------------------------------------
# We perform hmmix Steps 3, 4, and 5 on each ingroup sample.
# Loop through haplotype samples in ingroup_haps.txt.

# -----start loop-------
while read ingroup_hap; do
	echo "Ingroup_hap: "$ingroup_hap

	# hmmix Step 3: create ingroup
	echo "create_ingroup is starting"
	ingroup_file="obs.${ingroup_hap}.txt"
	hmmix create_ingroup -ind=$ingroup_hap -out=$ingroup_file -vcf=external_files/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf -out=obs -ancestral=external_files/homo_sapiens_ancestor_12.fa -weights=external_files/strickmask.bed -outgroup=outgroup.txt

	# hmmix Step 4: train HMM
	echo "train is starting"
	trained_file="trained_files/trained.${ingroup_hap}.haploid.json"
	hmmix train -obs=$ingroup_file -weights=external_files/strickmask.bed -mutrates=mutationrate.bed -out=$trained_file -haploid

	# hmmix Step 5: decode archaic states
	echo "decode is starting"
	decoded_file="decoded_files/decoded.${ingroup_hap}.autosome"
	hmmix decode -obs=$ingroup_file -weights=external_files/strickmask.bed -mutrates=mutationrate.bed -param=$trained_file -out=$decoded_file -haploid

	# tidy up
	mv $ingroup_file obs_files/
	
done <ingroup_haps.txt
# -----end loop-------

# ----------------------------------------------------
# * Run hmmix on each sample (SLURM example)
# ----------------------------------------------------
# We perform hmmix Steps 3, 4, and 5 on each ingroup sample.
# Submit SLURM array job for each haplotype sample in ingroup_haps.txt.

# !!!: to do this in practice:
# 		[ ] move below code to a separate file
#       [ ] uncomment one level
# 		[ ] set email address
# 		[ ] confirm upper limit of `--array` matches length of `ingroup_haps.txt`
# 		[ ] adjust the `source` and `module` commands to fit your installation 

# ---------vvvv SLURM batch submission file vvvvv------

# #SBATCH --mail-type=FAIL,END
# #SBATCH --mail-user=<your@email.here>

# #SBATCH -e arrayjob-%a.err
# #SBATCH -o arrayjob-%a.out

# #SBATCH -t 4:00:00
# #SBATCH --mem 4Gb

# #SBATCH --array=1-1835

# source ~/hmmix_env/bin/activate
# module load python bcftools vcftools

# ingroup_hap=`cat ingroup_haps.txt | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

# echo "Ingroup_hap: "$ingroup_hap

# # hmmix Step 3: create ingroup
# echo "create_ingroup is starting"
# ingroup_file="obs.${ingroup_hap}.txt"
# hmmix create_ingroup -ind=$ingroup_hap -out=$ingroup_file -vcf=external_files/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf -out=obs -ancestral=external_files/homo_sapiens_ancestor_12.fa -weights=external_files/strickmask.bed -outgroup=outgroup.txt

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
