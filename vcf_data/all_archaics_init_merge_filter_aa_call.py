### Dependencies ###
import gzip
import sys

## sys.argv[1] = unfiltered merged archaic vcf ##
## sys.argv[2] = chromosome ##
## sys.argv[3] = dataset prefix ##


# Define a function to qc a line from the unfiltered all archaics merged vcf.
def process_merged_line(line, chrom, anc_seq, dup_set):
    #########################
    ## Extract Information ##
    #########################
    # Split the line by tabs.
    spline = line.split()
    # Grab the refernce and alternative alleles.
    ref = spline[3]
    alt = spline[4]
    # Grab the position.
    pos = int(spline[1])
    # Determine the ancestral allele.
    anc_allele = anc_seq[pos-1]
    # Grab the archaic genotype info.
    den_gt = spline[-1]
    vin_gt = spline[-2]
    cha_gt = spline[-3]
    alt_gt = spline[-4]
    # Construct a list of all archaic genotype info.
    arc_gts = [alt_gt, cha_gt, vin_gt, den_gt]
    # Construct a list of genotype qualities.
    arc_gqs = [gt.split(':')[-1] for gt in arc_gts]
    # Intialize a list of qc'ed archaic genotypes.
    arc_qced_gts = []
    # For every archaic.
    for i, gq in enumerate(arc_gqs):
        # If the length of the genotype quality is one.
        if len(gq) == 1:
            # Set the genotype information to missing as the archaic failed qc.
            arc_qced_gts.append('./.:.:.:.:.:.:.:.')
        # Else-if the archaic has a genotype quality lower than 40.
        elif int(gq) < 40:
            # Set the genotype information to missing as the archaic failed qc.
            arc_qced_gts.append('./.:.:.:.:.:.:.:.')
        # Else, the archaic has a genotype quality of 40 or higher.
        else:
            # Keep the original genotype.
            arc_qced_gts.append(arc_gts[i])
    # Intialize a list to denote if an archaic passed qc (0 = missing, 1 = called).
    arc_qc = ['1', '1', '1', '1']
    # For every archaic.
    for i, qced_gt in enumerate(arc_qced_gts):
        # If the archaic has missing genotype information.
        if qced_gt[0] == '.':
            # Record the misisng archaic.
            arc_qc[i] = '0'
    
    
    #########################
    ## QC Failed Condtions ##
    #########################
    # If the current position has duplicated records.
    if pos in dup_set:
        # Return the qc information.
        return '\t'.join([str(chrom), str(pos), '0', 'dup_record'])+'\n', None, None
    # Else-if the length of refernce and alternative fields is not exactly 2,
    # ie if the site is a not mono- or biallelic snp.
    elif (len(ref) + len(alt)) != 2:
        # Return the qc information.
        return '\t'.join([str(chrom), str(pos), '1', 'multi_allelic_or_sv'])+'\n', None, None
    # Else-if all archaic genotypes are missing.
    elif all(gt[0] == '.' for gt in arc_gts):
        # Return the qc information.
        return '\t'.join([str(chrom), str(pos), '2', 'all_archaics_missing_gt'])+'\n', None, None
    # Else-if all archaics genotype quality is less than 40.
    elif all(gt[0] == '.' for gt in arc_qced_gts):
        # Return the qc information.
        return '\t'.join([str(chrom), str(pos), '3', 'all_archaics_low_gq'])+'\n', None, None
    # Else-if the ancestral allele is not known.
    elif (anc_allele == '.') | (anc_allele == '-') | (anc_allele == 'N'):
        # Return the qc information.
        return '\t'.join([str(chrom), str(pos), '4', 'no_anc_call'])+'\n', None, None
    # Else-if the ancestral allele does not match the refernce or alternate allele
    # and there is an allele call for alternate allele.
    elif (anc_allele.upper() != ref) & (anc_allele.upper() != alt) & (alt != '.'):
        # Return the qc information.
        return '\t'.join([str(chrom), str(pos), '5', 'not_snp_after_anc_call'])+'\n', None, None
    
    ###############################
    ## Invariant Site Conditions ##
    ###############################
    # Else-if the reference allele is the ancestral allele and there is no alternative allele.
    elif (anc_allele.upper() == ref) & (alt == '.'):
        # Return the all sites information.
        return None, '\t'.join([str(chrom), str(pos), arc_qc[0], arc_qc[1], arc_qc[2], arc_qc[3], '0', 'all_anc_ref_anc'])+'\n', None
    # Else-if the reference allele is not the ancestral allele and there is no alternative allele.
    elif (anc_allele.upper() != ref) & (alt == '.'):
        # Return the all sites information.
        return None, '\t'.join([str(chrom), str(pos), arc_qc[0], arc_qc[1], arc_qc[2], arc_qc[3], '1', 'all_der_alt_anc'])+'\n', None
    
    #############################
    ## Variant Site Conditions ##
    #############################
    # Else-if the reference allele is the ancestral allele and there is an alternative allele.
    elif (anc_allele.upper() == ref) & (alt != '.'):
        # Update the info field.
        spline[7] = f'AA={anc_allele}|||'
        # Update the archaic genotype information.
        spline[-4] = arc_qced_gts[0]
        spline[-3] = arc_qced_gts[1]
        spline[-2] = arc_qced_gts[2]
        spline[-1] = arc_qced_gts[3]
        # Append the ancestral genotype.
        spline.append('0|0:.:.:.:.:.:.:.')
        # Return the all sites informtion and the new vcf line.
        return None, '\t'.join([str(chrom), str(pos), arc_qc[0], arc_qc[1], arc_qc[2], arc_qc[3], '2', 'arc_var_ref_anc'])+'\n', '\t'.join(spline)+'\n'
    # Else-if the alternative allele is the ancestral allele.
    elif (anc_allele.upper() == alt):
        # Update the info field.
        spline[7] = f'AA={anc_allele}|||'
        # Update the archaic genotype information.
        spline[-4] = arc_qced_gts[0]
        spline[-3] = arc_qced_gts[1]
        spline[-2] = arc_qced_gts[2]
        spline[-1] = arc_qced_gts[3]
        # Append the ancestral genotype.
        spline.append('1|1:.:.:.:.:.:.:.')
        # Return the all sites informtion and the new vcf line.
        return None, '\t'.join([str(chrom), str(pos), arc_qc[0], arc_qc[1], arc_qc[2], arc_qc[3], '3', 'arc_var_alt_anc'])+'\n', '\t'.join(spline)+'\n'

# Define a function to filter the unfiltered all archaics merged vcf.
def filter_init_merge_vcf(vcf, chrom, prefix):
    # Inialize a set to store duplicated positions.
    dup_set = set()
    # Open the duplicated records regions file.
    with gzip.open(f'./dups/{prefix}_dup_sites_chr{chrom}.txt.gz', 'rt') as dup_data:
        # For every line in the duplicated records regions file.
        for line in dup_data:
            # Update the duplicated positions set.
            dup_set.add(int(line.split()[1]))
    # Intialize an empty string for storing the data.
    anc_seq = ''
    # Open the fasta file.
    with open(f'./epo_calls/homo_sapiens_ancestor_{chrom}.fa') as fa_data:
        # For every line in the fasta file.
        for line in fa_data:
            # If the line is a header.
            if line.startswith('>'):
                # Continue to the next line.
                continue
            # Else, the line contains nucleotides.
            else:
                # Append the line to the sequence string.
                anc_seq += line.strip().replace(' ','')
    # Intialize a list to store the the qc information.
    qc_file_lines = [] # Format: chrom, pos, id, qc.
    # Intialize a list to store the all sites information.
    all_sites_lines = [] # Format: chrom, pos, alt, cha, vin, den, id, qc.
    # Intialize a list to store the new vcf information.
    new_vcf_lines = []
    # Read the vcf file and intilialize the qc and all sites report.
    with gzip.open(vcf, 'rt') as vcf_data, \
         gzip.open(f'./bookkeeping/{prefix}_aa_calls_failed_qc_chr{chrom}.txt.gz', 'wt') as qc_file, \
         gzip.open(f'./bookkeeping/{prefix}_aa_calls_all_sites_chr{chrom}.txt.gz', 'wt') as sites_file:
        # Iterate through every line in the original vcf file.
        for line in vcf_data:
            # If the line contains meta information.
            if line.startswith('##'):
                # Append the line to the new vcf file list.
                new_vcf_lines.append(line)
            # Else-if, the line contains the header information
            elif line.startswith('#'):
                # Split the line by tabs.
                spline = line.split()
                # Append the ancestor column.
                spline.append('Ancestor')
                # Append the line to the new vcf file list.
                new_vcf_lines.append('\t'.join(spline)+'\n')
            # Else, the line contains genotype information.
            else:
                # Process the current line.
                qc_info, site_info, vcf_info = process_merged_line(line, chrom, anc_seq, dup_set)
                # If the current line failed qc.
                if qc_info:
                    # Append the qc file list.
                    qc_file_lines.append(qc_info)
                # Else-if the current line passed qc but is invariant.
                elif site_info and not vcf_info:
                    # Append the all sites list.
                    all_sites_lines.append(site_info)
                # Else, the current line passed qc and is variable.
                else:
                    # Append the all sites list.
                    all_sites_lines.append(site_info)
                    # Append the new vcf file list.
                    new_vcf_lines.append(vcf_info)
            # If there are 50_000 lines in the qc file list.
            if len(qc_file_lines) == 50_000:
                # Write the qc lines to the qc file in chunks to improve performance.
                qc_file.writelines(qc_file_lines)
                # Clear the written qc lines.
                qc_file_lines.clear()
            # If there are 50_000 lines in the all sites list.
            if len(all_sites_lines) == 50_000:
                # Write the qc lines to the qc file in chunks to improve performance.
                sites_file.writelines(all_sites_lines)
                # Clear the written all sites lines.
                all_sites_lines.clear()
            # If there are 50_000 lines in the new vcf file list.
            if len(new_vcf_lines) == 50_000:
                # Write the lines that passed qc to the new vcf file in chunks to improve performance.
                sys.stdout.writelines(new_vcf_lines)
                # Clear the written vcf lines.
                new_vcf_lines.clear()
        # If there are still qc lines to be written.
        if qc_file_lines:
            # Write the remaining qc lines.
            qc_file.writelines(qc_file_lines)
        # If there are still all sites lines to be written.
        if all_sites_lines:
            # Write the remaining all sites lines.
            sites_file.writelines(all_sites_lines)
        # If there are still new vcf lines to be written.
        if new_vcf_lines:
            # Write the remaining new vcf lines.
            sys.stdout.writelines(new_vcf_lines)
        
# Filter the unfiltered all archaics merged vcf.
filter_init_merge_vcf(vcf=sys.argv[1], chrom=sys.argv[2], prefix=sys.argv[3])