### Dependencies ###
import gzip
import sys
### sys.argv[1] = filtered vcf ###
### sys.argv[2] = chromosome ###

# Define function to filter the tgp joint vcf file.
def tgp_arc_eff_seq_len(vcf, chrom):
    """
    ###########################################################################
    INPUT
        vcf: A gzipped unfiltered VCF file where the last four indicies are the
             archaic.
        chrom: Chromosome for the VCF.
    ---------------------------------------------------------------------------
    OUTPUT: Filtered VCF file to standard out.
    ###########################################################################
    """
    # Load the fasta reference for the ancestral allele calls
    anc_fa = './epo_calls/homo_sapiens_ancestor_{0}.fa'.format(chrom)
    # Intialize an empty string for storing the data.
    anc_seq = ''
    # Open the fasta file.
    with open(anc_fa) as data:
        # For every line in data...
        for line in data:
            # If the line is a header...
            if line.startswith('>'):
                # Continue to the next line.
                continue
            # Else...
            else:
                # Append the line to the sequence string.
                anc_seq += line.strip().replace(' ','')
    # Intialize the filtered vcf file to be outputted to stdout.
    new_vcf = sys.stdout
    # Intialize a file to print the QC and all sites data report.
    qc_file = open('./vcf_bookkeeping/tgp_all_archaics_merged_all_sites_qc/chr{0}_qc_report.txt'.format(chrom), 'w')
    sites_file = open('./vcf_bookkeeping/tgp_all_archaics_merged_all_sites_qc/chr{0}_all_sites_report.txt'.format(chrom), 'w')
    # Open the original vcf file.
    with gzip.open(vcf, 'rt') as data:
        # Iterate through every line in the original vcf file.
        for line in data:
            # If the line is a meta info line...
            if line.startswith('##'):
                # Write it to the new vcf.
                new_vcf.write(line)
            # Else-if the line is the header line...
            elif line.startswith('#'):
                # Split the line by tabs.
                spline = line.split()
                # Append the ancestor column.
                spline.append('Ancestor')
                # Write it to the new vcf.
                new_vcf.write('\t'.join(spline)+'\n')
            # Else...
            else:
                # Split the line by tabs.
                spline = line.split()
                # Grab the refernce and alternative alleles.
                ref = spline[3]
                alt = spline[4]
                # Grab the position.
                pos = int(spline[1])
                # Grab the genotype info field.
                info = spline[7]
                # Grab the archaic site info.
                den_gt = spline[-1]
                vin_gt = spline[-2]
                cha_gt = spline[-3]
                alt_gt = spline[-4]
                # Construct a list of all archaic genotype calls.
                arc_gt_list = [alt_gt, cha_gt, vin_gt, den_gt]
                # Determine the ancestral allele.
                anc_allele = anc_seq[pos-1]
                # If the ancestral allele is not known.
                if ((anc_allele == '.') | (anc_allele == '-') | (anc_allele == 'N')):
                    # Print the QC report.
                    qc_file.write('\t'.join([str(chrom), str(pos), '3', 'no_anc_call'])+'\n')
                    # Continue to the next line.
                    continue
                # Else-if the ancestral allele does not match the refernce or alternate
                # allele and there is an allele call for alternate allele...
                elif ((anc_allele.upper() != ref) & (anc_allele.upper() != alt) & (alt != '.')):
                    # Print the QC report.
                    qc_file.write('\t'.join([str(chrom), str(pos), '4', 'not_snp_after_anc_call'])+'\n')
                    # Continue to the next line.
                    continue
                # Else...
                else:
                    # Intialize an all sites report list (chr, pos, alt, cha, vin, den).
                    # 0 = missing, 1 = called
                    report_list = [str(chrom), str(pos), '1', '1', '1', '1']
                    # For each archaic...
                    for idx in range(4):
                        # If the archaic has missing data.
                        if arc_gt_list[idx][0] == '.':
                            # Record the misisng archaic.
                            report_list[idx+2] = '0'
                    # Print the all sites summary.
                    sites_file.write('\t'.join(report_list)+'\n')
                    # If the ancestral allele is the refernce allele...
                    if (anc_allele.upper() == ref):
                        # Set ancestral genotype to be homozygous refernce.
                        anc_gt = '0|0:.:.:.:.:.:.:.'
                        # If the tgp has an allele call...
                        if ('AA' in info):
                            # Append the Ancestor column.
                            spline.append(anc_gt)
                            # Write it to the new vcf.
                            new_vcf.write('\t'.join(spline)+'\n')
                        # Else...
                        else:
                            # Append the info field.
                            spline[7] = spline[7]+';AA={0}|||'.format(anc_allele)
                            # Append the Ancestor column.
                            spline.append(anc_gt)
                            # Write it to the new vcf.
                            new_vcf.write('\t'.join(spline)+'\n')
                    # Else-if the ancestral allele is the alternative allele...
                    elif (anc_allele.upper() == alt):
                        # Set ancestral genotype to be homozygous alternative.
                        anc_gt = '1|1:.:.:.:.:.:.:.'
                        # If the tgp has an allele call...
                        if ('AA' in info):
                            # Append the Ancestor column.
                            spline.append(anc_gt)
                            # Write it to the new vcf.
                            new_vcf.write('\t'.join(spline)+'\n')
                        # Else...
                        else:
                            # Append the info field.
                            spline[7] = spline[7]+';AA={0}|||'.format(anc_allele)
                            # Append the Ancestor column.
                            spline.append(anc_gt)
                            # Write it to the new vcf.
                            new_vcf.write('\t'.join(spline)+'\n')
    # Close the QC and all sites report file.
    qc_file.close()
    sites_file.close()
    return

# Add the ancestral allele calls.
tgp_arc_eff_seq_len(vcf=str(sys.argv[1]), chrom=str(sys.argv[2]))