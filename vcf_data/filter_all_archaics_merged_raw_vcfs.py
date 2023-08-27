### Dependencies ###
import gzip
import sys
### sys.argv[1] = unfiltered vcf ###
### sys.argv[2] = chromosome ###


# Define function to filter the archaic joint vcf file.
def all_archaics_merged_raw_vcf(vcf, chrom):
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
    # Intialize the filtered vcf file to be outputted to stdout.
    new_vcf = sys.stdout
    # Intialize a file to print the QC report.
    qc_file = open('./vcf_bookkeeping/all_archaics_merged_raw_qc/chr{0}_qc_report.txt'.format(chrom), 'w')
    # Open the original vcf file.
    with gzip.open(vcf, 'rt') as data:
        # Iterate through every line in the original vcf file.
        for line in data:
            # If the line is a meta info line...
            if line.startswith('#'):
                # Write it to the new vcf.
                new_vcf.write(line)
            # Else...
            else:
                # Split the line by tabs.
                spline = line.split()
                # Grab the refernce and alternative alleles.
                ref = spline[3]
                alt = spline[4]
                # Grab the position.
                pos = int(spline[1])
                # Grab the archaic site info.
                den_gt = spline[-1]
                vin_gt = spline[-2]
                cha_gt = spline[-3]
                alt_gt = spline[-4]
                # Construct a list of all archaic genotype calls.
                arc_gt_list = [alt_gt, cha_gt, vin_gt, den_gt]
                # If the length of refernce and alternative fields is not
                # exactly 2, ie if the site is a not mono- or biallelic SNP...
                if ((len(ref) + len(alt)) != 2):
                    # Print the QC report.
                    qc_file.write('\t'.join([str(chrom), str(pos), '0', 'multi_allelic_or_sv'])+'\n')
                    # Continue to the next line.
                    continue
                # Else-if all archaic sites are missing...
                elif all(gt[0] == '.' for gt in arc_gt_list):
                    # Print the QC report.
                    qc_file.write('\t'.join([str(chrom), str(pos), '1', 'all_arcs_missing'])+'\n')
                    # Continue to the next line.
                    continue
                # Construct a list of GQ scores.
                arc_gq_list = [gt.split(':')[-1] for gt in arc_gt_list]
                # For every archaic...
                for idx in range(4):
                    # If the length of the GQ is one...
                    if (len(arc_gq_list[idx]) == 1):
                        # Set the genotype information to missing.
                        arc_gt_list[idx] = './.:.:.:.:.:.:.:.'
                    # Else-if the archaic has a GQ lower than 40...
                    elif (int(arc_gq_list[idx]) < 40):
                        # Set the genotype information to missing.
                        arc_gt_list[idx] = './.:.:.:.:.:.:.:.'
                # If after GQ QC all of the genotypes are missing...
                if all(gt[0] == '.' for gt in arc_gt_list):
                    # Print the QC report.
                    qc_file.write('\t'.join([str(chrom), str(pos), '2', 'all_arcs_fail_gq_qc'])+'\n')
                    # Continue to the next line.
                    continue
                # Else...
                else:
                    # Update the archaic site information.
                    spline[-4] = arc_gt_list[0]
                    spline[-3] = arc_gt_list[1]
                    spline[-2] = arc_gt_list[2]
                    spline[-1] = arc_gt_list[3]
                    # Write it to the new vcf.
                    new_vcf.write('\t'.join(spline)+'\n')
    # Close the QC report file.
    qc_file.close()
    return


# Filter the given VCF file.
all_archaics_merged_raw_vcf(vcf=str(sys.argv[1]), chrom=str(sys.argv[2]))