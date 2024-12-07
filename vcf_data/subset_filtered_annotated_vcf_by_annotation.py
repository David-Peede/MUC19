### Dependencies ###
import gzip
import sys

## sys.argv[1] = chromosome ##
## sys.argv[2] = dataset prefix ##
## sys.argv[3] = coding mutation type ##


# Define a function to extract coding sequence information from a annotated vcf
def process_annotated_line(line, chrom, suffix):
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
    # Grab the info field.
    info = spline[7]
    
    ###############################
    ## Synonmous Site Conditions ##
    ###############################
    # If the current position has a synonymous_variant annotation.
    if (suffix == 'syn') & ('synonymous_variant' in info):
        # Return the qc information.
        return '\t'.join([str(chrom), str(pos), 'synonymous_variant', str(info)])+'\n', line
    # Else-if the current position has a stop_retained_variant annotation.
    elif (suffix == 'syn') & ('stop_retained_variant' in info):
        # Return the qc information.
        return '\t'.join([str(chrom), str(pos), 'stop_retained_variant', str(info)])+'\n', line
    
    ###################################
    ## Non-synonmous Site Conditions ##
    ###################################
    # Else-if the current position has a missense_variant annotation.
    elif (suffix == 'non') & ('missense_variant' in info):
        # Return the qc information.
        return '\t'.join([str(chrom), str(pos), 'missense_variant', str(info)])+'\n', line
     # Else-if the current position has a start_lost annotation.
    elif (suffix == 'non') & ('start_lost' in info):
        # Return the qc information.
        return '\t'.join([str(chrom), str(pos), 'start_lost', str(info)])+'\n', line
     # Else-if the current position has a stop_lost annotation.
    elif (suffix == 'non') & ('stop_lost' in info):
        # Return the qc information.
        return '\t'.join([str(chrom), str(pos), 'stop_lost', str(info)])+'\n', line
     # Else-if the current position has a initiator_codon_variant annotation.
    elif (suffix == 'non') & ('initiator_codon_variant' in info):
        # Return the qc information.
        return '\t'.join([str(chrom), str(pos), 'initiator_codon_variant', str(info)])+'\n', line
     # Else-if the current position has a stop_gained annotation.
    elif (suffix == 'non') & ('stop_gained' in info):
        # Return the qc information.
        return '\t'.join([str(chrom), str(pos), 'stop_gained', str(info)])+'\n', line
    
    ################################
    ## Non-coding Site Conditions ##
    ################################
    else:
        return None, None
    

# Define a function to subset a filtered and annotated vcf by coding mutation type.
def subset_vcf_by_annotation(chrom, prefix, suffix):
    # Intialize the annoatated vcf file.
    vcf = f'./filtered_merge_ann/{prefix}_all_annotations_chr{chrom}.vcf.gz'
    # Intialize a list to store the the qc information.
    qc_file_lines = [] # Format: chrom, pos, id, info.
    # Intialize a list to store the new vcf information.
    new_vcf_lines = []
    # Read the vcf file and intilialize the qc sites report.
    with gzip.open(vcf, 'rt') as vcf_data, \
         gzip.open(f'./ann_summary/{prefix}_{suffix}_mut_sites_chr{chrom}.txt.gz', 'wt') as qc_file:
        # Iterate through every line in the original vcf file.
        for line in vcf_data:
            # If the line contains meta information.
            if line.startswith('#'):
                # Append the line to the new vcf file list.
                new_vcf_lines.append(line)
            # Else, the line contains genotype information.
            else:
                # Process the current line.
                qc_info, vcf_info = process_annotated_line(line, chrom, suffix)
                # If the current line containes coding information.
                if qc_info:
                    # Append the qc file list.
                    qc_file_lines.append(qc_info)
                    # Append the new vcf file list.
                    new_vcf_lines.append(vcf_info)
            # If there are 50_000 lines in the qc file list.
            if len(qc_file_lines) == 50_000:
                # Write the qc lines to the qc file in chunks to improve performance.
                qc_file.writelines(qc_file_lines)
                # Clear the written qc lines.
                qc_file_lines.clear()
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
        # If there are still new vcf lines to be written.
        if new_vcf_lines:
            # Write the remaining new vcf lines.
            sys.stdout.writelines(new_vcf_lines)
        
# Filter the unfiltered tgp/sgdp phased archaic merged vcf.
subset_vcf_by_annotation(chrom=sys.argv[1], prefix=sys.argv[2], suffix=sys.argv[3])