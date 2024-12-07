### Dependencies ###
import gzip
import sys

## sys.argv[1] = vcf ##
## sys.argv[2] = chromosome ##
## sys.argv[3] = dataset prefix ##


# Define a function to identify what duplicated entires from a vcf.
def identify_dups(vcf, chrom, prefix):
    # Intialize a list to store the regions file information.
    dup_file_lines = []
    # Intialize a list to store the duplicated records vcf information.
    dup_vcf_lines = []
    # Intilailize the previous line and position.
    prev_line = None
    prev_pos = None
    # Inialize a set to store duplicated positions.
    dup_set = set()
    # Read the vcf file and intilialize the duplicated sites regions file.
    with gzip.open(vcf, 'rt') as data, \
         gzip.open(f'./dups/{prefix}_dup_sites_chr{chrom}.txt.gz', 'wt') as dup_file:
        # Iterate through every line in the original vcf file.
        for line in data:
            # If the line contains meta or header information.
            if line.startswith('#'):
                # Write the line to the duplicated sites vcf list.
                dup_vcf_lines.append(line)
             # Else, the line contains genotype information.
            else:
                # Grab the current poistion.
                pos = line.split()[1]
                # If the current position is already known to be duplicated.
                if pos in dup_set:
                    # Append the current line to the duplicated sites vcf file list.
                    dup_vcf_lines.append(line)
                # Else-if the current position is the same as the previous position.
                elif pos == prev_pos:
                    # Add the position to the duplicated set.
                    dup_set.add(pos)
                    # Append the duplicated position to the regions file list.
                    dup_file_lines.append(chrom + '\t' + pos + '\n')
                    # Append the previous line to the duplicated sites vcf file list.
                    dup_vcf_lines.append(prev_line)
                    # Append the current line to the duplicated sites vcf file list.
                    dup_vcf_lines.append(line)
                # Else, the current position is unique.
                else:
                    # Update the previous line and position.
                    prev_line = line
                    prev_pos = pos
            # If there are 2_500 lines in the duplicated records regions file list.
            if len(dup_file_lines) == 2_500:
                # Write the duplicated sites information to the duplicated regions file
                # in chunks to improve performance.
                dup_file.writelines(dup_file_lines)
                # Clear the written lines.
                dup_file_lines.clear()
            # If there are 2_500 lines in the duplicated records vcf file list.
            if len(dup_vcf_lines) == 2_500:
                # Write the duplicted vcf lines to the vcf file in chunks to improve performance.
                sys.stdout.writelines(dup_vcf_lines)
                # Clear the written lines.
                dup_vcf_lines.clear()
        # If there are still duplicated records regions file lines to be written.
        if dup_file_lines:
            # Write the remaining lines.
            dup_file.writelines(dup_file_lines)
        # If there are still duplicated records vcf lines to be written.
        if dup_vcf_lines:
            # Write the remaining lines.
            sys.stdout.writelines(dup_vcf_lines)


# Generate the duplicated sites vcf and regions file.
identify_dups(vcf=sys.argv[1], chrom=sys.argv[2], prefix=sys.argv[3])