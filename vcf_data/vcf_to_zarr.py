### Dependencies ###
import allel
import numcodecs
import numpy as np
import sys
import zarr

### sys.argv[1] = vcf file prefix ###
### sys.argv[2] = zarr file prefix ###
### sys.argv[3] = chromosome number ###


# Define file paths.
vcf_path = '{0}.vcf.gz'.format(str(sys.argv[1]))
zarr_path = './zarr_data/{0}.zarr'.format(str(sys.argv[2]))

# Convert the vcf file to a zarr array.
allel.vcf_to_zarr(vcf_path, zarr_path, group=str(sys.argv[3]), fields=['GT', 'POS'], log=sys.stdout, overwrite=True)