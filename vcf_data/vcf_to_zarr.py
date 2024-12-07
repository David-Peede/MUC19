### Dependencies ###
import allel
import numcodecs
import numpy as np
import sys
import zarr

### sys.argv[1] = vcf file ###
### sys.argv[2] = zarr file prefix ###
### sys.argv[3] = chromosome number ###


# Convert the vcf file to a zarr array.
allel.vcf_to_zarr(
    str(sys.argv[1]), f'../zarr_data/{sys.argv[2]}.zarr', group=str(sys.argv[3]), 
    fields=['GT', 'POS'], log=sys.stdout, overwrite=True,
)