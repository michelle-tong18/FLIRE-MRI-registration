# FLIRE-MRI-registration

## Description
Fast Longituinal Image REgistration (FLIRE): a longitudinal registration algorithm to align T1-weighted Breast MRI acquired at multiple timepoints.
FLIRE registration consists of (1) preprocessing, (2) an intial linear registration, and (3) a nonlinear deformable registration. The algorithm uses 
signal intensity matching and Gaussian smoothing to calculate three displacement fields. Each displacement field corresponds to movement in one 
dimension (dx, dy, dz) for each voxel. Gaussian smoothing heavily blurs both the moving and target images by convolution with an isotropic Gaussian 
kernel of a specified standard deviation. Smoothing is reduced in each successive iteration so that smaller differences in signal intensity become 
distinguishable. This technique aligns the most prominent signal intensity features before tuning finer details, which reduces the likelihood of 
converging to a suboptimal solution. The FLIRE deformable registration algorithm was adapted from a non-linear registration technique for brain by 
Holland and Dale [[1](https://doi.org/10.1016/j.media.2011.02.005)]. 

## Getting Started

### Dependencies
The scripts in this repository have been tested with Matlab R2020a running in Rocky Linux release 8.8 only.

If not using linux os, download the Matlab scripts and related dependencies for `movepixels.m`, `movepixels_3d_double.c`, and `movepixels_3d_single.c`
from [Mathworks File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/21451-multimodality-non-rigid-demon-algorithm-image-registration?s_tid=srchtitle).


### Installing
1. Set up a github API Key
2. Clone the repo
```
git clone git@github.com:michelle-tong18/FLIRE-MRI-registration.git
```

## Running the algorithm

### Preprocess longitudinal MRI data. 
This portion of the code will read in dicom files to extract image volumes as save in a standardized format for downstream processing.
This script is able to save multiple files from a breast MRI exam. However for this project we are mainly interested 
in T1-weighted breast MR images without fat saturation.
1. Data will be converted from dicoms to .mat files. 
2. (Optional) Artifact correction will be performed on diffusion data.
* Inputs: The expected directory structure is as follows: 
    * `/.../patient_id/date/scan_series/*.dcm`
* Relevant Outputs: The output directory structure is as follows: 
    * `/.../patient_id/date/T1_noFS_PRE.mat`

Modify the script **`run_preprocessing.m`** to define project specific variables. Then run to preprocess data.

### Register longitudinal MRI data.
This portion of the code registeres follow-up T1-weighted breast MRI to the baseline timepoint. 
1. Preprocessing will normalize input images by the 98th percentile signal intensity.
2. An intial registration will be performed which is the weighted sum of three affine registrations.
3. A nonlinear registration will be performed.
* Inputs: T1-weighted Breast MRI. The expected directory structure is as follows: 
    * `/.../patient_id/date/T1_noFS_PRE.mat`
* Outputs: Saved displacement field, registered T1-weighted images, and verbose data used during the registration. The output directory structure is as follows:
    * `/.../patient_id/method_name/displacement_field.mat`      
    * `/.../patient_id/method_name/T1_noFS_PRE_reg.mat`
    * `/.../patient_id/method_name/datastructs_method_name.mat`

Modify the script **`run_FLIRE_registration.m`** to define project specific variables. Then run to register data.

## Contributors
* Michelle Tong - mwtong@ucsd.edu - preprocessing and registration code
* Anders Dale - registration and legacy code
* Lauren Fang - preprocessing code
* Nathan White - legacy code
* Hon Yu - consolidating code

## Acknowledgements
* This work was supported by the California Breast Cancer Research Program Early Career Award [25IB-0056]; 
GE Healthcare; NIH [R37 CA249659]; and Kruger Wyeth Research fund.
* Code snippets from various sources including but not limited to:
    * [Mathworks File Exchange - movepixels](https://www.mathworks.com/matlabcentral/fileexchange/21451-multimodality-non-rigid-demon-algorithm-image-registration?s_tid=srchtitle)
    * [Mathworks File Exchange - imshow3D](https://www.mathworks.com/matlabcentral/fileexchange/41334-imshow3d)
    * [Mathworks File Exchange - imshow3Dfull](https://www.mathworks.com/matlabcentral/fileexchange/47463-imshow3dfull)
    * Distortion correction code for diffusion MRI

## License
Distributed under the MIT License. See LICENSE.txt for more information.

## Citations
* Tong MW, Yu HJ, Andreassen MMS, Loubrie S, Rodriguez-Soto AE, Seibert TM, Rakow-Penner R, Dale M. Longitudinal registration of T1-weighted breast MRI: A registration algorithm (FLIRE) and clinical application.
MRI 2024;113.doi: [10.1016/j.mri.2024.110222](https://doi.org/10.1016/j.mri.2024.110222)
* Holland D, Dale A. Alzheimer's 
Disease Neuroimaging Initiative.  Nonlinear registration of longitudinal images and measurement of change in regions of interest. 
Med Imag Anal 2011;15(4):489-97.doi: [10.1016/j.media.2011.02.005](https://doi.org/10.1016/j.media.2011.02.005)
