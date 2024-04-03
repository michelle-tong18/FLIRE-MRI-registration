%% WRAPPER SCRIPT
% MT Last update: 2022

% This portion of the code registeres follow-up T1-weighted breast MRI to the baseline timepoint. 
% 1. Preprocessing will normalize input images by the 98th percentile signal intensity.
% 2. An intial registration will be performed which is the weighted sum of three affine registrations.
% 3. A nonlinear registration will be performed.

% Inputs: T1-weighted Breast MRI. The expected directory structure is as follows: 
%     - `/.../patient_id/date/T1_noFS_PRE.mat`
% Outputs: Saved displacement field, registered T1-weighted images, and verbose data used during the registration. The output directory structure is as follows:
%     - `/.../patient_id/method_name/displacement_field.mat`      
%     - `/.../patient_id/method_name/T1_noFS_PRE_reg.mat`
%     - `/.../patient_id/method_name/datastructs_method_name.mat`

%% Part 0: Set Up Matlab ---------------------------------------------------
% Set the path to where all the Matlab dependencies are located
addpath(genpath(pwd));

%% Part 1: Identify Data ---------------------------------------------------
% Part 1.1: Initialize Data
% Specify path to directory containing preprocessed MR images that are saved as .mat files.
% This directory contains `/dirname_proc/patient_id/date/T1_noFS_PRE.mat`
dirname_proc = '/directory_to_where_the_source_images_are_located';
% Specify the directory to save outputs.
% Outputs will be saved in /outdir/patient_id/method_name/*.mat`
outdir = '/directory_to_where_the_output_should_be_directed';

dirlist = listdir(pthRawData,'dirs','cell');
% dirlist = {...
%     'patiend_id1',...
%     'patient_id2', ... 
%     };

veclist=[1:length(dirlist)]; %idx of patients in dirlist for analysis
dirlist(veclist)

forceflag = true;          %Options: true, false, 0, 1
method = 'FLIRE_for_github';


%% Part 2: Registration
%for each patient
for veci = 1:length(veclist)
    diri = veclist(veci);
    status = FLIRE_for_T1_Breast_MRI(dirlist,dirname_proc,outdir,method,diri,forceflag)
end
