%% WRAPPER SCRIPT
% MT & LKF Last update: 9/1/2020
% Primary script to run - sorts through patients/dates and calls main.m to
% process scan folders

% This portion of the code will read in dicom files to extract image volumes as save 
% in a standardized format for downstream processing. This script is able to save 
% multiple files from a breast MRI exam. However for this project we are mainly 
% interested in T1-weighted breast MR images without fat saturation.

% 1. Data will be converted from dicoms to .mat files. 
% 2. (Optional) Artifact correction will be performed on diffusion data.

% Inputs: The expected directory structure is as follows: 
%     - `/.../patient_id/date/scan_series/*.dcm`
% Relevant Outputs: The output directory structure is as follows: 
%     - `/.../patient_id/date/T1_noFS_PRE.mat`

%% add path to relevant scripts
addpath(genpath(pwd));

%% initialize variables 
global force_save
force_save = false; % 1 (true) to reprocess all data
towritedcm = false; % 1 (true) to write dicoms for *RSI_uw_nc_eddy_averaged data

%% PART I - Select project folder (i.e. "Breast") / UCSD Data
% Specify path to directory containing MR images saved as dicom files.
% This directory contains `/pthRawData/patient_id/date/scan_series/*.dcm`
pthRawData = '/directory_to_where_the_source_images_are_located';
% Specify the directory to save MR images as .mat files.
% Outputs will be saved in `/pthOut/patient_id/date/*.mat`
% Longitudinal registration uses the files `/pthOut/patient_id/date/T1_noFS_PRE.mat`
pthOut = '/directory_to_where_the_output_should_be_directed';

% list of patient directory names as a struct with 'name' field
listPatients = listdir(pthRawData,'dirs','struct');
% patientList = {...
%     'patiend_id1',...
%     'patient_id2', ... 
%     }; 
% listPatients = cell2struct(patientList,'name',1);

if ~(exist(pthOut)); mkdir(pthOut); end
pthQA = fullfile(pthOut,'QA');
if ~(exist(pthQA)); mkdir(pthQA); end
pthQAMask = fullfile(pthQA,'mask');
if ~(exist(pthQAMask)); mkdir(pthQAMask); end
pthQADistCorr = fullfile(pthQA,'dist_corr');
if ~(exist(pthQADistCorr)); mkdir(pthQADistCorr); end

%% PART II - Sort Files / UCSD Data
for i = 1:length(listPatients)
    pthRawPatient = fullfile(pthRawData,listPatients(i).name);
    pthOutPatient = fullfile(pthOut,listPatients(i).name);
    % create patient folder with identical name in proc2
    if ~exist(pthOutPatient); mkdir(pthOutPatient); end
    
    listDates = listdir(pthRawPatient,'dirs','struct');
    for j = 1:length(listDates)
        pthRawDate = fullfile(pthRawPatient,listDates(j).name);
        pthOutDate = fullfile(pthOutPatient,listDates(j).name);
        % create date folder with identical name in proc2
        if ~exist(pthOutDate); mkdir(pthOutDate); end
        
        %make subfolders
        subfolders = {fullfile(pthOutDate,'Ax_DWI_ASSET'),...
            fullfile(pthOutDate,'fullFOV_RSI'),...
            fullfile(pthOutDate,'reducedFOV_RSI'),...
            fullfile(pthOutDate,'Ax_DWI_ASSET','dist_corr'),...
            fullfile(pthOutDate,'fullFOV_RSI','dist_corr'),...
            fullfile(pthOutDate,'reducedFOV_RSI','dist_corr')};
        for k = 1:length(subfolders)
            if ~(exist(subfolders{k})); mkdir(subfolders{k}); end
        end
                
        listScans = listdir(fullfile(pthRawPatient,listDates(j).name),'dirs','struct');
        % run main on listScan, pthOutDate
        fprintf('Patient %d, Date %d',i,j)
        preprocess_dicom_data(pthOutDate, pthRawDate, listScans, towritedcm);
        
        %removes empty subfolders, do not add 's' flag
        for k = length(subfolders):-1:1
            status = rmdir(subfolders{k});
        end
        status = rmdir(pthOutDate);
    end
    status = rmdir(pthOutPatient);
end
beep
disp('Completed.');
