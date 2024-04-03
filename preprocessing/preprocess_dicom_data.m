%% MAIN FUNCTION TO PREPROCESS SCAN SERIES
% MT & LKF Last update: 8/23/2020
% input:   pthOutDate - path to proc date folder
%          pthRawDate - path to raw date folder
%          listScans - stuct series (scan folder) list
%          towritedcm (optional) - whether to write corrected images to dicom (optional - assumes false)
% output:  saves a .mat file with documents that could not save in date folder
%
% process: call functions to read dicoms and perform distortion correction,
% each function processes all the scan series in the given date folder
% before returning to main

function listErrors = preprocess_dicom_data(pthOutDate, pthRawDate, listScans, towritedcm)
if ~exist('towritedcm','var')
    towritedcm = false;
end
listErrors = [];

%% PART 1: Read dicoms to *.mgz and *.mat for selected scan types
listErrors = breast_read_dicoms(pthOutDate, pthRawDate, listScans);

%% PART 2: Diffusion Artifact Correction

% Note: for folders AX_MMIL_MULTISHELL and MMIL_MULTISHELL_DIFFUSION
listRSI = struct('pthRSI',{fullfile(pthOutDate,'fullFOV_RSI'), fullfile(pthOutDate,'reducedFOV_RSI')});
% listRSI = struct('pthRSI', {fullfile(pthOutDate)});

for ii = 1:length(listRSI)
    pthOutRSI = listRSI(ii).pthRSI;
    % part 2.1: Distortion Correction
    % initialize parameters kvec and lvec
    kvec = 30:15:60; lvec = 500:500:1500; %init params based on breast
    % kvec = 30:5:40; lvec = 1000:500:2500; %init params based on prostrate
    % kvec = 30:15:60; lvec = 500:500:2500; %init params based on cervix
    % kvec =45; lvec=1000; %testing params to run code fast
    listErrors = [listErrors, AC_dist_corr(pthOutRSI,kvec,lvec)];

    % part 2.2/3/4: Noise Correction, Eddy Current Correction, Average across diffusion directions
    % --- OPTION A: to perform each step independently (input vol loaded from file)
    listErrors = [listErrors, AC_noise(pthOutRSI)];
    listErrors = [listErrors, AC_eddy(pthOutRSI)];
    listErrors = [listErrors, AC_average(pthOutRSI)];

    % --- OPTION B: to create database from scratch (runs faster because there is less file loading)
%     listErrors = [listErrors, AC_noise_eddy_average(pthOutRSI)];     
%     listErrors = [listErrors, AC_diff_uw_nc(pthOutRSI)];
%     listErrors = [listErrors, AC_diff_nc_eddy(pthOutRSI)];

    % part 2.4: Write back to a dicom files
    %if (towritedcm)
        % resave as a dicom
    %end

    % part 2.5: Normalize
    listErrors = [listErrors, AC_norm(pthOutRSI)];
end

%% PART 3: Save Errors mat
if ~(isempty(listErrors)); 
    if exist(fullfile(pthOutDate,'errors.mat'))
        errorsPrev = load(fullfile(pthOutDate,'errors.mat'));
        errorsPrev = errorsPrev.listErrors;
        listErrors = [errorsPrev, listErrors];
    end
    listErrors = sort_errorList(listErrors);
    save(fullfile(pthOutDate,'errors.mat'),'listErrors'); 
    fprintf('Saved: \n\t File: errors.mat\n\t Path: %s\n\n', pthOutDate);
end
end
