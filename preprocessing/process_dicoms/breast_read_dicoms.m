%% MAIN FUNCTION TO READ DICOMS
% MT & LKF Last update: 9/11/2020

% input:  pthOutDate - path to proc folder
%         pthRawDate - path to raw date folder
%         listScans - struct of scan folders in date folder
% output: .mgz and .mat files

%-----Summary-----
% folder                    --> file name (.mat and .mgz unless specified)
% AX_T2_FSE                 --> T2_FS_PRE       (T2 fat suppression)
% AX_DYN_PRE                --> T1_FS_PRE       (T1 fat suppression)
% AX_NFS_PRE                --> T1_noFS_PRE     (T1 no fat suppression)
% AX_MMIL_MULTISHELL        --> reducedFOV_RSI  (reduced field of view - does not include b0rvs)
%                    
%--> reducedFOV_RSI_b0_fwd.mgz,reducedFOV_RSI_b0_rvs
%                           --> reducedFOV_bVals.mat, reducedFOV_qmat.mat
% MMIL_MULTISHELL_DIFFUSION --> fullFOV_RSI     (normal field of view - does not include b0rvs)
%                           --> fullFOV_RSI_b0_fwd.mgz,fullFOV_RSI_b0_rvs
%                           --> fullFOV_bVals.mat, fullFOV_qmat.mat
% Ph#_AX_DYN for #1-5       --> DCE.mat         (phase encoding directions)
%                           --> DCE_Ph2_SUB.mgz
%----- end of summary -----

%special notes
%   if RSI.mgz and .mat are already saved, code will not save b0.mgz or bVals.mat
%   unless force_save == 1 -> this will force run all code

function listErrors = breast_read_dicoms(pthOutDate, pthRawDate, listScans)
listErrors = []; % for Errors struct
global force_save

% for each scan folder ...
for i = 1:numel(listScans)
    scanName = listScans(i).name;
    % determine the full folder path to the raw dicoms
    pthRawScan = fullfile(pthRawDate, scanName);
    
    %----------------------------------------------------------------------
    % PART 1A: determine whether to process folder from the first image
    % list all MRDC and DCM images
    listImgs = [dir(fullfile(pthRawScan,'*.dcm'))];
    if (isempty(listImgs))
        continue
    end
    
    % determine the scan type from the first image
    try
        info = dicominfo(fullfile(pthRawScan,listImgs(1).name));
        % determine the scan description, replace spaces in name with '_'
        description = strrep(info.SeriesDescription,' ','_');
    catch
        description = 'unknown';
    end
    % determine whether to process folder and set name to save file
    switch description
        case 'Ax_T2_FSE_FS';                      saveName='T2_FS_PRE'; pthOutSave = pthOutDate;
        case 'AX_DYN_PRE_TEST_2.4_326_62_BW83';   saveName='T1_FS_PRE'; pthOutSave = pthOutDate;
        case 'AX_NFS_PRE_2.4_400_BW83_62_FA_15';  saveName='T1_noFS_PRE'; pthOutSave = pthOutDate;
        case 'Ph1/AX_DYN_MP_326_2.4_62_NO_DELAY'; saveName='DCE'; pthOutSave = pthOutDate;
        case 'MMIL_MULTISHELL_DIFFUSION';         saveName='fullFOV_RSI'; pthOutSave = fullfile(pthOutDate,'fullFOV_RSI');
        case 'AX_MMIL_MULTISHELL';                saveName='reducedFOV_RSI'; pthOutSave = fullfile(pthOutDate,'reducedFOV_RSI');
        case 'Ax_Focus_DWI';                      saveName='RSI'; pthOutSave = pthOutDate;
        case 'Ax_DWI_ASSET';                      saveName='Ax_DWI_ASSET'; pthOutSave = fullfile(pthOutDate,'Ax_DWI_ASSET');
        otherwise;                                continue;
    end
    %----------------------------------------------------------------------
    % PART 1B: check if images are already saved
    % check if images were already saved to determine whether to read dicoms
    if ~exist(fullfile(pthOutSave,[saveName '.mgz'])) && (isunix || ismac)
    elseif ~exist(fullfile(pthOutSave,[saveName '.mat']))
    elseif force_save;
    else % images already created, don't save
        continue; % use continue to avoid errors saving multiple times
    end
    %----------------------------------------------------------------------
    % PART 1C: Check for special cases
    %condition to disgard reduced FOV data
    if strcmp(saveName,'reducedFOV_RSI') && (info.PercentPhaseFieldOfView > 70)
        msg = sprintf('Error: Not saving because PercentPhaseFieldOfView is %d (greater than 70).',info.PercentPhaseFieldOfView);
        listErrors = add_error(listErrors, saveName, pthRawScan, msg);
        continue;
    end
    % check for incorrect number of images within series
    if (length(listImgs) ~= info.ImagesInAcquisition) && ~(strcmp(saveName,'DCE'))
        % determine a new list of a subset of images in folder
        listImgsCells = sort({listImgs(:).name})';
        listImgs = [];
        
        % check last to first acquired series for complete scan series
        % handles multiple volumes and same volume repeated within series:
        % 'IM-1234-0001.dcm', 'IM-5678-0001.dcm', 'IM-1234-0001-0001.dcm', 'IM-1234-0001-0002.dcm'
        listUniqSeries = unique(cellfun(@(x) strcat(x(1:8),'*',x(13:end)), listImgsCells, 'UniformOutput', 0),'stable');
        for j=numel(listUniqSeries):-1:1
            uniqueSeries = listUniqSeries{j};
            listImgs = dir(fullfile(pthRawScan,uniqueSeries));
            % if images have a prefix but no suffix, remove images with suffix
            if length(uniqueSeries)==13; % i.e. 'IM-1234-*.dcm' for 'IM-1234-5678.dcm'
                listImgs(cellfun('length',{listImgs.name})>16)=[];
            end
            
            % check for the proper number of images
            if (length(listImgs) == info.ImagesInAcquisition)
                %warning
                switch j
                    case 1;    str = 'st';
                    case 2;    str = 'nd';
                    case 3;    str = 'rd';
                    otherwise; str = 'th';
                end
                msg = sprintf('Warning: Reading scan acquired %d%s of %d in folder, %s',j,str,numel(listUniqSeries),uniqueSeries);
                listErrors = add_error(listErrors, saveName, pthRawScan, msg);
                break
            else
                listImgs = []; % reset listImgs for next iteration
            end
        end
        % error
        if isempty(listImgs)
            msg = sprintf('Error: Cannot read dicoms. Incorrect number of images in all %d series acquired.',numel(listUniqSeries));
            listErrors = add_error(listErrors, saveName, pthRawScan, msg);
            continue
        end
    end
    
    %----------------------------------------------------------------------
    % PART 2: For the desired series convert images into a .mgz and .mat
    %read dicoms and create vol variable
    fprintf('Reading images: \n\t Folder: %s \n\t Path: %s\n\n', scanName, pthRawScan);
    pthOriginal = pwd;
    %----------------------------------------------------------------------
    %PART 2A: create volume for data type
    %OPTION A: read diffusion data
    if strcmp(saveName,'reducedFOV_RSI') || strcmp(saveName,'fullFOV_RSI')
        % move to folder with the dicom imgs
        cd(pthRawScan);
        [rsidat,M,~,bValsReadDcm,gwarpInfoRsi] = ReadDicomDiffusionData_new({listImgs(:).name});
        vol = ctx_mgh2ctx(rsidat,M); % volume output is transposed/sideways
        cd(pthOriginal)
        vol_rvs_imgs = vol.imgs(:,:,:,1);
        vol_fwd_imgs = squeeze(vol.imgs(:,:,:,2));
        vol.imgs = vol.imgs(:,:,:,2:end);   %remove b0 rvs image from RSI.mgz/.mat
    elseif strcmp(saveName,'Ax_DWI_ASSET')
        % first reverse b0
        try
            cd(fullfile(fileparts(pthRawScan), 'Ax_DWI_ASSET_rev'))
            listImgs_rev = listdir(fullfile(fileparts(pthRawScan), 'Ax_DWI_ASSET_rev'), 'files', 'struct');
            [dwidat,M] = ReadDicomDiffusionData_new({listImgs_rev(:).name});
            vol = ctx_mgh2ctx(dwidat,M); % volume output is transposed/sideways
            vol_rvs_imgs = flip(vol.imgs(:, :, :, 1), 1);
        catch
            vol_rvs_imgs = [];
        end
        % then forward volumes
        cd(pthRawScan);
        [dwidat,M,~,bValsReadDcm,gwarpInfoRsi] = ReadDicomDiffusionData_new({listImgs(:).name});
        vol = ctx_mgh2ctx(dwidat,M); % volume output is transposed/sideways
        vol.imgs = vol.imgs(:, :, :, [4 1 2 3]); % b0 first
        vol_fwd_imgs = squeeze(vol.imgs(:,:,:,1));
        %OPTION B: read phase encoding directions
    elseif strcmp(saveName,'DCE')
        [vol,listErrors] = Read_DCE(pthRawDate,listErrors);
        if isempty(vol.imgs)
            continue
        end
        %OPTION C: read T2 Data
    elseif strcmp(saveName,'T2_FS_PRE')
        try
            % move to folder with the dicom imgs
            cd(pthRawScan);
            [vol,gradwarpinfo] = Read_DICOM_2D_Directory({listImgs(:).name});
            cd(pthOriginal)
            % if multiple volumes in vol cell array, warning and documentation
            if length(vol) > 1
                msg = 'Warning: Multiple volumes after reading T2 dicoms, only saving the first volume';
                listErrors = add_error(listErrors, saveName, pthRawScan, msg);
            end
            vol = vol{1};
        catch
            msg = strcat('Error: Could not run Read_DICOM_2D_Directory.m to create vol of dicom images. ',...
                'Run code and error message will be displayed in the command window.');
            listErrors = add_error(listErrors, saveName, pthRawScan, msg);
            cd(pthOriginal)
            vol.imgs = [];
            continue
        end
        %OPTION D: read T1 Data, saveName = 'T1_FS_PRE', 'T1_noFS_PRE'
    else
        try
            % move to folder with the dicom imgs
            cd(pthRawScan);
            [vol,gradwarpinfo] = QD_Read_DICOM_3D_Directory({listImgs(:).name});
            cd(pthOriginal)
        catch
            msg = strcat('Error: Could not run QD_Read_DICOM_3D_Directory.m to create vol of dicom images. ',...
                'Run code and error message will be displayed in the command window.');
            listErrors = add_error(listErrors, saveName, pthRawScan, msg);
            cd(pthOriginal)
            vol.imgs = [];
            continue
        end
    end
    
    % if unable to read the volume, do not save and document error
    if isempty(vol.imgs)
        msg = 'Error: Unable to read dicom images to create volume. Run code and error message will be displayed in the command window.';
        listErrors = add_error(listErrors, saveName, pthRawScan, msg);
        continue
    end
    
    %----------------------------------------------------------------------
    %PART 2B: Save special cases
    %save bVals, b0 fwd, and b0 rvs
    if strcmp(saveName,'reducedFOV_RSI') || strcmp(saveName,'fullFOV_RSI')
        pthOutDistCorr = fullfile(pthOutSave,'dist_corr');
        fovType = saveName(1:(end-3));
        
        % save bVals and qmat
        [bVals,qmat] = get_bVals_qmat(bValsReadDcm,info);
        saveNamebVals = strcat(fovType,'bVals.mat');
        saveNameQmat = strcat(fovType,'qmat.mat');
        save(fullfile(pthOutSave,saveNamebVals),'bVals');
        fprintf('Saved: \n\t File: %s \n\t Path: %s\n\n', saveNamebVals, pthOutSave);
        save(fullfile(pthOutSave,saveNameQmat),'qmat');
        fprintf('Saved: \n\t File: %s \n\t Path: %s\n\n', saveNameQmat, pthOutSave);
        
        %save fwd as .mgz and .mat
        saveNameFwd = strcat(fovType,'b0_fwd');
        vol_fwd = vol;
        vol_fwd.imgs = vol_fwd_imgs;
        isSaved = save_vol2mgz (vol_fwd, saveNameFwd, pthOutDistCorr);
        if(~isSaved); listErrors = add_error(listErrors, [saveNameFwd,'.mgz'], pthRawScan); end
        
        vol_fwd.imgs = permute(vol_fwd.imgs,[2 1 3]);
        save(fullfile(pthOutDistCorr,[saveNameFwd '.mat']),'vol_fwd');
        fprintf('Saved: \n\t File: %s.mat \n\t Path: %s\n\n', saveNameFwd, pthOutDistCorr);
        
        %save rvs as .mgz and .mat
        saveNameRvs = strcat(fovType,'b0_rvs');
        vol_rvs = vol;
        vol_rvs.imgs = vol_rvs_imgs;
        nSlices = size(vol_rvs.imgs,3);
        %LR filp img across vertical axis, since ctx vol is transposed flip up/down
        if strcmp(info.InPlanePhaseEncodingDirection,'ROW')
            for idxSlice = 1:nSlices; vol_rvs.imgs(:,:,idxSlice) = flip(vol_rvs.imgs(:,:,idxSlice),1); end
            % if code crashes (run on before R2013b), use flipud(vol_rvs.imgs(:,:,idxSlice));
            %AP flip img across horizontal axis, since ctx vol is transposed flip left/right
        elseif strcmp(info.InPlanePhaseEncodingDirection,'COL')
            for idxSlice = 1:nSlices; vol_rvs.imgs(:,:,idxSlice) = flip(vol_rvs.imgs(:,:,idxSlice),2); end
            % if code crashes (run on before R2013b), use fliplr(vol_rvs.imgs(:,:,idxSlice));
        end
        
        %save as .mgz
        isSaved = save_vol2mgz (vol_rvs, saveNameRvs, pthOutDistCorr);
        if(~isSaved); listErrors = add_error(listErrors,[saveNameRvs,'.mgz'],pthRawScan); end
        
        %save as .mat, note b0rvs is not saved in RSI.mat or .mgz
        vol_rvs.imgs = permute(vol_rvs.imgs,[2 1 3]);
        save(fullfile(pthOutDistCorr,[saveNameRvs '.mat']),'vol_rvs');
        fprintf('Saved: \n\t File: %s.mat \n\t Path: %s\n\n', saveNameRvs, pthOutDistCorr);
        
    elseif strcmp(saveName,'Ax_DWI_ASSET')
        pthOutDistCorr = fullfile(pthOutSave,'dist_corr');
        fovType = strcat(saveName, '_');
        
        % save bVals and qmat (find how to process qmat)
        bVals = get_bVals_iSPY({listImgs(:).name});
        saveNamebVals = strcat(fovType,'bVals.mat');
        save(fullfile(pthOutSave,saveNamebVals),'bVals');
        fprintf('Saved: \n\t File: %s \n\t Path: %s\n\n', saveNamebVals, pthOutSave);
        
        %save fwd as .mgz and .mat
        saveNameFwd = strcat(fovType,'b0_fwd');
        vol_fwd = vol;
        vol_fwd.imgs = vol_fwd_imgs;
        isSaved = save_vol2mgz (vol_fwd, saveNameFwd, pthOutDistCorr);
        if(~isSaved); listErrors = add_error(listErrors, [saveNameFwd,'.mgz'], pthRawScan); end
        
        vol_fwd.imgs = permute(vol_fwd.imgs,[2 1 3]);
        save(fullfile(pthOutDistCorr,[saveNameFwd '.mat']),'vol_fwd');
        fprintf('Saved: \n\t File: %s.mat \n\t Path: %s\n\n', saveNameFwd, pthOutDistCorr);
        
        %save rvs as .mgz and .mat
        saveNameRvs = strcat(fovType,'b0_rvs');
        vol_rvs = vol;
        vol_rvs.imgs = vol_rvs_imgs;
        nSlices = size(vol_rvs.imgs,3);
        %LR filp img across vertical axis, since ctx vol is transposed flip up/down
        if strcmp(info.InPlanePhaseEncodingDirection,'ROW')
            for idxSlice = 1:nSlices; vol_rvs.imgs(:,:,idxSlice) = flip(vol_rvs.imgs(:,:,idxSlice),1); end
            % if code crashes (run on before R2013b), use flipud(vol_rvs.imgs(:,:,idxSlice));
            %AP flip img across horizontal axis, since ctx vol is transposed flip left/right
        elseif strcmp(info.InPlanePhaseEncodingDirection,'COL')
            for idxSlice = 1:nSlices; vol_rvs.imgs(:,:,idxSlice) = flip(vol_rvs.imgs(:,:,idxSlice),2); end
            % if code crashes (run on before R2013b), use fliplr(vol_rvs.imgs(:,:,idxSlice));
        end
        
        %save as .mgz
        isSaved = save_vol2mgz (vol_rvs, saveNameRvs, pthOutDistCorr);
        if(~isSaved); listErrors = add_error(listErrors,[saveNameRvs,'.mgz'],pthRawScan); end
        
        %save as .mat, note b0rvs is not saved in RSI.mat or .mgz
        vol_rvs.imgs = permute(vol_rvs.imgs,[2 1 3]);
        save(fullfile(pthOutDistCorr,[saveNameRvs '.mat']),'vol_rvs');
        fprintf('Saved: \n\t File: %s.mat \n\t Path: %s\n\n', saveNameRvs, pthOutDistCorr);
        
        %save DCE_Ph2_SUB.mgz
    elseif strcmp(saveName,'DCE') && exist(fullfile(pthOutDate,'T1_FS_PRE.mgz'))
        vol_T1_FS_Pre = QD_ctx_load_mgh(fullfile(pthOutDate,'T1_FS_PRE.mgz'));
        listErrors = save_DCE_Ph2_SUB(vol, vol_T1_FS_Pre, pthOutSave, listErrors);
    elseif strcmp(saveName,'T1_FS_PRE') && exist(fullfile(pthOutDate,'DCE.mgz'))
        vol_DCE = QD_ctx_load_mgh(fullfile(pthOutDate,'DCE.mgz'));
        listErrors = save_DCE_Ph2_SUB(vol_DCE, vol, pthOutSave, listErrors);
    end
    
    %----------------------------------------------------------------------
    % PART 2C: Save vol
    % save images as .mgz
    isSaved = save_vol2mgz (vol, saveName, pthOutSave);
    if(~isSaved); listErrors = add_error(listErrors, [saveName '.mgz'], pthRawScan); end
    
    % add orientation tag (AP or LR) to vol and save in .mat
    if strcmp(info.InPlanePhaseEncodingDirection,'ROW')
        vol.pedir = 'LR';
    elseif strcmp(info.InPlanePhaseEncodingDirection,'COL')
        vol.pedir = 'AP';
    else
        vol.pedir = 'not specified';
    end
    
    % save images as .mat
    vol.imgs = permute(vol.imgs,[2 1 3 4]);
    save(fullfile(pthOutSave,[saveName '.mat']),'vol');
    fprintf('Saved: \n\t File: %s.mat \n\t Path: %s\n\n', saveName, pthOutSave);
end

end

function [vol, listErrors] = Read_DCE(pthRawDate,listErrors)
vol.imgs = []; scanName = {};

% check if folders exist
for nPhase = 1:5
    % set folder name
    scanName{nPhase} = sprintf('Ph%d-AX_DYN_MP_326_2.4_62_NO_DELAY',nPhase);
    pthRawScan = fullfile(pthRawDate,scanName{nPhase});
    % to account for cases with series number Ana added
    if ~exist(pthRawScan)
        try
            listScan = [dir(sprintf('%s/Ph%d-AX_DYN_MP_326_2.4_62_NO_DELAY_*',pthRawDate,nPhase))];
            scanName{nPhase} = listScan(1).name;
            pthRawScan = fullfile(pthRawDate,scanName{nPhase});
        end
    end
    
    % check if folder exists to add to vol
    if ~exist(pthRawScan, 'dir')
        msg = sprintf('Error: Unable to create DCE because folder %s does not exist.',scanName{nPhase});
        warning(sprintf('%s\n\t Folder: %s\n\t Path: %s\n\n',msg,scanName{nPhase},pthRawScan));
        listErrors = add_error(listErrors, 'DCE', pthRawScan, msg);
        vol.imgs = [];
        return
    end
end

% check number of images for series number
listImgs = [dir(fullfile(pthRawDate,scanName{1},'*.dcm'))];
listImgsCells = sort({listImgs(:).name})';
% check last to first acquired series for complete scan series
% handles multiple volumes and same volume repeated within series:
% 'IM-1234-0001.dcm', 'IM-5678-0001.dcm', 'IM-1234-0001-0001.dcm', 'IM-1234-0001-0002.dcm'
listUniqSeries = unique(cellfun(@(x) strcat(x(1:8),'*',x(13:end)), listImgsCells, 'UniformOutput', 0),'stable');

for j=numel(listUniqSeries):-1:1
    uniqueSeries = listUniqSeries{j};
    
    for nPhase = 1:5
        % set folder name
        pthRawScan = fullfile(pthRawDate,scanName{nPhase});
        
        listImgs = [];
        listImgs = dir(fullfile(pthRawScan,uniqueSeries));
        
        % check number of images
        if nPhase == 1
            nImgs = length(listImgs);
        elseif length(listImgs) ~= nImgs
            vol.imgs = [];
            if j==1
                msg = sprintf('Error: Cannot read dicoms in all %d series acquired. Incorrect number of images in series %s.',numel(listUniqSeries),uniqueSeries);
                listErrors = add_error(listErrors, 'DCE', pthRawScan, msg);
                return
            else
                break
            end
        end
        
        % read and create volume structure
        try
            pthOriginal = pwd;
            cd(pthRawScan)
            [volPhase,gradwarpinfo] = QD_Read_DICOM_3D_Directory({listImgs(1:296).name});
            cd(pthOriginal)
        catch
            cd(pthOriginal)
            vol.imgs = [];
            if j==1
                msg = sprintf(strcat('Error: Cannot read dicoms in all %d series acquired. Could not run QD_Read_DICOM_3D_Directory.m for series %s to create vol of dicom images. ',...
                    'Run code and error message will be displayed in the command window.'),numel(listUniqSeries),uniqueSeries);
                listErrors = add_error(listErrors, 'DCE', pthRawScan, msg);
                return
            else
                break
            end
        end
        
        if nPhase ==1
            vol = volPhase;
        else
            vol.imgs(:,:,:,nPhase) = volPhase.imgs;
        end
        
        seriesNum = sprintf('%04d',str2num(uniqueSeries(4:7))+1); %MT Updated to have 4 chars with leading 0s
        uniqueSeries = strcat(uniqueSeries(1:3),seriesNum,uniqueSeries(8:end));
    end
    
    
    if ~isempty(vol.imgs) && (numel(listUniqSeries) ~= 1)
        switch j
            case 1;    str = 'st';
            case 2;    str = 'nd';
            case 3;    str = 'rd';
            otherwise; str = 'th';
        end
        
        %warning
        msg = sprintf('Warning: Reading scan acquired %d%s of %d in folder, %s',j,str,numel(listUniqSeries),uniqueSeries);
        listErrors = add_error(listErrors, 'DCE', pthRawScan, msg);
        return
    end
end

end

function listErrors = save_DCE_Ph2_SUB(vol_DCE, vol_T1_FS_Pre, pthOutSave, listErrors)
vol = vol_DCE;
vol.imgs = vol_DCE.imgs(:,:,:,2) - vol_T1_FS_Pre.imgs;

%save as .mgz
isSaved = save_vol2mgz (vol, 'DCE_Ph2_SUB', pthOutSave);
if(~isSaved); listErrors = add_error(listErrors,'DCE_Ph2_SUB.mgz',pthOutSave); end

%save as .mat
vol.imgs = permute(vol.imgs,[2 1 3 4]);
save(fullfile(pthOutSave,'DCE_Ph2_SUB.mat'),'vol');
fprintf('Saved: \n\t File: DCE_Ph2_SUB.mat \n\t Path: %s\n\n', pthOutSave);
end