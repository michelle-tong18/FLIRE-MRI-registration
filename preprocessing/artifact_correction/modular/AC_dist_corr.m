%% Distortion Correction
% MT & LKF Last update: 11/16/2020
%% UPDATED SOME LINES FOR mTE PROJECT, SEE COMMENTS

% input:  pthOut        path to folder with full/reduced_RSI.mgz
%         needs files   full/reduced_b0_fwd.mgz in dist_corr folder
%                       full/reduced_b0_rvs.mgz in dist_corr folder
%                       full/reduced_RSI.mgz
% output: creates files full/reduced_volmask.mgz
%                       patient_date_volmask.jpg in proc/QA/mask folder
%                       unwarpB0 files in dist_corr folder
%                       full/reducedFOV_RSI_uw.mgz and .mat
% returns: listErrors   struct of files that were not able to save

function listErrors = AC_dist_corr(pthOut,kvec,lvec)
    listErrors = [];
    global force_save

    %list all RSI .mgz files
    listFiles = [dir(fullfile(pthOut,'*RSI.mgz'))];

    for i = 1:length(listFiles)
        [~,fileName] = fileparts(listFiles(i).name); %fileName does not have extension (.mgz)
        saveName = [fileName,'_uw'];
%         % EDITED FOR mTE PROJECT
%         [~,fovIdx] = regexp(fileName,'FOV_'); 
%         fovType = fileName(1:fovIdx(1));
        fovType = fileName(1:end-3);
        
        %skip if output already exists
        if exist(fullfile(pthOut,[saveName,'.mgz'])) && ...
                exist(fullfile(pthOut,[saveName,'.mat'])) && ~force_save
            continue
        end

        % determine phase encoding direction
        vol_raw = load( fullfile(pthOut, [fileName,'.mat']) );
        vol_raw = vol_raw.vol;
        if isfield(vol_raw, 'pedir') 
            pe_dim = vol_raw.pedir;
        else
            pe_dim = 'AP';
            msg = 'Warning: Cannot determine phase encoding direction for distortion correction, assuming A/P';
            listErrors = add_error(listErrors, saveName, pthOut, msg);
        end
        
        %% Part 1: make mask
        pthDistCorr = fullfile(pthOut,'dist_corr');        
        saveNameMask = strcat(fovType,'volmask');
        nameb0Fwd = strcat(fovType,'b0_fwd');
        nameb0Rvs = strcat(fovType,'b0_rvs');
        vol_fwd = ctx_load_mgh( fullfile(pthDistCorr,[nameb0Fwd,'.mgz']) );
        vol_rvs = ctx_load_mgh( fullfile(pthDistCorr,[nameb0Rvs,'.mgz']) );
        vol_f = vol_fwd.imgs;
        vol_r = vol_rvs.imgs;
        
        if any(strfind(pthOut,'reduced')); 
            percent = 0.85;
        else
            percent = 0.7;
        end
        
        mu = (vol_f+vol_r)/2;         % avg SI/voxel
        sig = mean(mu(find(mu>0)));   % avg SI
        volmask = zeros(size(mu));    % initialize mask
        volmask(find(mu>percent*sig)) = 1; % mask = 1 where SI is greater than 70% avg SI
        volmask = double(real(smooth3d(single(volmask),15,15,15))); % Smooth and dilate mask, moving avg filter
        thr = graythresh(volmask);    % determines global thr for mask
        volmask = 1.0*(volmask>thr);
        
        % to get correct idxs, volmask must the the same orientation as dist corr images
        if (strcmp(pe_dim,'LR'))
            % transpose mask when saving (but not idxs)
            volmask = permute(volmask, [2 1 3]);
            idxs = find(volmask);
            volmask = permute(volmask, [2 1 3]);
        else % for AP and unspecified phase encoding direction
            idxs = find(volmask);
        end
        
        % save mask as .mgz
        mask = vol_fwd;
        mask.imgs = volmask;
        ctx_save_mgh(mask, fullfile(pthOut, [saveNameMask '.mgz']));
        
        % save masks as .mat - dont need to save b/c saved in RSI_uw.mat struct
        volmask = permute(volmask,[2 1 3]);
        mask.imgs = volmask;
%        save(fullfile(pthOut,[saveNameMask '.mat']),'mask');
%        fprintf('Saved: \n\t File: %s.mat \n\t Path: %s\n\n', saveNameMask, pthOut);

        % save .jpg of slices of the mask volume
        % EDITED FOR mTE PROJECT
%         [pthmTE,FOV] = fileparts(pthOut);
%         [pthOutDate,~] = fileparts(pthmTE);
        [pthOutDate,~] = fileparts(pthOut);
        [pthPatient,date] = fileparts(pthOutDate);
        [pthProc,patient] = fileparts(pthPatient);
        pthMask = fullfile(pthProc,'QA','mask');
        saveNameMaskImg = [patient,'_',date,'_',saveNameMask];

        nSlices = size(volmask,3);
        if nSlices > 2
            idxHalve = round(nSlices/2);
            idxFourth = round(nSlices/4);
            fig = figure();
            subplot(1,3,1); imshow(volmask(:,:,idxFourth)); xlabel(sprintf('Slice %d',idxFourth)); axis('square'); set(gca,'XTick',[], 'YTick', [])
            subplot(1,3,2); imshow(volmask(:,:,idxHalve)); xlabel(sprintf('Slice %d',idxHalve)); axis('square'); set(gca,'XTick',[], 'YTick', [])
            subplot(1,3,3); imshow(volmask(:,:,idxFourth*3)); xlabel(sprintf('Slice %d',idxFourth*3)); axis('square'); set(gca,'XTick',[], 'YTick', [])
            subplot(1,3,2); title(saveNameMaskImg,'Interpreter','none')
            saveas(fig,fullfile(pthMask,[saveNameMaskImg,'.jpg']));
            if ishandle(fig); close(fig); end
        else
            %error unable to save .jpg because not enough slices for img
            msg = 'Error: Unable to create mask image because less than 3 slices in vol';
            listErrors = add_error(listErrors, [saveNameMaskImg,'.jpg'], pthMask, msg);            
        end

        %% Part 2: Apply distortion correction
        if (strcmp(pe_dim,'LR'))
            % transpose and save RSI.mgz, b0_fwd.mgz, b0_rvs.mgz
            nameb0FwdOrig = nameb0Fwd;
            nameb0RvsOrig = nameb0Rvs;
            fileNameOrig = fileName;
            saveNameOrig = saveName;
            
            nameb0Fwd = strcat(nameb0Fwd,'_transposed');
            nameb0Rvs = strcat(nameb0Rvs,'_transposed');
            fileName = strcat(fileName,'_transposed');
            saveName = strcat(saveName,'_transposed');
            
            vol_fwd.imgs = permute(vol_f,[2 1 3]);
            ctx_save_mgh(vol_fwd, fullfile(pthDistCorr, [nameb0Fwd '.mgz']));
            fprintf('Saved: \n\t File: %s.mgz\n\t Path: %s\n\n', nameb0Fwd, pthDistCorr);
            
            vol_rvs.imgs = permute(vol_r,[2 1 3]);
            ctx_save_mgh(vol_rvs, fullfile(pthDistCorr, [nameb0Rvs '.mgz']));
            fprintf('Saved: \n\t File: %s.mgz\n\t Path: %s\n\n', nameb0Rvs, pthDistCorr);
            
            % imgs have correct orientation vol_raw struct loaded as .mat 
            ctx_save_mgh(vol_raw, fullfile(pthOut, [fileName '.mgz']));
            fprintf('Saved: \n\t File: %s.mgz\n\t Path: %s\n\n', fileName, pthOut);
        end
        
        %% Start Distortion Correction ------------------------------------------------------
        % determine which kernel and lamba maximize correlation
        corrKL = zeros(length(kvec),length(lvec)); maxCorr = -10;
        for k = 1:length(kvec);
            kernelWidthMax = kvec(k); %Pixel width of Gaussian kernel for blurring; an odd number >= 3, or 0
            for l = 1:length(lvec);
                lambda2 = lvec(l);    %Coefficiennt of gradient in cost fn. > 0
                
                % apply distortion correction to b0f and b0r, saves distortion corrected imgs
                cmd = sprintf('/bin/bash -c ''preprocessing/artifact_correction/distortion_correction/unwarpB0Dev_v1/epic -f %s -r %s -od %s -kernelWidthMax %d -lambda2 %d -xs %d -ys %d -zs %d -resample -nchunksZ %u''', ...
                    fullfile(pthDistCorr,[nameb0Fwd,'.mgz']),fullfile(pthDistCorr,[nameb0Rvs,'.mgz']),...
                    pthDistCorr,kernelWidthMax,lambda2,2,2,2,1); %last inputs, (..., voxel_resampling x 3, nchunksZ);
                system(cmd);
        
                % load distortion corrected data
                volF_uw = ctx_load_mgh(fullfile(pthDistCorr,sprintf('%s.mgz',nameb0FwdOrig)));
                volR_uw = ctx_load_mgh(fullfile(pthDistCorr,sprintf('%s.mgz',nameb0RvsOrig)));

                corrKL(k,l) = corr(volF_uw.imgs(idxs),volR_uw.imgs(idxs))
                fprintf('%s -- %s.m:\n    DISCO Distortion correlation - %g.\n',pthDistCorr,fileName,corrKL(k,l));

                % find K and L with the most correlation
                if (corrKL(k,l) > maxCorr)
                    kernelWidthMax_MaxCorr = kernelWidthMax;
                    lambda2_MaxCorr = lambda2;
                    maxCorr = corrKL(k,l);
                    if (maxCorr > 0.85)
                        break
                    end
                end
            end
        end

        % apply distortion correction using max L and max K
        cmd = sprintf('/bin/bash -c ''preprocessing/artifact_correction/distortion_correction/unwarpB0Dev_v1/epic -f %s -r %s -od %s -kernelWidthMax %d -lambda2 %d -xs %d -ys %d -zs %d -resample -nchunksZ %u''', ...
            fullfile(pthDistCorr,[nameb0Fwd,'.mgz']),fullfile(pthDistCorr,[nameb0Rvs,'.mgz']),...
            pthDistCorr,kernelWidthMax_MaxCorr,lambda2_MaxCorr,2,2,2,1); %last inputs, (..., voxel_resampling x 3, nchunksZ);
        system(cmd);
        saveNameD = sprintf('d_L2_%d_K%d.mgz',lambda2_MaxCorr,kernelWidthMax_MaxCorr);
        movefile(fullfile(pthDistCorr,'d.mgz'),fullfile(pthDistCorr,saveNameD));    
        
        % Apply distortion correction to complete RSI data
        cmd = sprintf('/bin/bash -c ''preprocessing/artifact_correction/distortion_correction/applyUnwarpB0_mod/applyUnwarpB0 -f %s -r %s -d %s -fo %s -ro %s -od %s''', ...
            fullfile(pthOut,[fileName,'.mgz']),fullfile(pthDistCorr,[nameb0Rvs,'.mgz']),fullfile(pthDistCorr,saveNameD),...
            [saveName,'.mgz'],[saveName,'_rB0uw.mgz'],pthOut); 
        system(cmd);
        
        %% End Distortion Correction ------------------------------------------------------
        % rename files with FOV prefix
        if (~isempty(fovType)) || (strcmp(pe_dim,'LR'))
            fprintf('Renaming distortion corrected .mgz files and changing orientation if necessary. \n\t Files: *.mgz\n\t Path: %s\n\n', pthDistCorr);
            listDCFiles = dir(fullfile(pthDistCorr,'*.mgz'));
            listDCFilesNew = {listDCFiles(:).name};
            listDCFilesNew(ismember(listDCFilesNew, {'.', '..','desktop.ini'})) = [];
            listDCFilesNew(~cellfun('isempty', strfind(listDCFilesNew,'_b0_fwd') )) = [];
            listDCFilesNew(~cellfun('isempty', strfind(listDCFilesNew,'_b0_rvs'))) = [];
            for ii = 1:length(listDCFilesNew)
                fileNameIn = listDCFilesNew{ii};
                
                if (strcmp(pe_dim,'LR'))
                    [vol_tmp M_tmp] = QD_load_mgh( fullfile(pthDistCorr,fileNameIn) );
                    vol_tmp = permute(vol_tmp,[2 1 3 4]);
                    QD_save_mgh(vol_tmp,fullfile(pthDistCorr,fileNameIn),M_tmp);
                    fprintf('\t corrected orientation of --- %s\n', fileNameIn);
                end
                
                % rename as with FOV and transposed
                fileNameOut = strcat(fovType,fileNameIn);
                movefile(fullfile(pthDistCorr,fileNameIn),fullfile(pthDistCorr,fileNameOut));
            end
        end
        
        % Load distortion corrected vol
        try
            vol = ctx_load_mgh(fullfile(pthOut,[saveName, '.mgz']));
            vol_uw_r = ctx_load_mgh(fullfile(pthOut,[saveName,'_rB0uw.mgz']));
        catch
            msg = 'Error: Epic distortion function was unable to create DWI vol';
            listErrors = add_error(listErrors, saveName, pthOut, msg);
            continue
        end
        
        % save correct orientation of files in dist_corr and remove transposed files
        if (strcmp(pe_dim,'LR'))
            fprintf('\t Path: %s\n', pthOut);
            listDCNew = {[saveName,'.mgz'],[saveName,'_rB0uw.mgz']};
            listDCOut = {[saveNameOrig,'.mgz'],[saveNameOrig,'_rB0uw.mgz']};
            for jj = 1:length(listDCNew)
                fileNameIn = listDCNew{jj};
                fileNameOut = listDCOut{jj};
                
                [vol_tmp M_tmp] = QD_load_mgh( fullfile(pthOut,fileNameIn) );
                vol_tmp = permute(vol_tmp,[2 1 3 4]);
                QD_save_mgh(vol_tmp,fullfile(pthOut,fileNameIn),M_tmp);
                fprintf('\t corrected orientation of --- %s\n', fileNameIn);
                movefile(fullfile(pthOut,fileNameIn),fullfile(pthOut,fileNameOut));
            end
           
            % remove transposed imgs created for dist corr
            delete(fullfile(pthDistCorr,[nameb0Fwd,'.mgz']));
            delete(fullfile(pthDistCorr,[nameb0Rvs,'.mgz']));
            delete(fullfile(pthOut,[fileName,'.mgz']));
            
            %update for save .mat
            saveName = saveNameOrig;
        else
            % update to save as a .mat
            vol.imgs = permute (vol.imgs, [2 1 3 4]);
            vol_uw_r.imgs = permute(vol_uw_r.imgs, [2 1 3 4]);
        end
        
        % save as a .mat (note 2 lines of code above)
        vol.mask = volmask;
        vol.maskThreshold = sprintf('%.2f * average signal intensity',percent);
        vol.kernelWidthMax = kernelWidthMax_MaxCorr;
        vol.lambda2 = lambda2_MaxCorr;
        vol.maxCorr = maxCorr;
        vol.kvec = kvec;
        vol.L2vec = lvec;
        vol.corrKL = corrKL;
        save(fullfile(pthOut,[saveName, '.mat']),'vol');  % save in postprocessing folder
        fprintf('Saved: \n\t File: %s.mat\n\t Path: %s\n\n', saveName, pthOut);
%            imshow4D(vol.imgs) %to check imgs

        % save rvs as .mat
        vol_uw = vol;
        clear vol
        vol = vol_uw_r;
        save(fullfile(pthOut,[saveName,'_rB0uw.mat']),'vol');  % save in postprocessing folder
        fprintf('Saved: \n\t File: %s_rB0uw.mat\n\t Path: %s\n\n', saveName, pthOut);
   
        %% save QA for distortion correction
        save_QA_dist_corr(vol_uw, vol_uw_r,pthOutDate,fovType);
    end
end

function save_QA_dist_corr(vol_uw, vol_uw_r,pthOutDate,fovType)
        
        % save .jpg of slices of the mask volume
        % EDITED for mTE PROJECT
        %[pthmTE,FOV] = fileparts(pthOut);
        %[pthOutDate,~] = fileparts(pthmTE);
        [pthPatient,date] = fileparts(pthOutDate);
        [pthProc,patient] = fileparts(pthPatient);
        pthQADC = fullfile(pthProc,'QA','dist_corr');
        saveNameQADCImg = [patient,'_',date,'_',fovType,'dist_corr'];
        imgTitle = sprintf('%s\nL2: %d, K: %d', saveNameQADCImg, vol_uw.lambda2, vol_uw.kernelWidthMax);
        
        volf = vol_uw.imgs(:,:,:,1);
        volr = vol_uw_r.imgs;
        
        fig = figure();
        nSlices = size(volf,3);
        idxHalve = round(nSlices/2);
        idxFourth = round(nSlices/4);
%         SI = colvec([...
%             volf(:,:,idxFourth)...
%             volf(:,:,idxHalve)...
%             volf(:,:,idxFourth*3)...
%             volr(:,:,idxFourth)...
%             volr(:,:,idxHalve)...
%             volr(:,:,idxFourth*3)]);
%         clim = [min(SI) max(SI)];
        
%         subplot(2,3,1); imagesc(volf(:,:,idxFourth)); xlabel(sprintf('Slice %d',idxFourth)); axis('square'); set(gca,'XTick',[], 'YTick', [])
%         subplot(2,3,2); imagesc(volf(:,:,idxHalve)); xlabel(sprintf('Slice %d',idxHalve)); axis('square'); set(gca,'XTick',[], 'YTick', [])
%         subplot(2,3,3); imagesc(volf(:,:,idxFourth*3)); xlabel(sprintf('Slice %d',idxFourth*3)); axis('square'); set(gca,'XTick',[], 'YTick', [])
%         subplot(2,3,4); imagesc(volr(:,:,idxFourth)); xlabel(sprintf('Slice %d',idxFourth)); axis('square'); set(gca,'XTick',[], 'YTick', [])
%         subplot(2,3,5); imagesc(volr(:,:,idxHalve)); xlabel(sprintf('Slice %d',idxHalve)); axis('square'); set(gca,'XTick',[], 'YTick', [])
%         subplot(2,3,6); imagesc(volr(:,:,idxFourth*3)); xlabel(sprintf('Slice %d',idxFourth*3)); axis('square'); set(gca,'XTick',[], 'YTick', [])
% 
%         sgtitle(imgTitle,'Interpreter','none')
%         subplot(2,3,2); title([fovType, 'RSI_uw_fB0'],'Interpreter','none');
%         subplot(2,3,5); title([fovType, 'RSI_uw_rB0'],'Interpreter','none');
 %       colorbar('Position',[0.92 0.12 0.025 0.75])
        
        saveas(fig,fullfile(pthQADC,[saveNameQADCImg,'.jpg']));
        fprintf('Saved: \n\t File: %s.jpg\n\t Path: %s\n\n', saveNameQADCImg, pthQADC);
   
        if ishandle(fig); close(fig); end
end
