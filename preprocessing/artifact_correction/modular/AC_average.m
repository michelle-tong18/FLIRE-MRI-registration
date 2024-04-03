%% Average
% MT & LKF Last update: 8/12/2020

% input:  savePth       path to folder with needed files
%         needs files   full/reducedFOV_RSI_uw_nc_eddy.mgz
%                       full/reducedFOV_bVals.mat
% output: creates files full/reducedFOV_RSI_uw_nc_eddy_averaged.mgz
% returns: listErrors   struct of files that were not able to save

function listErrors = AC_average(savePth)    
    listErrors = [];
    global force_save
    
    %list all RSI noise corrected .mgz files
    listFiles = [dir(fullfile(savePth,'*RSI_uw_eddy.mgz'))];

    %list all RSI noise un-corrected .mgz, eddy current corrected, files for ADC and kurtosis
    %mapping 
%     listFiles = [dir(fullfile(savePth,'*RSI_uw_eddy.mgz'))];  
    
    for i = 1:length(listFiles)
        [~,fileName] = fileparts(listFiles(i).name); %fileName does not have extension (.mgz)
        [vol_in, M] = QD_load_mgh( fullfile(savePth, [fileName,'.mgz']) );
        saveName = [fileName,'_averaged'];
       [~,fovIdx] = regexp(fileName,'FOV_');
        fovType = fileName(1:fovIdx(1));
%         fovType = 'reducedFOV_'; %BIDMC data      
        
        %skip if output already exists
%        if exist(fullfile(savePth,[saveName,'.mgz'])) && ...
%                exist(fullfile(savePth,[saveName,'.mat'])) && force_save==0
%            continue
%        end
        
        try
%             bvals = load( fullfile(savePth,'BIDMC_bvals.mat') );
%             bvals = bvals.bvals;
            bvals = load(fullfile(savePth,[fovType 'bVals.mat']));
            bvals = bvals.bVals;
            bvals = bvals(2:end); % remove b=0 rvs
            ubvals = unique(bvals);
        catch
            msg = 'Error: Unable to apply averaging because unable to load bVals';
            listErrors = add_error(listErrors, saveName, savePth, msg);
            continue
        end
        
        %% Begin Averaging ----------------------------------------------------
        dims = size(vol_in);
        vol_avg = zeros(dims(1),dims(2),dims(3),length(ubvals));
        for b = 1:length(ubvals);
		    b_inds = find(bvals==ubvals(b));
		    vol_avg(:,:,:,b) = mean(vol_in(:,:,:,b_inds), 4 );
        end
        %% End Averaging ----------------------------------------------------
        
        % save as a .mgz file
        vol = ctx_mgh2ctx(vol_avg,M);
        ctx_save_mgh(vol, fullfile(savePth, [saveName '.mgz']));
        fprintf('Saved: \n\t File: %s.mgz\n\t Path: %s\n\n', saveName, savePth);
        
        % save as a .mat file
        vol.imgs = permute(vol.imgs,[2 1 3 4]);
        save(fullfile(savePth,[saveName,'.mat']),'vol');
        fprintf('Saved: \n\t File: %s.mat\n\t Path: %s\n\n', saveName, savePth);
        %imshow4D(vol.imgs);
        
        % % Resample to resDCE
%         if exist(fullfile(fileparts(savePth),'DCE.mgz'))
%             vol_DCE = ctx_load_mgh(fullfile(fileparts(savePth),'DCE.mgz'));
%             vol_resDCE = vol_resample(vol,vol_DCE,eye(4),1);
%             ctx_save_mgh(vol_resDCE, fullfile(savePth,[saveName '_resDCE.mgz']));
%             vol_resDCE.imgs = permute(vol_resDCE.imgs,[2 1 3]);
%             save(fullfile(savePth,[saveName '_resDCE.mat']),'vol_resDCE')
%             fprintf('\t File: %s_resDCE.mgz and .mat\n', saveName);
%         end

    end
end