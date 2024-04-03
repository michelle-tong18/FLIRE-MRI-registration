%% Subtract Noise-UW to see changes in noise correction 
% MT & LKF Last Update: 8/16/2020

% input: pthSave        path to folder with needed files
%        needs files    full/reducedFOV_RSI_uw.mgz
%                       full/reducedFOV_RSI_uw_nc.mgz
% output: creates files full/reducedFOV_RSI_nc_diff.mgz/mat
% returns: listErrors   struct of files that were not able to save

function listErrors = AC_diff_uw_nc(savePth)
    listErrors = [];
    global force_save

    listFiles = [dir(fullfile(savePth,'*RSI_uw.mgz'))];
    for ii = 1:length(listFiles)
        [~,fileName] = fileparts(listFiles(ii).name);
        saveName = [fileName,'_nc_diff'];
        
        %skip if output already exists
        if exist(fullfile(savePth,[saveName,'.mgz'])) && ...
                exist(fullfile(savePth,[saveName,'.mat'])) && ~force_save
            return
        end
        
        try
            vol_uw = QD_ctx_load_mgh(fullfile(savePth, [fileName,'.mgz']));
            vol_nc = QD_ctx_load_mgh(fullfile(savePth, [fileName,'_nc.mgz']));
        catch
            msg = 'Error: Unable to load uw and nc volumes for image subtraction';
            listErrors = add_error(listErrors, saveName, savePth, msg);
            continue
        end
        vol = vol_nc;
        vol.imgs = vol_nc.imgs - vol_uw.imgs;
        
        % save as .mgz file
        QD_ctx_save_mgh(vol, fullfile(savePth, [saveName '.mgz']));
        fprintf('Saved: \n\t File: %s.mgz\n\t Path: %s\n\n', saveName, savePth);
        %isSaved = save_vol2mgz(vol,saveName,savePth);
        
        % save as .mat file
        vol.imgs = permute(vol.imgs,[2 1 3 4]);
        save(fullfile(savePth,[saveName,'.mat']),'vol');
        fprintf('Saved: \n\t File: %s.mat\n\t Path: %s\n\n', saveName, savePth);
    end
end