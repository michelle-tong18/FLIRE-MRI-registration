%% Subtract Eddy-Noise to see changes in vol Correction
% MT & LKF Last update: 8/16/2020

% input:  savePth       path to folder with needed files
%         needs files   full/reducedFOV_RSI_uw_nc.mgz
%                       full/reducedFOV_RSI_uw_nc_eddy.mgz
% output: creates files full/reducedFOV_RSI_eddy_diff.mgz/mat
% returns: listErrors   struct of files that were not able to save

function listErrors = AC_diff_nc_eddy(savePth)    
    listErrors = [];
    global force_save

    listFiles = [dir(fullfile(savePth,'*RSI_uw_nc.mgz'))];
    for ii = 1:length(listFiles)
        [~,fileName] = fileparts(listFiles(ii).name); %fileName does not have extension (.mgz)
        saveName = [fileName,'_eddy_diff'];
        
        %skip if output already exists
        if exist(fullfile(savePth,[saveName,'.mgz'])) && ...
                exist(fullfile(savePth,[saveName,'.mat'])) && ~force_save
            return
        end
        
        try
            vol_nc = ctx_load_mgh( fullfile(savePth, [fileName,'.mgz']) );
            vol_eddy = ctx_load_mgh( fullfile(savePth, [fileName,'_eddy.mgz']) );
        catch
            msg = 'Error: Unable to load nc and eddy corrected volumes for image substraction';
            listErrors = add_error(listErrors, saveName, savePth, msg);
            continue
        end
        vol = vol_eddy;
        vol.imgs = vol_eddy.imgs - vol_nc.imgs;
        
        % save as a .mgz file
        ctx_save_mgh(vol, fullfile(savePth, [saveName '.mgz']));
        fprintf('Saved: \n\t File: %s.mgz\n\t Path: %s\n\n', saveName, savePth);
        %isSaved = save_vol2mgz (vol, saveName, savePth);
        
        % save as a .mat file
        vol.imgs = permute(vol.imgs,[2 1 3 4]);
        save(fullfile(savePth,[saveName,'.mat']),'vol');
        fprintf('Saved: \n\t File: %s.mat\n\t Path: %s\n\n', saveName, savePth);
    end
end