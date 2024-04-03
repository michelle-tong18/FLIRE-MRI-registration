%% Normalize
% MT & LKF Last update: 8/18/2020

% input:  savePth       path to folder with needed files
%         needs files   full/reducedFOV_RSI_uw_nc_eddy_averaged.mgz
% output: creates files full/reducedFOV_RSI_uw_nc_eddy_averaged_norm.mgz
% returns: listErrors   struct of files that were not able to save

function listErrors = AC_norm(savePth)    
    listErrors = [];
    global force_save
    
    %list all RSI noise corrected .mgz files
    listFiles = [dir(fullfile(savePth,'*RSI_uw_nc_eddy_averaged.mgz'))];
    
    for i = 1:length(listFiles)
        [~,fileName] = fileparts(listFiles(i).name); %fileName does not have extension (.mgz)
        [vol_in, M] = QD_load_mgh( fullfile(savePth, [fileName,'.mgz']) );
        saveName = [fileName,'_norm'];
        
        %skip if output already exists
        if exist(fullfile(savePth,[saveName,'.mgz'])) && ...
                exist(fullfile(savePth,[saveName,'.mat'])) && ~force_save
            continue
        end
        
        %% Begin Normalization ----------------------------------------------
        vol_norm = MF_normalize(vol_in);
        %% End Normalization ------------------------------------------------
        
        % save as a .mgz file
        vol = ctx_mgh2ctx(vol_norm,M);
        ctx_save_mgh(vol, fullfile(savePth, [saveName '.mgz']));
        fprintf('Saved: \n\t File: %s.mgz\n\t Path: %s\n\n', saveName, savePth);
        
        % save as a .mat file
        vol.imgs = permute(vol.imgs,[2 1 3 4]);
        save(fullfile(savePth,[saveName,'.mat']),'vol');
        fprintf('Saved: \n\t File: %s.mat\n\t Path: %s\n\n', saveName, savePth);
%        imshow4D(vol.imgs);
        
    end
end