%% Apply Eddy Correction
% MT & LKF Last update: 8/23/2020

% input:  savePth       path to folder with needed files
%         needs files   full/reducedFOV_RSI_uw_nc.mgz
%                       full/reducedFOV_qmat.mat
%                       full/reducedFOV_bVals.mat
% output: creates files full/reducedFOV_RSI_uw_nc_eddy.mgz
% returns: listErrors   struct of files that were not able to save

function listErrors = AC_eddy(savePth,bVals,qmat)    
    listErrors = [];
    global force_save
    
    %list all RSI noise corrected .mgz files
    listFiles = [dir(fullfile(savePth,'*RSI_uw.mgz'))];    
    
    %list all RSI noise un-corrected .mgz files for ADC and kurtosis
    %mapping
%     listFiles = [dir(fullfile(savePth,'*RSI_uw.mgz'))];  
    
    for i = 1:length(listFiles)
        [~,fileName] = fileparts(listFiles(i).name); %fileName does not have extension (.mgz)
        [vol_in, M] = QD_load_mgh( fullfile(savePth, [fileName,'.mgz']) );
        saveName = [fileName,'_eddy'];
        [~,fovIdx] = regexp(fileName,'FOV_');
        fovType = fileName(1:fovIdx(1));
%         fovType = 'fullFOV_'; %ucsd data
        
        %skip if output already exists
%         if exist(fullfile(savePth,[saveName,'.mgz'])) && ...
%                 exist(fullfile(savePth,[saveName,'.mat'])) && ~force_save
%             continue
%         end
        
        try
            bvals = load( fullfile(savePth,[fovType,'bVals.mat']) ); %UCSD
            bvals = bvals.bVals;
            bvals = bvals(2:end);
%             bvals_qmat = load(fullfile(savePth,'BIDMC_bvals.mat')); %BIDMC data
%             bvals = bvals_qmat.bvals;
            ubvals = unique(bvals);
            
            qmat = load( fullfile(savePth,[fovType,'qmat.mat']) ); qmat = qmat.qmat;
%             qmat = bvals_qmat.qmat;
            qmat = qmat(2:end,:);
        catch
            msg = 'Error: Unable to apply eddy correction because unable to load bVals or qmat';
            listErrors = add_error(listErrors, saveName, savePth, msg);
            continue
        end
        
        % determine phase encoding direction
        vol_raw = load( fullfile(savePth, [saveName(1:end-8),'.mat']) );
%         %load ??FOV_RSI.mat UCSD
%         vol_raw = load( fullfile(savePth, [saveName,'.mat']) ); %BIDMC 
        vol_raw = vol_raw.vol;
        if isfield(vol_raw, 'pedir')
            pe_dim = vol_raw.pedir;
        else
            pe_dim = 'AP';
            msg = 'Warning: Cannot determine phase encoding direction for eddy current correction, assuming A/P';
            listErrors = add_error(listErrors, saveName, savePth, msg);
        end
        
        switch pe_dim
            case 'AP'; pe_num = 1;
            case 'LR'; pe_num = 2;
        end
               
        %% Begin Eddy Correction ----------------------------------------------------       
        [vol_eddy M_eddy]=QD_Eddy(vol_in,qmat,bvals,pe_num);

        %% End Eddy Correction ------------------------------------------------------        
        
        % save as a .mgz file
        vol = ctx_mgh2ctx(vol_eddy,M);
        ctx_save_mgh(vol, fullfile(savePth, [saveName '.mgz']));
        fprintf('Saved: \n\t File: %s.mgz\n\t Path: %s\n\n', saveName, savePth);
        
        % save as a .mat file
        vol.imgs = permute(vol.imgs,[2 1 3 4]);
        vol.MEddy = M_eddy;
        save(fullfile(savePth,[saveName,'.mat']),'vol');
        fprintf('Saved: \n\t File: %s.mat\n\t Path: %s\n\n', saveName, savePth);
%        imshow4D(vol.imgs)
        
    end
end