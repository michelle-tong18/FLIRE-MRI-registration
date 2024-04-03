%% Apply Noise Correction
% MT & LKF Last update: 8/23/2020

% input:  savePth       path to folder with needed files
%         needs files   full/reducedFOV_RSI_uw.mgz
%                       full/reducedFOV_bVals.mat
% output: creates files full/reducedFOV_RSI_uw_nc.mgz
% returns: listErrors   struct of files that were not able to save

function listErrors = AC_noise(savePth)    
    listErrors = [];
    global force_save
    
    %list all RSI noise corrected .mgz files
    listFiles = [dir(fullfile(savePth,'*RSI_uw.mgz'))];   
    for i = 1:length(listFiles)
        [~,fileName] = fileparts(listFiles(i).name); %fileName does not have extension (.mgz)
        [dwi_vol, M] = QD_load_mgh( fullfile(savePth, [fileName,'.mgz']) );
%        dwi_vol = dwi_vol(:,:,:,2:end); %already removed b0rvs img
        saveName = [fileName,'_nc'];
        [~,fovIdx] = regexp(fileName,'FOV_'); 
        fovType = fileName(1:fovIdx(1));
        
        %skip if output already exists
        if exist(fullfile(savePth,[saveName,'.mgz'])) && ...
                exist(fullfile(savePth,[saveName,'.mat'])) && ~force_save 
            continue
        end
        
        try
            bvals = load( fullfile(savePth,[fovType,'bVals.mat']) );
            bvals = bvals.bVals;
            bvals = bvals(2:end);
            ubvals = unique(bvals);
        catch
            msg = 'Error: Unable to apply noise correction because unable to load bVals';
            listErrors = add_error(listErrors, saveName, savePth, msg);
            continue
        end
        
        %% Start noise correction --------------------------------------------------
		dims = size(dwi_vol);
        %create background mask - 1s where signal intensity is 0 for all bValues
		maskvol_zero = not(max(dwi_vol,[],4)>0); 
		%create edge mask - 1s along each edge of 3D vol cube
        maskvol_edge = false(size(maskvol_zero)); maskvol_edge([1 end],:,:) = true; maskvol_edge(:,[1 end],:) = true; maskvol_edge(:,:,[1 end]) = true;
		%computes the Euclidean distance transform of the binary image BW
        maskvol_nz = bwdist(maskvol_zero|maskvol_edge)>2; ivec_nz = find(maskvol_nz); 
		%reshape 4D DWI vol into a 2D vol [M x N] where M is for each volume and N is for each bVal
        tmp = reshape(dwi_vol,[prod(dims(1:3)) dims(4)]);
        %determine bin size for histogram, 10000 evenly spaced bins from 0 to the largest number in tmp(ivec_nz,:)
        hv = linspace(0,max(colvec(tmp(ivec_nz,:))),10000);
		hc = hist(colvec(tmp(ivec_nz,bvals==ubvals(end))),hv);
		hc_sm = smoothdata(hc,'gaussian',20);
		[mvs mis] = findpeaks(hc_sm,'MinPeakProminence',0.1*max(hc_sm));
            %specifies how much the peak must stand out to be identified
		noiseval = hv(mis(1)); ivec = find(hv<=2*noiseval);
		[muHat,sigmaHat,muCI,sigmaCI] = normfit(hv(ivec),0.05,[],hc(ivec));

		vol_nc = dwi_vol - muHat;
        %% End noise correction -----------------------------------------------------
        % note: non-negativity is imposed on the images later during fitting to avoid SI bias
        
        % save as a .mgz file      
        vol = ctx_mgh2ctx(vol_nc,M);
        ctx_save_mgh(vol, fullfile(savePth, [saveName '.mgz']));
        fprintf('Saved: \n\t File: %s.mgz\n\t Path: %s\n\n', saveName, savePth);
        
        % save as a .mat file
        vol.imgs = permute(vol.imgs,[2 1 3 4]);
        save(fullfile(savePth,[saveName,'.mat']),'vol');
        fprintf('Saved: \n\t File: %s.mat\n\t Path: %s\n\n', saveName, savePth);
%        imshow4D(vol.imgs)

    end
end