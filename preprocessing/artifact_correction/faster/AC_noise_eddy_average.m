%% Apply Noise & Eddy Correction
% MT & LKF Last update: 11/16/2020
%% UPDATED SOME LINES FOR mTE PROJECT, SEE COMMENTS

% input:  savePth       path to folder with needed files
%         needs files   full/reducedFOV_RSI_uw.mgz
%                       full/reducedFOV_bVals.mat
% output: creates files full/reducedFOV_RSI_uw_nc.mgz
%                       full/reducedFOV_RSI_uw_nc_eddy.mgz
%                       full/reducedFOV_RSI_uw_nc_eddy_averaged.mgz
% returns: listErrors   struct of files that were not able to save
% function: if all 3 output files do not exist, all 3 files are created

function listErrors = AC_noise_eddy_average(savePth)
    listErrors = [];
    global force_save
    
    %list all RSI distortion corrected .mgz files
    listFiles = [dir(fullfile(savePth,'*RSI_uw.mgz'))];   
    for i = 1:length(listFiles)
        [~,fileName] = fileparts(listFiles(i).name); %fileName does not have extension (.mgz)
        % EDITED FOR mTE PROJECT
%         [~,fovIdx] = regexp(fileName,'FOV_'); 
%         fovType = fileName(1:fovIdx(1));
        fovType = fileName(1:end-6);
        
        [vol_dwi, M] = ctx_load_mgh( fullfile(savePth, [fileName,'.mgz']) );
%        vol_dwi = vol_dwi(:,:,:,2:end); %already removed b0rvs img
        
        %skip if output already exists
        if exist(fullfile(savePth,[fileName,'_nc.mgz'])) && ...
                exist(fullfile(savePth,[fileName,'_nc_eddy.mgz'])) && ...
                exist(fullfile(savePth,[fileName,'_nc_eddy_averaged.mgz'])) && ...
                force_save==0
            continue
        end
        
        %load bVals and qmat
        try
            bvals = load( fullfile(savePth,[fovType,'bVals.mat']) );
            bvals = bvals.bVals;
            bvals = bvals(2:end);
            ubvals = unique(bvals);
            qmat = load( fullfile(savePth,[fovType,'qmat.mat']) );
            qmat = qmat.qmat;
            qmat = qmat(2:end,:);
        catch
            msg = 'Error: Unable to apply noise correction because unable to load bVals and/or qmat';
            listErrors = add_error(listErrors, [fileName, '_nc.mgz'], savePth, msg);
            continue
        end
        
        %% Start noise correction --------------------------------------------------
		saveName = [fileName,'_nc'];
        
        dims = size(vol_dwi);
        %create background mask - 1s where signal intensity is 0 for all bValues
		maskvol_zero = not(max(vol_dwi,[],4)>0); 
		%create edge mask - 1s along each edge of 3D vol cube
        maskvol_edge = false(size(maskvol_zero)); maskvol_edge([1 end],:,:) = true; maskvol_edge(:,[1 end],:) = true; maskvol_edge(:,:,[1 end]) = true;
		%computes the Euclidean distance transform of the binary image BW
        maskvol_nz = bwdist(maskvol_zero|maskvol_edge)>2; ivec_nz = find(maskvol_nz); 
		%reshape 4D DWI vol into a 2D vol [M x N] where M is for each volume and N is for each bVal
        tmp = reshape(vol_dwi,[prod(dims(1:3)) dims(4)]);
        %determine bin size for histogram, 10000 evenly spaced bins from 0 to the largest number in tmp(ivec_nz,:)
        hv = linspace(0,max(colvec(tmp(ivec_nz,:))),10000);
		hc = hist(colvec(tmp(ivec_nz,bvals==ubvals(end))),hv);
		hc_sm = smoothdata(hc,'gaussian',20);
		[mvs mis] = findpeaks(hc_sm,'MinPeakProminence',0.1*max(hc_sm));
            %specifies how much the peak must stand out to be identified
		noiseval = hv(mis(1)); ivec = find(hv<=2*noiseval);
		[muHat,sigmaHat,muCI,sigmaCI] = normfit(hv(ivec),0.05,[],hc(ivec));

		vol_nc = vol_dwi - muHat;
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
        
        %% Begin Eddy Correction ----------------------------------------------------
        saveName = [saveName,'_eddy'];
        clear vol
        
        % determine phase encoding direction
        vol_raw = load( fullfile(savePth, [saveName(1:end-11),'.mat']) ); %load ??FOV_RSI.mat
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
        
        [vol_eddy M_eddy]=QD_Eddy(vol_nc,qmat,bvals,pe_num);
        
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
%        imshow4D(vol.imgs);

        %% Begin Averaging ----------------------------------------------------
        saveName = [saveName,'_averaged'];
        clear vol
        
        dims = size(vol_eddy);
        vol_avg = zeros(dims(1),dims(2),dims(3),length(ubvals));
        for b = 1:length(ubvals);
		    b_inds = find(bvals==ubvals(b));
		    vol_avg(:,:,:,b) = mean(vol_eddy(:,:,:,b_inds), 4 );
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
%        imshow4D(vol.imgs);

    end
end