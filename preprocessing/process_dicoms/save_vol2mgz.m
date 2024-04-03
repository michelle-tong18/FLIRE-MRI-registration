%% Script to save vol to mgz - can handle spaces in path and windows environment
% MT & LKF Last update: 8/23/2020
%
% output: saves vol as saveName.mgz in savePth

function isSaved = save_vol2mgz (vol, saveName, savePth)
    if (ismac || isunix) %if a mac or unix platform ...
        % use for files with a space in path, i.e. 'Google Drive', 'One Drive'
%        QD_ctx_save_mgh(vol, fullfile('~', [saveName '.mgz'])) % function contains a unix command - must be run in appropriate environment (Mac/Linux OS)
%        movefile(fullfile('~',[saveName '.mgz']), fullfile(savePth,[saveName '.mgz']))
        % use for files without a space in path
        ctx_save_mgh(vol, fullfile(savePth, [saveName '.mgz'])) % function contains a unix command - must be run in appropriate environment (Mac/Linux OS)
        fprintf('Saved: \n\t File: %s.mgz\n\t Path: %s\n\n', saveName, savePth);
        isSaved = 1;
    else
        msg = 'CANNOT SAVE - Operating system/platform not supported.';
        warning(sprintf('Error: %s\n\t File: %s.mgz\n\t Path: %s\n\n',msg,saveName,savePth));
        isSaved = 0;
    end
end