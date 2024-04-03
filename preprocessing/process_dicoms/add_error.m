%% Add errors to struct
% MT & LKF Last update: 8/23/2020

function listErrors = add_error(listErrors, saveName, pthToError, msg)
    if nargin < 4
        msg = 'Error:';
    end
    warning(sprintf('%s\n\t File: %s\n\t Path: %s\n\n',msg,saveName,pthToError));
    
    idx = length(listErrors) + 1;
    listErrors(idx).date = datestr(now, 'yyyymmdd_HH:MM');
    listErrors(idx).file = saveName;
    listErrors(idx).path = pthToError;
    listErrors(idx).reason = msg;
end