%% Sort errors struct
% MT & LKF Last update: 8/23/2020

function listErrors = sort_errorList(listErrorsIn)
    % sort
    errorTable = struct2table(listErrorsIn);
    errorTable = sortrows(errorTable,{'path','file','date'},{'ascend','ascend','ascend'});
    listErrors = table2struct(errorTable);
    
    % match orientation of input, listErrors may be 1xM or Mx1
    if size(listErrorsIn,1) == 1
        listErrors = listErrors'; %transposes to 1xM
    end
end