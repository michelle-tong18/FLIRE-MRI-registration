function bVals = get_bVals_iSPY(listImages)

for i = 1:length(listImages)
    try
        info = dicominfo(listImages{i});
        bVals(i) = info.DiffusionBValue;
    catch
        bVals(i) = 0;
    end
    bVals = sort(unique(bVals));
end