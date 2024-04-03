function [bVals,qmat] = get_bVals_qmat(bValsIn,info,nb0vols)
tensordat_ID = ['tensor' mat2str(info.Private_0019_10b6) '.dat'];
% tensordat_ID = ['tensor97.dat'];
nvols = info.ImagesInAcquisition/info.Private_0021_104f(1);
if nargin < 3; nb0vols = 2; end
qmat = QD_read_GE_tensor_dat(nvols-nb0vols,nb0vols,tensordat_ID);
bmax = max(bValsIn);
bVals = bmax*sum(qmat.^2,2);
for i=1:length(bVals)
    if bVals(i) > 950
        bVals(i) = roundn(bVals(i), 2);
    elseif bVals(i) > 100
        bVals(i) = round(bVals(i),-1);
     elseif bVals(i) < 100
        bVals(i) = roundn(bVals(i), 0);
    end
end

end