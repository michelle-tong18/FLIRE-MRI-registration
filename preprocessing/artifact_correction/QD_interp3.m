function vol_res = QD_interp3(vol,dvol1,dvol2,dvol3,method)

if ~exist('method','var') | isempty(method), method = 'linear'; end
[vdim1 vdim2 vdim3] = size(vol);
[indvol2,indvol1,indvol3] = meshgrid(1:vdim2,1:vdim1,1:vdim3);
vol_res = interp3(vol,indvol2+dvol2,indvol1+dvol1,indvol3+dvol3,method);
vol_res(find(isnan(vol_res)))=0;
