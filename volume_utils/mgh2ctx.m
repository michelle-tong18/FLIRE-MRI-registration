function vol_ctx = mgh2ctx(vol,M)
% function vol_ctx = QD_mgh2ctx(vol,M)
% 
% QD_mgh2ctx.m - takes 1-based array of voxel intensities and 1-based M_vox2_ras matrix
%                and produces ctx 1-based volume structure
%
% Input:
%   vol: nd-array of 1-based voxel intensities.
%        Note, only first volume is used for geometry information
%     M: 1-based M_vox2ras matrix
%
% last mod: 1/7/2011 by N White

[width,height,depth] = size(vol) ;

xsize = norm(M(:,1));
ysize = norm(M(:,2));
zsize = norm(M(:,3));
x_r = M(1,1)/xsize;
x_a = M(2,1)/xsize;
x_s = M(3,1)/xsize;
y_r = M(1,2)/ysize;
y_a = M(2,2)/ysize;
y_s = M(3,2)/ysize;
z_r = M(1,3)/zsize;
z_a = M(2,3)/zsize;
z_s = M(3,3)/zsize;
ci = (width)/2 ; cj = (height)/2 ; ck = (depth)/2 ;
c_r = M(1,4) + (M(1,1)*ci + M(1,2)*cj + M(1,3)*ck) ;
c_a = M(2,4) + (M(2,1)*ci + M(2,2)*cj + M(2,3)*ck);
c_s = M(3,4) + (M(3,1)*ci + M(3,2)*cj + M(3,3)*ck);

vol_ctx.imgs = double(vol);
vol_ctx.Mvxl2lph = M_RAS_TO_LPH * M;
vol_ctx.dimc = size(vol,2);
vol_ctx.dimr = size(vol,1);
vol_ctx.dimd = size(vol,3);
vol_ctx.vx = xsize;
vol_ctx.vy = ysize;
vol_ctx.vz = zsize;
vol_ctx.lphcent=[-c_r;-c_a;c_s]; % 1-based
vol_ctx.DirCol=[-x_r;-x_a;x_s];
vol_ctx.DirRow=[-y_r;-y_a;y_s];
vol_ctx.DirDep=[-z_r;-z_a;z_s];

% save ctx using ultimately fs_save_mgh it assumes Mvxl2lph is 1-based