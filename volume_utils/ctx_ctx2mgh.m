function [vol,M] = ctx_ctx2mgh(vol_ctx)
% function [vol M] = QD_ctx2mgh(vol_ctx)
% 
% QD_ctx2mgh.m - converts 1-based ctx vol structure (double) to 1-based (single) nd-array of voxel intensites
%                and 1-based M_vox2ras matrix suitable for saving mgh file.

vol = single(vol_ctx.imgs);
M = M_LPH_TO_RAS * vol_ctx.Mvxl2lph;
