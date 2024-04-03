function qmat = DTI_read_GE_tensor_dat(ndirs,nb0,fname)
%function qmat = DTI_read_GE_tensor_dat(ndirs,[nb0],[fname])
%
% Purpose: read diffusion vectors from GE tensor.dat file
%
% Input:
%  ndirs: number of diffusion directions
%  nb0 (optional): number of b0 images
%   default: 1
%  fname (optional): file name containing diffusion vectors
%   default: 'tensor.dat'
%
% Output:
%  qmat: matrix of diffusion direction vectors
%   size will be ndirs+nb0 x 3
%   (0 0 0 will be prepended as first nb0 rows)
%

qmat = [];

if nargin<1
  help(mfilename);
  return;
end;

min_ndirs = 1;
max_ndirs = 150;

if ~exist('nb0','var') nb0=[]; end;
if isempty(nb0) nb0=1; end;
if ~exist('fname','var'), fname=[]; end;
if isempty(fname) fname='tensor.dat'; end;

if ndirs>max_ndirs
  fprintf('%s: error: ndirs must be less than %d\n',mfilename,max_ndirs);
  return;
end;
if ndirs<min_ndirs
  fprintf('%s: error: ndirs must be greater than %d\n',mfilename,min_ndirs);
  return;
end;

fid=fopen(fname,'rt');
if fid==-1
  fprintf('%s: error reading %s\n',mfilename,fname);
  return;
end;

qmat = zeros(ndirs+nb0,3);
qmat(1:nb0,:) = zeros(nb0,3);
N=min_ndirs;
while (~feof(fid))
  temp=fgetl(fid);
  if strcmp(temp(1),'#'), continue; end;
  N = sscanf(temp,'%d');
  if length(N)==1 && N==ndirs, break; end;
end

for i=1:ndirs
  temp=fgetl(fid);
  A = sscanf(temp,'%f');
  if length(A)~=3
    fprintf('%s: error reading direction vector\n',mfilename);
    return;
  end;
  qmat(nb0+i,:)=A';
  
end;

fclose(fid);

