function QD_save_mgh(vol, fname, M, mr_parms)
%function QD_save_mgh(vol, fname, M, mr_parms);
%
% QD_save_mgh.m - saves matlab 1-based array (vol) and M_vox_2_ras matrix (M) as .mgh (.mgz) file
% 
% Input:
%    vol: 3d or 4d array of voxel intensities
%    M:   is the 4x4 vox2ras transform such that
%         y(i1,i2,i3), xyz = M*[i1 i2 i3 1] where the
%         indicies are 1-based. M is converted to
%         0-based before saving.
%
%    mr_parms: = [tr flipangle te ti fov]
%
% early mod:  12/13/10 by C Roddey
% last mod:  1/10/2011 by N White
%

if(nargin < 2 | nargin > 4)
  help(mfilename);
  return;
end

if ~exist('M','var') | isempty(M), M=eye(4); end;
if ~exist('mr_parms','var') | isempty(mr_parms), mr_parms = [0 0 0 0]; end
if(length(mr_parms) < 4)
  error(sprintf('mr_parms length = %d, must be 4 or 5\n', ...
	  length(mr_parms)));
end

if ((length(fname) >=4 && strcmpi(fname((length(fname)-3):length(fname)), '.MGZ')) | ...
		(length(fname) >=3 && strcmpi(fname((length(fname)-2):length(fname)), '.GZ')))
  zipflag = 1;
  orig_fname = fname;
  % rand('twister',sum(100*clock));
  rng('default');
  tmp_stem = num2str(round(rand(1)*10000000));
  if verLessThan('matlab', '7.1.1')
      [out_path,out_stem,out_ext,out_ver] = fileparts(fname);
  else
      [out_path,out_stem,out_ext] = fileparts(fname); % R2010b removed fourth output argument
  end

  if isempty(out_path), out_path = './'; end;
  fname = sprintf('%s/%s_%s.mgh',out_path,out_stem,tmp_stem);
else
  zipflag = 0;
end;

% These dont appear to be used %
MRI_UCHAR =  0 ;
MRI_INT =    1 ;
MRI_LONG =   2 ;
MRI_FLOAT =  3 ;
MRI_SHORT =  4 ;
MRI_BITMAP = 5 ;
MRI_TENSOR = 6 ;

fid = fopen(fname, 'wb', 'b') ;
if(fid == -1)
  error(sprintf('could not open %s for writing\n',fname));
end

[ndim1,ndim2,ndim3,frames] = size(vol) ;
fwrite(fid, 1, 'int') ;		% magic #
fwrite(fid, ndim1, 'int') ; 
fwrite(fid, ndim2, 'int') ; 
fwrite(fid, ndim3, 'int') ; 
fwrite(fid, frames, 'int') ;	% # of frames
if(ndims(vol) == 5)
  is_tensor = 1 ;
  fwrite(fid, MRI_TENSOR, 'int') ; % type = MRI_TENSOR
else
  is_tensor = 0 ;
  fwrite(fid, MRI_FLOAT, 'int') ;  % type = MRI_FLOAT
end

fwrite(fid, 1, 'int') ;          % dof (not used)
dof = fread(fid, 1, 'int') ; 

UNUSED_SPACE_SIZE= 256;
USED_SPACE_SIZE = (3*4+4*3*4);  % space for ras transform

M = M*[1 0 0 1; 0 1 0 1; 0 0 1 1; 0 0 0 1]; % Convert from 1-based to 0-based indexing

MdcD = M(1:3,1:3);
delta = sqrt(sum(MdcD.^2));

Mdc = MdcD./repmat(delta,[3 1]);
Pcrs_c = [ndim1/2 ndim2/2 ndim3/2 1]'; %'
Pxyz_c = M*Pcrs_c;
Pxyz_c = Pxyz_c(1:3);

fwrite(fid, 1,      'short') ;       % ras_good_flag = 1
fwrite(fid, delta,  'float32') ; 
fwrite(fid, Mdc,    'float32') ; 
fwrite(fid, Pxyz_c, 'float32') ; 

unused_space_size = UNUSED_SPACE_SIZE-2 ;
unused_space_size = unused_space_size - USED_SPACE_SIZE ;
fwrite(fid, zeros(unused_space_size,1), 'char') ;

tmp = find(isnan(vol));
if ~isempty(tmp)
  vol(tmp) = 0;
end;
fwrite(fid,vol,'float32');

fwrite(fid, mr_parms, 'float32') ; 
fclose(fid) ;

if zipflag
  cmd = sprintf('gzip %s ; mv %s.gz %s', fname, fname, orig_fname);
	[status,result] = unix(cmd);
  if status
    disp(result);
    error('cmd %s failed',cmd);
  end;
end	

return;

