function [vol,M] = read_dicom_4dvol(fnames,nrow,ncol,nslices,nreps,...
    convert_flag,tmpdir);
%function [vol,M] = read_dicom_4dvol(fnames,nrow,ncol,nslices,nreps,...
%  [convert_flag],[tmpdir]);
%
% Purpose: read 4d volume (3d x nreps) from dicoms
%
% Required Input:
%   fnames: cell array of dicom file names (must be pre-sorted)
%   nrow: number of rows
%   ncol: number of columns
%   nslices: number of slices
%   nreps: number of repetitions
%
% Optional Input:
%   convert_flag: [0|1|2] use mri_convert to convert dicoms
%     0: do not use mri_convert
%     1: use mri_convert for Seimens mosaic dicom format
%     2: use mri_convert for Philips dicom format
%     {default = 0}
%   tmpdir: temporary directory for loading mosaics
%     {default = '/tmp'}
%
% Created:  10/31/06 Don Hagler
% Last Mod: 12/13/10 Josh Kuperman
% Last Mod: 09/11/12 Don Hagler
% Last Mod: 10/02/23 Ana Rodriguez-Soto (added waitbar)

% based on read_dicomvol.m

%% todo: load mosaics directly, without using mri_convert?

M = eye(4);
vol=[];
nfiles = length(fnames);

if ~exist('convert_flag','var') || isempty(convert_flag), convert_flag=0; end;
if ~exist('tmpdir','var'), tmpdir=[]; end;

if convert_flag % use mri_convert
    fname = char(fnames{1});
    if ~exist(fname,'file')
        fprintf('%s: ERROR: file %s not found\n',mfilename,fname);
        return;
    end;
    if ~isempty(tmpdir)
        [tmp,tmp_fstem] = fileparts(tempname);
        fname_tmp = sprintf('%s/%s.mgh',tmpdir,tmp_fstem);
    else
        fname_tmp = sprintf('%s.mgh',tempname);
    end;
    if convert_flag == 1 % Siemens mosaic dicoms
        if nfiles ~= nreps
            fprintf('%s: ERROR: nfiles (%d) does not match nreps (%d)\n',...
                mfilename,nfiles,nreps);
            return;
        end;
        cmd = sprintf('mri_convert -it siemens_dicom %s %s',fname,fname_tmp);
        [status,result] = unix(cmd);
        if status
            fprintf('%s: ERROR: failed to convert mosaic dicoms:\n',mfilename);
            disp(result);
            return;
        end;
        [vol,M] = fs_load_mgh(fname_tmp);
    else % Philips dicoms
        cmd = sprintf('mri_convert %s %s',fname,fname_tmp);
        [status,result] = unix(cmd);
        if status
            fprintf('%s: ERROR: failed to convert Philips dicoms:\n',mfilename);
            disp(result);
            return;
        end;
        [vol,M] = fs_load_mgh(fname_tmp);
    end;
    delete(fname_tmp);
else
    if nfiles > nslices*nreps
        fprintf('%s: WARNING: nfiles (%d) is greater than nslices*nreps (%dx%d = %d) -- will use first %d reps only\n',...
            mfilename,nfiles,nslices,nreps,nslices*nreps,nreps);
    elseif nfiles ~= nslices*nreps
        fprintf('%s: ERROR: nfiles (%d) does not match nslices*nreps (%dx%d)\n',...
            mfilename,nfiles,nslices,nreps);
        return;
    end;
    M = read_dicom_M({fnames{1:nslices}});
    if isempty(M)
        vol=[];
        return;
    end;
    M = M_LPH_TO_RAS*M; % Convert from DICOM LPH to RAS coordinates
    dcminfo = dicominfo(char(fnames{1}));
    nrow = dcminfo.Rows;
    ncol = dcminfo.Columns;
    vol = zeros(nrow,ncol,nslices,nreps, 'single');
    n=1;

    f = waitbar(0, 'Loading data...');
    for r=1:nreps
        for z=1:nslices
            waitbar((((r-1)*nreps)+z)/(nreps*nslices), f, 'Loading data...')
            fname = char(fnames{n});
            if ~exist(fname,'file')
                fprintf('%s: ERROR: file %s not found\n',mfilename,fname);
                vol=[];
                return;
            end;
            im = dicomread(fname);
            if(isempty(im))
                fprintf('%s: ERROR: could not load pixel data from %s\n',mfilename,fname);
                vol=[];
                return;
            end
            vol(:,:,z,r) = im';
            n=n+1;
        end
    end;
    close(f)
end;

