function [vol,Mvxl2ras,qmat,bvals,gwinfo, info] = ReadDicomDiffusionData_new(indir,ignoreGradInfo,varargin)
% ReadDicomDiffusionData  Read DICOM diffusion data


if ~exist('ignoreGradInfo','var')
    ignoreGradInfo = false;
end

if iscell(indir)
    file_list = indir;
else
    file_list = recursive_dir(indir);
end

totfiles = length(file_list);
fnames = {};
instancenumber = [];
qmat = [];
bvals = [];
iter = 1;
fprintf('%s -- %s.m:    Reading DICOM headers...',datestr(now),mfilename);
for i = 1:totfiles
    fname = char(file_list(i));
    try
        info = dicominfo(fname);
        info = fix_impax_dcm_tags(info);
        instancenumber(iter) = info.InstanceNumber;
        if ~ignoreGradInfo
            
          if strcmp(class(info.Private_0043_1039), 'char') == 1
             bval_num = str2num(info.Private_0043_1039);
             bvals = [bvals; bval_num(1)];
             qvecs = [str2num(info.Private_0019_10bb) str2num(info.Private_0019_10bc) str2num(info.Private_0019_10bd)];
             qmat = [qmat; qvecs];
          else
             bval = info.Private_0043_1039;
             bvals = [bvals; bval(1)];
             qvecs = [info.Private_0019_10bb info.Private_0019_10bc info.Private_0019_10bd];
             qmat = [qmat; qvecs];
          end

        end
        fnames{iter} = fname;
        iter = iter+1;
    catch ME
        warning(sprintf('%s - %s\n',ME.message,fname))
        continue
	% keyboard
    end
end
fprintf('\n');

if isempty(fnames);error('%s - Cant find any DICOM files in %s\n',mfilename,indir);end
fprintf('%s -- %s.m:    Loading DICOMs...\n',datestr(now),mfilename);

[tmp,sortindx] = sort(instancenumber);
if ~ignoreGradInfo
    qmat_sort = qmat(sortindx,:);
    bvals_sort = bvals(sortindx);
end
nr = info.Rows;
nc = info.Columns;
ns = info.Private_0021_104f;

if length(ns) > 1
    fprintf('Warning: Expecting single integer value for private tag 0021_104\n');
    fprintf('Warning: Using first index for number of slices = %d\n',ns(1));
    ns = double(ns(1));
end

nreps = info.ImagesInAcquisition/ns;

% TODO: Cross reference with QD_Read_DICOM_3D_Directory gwinfo calculation 
% make sure isocenter flag is computed properly

%gwinfo = ctx_get_gradwarpinfo(info);

[gwinfo,errmsg] = get_gradwarpinfo(info);
if isfield(gwinfo, 'ambiguousgwtype')
    if gwinfo.ambiguousgwtype == 1
        gwinfo.gwtype = ctx_get_gwtype(vol, 'isoctrflag', gwinfo.isoctrflag, 'gwtypelist', [2 3 7 8]);
    end
end

if ~ignoreGradInfo
    qmat = qmat_sort(1:ns:end,:); % single diffusion direction per volume is standard.
    bvals = bvals_sort(1:ns:end);   % single b-value per volume as well.
end

if isempty(varargin)
    [vol,Mvxl2ras] = read_dicom_4dvol(fnames(sortindx),nr,nc,ns,nreps); % based on read_dicomvol.m
else
    nr = varargin{1}(1);
    nc = varargin{1}(2);
    ns = varargin{1}(3);
    nreps = varargin{1}(4);
    
    [vol,Mvxl2ras] = read_dicom_4dvol(fnames(sortindx),nr,nc,ns,nreps); % based on read_dicomvol.m
end

end
