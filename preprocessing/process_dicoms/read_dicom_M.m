function M = read_dicom_M(fnames);
%function M = read_dicom_M(fnames);
%
%

M = eye(4);
if ~iscell(fnames), fnames = {fnames}; end;
nfiles = length(fnames);

fname=char(fnames{1});
if ~exist(fname,'file')
  fprintf('%s: ERROR: file %s not found\n',mfilename,fname);
  M=[];
  return;
end;
dcminfo_1 = dicominfo(fname);

fname=char(fnames{end});
if ~exist(fname,'file')
  fprintf('%s: ERROR: file %s not found\n',mfilename,fname);
  M=[];
  return;
end;
dcminfo_end = dicominfo(fname);
M(1:3,1:2) = [dcminfo_1.PixelSpacing(1)*dcminfo_1.ImageOrientationPatient(1:3) dcminfo_1.PixelSpacing(2)*dcminfo_1.ImageOrientationPatient(4:6)];
M(1:3,3) = (dcminfo_end.ImagePositionPatient-dcminfo_1.ImagePositionPatient)/(nfiles-1);
M(1:3,4) = dcminfo_1.ImagePositionPatient-M(1:3,:)*[1 1 1 1]'; % Adjust for Matlab 1-based indexing


