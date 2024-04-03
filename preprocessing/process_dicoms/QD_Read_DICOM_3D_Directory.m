function [vol,gradwarpinfo] = QD_Read_DICOM_3D_Directory(indir)
%

if iscell(indir)
    file_list = indir;
else
    file_list = recursive_dir(indir);
end

if isempty(file_list)
    fprintf('%s -- %s.m:    ERROR: No files found in %s\n', datestr(now),mfilename,indir)
    return
end

totfiles = length(file_list);
fcntr = 0;
for fnum = 1:totfiles
  fname = char(file_list(fnum));
  try
      tmp = dicominfo(fname);
      fcntr = fcntr+1;
      dcminfo(fcntr) = tmp;
  catch
      fprintf(1,'Cannot read DICOM file %s\n',fname);
      %      return;
  end
  if fcntr == 1
     StudyInstanceUID = dcminfo(fcntr).StudyInstanceUID;
     SeriesInstanceUID = dcminfo(fcntr).SeriesInstanceUID;
  end
  if ~strcmp('3D',dcminfo(fcntr).MRAcquisitionType)
      fprintf(1,'MRAcquisitionType not 3D\n');
      %return;
  end
 if ~strcmp(StudyInstanceUID,dcminfo(fcntr).StudyInstanceUID) | ~strcmp(SeriesInstanceUID,dcminfo(fcntr).SeriesInstanceUID)
     fprintf(1,'Multiple StudyInstanceUIDs or SeriesInstanceUIDs in directory\n');
     return;
 end
end

for fcntr = 1:length(dcminfo)
  dcminfo(fcntr).filename = dcminfo(fcntr).Filename;
  dcminfo(fcntr).filenamelen = length(dcminfo(fcntr).filename);
end
[gradwarpinfo,errmsg] = get_gradwarpinfo(dcminfo);
        
[tmp,sortindx] = sort([dcminfo.InstanceNumber]);
dcminfo = dcminfo(sortindx);
fnamelist = {};
for i = 1:length(dcminfo);fnamelist{i} = dcminfo(i).Filename;end

[vol,M] = read_dicomvol(fnamelist);
vol = ctx_mgh2ctx(vol,M);

if isfield(gradwarpinfo, 'ambiguousgwtype')
    if gradwarpinfo.ambiguousgwtype == 1
        gradwarpinfo.gwtype = ctx_get_gwtype(vol,'isoctrflag', gradwarpinfo.isoctrflag,'gwtypelist',[2 3 7 8]);
    end
end
