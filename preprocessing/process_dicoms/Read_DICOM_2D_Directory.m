function [vols,gradwarpinfo] = Read_DICOM_2D_Directory(indir)
%

%warnstatussave = warning;
warning('off');
fcntr = 0;
if iscell(indir)
    file_list = indir;
else
    file_list = recursive_dir(indir);
end
totfiles = length(file_list);
for fnum = 1:totfiles

    fname = char(file_list(fnum));
    try
        tmp = dicominfo(fname);
        fcntr = fcntr+1;
        if fcntr > 1 & length(fieldnames(dcminfo(1))) ~= length(fieldnames(tmp)),
          fprintf('%s -- %s.m:     DICOM file header does not match other images in series %s\n',datestr(now),mfilename,fname);
          % copy the first entry
          e = dcminfo(1);
          % only replace what we can replace
          names1 = fieldnames(e);
          names2 = fieldnames(tmp);
          for key=1:length(names1),
            if sum(ismember(names1{key},names2)), % if we have that key in the structure
               e.(names1{key}) = tmp.(names1{key});
            end
          end;
          dcminfo(fcntr) = e;
        else
           dcminfo(fcntr) = tmp;
        end
    catch
        fprintf('%s -- %s.m:     Cannot read DICOM file %s\n',datestr(now),mfilename,fname);
    end
    if fcntr == 1
       StudyInstanceUID = dcminfo(fcntr).StudyInstanceUID;
       SeriesInstanceUID = dcminfo(fcntr).SeriesInstanceUID;
    end

    if exist('dcminfo','var')
       if ~strcmp(StudyInstanceUID,dcminfo(fcntr).StudyInstanceUID) | ~strcmp(SeriesInstanceUID,dcminfo(fcntr).SeriesInstanceUID)
           fprintf('%s -- %s.m:     Not the same StudyInstanceUID or SeriesInstanceUID\n',datestr(now),mfilename);
           return;
       end
    end
end

for fcntr = 1:length(dcminfo)
    dcminfo(fcntr).filename = dcminfo(fcntr).Filename;
    dcminfo(fcntr).filenamelen = length(dcminfo(fcntr).filename);
end

% change this when ctx_get_gradwarpinfo.m is updated in MMPS 
[gradwarpinfo,errmsg] = QD_get_gradwarpinfo(dcminfo);

tmp = [dcminfo.InstanceNumber];
[tmp,sortindx] = sort(tmp);
dcminfo = dcminfo(sortindx);

% Should add loop over multiple averages, pers, etc.
% Should make surre all images are "legal" MR type (not secondary capture - ets on Bydder) {dcminfo(:).ImageType}
[loc_list,A,B] = unique(cell2mat({dcminfo(:).SliceLocation})); % May not work for Siemens -- should use ImagePositionPatient?
nvols = length(B)/length(loc_list);
vols = cell(1,nvols);
for volindx = 1:nvols
  offset = (volindx-1)*length(loc_list); imlist = offset+[1:length(loc_list)]; % find(cell2mat({dcminfo(:).SliceLocation})==loc_list(locindx));
  ctx_vol = ctx_read_dicomfiles(dcminfo(imlist));
  vols{volindx} = ctx_vol;
end

% warning only returns gradwarp info for first volume.  needs to be fixed.
if isfield(gradwarpinfo, 'ambiguousgwtype')
    if gradwarpinfo.ambiguousgwtype == 1
        gradwarpinfo.gwtype = ctx_get_gwtype(vols{1},'isoctrflag', gradwarpinfo.isoctrflag,'gwtypelist',[2 3 7 8]);
    end
end

warning('on');

