function [gradwarpinfo,errmsg] = get_gradwarpinfo(dcminfo)
%

errmsg='';

Manufacturer = dcminfo(1).Manufacturer;
FieldStrength = dcminfo(1).MagneticFieldStrength;
if isfield(dcminfo(1),'ManufacturersModelName')
    ManufacturersModelName = dcminfo(1).ManufacturersModelName;
else
    ManufacturersModelName = dcminfo(1).ManufacturerModelName;
end
switch lower(Manufacturer)
  case 'siemens'
    gradwarpinfo.unwarpflag = 0; % full 3D unwarping
    %gradwarpinfo.isoctrflag = 0; % assume not using isocenter scanning (NB: should check DICOM header info on new systems)
    SiemensCsaParse_ReadDicomTag(['Private_0029_1010']);
    SiemensCsaParse_ReadDicomTag(['Private_0029_1020']);
    gradwarpinfo.isoctrflag = dcminfo.csa.Isocentered;
    
    switch lower(ManufacturersModelName)
      case {'sonata','trio','sonatavision'}
        gradwarpinfo.gwtype = 0;
      case {'allegra'}
        gradwarpinfo.gwtype = 1;
      case {'avanto','triotim'}
        gradwarpinfo.gwtype = 4;
      case {'espree','axxess'}
        gradwarpinfo.gwtype = 5;
      case {'symphony','symphonyvision','symphonytim'} % Quantum
        gradwarpinfo.gwtype = 6;
      case {'skyra'} % skyra
        gradwarpinfo.gwtype = 11;
      case {'connectome'} %% todo: check this
        gradwarpinfo.gwtype = 12;
          %gradwarpinfo.isoctrflag = 1;
      case {'prisma','prisma_fit'} % PRISMA
        gradwarpinfo.gwtype = 13;
      otherwise
        errmsg=sprintf('%s: Unknown gradient model %s %s\n',mfilename,Manufacturer,ManufacturersModelName);
        return;
    end
  case 'ge medical systems'
    gradwarpinfo.unwarpflag = 1; % through-plane only (should check NOGRAD CV)
    gradwarpinfo.isoctrflag = 1; % assume using isocenter scanning (should check DICOM tag?)
    if ismember(deblank(lower(ManufacturersModelName)),{'discovery mr450','discovery mr750'})
        gradwarpinfo.gwtype = 9;
        return;
    end
    if ismember(deblank(lower(ManufacturersModelName)),{'discovery mr750w','signa pet/mr','signa architect'})
        gradwarpinfo.gwtype = 10;
        return;
    end
    if ismember(deblank(lower(ManufacturersModelName)),{'signa creator'}) && FieldStrength == 3
        gradwarpinfo.gwtype = 9;
        return;
    end
    if ismember(deblank(lower(ManufacturersModelName)),{'signa premier'})
        gradwarpinfo.gwtype = 14;
        return;
    end
    if isfield(dcminfo(1),'Private_0043_106f')
        tmp = dcminfo(1).Private_0043_106f;
        if length(tmp)>=2 & tmp(1)>=48 & tmp(1)<=57 & tmp(2)==92 % Does it look like ASCII nums delimited by \?
            tmp = tmp(1:2:end)-48;
        end
        if length(tmp) >=4
            key = tmp(4);
        else
            key = 0;
        end
        if key == 1
            gradwarpinfo.gwtype = 7; % Whole
        elseif key == 2;
            gradwarpinfo.gwtype = 8; % Zoom
        elseif key == 0;
            %        gradwarpinfo.gwtype = 2; % BRM mode??? NB: Needs to be checked!
            gradwarpinfo.ambiguousgwtype = 1; % Ambiguous
            errmsg=sprintf('%s: Ambiguous gradient info for %s %s\n',mfilename,Manufacturer,ManufacturersModelName);
            return;
        else
            gradwarpinfo.ambiguousgwtype = 1; % Ambiguous
            errmsg=sprintf('%s: TwinSpeed mode = %d unknown for %s %s\n',mfilename,key,Manufacturer,ManufacturersModelName);
            return;
        end
        return
    end
    switch lower(ManufacturersModelName)
      case {'brm'} % Check on actual name
        gradwarpinfo.gwtype = 2;
      case {'crm'} % Check on actual name
        gradwarpinfo.gwtype = 3;
      case {'signa excite','signa hdx'}
        gradwarpinfo.ambiguousgwtype = 1; % Ambiguous
        errmsg=sprintf('%s: Missing DICOM tag Private_0043_106f for %s %s\n',mfilename,Manufacturer,ManufacturersModelName);
        return;
      case 'genesis_signa'
        %        gradwarpinfo.gwtype = 2; % BRM mode??? NB: Needs to be checked!
        %        GET INFO FROM TEXTFILE FOR SYSTEM           XXX
        gradwarpinfo.ambiguousgwtype = 1; % Ambiguous
        errmsg=sprintf('%s: Ambiguous gradient info for %s %s\n',mfilename,Manufacturer,ManufacturersModelName);
        return;
      otherwise
        gradwarpinfo.ambiguousgwtype = 1; % Ambiguous
        errmsg=sprintf('%s: Unknown gradient model %s %s\n',mfilename,Manufacturer,ManufacturersModelName);
        return;
    end
  case 'philips medical systems'
    gradwarpinfo = []; % Default to no gradwarp for Philips scanners
  otherwise
    gradwarpinfo = []; % Default to no gradwarp for unknown scanner mfgrs
end
function SiemensCsaParse_ReadDicomTag(strTag)
    currdx=0;
    
    if ~strcmp(char(private_read(4)),'SV10') || ~all(private_read(4)==[4 3 2 1])
        error('Unsupported CSA block format');
    end
    
    % This parsing code is translated from gdcm (http://gdcm.sf.net/)
    numElements = double(private_readuint32(1));
    
    % Sanity check
    if private_readuint32(1)~=77
        error('Unsupported CSA block format');
    end
    
    for tagdx=1:numElements
        tagName = private_readstring(64);
        
        % Fix up tagName
        tagName(tagName == '-') = [];
        
        vm = private_readuint32(1);
        vr = private_readstring(4);
        syngodt = private_readuint32(1);
        nitems = double(private_readuint32(1));
        
        checkbit = private_readuint32(1);
        
        if checkbit ~= 77 && checkbit ~= 205
            error('Unsupported CSA block format');
        end
        
        data = {};
        for itemdx=1:nitems
            header = double(private_readuint32(4));
            
            if (header(3) ~= 77 && header(3) ~= 205) || ...
                    (header(1) ~= header(2)) || ...
                    (header(1) ~= header(4))
                error('Unsupported CSA block format');
            end
            
            data{itemdx} = private_readstring(header(1));
            
            % Dump junk up to DWORD boundary
            private_read(mod(mod(4-header(1),4),4));
        end
        
        % Store this in the csa structure
        switch vr
            case {'CS', 'LO', 'LT', 'SH', 'SS', 'UI', 'UT', 'UN'} % Strings and unknown byte string
                if numel(data) < vm
                    % Pad if necessary. Siemens CSA format omits null strings.
                    data{vm} = '';
                end
                
                if vm == 1
                    dcminfo.csa.(tagName) = data{1};
                else
                    dcminfo.csa.(tagName) = data(1:vm);
                end
            case {'DS', 'FD', 'FL', 'IS', 'SL', 'ST', 'UL', 'US'} % Numbers
                dataNumeric = arrayfun(@str2double,data);
                
                if numel(dataNumeric) < vm
                    % Zero pad if necessary. Siemens CSA format omits zeros.
                    dataNumeric(vm) = 0;
                end
                
                dcminfo.csa.(tagName) = dataNumeric(1:vm);
            otherwise
                warning('RodgersSpectroTools:UnknownVrType','Unknown VR type: "%s".',vr)
        end
    end
        
        
    %% Helper functions to simulate file I/O
    function [out] = private_read(numBytes)
        out = dcminfo.(strTag)(currdx+(1:numBytes)).';
        currdx=currdx+numBytes;
    end
    
    function [out] = private_readuint32(num)
        out=typecast(private_read(4*num),'uint32');
    end
    
    function [out] = private_readstring(maxchar)
        out = reshape(char(private_read(maxchar)),1,[]);
        terminator = find(out==0,1);
        if numel(terminator)>0
            out=out(1:(terminator-1));
        end
    end
    
end

end
