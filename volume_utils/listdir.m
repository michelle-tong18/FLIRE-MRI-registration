function [outList] = listdir(inputdir,type,format)

% Lists the contents of a directory
%
%   Usage:
%
%   [outlist] = listdir(inputdir,type)
% 
%   Inputs:
%       inputdir = path to input directory
%       type = 
%           'files' to create a list of files in <inputdir>
%           'dirs' to create a list of directories in <inputdir>
%       format = (optional, assumes cell)
%           'cell' to output cell array
%           'struct' to output struct
%
%   Output: list of directory contents, based on type
%
%   Written by Andrew S Bock Feb 2014
%   Updated by Michelle Tong Aug 2020 to include struct output

if ~exist('format','var'); format = 'cell'; end

D = dir(inputdir);
% remove hidden files from the list
isHidden = strncmp('.',{D(:).name}',1); 
D(isHidden) = [];
% Find directories 
dirList = cell2mat({D(:).isdir}'); % '1' if dir, '0' if not
% Sort by type desired (directory list vs file list)
if strcmp(type,'files')
    D(dirList) = [];
    if strcmp(format,'cell')
        outList = {D(:).name};
    else
        outList = D;
    end
elseif strcmp(type,'dirs')
    D(dirList==0) = [];
    if strcmp(format,'cell')
        outList = {D(:).name};
    else
        outList = D;
    end
else
    disp('missing "type" input')
end
