function [y,ypd] = Get_filenames(directory, exp2filter, extra)
% function [y,ypd] = Get_filenames(directory, exp2filter, extra)
%
% 1. Description:
%   Get file names using a filter specified by exp2filter. The y-variable 
%   corresponds to a cell-array
% 
%   directory     - Directory to be checked
%   exp2filter    - expression to filter file names (e.g. 'sb*.wav', i.e. all
%                   the files starting with 'sb' having wav extension)
%
%   y             - cell array
% 
%   Last update (on 27/02/2020): when no filter is specified the current '.' and
%                   parent directory are ignored.
% 
% 2.1 Example Unix-based:
%   directory        = '~/Documenten/fda_eval_VlMatrix/wav/';
%   extra.bExtension = 0; % To delete extension
%   file_orig        = Get_filenames(directory,['*.wav'],extra);
% 
% 2.2 Example cross-platform:
%   directory = uigetdir('Select a directory with wav-files');
%   file_orig = Get_filenames(directory,['*.wav']);
% 
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium, 2014-2016
% Created in     : 2013-2014
% Last update on : 30/07/2014
% Last use on    : 20/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    extra.bExtension = 1;
end

if nargin < 2 
    exp2filter = '';
end

ytmp = dir([directory filesep exp2filter]);

if isempty(ytmp)
    y = [];
end

id2remove = [];
for i = 1:length(ytmp)
    
    if strcmp(ytmp(i).name,'.') || strcmp(ytmp(i).name,'..') % current and parent directories...
        id2remove(end+1) = i;
    end
    
    if extra.bExtension 
        y{i}   = ytmp(i).name;
    else
        str_tmp = strsplit( ytmp(i).name,'.');
        
        y{i} = [];
        for k = 1:length(str_tmp)-1 % the last part after the last point is removed
            y{i} = [y{i} str_tmp{k}];
        end
        % y{i} = Delete_extension( ytmp(i).name, exp2filter);
    end
    ypd{i} = [directory filesep y{i}];
end

if ~isempty(id2remove)
    y(id2remove) = [];
    ypd(id2remove) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
