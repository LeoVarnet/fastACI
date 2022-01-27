function [subjectname,dirname_new] = Get_subjectname_from_dirname(dirname,new_subjectname)
% function [subjectname,dirname_new] = Get_subjectname_from_dirname(dirname,new_subjectname)
%
% 1. Description: 
%
% 2. Example:
%       subj = Get_subjectname_from_dirname(cfg_game.dir_noise);
%
% Programmed by Alejandro Osses, ENS, France 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirname_new = [];
  
filesep_here = dirname(end);
switch filesep_here
    case {'/','\'}
        % Nothing to do
    otherwise
        error('A file separator string is expected as last character...')
end
str = strsplit(dirname(1:end-1),filesep_here);
subjectname = str{end-1};

if nargin >= 2
    dirname_new = fileparts(fileparts(dirname(1:end-1))); % two levels up
    dirname_new = [dirname_new filesep_here new_subjectname filesep_here str{end} filesep_here];
end
    
disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end