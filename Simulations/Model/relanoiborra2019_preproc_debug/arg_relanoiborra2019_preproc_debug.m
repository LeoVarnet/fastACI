function definput = arg_relanoiborra2019_preproc_debug(definput)
% function definput = arg_relanoiborra2019_preproc_debug(definput)
%
% This function is temporally named and located in this folder. These are
%     extra settings that are not yet available within AMT.
%
% Author: Alejandro Osses

definput.flags.afb = {'erbspacebw','erbspace'};
definput.keyvals.basef = 8000; % one of the fc's will be exactly basef
definput.keyvals.subject='NH';

% relanoiborra2019: basef=8000 Hz, to match this frequency
%                   bwmul=0.5 because their filter design uses 60 bands
%                      spaced at 0.5 ERBN between 100 and 8000 Hz:
definput.groups.drnl_relanoiborra2019_preproc = {'outerear','middleear', ...
    'hearing_profile','NH', ...
    'internalnoise', ...
    'bwmul',0.5};

% Settings to be used by Alejandro
definput.groups.afb_relanoiborra2019_preproc = {'outerear','middleear', ...
    'basef',[], ...
    'subject','NH', ... % my old local name: hearing_profile 
    'no_internalnoise', ...
    'bwmul',1}; 
