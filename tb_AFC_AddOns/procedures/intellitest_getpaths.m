function misc = intellitest_getpaths(type) 
% function misc = intellitest_getpaths(type) 
%   
%   Get local paths that need to be used for running the Flemish Matrix 
%   material. If 'type' is specified, then 'misc' will be a character 
%   containing the
%   that specific path. See example 2 (below).
% 
% % Example 1: to get all important directories
%       misc = vlmatrix_getpaths;
%
% % Example 2: to get all important directories
%       misc = vlmatrix_getpaths('VlMatrix');
%
% Programmed by Alejandro Osses (ale.a.osses@gmail.com), Hearing Technology,
% WAVES, UGent, Belgium, 2018-2019
% Created on    : 04/10/2019
% Last update on: 04/10/2019
% Last use on   : 04/10/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global global_vars

f = mfilename('fullpath');
try
    % This will be the folder, if fastACI_baseath.m exists:
    dir_MATLAB = fastACI_basepath;
catch
    dir_MATLAB = [fileparts(fileparts(fileparts(f(1:end-1)))) filesep]; % 3 levels up
end

misc.material_location = [dir_MATLAB 'Stimuli' filesep 'Intellitest' filesep];
if ~exist(misc.material_location,'dir')
    if isunix
        misc.material_location  = '/home/alejandro/Documents/Databases/data/speech-materials/french/Intellitest/VCVCVs/'; 
    else
       % error('Leo: Put the exact location of your VCVCVs folder here...')
        misc.material_location  = 'C:\Users\LSP005\ownCloud\Professionnel\Project_fastACI\Intellitest\VCVCVs\'; 
    end
end

% misc.tb_AFC      = [dir_MATLAB 'tb_AFC' filesep]; % Path to AFC toolbox
misc.tb_AFC_AddOns = [dir_MATLAB 'tb_AFC_AddOns' filesep]; 
str = strsplit(mfilename,'_');
str = str{1};
misc.TestLists    = [misc.tb_AFC_AddOns 'procedures'  filesep str filesep 'TestLists' filesep];
if isfield(global_vars,'result_path')
    result_path = global_vars.result_path;
else
    result_path = [misc.tb_AFC_AddOns 'experiments' filesep 'interimresults' filesep]; % default result path
end
misc.result_path  = result_path;
misc.control_path = misc.result_path;

if nargin==1
    if (~isfield(misc,type))
        error('Invalid type');
    end
    misc=misc.(type);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
