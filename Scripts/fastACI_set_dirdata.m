function fastACI_set_dirdata
% function fastACI_set_dirdata
%
% 1. Description:
%   This script generates the a script named fastACI_dir_data.m which will
%   contain the location of the dir_data folder that will be used by fastACI
%   to store all the experimental data including the stimuli and it is also
%   the default location for the post-processed data. 
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

text2show = 'New ''fastACI_data'' folder: Please indicate where you plan to store all the extra toolbox binary data';
fprintf('\t%s\n',text2show);

target_path = uigetdir(fileparts(which(mfilename)),text2show);
target_path = [target_path filesep 'fastACI_data' filesep];

p = [];
p.script_name = 'fastACI_dir_data';
p.target_path = target_path;

if ispc
    p.target_path = Convert_to_double_filesep(p.target_path);
end

text_to_write = readfile_replace('fastACIsim_path_template.txt',p);

dir_file = [fastACI_basepath 'Local' filesep];
fname_file = [dir_file p.script_name '.m'];

if exist(fname_file,'file')
    fprintf('----------------------------------------------------------------------------\n')
    fprintf('file %s exists, \npress any key to continue (will overwrite) or press ctrl+C to cancel \n',fname_file);
    fprintf('----------------------------------------------------------------------------\n')
    pause
end

if ~exist(target_path,'dir')
    mkdir(target_path);
end

fid = fopen(fname_file, 'w');
fwrite(fid, text_to_write);
fclose(fid);
%%%
