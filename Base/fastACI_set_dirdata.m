function fastACI_set_dirdata
% function fastACI_set_dirdata
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target_path = uigetdir(fileparts(which(mfilename)),'Please indicate the folder where you plan to store all your fastACI data');
target_path = [target_path filesep 'fastACI' filesep];

p = [];
p.script_name = 'fastACI_dirdata';
p.target_path = target_path;

text_to_write = readfile_replace('fastACIsim_path_template.txt',p);

dir_file = [fastACI_basepath 'Utility' filesep];
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
