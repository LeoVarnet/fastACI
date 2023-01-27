function fastACI_set_amtoolbox
% function fastACI_set_amtoolbox
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

text2show = 'Please indicate the folder where AMT 1.0 is located';
fprintf('\t%s\n',text2show);

target_path = uigetdir(fileparts(which(mfilename)),text2show);
target_path = [target_path filesep];

p = [];
p.script_name = 'fastACI_dir_amtoolbox';
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
