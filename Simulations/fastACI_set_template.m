function fastACI_set_template(location, cfg_game)
% function fastACI_set_template(location, cfg_game)
%
% This file is called from aci_detect.m, and stores a template that has just
% been generated.
% You can also use this script to overwrite manually a previously existing
% template, for this you just need to run:
%     fastACI_set_template(location);
% where 'location' is the full path to your stored template.
%
% location = '/home/alejandro/Documents/Databases/data/fastACI_data/speechACI_Logatome-abda-S43M/osses2022a/';
% fastACI_set_template(location);
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error('Temporally disabled...')
p = [];

%%%
if nargin == 2
    p.modelname = cfg_game.Subject_ID;
else
    
    if strcmp(location(end),filesep)
        Nf = length(location)-1;
    else
        Nf = length(location);
    end
    str = strsplit(location(1:Nf),filesep);
    p.modelname = str{end};
    
    location = [fileparts(location(1:Nf)) filesep];
    p.target_path = ['''' location ''''];
end
%%%

if nargin < 2
    % Nothing to do
else
    idxs = strfind(location,cfg_game.experiment);
    idxfs = idxs + length(cfg_game.experiment);
    p.target_path = sprintf('[''%s'' experiment_full ''%s'' experiment_full ''%s''];',location(1:idxs(1)-1),location(idxfs(1):idxs(2)-1), location(idxfs(2):end));
    p.modelname = cfg_game.Subject_ID;
end

if ispc
    p.target_path = Convert_to_double_filesep(p.target_path);
end

text_to_write = readfile_replace('fastACI_location_template.txt',p);

dir_file = [fastACI_basepath 'Utility' filesep];
fname_file = [dir_file 'fastACI_file_template.m'];

fid = fopen(fname_file, 'w');
fwrite(fid, text_to_write);
fclose(fid);
%%%
