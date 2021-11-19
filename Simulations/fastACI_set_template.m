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
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = [];
if nargin < 2
    p.target_path = ['''' location ''''];
else
    idxs = strfind(location,cfg_game.experiment);
    idxfs = idxs + length(cfg_game.experiment);
    % p.targetpath = [location(1:idxs(1)-1) cfg_game.experiment_full location(idxfs(1):idxs(2)-1) cfg_game.experiment_full location(idxfs(2):end)];
    p.target_path = sprintf('[''%s'' experiment_full ''%s'' experiment_full ''%s''];',location(1:idxs(1)-1),location(idxfs(1):idxs(2)-1), location(idxfs(2):end));
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
