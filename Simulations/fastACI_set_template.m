function fastACI_set_template(location)
% function fastACI_set_template(location)
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = [];
p.target_path = location;

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
