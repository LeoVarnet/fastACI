function filepath = Convert_to_double_filesep(filepath)
% function filepath = readfile_replace(filepath)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idxs = strfind(filepath,filesep);

for i = length(idxs):-1:1
    filepath = [filepath(1:idxs(i)-1) filesep filepath(idxs(i):end)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
