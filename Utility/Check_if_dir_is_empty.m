function bEmpty = Check_if_dir_is_empty(dir2check,filter2use)
% function bEmpty = Check_if_dir_is_empty(dir2check,filter2use)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

files = Get_filenames(dir2check,filter2use);
bEmpty = isempty(files);