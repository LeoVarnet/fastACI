function Check_readme(dir2check)
% function Check_readme(dir2check)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(dir2check(end),filesep)
    dir_name = dir2check(1:end-1);
else
    dir_name = dir2check;
end
dir_name = strsplit(dir_name,filesep);
dir_name = dir_name{end};
fname = [dir2check 'Readme_' dir_name];
if exist([fname '.m'],'file')
    eval(['run(''' fname ''');'])
end