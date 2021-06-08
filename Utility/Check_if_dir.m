function Check_if_dir(dir2check,number_lvl_up)
% function Check_if_dir(dir2check,number_lvl_up)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirs2check{1} = dir2check;

if strcmp(dir2check(end),filesep)
    dir2check = dir2check(1:end-1);
end

for i = 2:number_lvl_up
    dir2check = fileparts(dir2check); % goes one level up
    dirs2check{i} = [dir2check filesep];
end

for i = number_lvl_up:-1:1
    
    if ~exist(dirs2check{i},'dir')
        fprintf('\tCreating folder %s\n',dirs2check{i});
        
        mkdir(dirs2check{i});
    end
    
end
