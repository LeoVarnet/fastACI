function crosspred = Check_if_crosspred(fnameACI,file_crosspred)
% function crosspred = Check_if_crosspred(fnameACI,file_crosspred)
%
% Used in l20211217_crossprediction_version_Alejandro_v2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
crosspred = [];

% info_toolbox = Get_toolbox_info(script_str);
if ischar(file_crosspred)
    file_crosspred = {file_crosspred};
end
N_crosspred = length(file_crosspred);
idx = nan([1 N_crosspred]);

[dir_where,fname] = fileparts(fnameACI);
dir_where = [dir_where filesep fname filesep]; 
if exist(dir_where,'dir')
    fname = [dir_where 'Crosspred.mat'];
    if exist(fname,'file')
        var = load(fname);
        N_existing = length(var.crosspred);
        for i = 1:N_crosspred % each requested file_crosspred
            for j = 1:N_existing
                if strcmp(file_crosspred{i},var.crosspred(j).ACI_crosspred)
                    idx(i) = j;
                end
            end
        end
        
        if sum(isnan(idx))==0 
            % then all files were found
            crosspred = var.crosspred(idx);
        end
    end
end