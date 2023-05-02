function List_experiments = Doc_list_Experiments(dir2look)
% function List_experiments = Doc_list_Experiments(dir2look)
%
% 1. Description: This script lists all the fastACI experiments that are
%      available locally, under the folder fastACI/Experiments/
%
% 2. Example: List of experiments in the default fastACI folder:
%       Doc_list_Experiments;
%  
% Author: Alejandro Osses 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    dir2look = [fastACI_basepath 'Experiments' filesep]; 
end
dirs = Get_filenames(dir2look,'*');

look4compulsory = {'cfg.m','user.m','set.m','init.m'};

List_experiments = [];

N_dirs = length(dirs);
for i = 1:N_dirs
    % Within each folder
    dir2look_here = [dir2look dirs{i} filesep];
    files = Get_filenames(dir2look_here,'*.m'); % looks for m-files
    experiment_candidates = il_get_experiment_candidates(files);
    for i_can = length(experiment_candidates):-1:1
        candidate = experiment_candidates{i_can};
        bMatch = 1; % will be changed to 0 if one of the compulsory files
                    % is not found
        for i_look = 1:length(look4compulsory)
            file2search = [candidate '_' look4compulsory{i_look}];
            if ~exist([dir2look_here file2search],'file')
                bMatch = 0;
            end
        end
        if bMatch == 1
            List_experiments{end+1,1} = experiment_candidates{i_can};
        else
            experiment_candidates(i_can) = []; % candidate removed
        end
    end
end
if nargout == 0
    fprintf('%s.m: %.0f fastACI experiments were found on disk (in %s):\n',mfilename,length(experiment_candidates),dir2look);
    if ~isempty(List_experiments)
        Show_cell(List_experiments);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out_labels = il_get_experiment_candidates(files)

out_labels = [];
for i = 1:length(files)
    str_tmp = strsplit(files{i},'_');
    
    out_labels{i} = str_tmp{1};
    for j = 2:length(str_tmp)-1
        out_labels{i} = [out_labels{i} '_' str_tmp{j}];
    end
end
out_labels = unique(out_labels);
    
