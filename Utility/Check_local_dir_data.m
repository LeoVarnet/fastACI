function [dir_results,dir_results_completed,dir_data] = Check_local_dir_data(experiment,Subject_ID)
% function [dir_results,dir_results_completed,dir_data] = Check_local_dir_data(experiment,Subject_ID
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    % New location since 8/06/2021:
    dir_data = fastACI_paths('dir_data');
    dir_results           = [dir_data experiment filesep Subject_ID filesep 'Results'               filesep];
    dir_results_completed = [dir_data experiment filesep Subject_ID filesep 'Results_past_sessions' filesep];
    
    if ~exist(dir_results,'dir')
        % if it does not exist
        dir_lvl = 3; % checks dir_results and two levels up (after dir_data)
        Check_if_dir(dir_results,dir_lvl);
    end
    if ~exist(dir_results_completed,'dir')
        mkdir(dir_results_completed);
    end
catch
    warning('It seems that you don''t have the expected fastACI folder structure...')
    % Previous location, in the same folder as the fastACI toolbox:
    dir_main = fileparts(which(mfilename)); % path will be the folder where this file is located...
    dir_main = [dir_main filesep]; % adds separator / if unix or \ if windows
    
    dir_results = [dir_main 'Interim_results' filesep];
    dir_results_completed = [dir_results Subject_ID filesep];
end
