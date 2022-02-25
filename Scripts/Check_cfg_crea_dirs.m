function cfg_crea = Check_cfg_crea_dirs(cfg_crea, dir_new)
% function cfg_crea = Check_cfg_crea_dirs(cfg_crea, dir_new)
%
% Description:
%     This script checks whether all the folders specified in cfg_crea exist
%     locally. If a folder does not exist, assumes that cfg_crea was created
%     using another computer, and replaces the data folder in cfg_cre by the
%     folder fastACI_paths('dir_data').
%
%     Specify dir_new if you want to update the 'main dir' in a creation file.
%     
% Tested for the first time from g20211129_getting_waveforms_SLV.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bUpdate = 0;
if ~isfield(cfg_crea,'dir_data_experiment')
    bUpdate = 1;
elseif ~exist(cfg_crea.dir_data_experiment,'dir')
    bUpdate = 1;
end

dir_data = fastACI_paths('dir_data');
if nargin >= 2
    dir_data = dir_new;
    bUpdate = 1;
end

if ~isfield(cfg_crea,'experiment_full')
    cfg_crea.experiment_full = cfg_crea.experiment;
end

if bUpdate
    dir_data_experiment = [dir_data cfg_crea.experiment_full filesep]; % local directory
    if ~isfield(cfg_crea,'dir_data_experiment')
        cfg_crea.dir_data_experiment = dir_data_experiment;
    end
    filesep_orig = cfg_crea.dir_noise(end); % but it reads the file separator
                                            % from a previously existing directory
    str_dir_data = strsplit(cfg_crea.dir_data_experiment(1:end-1),cfg_crea.dir_data_experiment(end));
    str_dir = strsplit(cfg_crea.dir_target(1:end-1),filesep_orig);
    
    if isfield(cfg_crea,'Subject_ID')
        dir_target = [dir_data_experiment cfg_crea.Subject_ID filesep str_dir{end} filesep];
    end
    if ~exist(dir_target,'dir')
        % If no Subject_ID, then this savegame_file is either very old, or
        %   the subject uses the waveforms of another participant
        dir_target = dir_data_experiment;
        idxi = length(str_dir_data)+1;
        idxf = length(str_dir);
        for i = idxi:idxf
            dir_target = [dir_target str_dir{i} filesep];
        end
    end
    
    str_dir = strsplit(cfg_crea.dir_noise(1:end-1),filesep_orig);
    
    dir_noise = dir_data_experiment;
    if isfield(cfg_crea,'Subject_ID')
        dir_noise = [dir_noise cfg_crea.Subject_ID filesep str_dir{end} filesep];
    end
    if ~exist(dir_noise,'dir')
        dir_noise = dir_data_experiment; % we start the re-assignment again
        idxi = length(str_dir_data)+1;
        idxf = length(str_dir);
        for i = idxi:idxf
            dir_noise = [dir_noise str_dir{i} filesep];
        end
    end
    
    if ~exist(dir_noise,'dir')
        % This means that all our efforts so far have not succeeded
        msg_here = ['Please indicate where the noise folder (' str_dir{end} ') is...'];
        fprintf('%s\n',msg_here);
        dir_noise = uigetdir(dir_data_experiment,msg_here);
        dir_noise = [dir_noise filesep];
    end
    if ~exist(dir_target,'dir')
        % This means that all our efforts so far have not succeeded
        msg_here = 'Please indicate where the speech (target) folder is...';
        fprintf('%s\n',msg_here);
        dir_target = uigetdir(dir_data_experiment,msg_here);
        dir_target = [dir_target filesep];
    end
    
    cfg_crea.dir_data_experiment_orig = cfg_crea.dir_data_experiment;
    cfg_crea.dir_data_experiment      = dir_data_experiment;
    
    cfg_crea.dir_target_orig          = cfg_crea.dir_target;
    cfg_crea.dir_target               = dir_target;
    
    cfg_crea.dir_noise_orig           = cfg_crea.dir_noise;
    cfg_crea.dir_noise                = dir_noise;
end