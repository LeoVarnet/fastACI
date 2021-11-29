function cfg_crea = Check_cfg_crea_dirs(cfg_crea)
% function cfg_crea = Check_cfg_crea_dirs(cfg_crea)
%
% Tested for the first time from g20211129_getting_waveforms_SLV.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bUpdate = 0;
if ~exist(cfg_crea.dir_data_experiment,'dir')
    bUpdate = 1;
end

if bUpdate
    dir_data_experiment = [fastACI_paths('dir_data') cfg_crea.experiment_full filesep];
    filesep_orig = cfg_crea.dir_data_experiment(end);
    
    str_dir_data = strsplit(cfg_crea.dir_data_experiment(1:end-1),filesep_orig);
    str_dir = strsplit(cfg_crea.dir_target(1:end-1),filesep_orig);
    
    dir_target = dir_data_experiment;
    for i = length(str_dir_data)+1:length(str_dir)
        dir_target = [dir_target str_dir{i} filesep];
    end
    
    str_dir = strsplit(cfg_crea.dir_noise(1:end-1),filesep_orig);
    
    dir_noise = dir_data_experiment;
    for i = length(str_dir_data)+1:length(str_dir)
        dir_noise = [dir_noise str_dir{i} filesep];
    end
    
    cfg_crea.dir_data_experiment_orig = cfg_crea.dir_data_experiment;
    cfg_crea.dir_data_experiment      = dir_data_experiment;
    
    cfg_crea.dir_target_orig          = cfg_crea.dir_target;
    cfg_crea.dir_target               = dir_target;
    
    cfg_crea.dir_noise_orig           = cfg_crea.dir_noise;
    cfg_crea.dir_noise                = dir_noise;
end