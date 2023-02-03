function [cfg_crea, bReproducible, bCopied] = fastACI_experiment_init_from_cfg_crea(cfg_crea_file,dir_data)
% function [cfg_crea, bReproducible, bCopied] = fastACI_experiment_init_from_cfg_crea(cfg_crea_file,experiment_full,Subject_ID, Condition)
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    dir_data = fastACI_dir_data;
end
load(cfg_crea_file,'cfg_crea');

bCopied = 0;
bReproducible = 1;

if ~isfield(cfg_crea,'seeds_order')
    fprintf('No field seeds_order... No reproducible waveforms...\n');
    bReproducible = 0;
end

if ~isfield(cfg_crea,'seeds_order_method')
    fprintf('No field seeds_order... No reproducible waveforms...\n');
    bReproducible = 0;
end

if bReproducible
    fprintf('The waveforms seem to be reproducible...\n')
end

N_created = 0;
if bReproducible
    cfg_crea = Check_cfg_crea_dirs(cfg_crea, dir_data, bReproducible);
    if ~exist(cfg_crea.dir_data_experiment,'dir')
        mkdir(cfg_crea.dir_data_experiment);
    end
    bGenerate_noises = Check_if_dir_is_empty(cfg_crea.dir_noise,'*.wav');
    
    if bGenerate_noises
        %%% Checking if an *.init file is found on disc:
        script_name = sprintf('%s_init',cfg_crea.experiment);
        if exist([script_name '.m'],'file')
            fprintf('\tScript %s.m found on disc...\n',script_name);
            exp2eval = sprintf('%s(cfg_crea);',script_name); % init script without output arguments...
            eval(exp2eval);
        end
    
        % dir_noise_orig = [cfg_crea.dir_noise 'Original' filesep];
        files = Get_filenames(cfg_crea.dir_noise,'*.wav');
        N_created = length(files);
    end
end

if bReproducible
    bCopied = Copy_crea_or_savegame_to_subject_dir(cfg_crea,cfg_crea_file,'cfg_crea');
    
    if bCopied
        dir_name_results = '1-experimental_results';
        dir_save = [fileparts(fileparts(cfg_crea_file(1:end-1))) filesep dir_name_results filesep];
        if exist(dir_save,'dir')
            file_crea = Get_filenames(dir_save,['savegame*' cfg_crea.Condition '.mat']);
            if length(file_crea) == 1
                load([dir_save file_crea{1}],'cfg_game')
            end
            cfg_game.dir_target = cfg_crea.dir_target;
            cfg_game.dir_noise = cfg_crea.dir_noise;
            if isfield(cfg_game,'dir_data_experiment')
                cfg_game.dir_data_experiment = [fileparts(fileparts(cfg_game.dir_noise(1:end-1))) filesep];
            end
            cfg_crea_file = [dir_save file_crea{1}];
            bCopied = Copy_crea_or_savegame_to_subject_dir(cfg_game,cfg_crea_file,'cfg_game');
        else
            fprintf('%s.m: No ''%s'' was found. Probably the data for this experiment have not been collected yet\n',upper(mfilename),dir_name_results);
        end
    end
end