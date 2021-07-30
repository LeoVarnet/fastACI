function [cfg_inout,s_current,bGenerate_stimuli] = Check_seeds_and_initialise(cfg_inout)
% function [cfg_inout,s_current,bGenerate_stimuli] = Check_seeds_and_initialise(cfg_inout)
%
% 1. Description:
%       It assigns the following fields:
%       - cfg_inout.seeds_order
%       - cfg_inout.seeds_order_method
%       - cfg_inout.stim_order = randperm(cfg_inout.N); 
%       - cfg_inout.bRove_level = var.cfg_crea.bRove_level;
%       - cfg_inout.Rove_level  = var.cfg_crea.Rove_level;
%       - cfg_inout.Rove_range  = var.cfg_crea.Rove_range;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(cfg_inout,'seeds_order') && ~isfield(cfg_inout,'stim_order')
    bDo_the_check = 1;
else
    % In this case, the participant already exists, so the current script is skipped
    bDo_the_check = 0;
    bGenerate_stimuli = 1; % overwritten to 0 if dir_noise exists
    if exist(cfg_inout.dir_noise,'dir')
        bGenerate_stimuli = 0;
    end
end

%%%
s_current = rng; % gets current seed
%%%
if bDo_the_check
    list_other_subjects = Get_filenames(cfg_inout.dir_data_experiment,'*');
    for i = length(list_other_subjects):-1:1
        crea_extern = [cfg_inout.dir_data_experiment list_other_subjects{i} filesep 'Results' filesep];
        files = Get_filenames(crea_extern,['cfgcrea*' cfg_inout.Condition '.mat']);

        if isempty(files)
            list_other_subjects(i) = [];
        end
    end

    if ~isempty(list_other_subjects) % ~isempty(files)
        % In case there are completed participants found on disk:
        disp('Do you want to copy the seeds from any of the following subjects?')
        Show_cell(list_other_subjects);
        bInput = input('Type the number of the subject you want to copy (default = 0, new seed): ');
    else
        % Then seeds will be newly generated:
        bInput = 0;
    end

    if bInput == 0
        % Initialises the seed numbers on the fly:
        seed_number_max = 4*cfg_inout.N; % arbitrary number, seeds will range from 0 to 4*N
        method = 'perm';
        type_seed = 'shuffle';

        try
            rng(type_seed); % seed based on the computer's clock
        catch me
            rng('default');
            rng(type_seed); % seed based on the computer's clock
        end

        switch method
            case 'unif'
                error('Not tested recently')
            case 'perm'
                numbers = randperm(seed_number_max);
                seed_numbers = numbers(1:cfg_inout.N); % takes only the first 'N' numbers
        end

        cfg_inout.seeds_order = seed_numbers; % to be used sequentially
        cfg_inout.seeds_order_method = method;

        if cfg_inout.randorder == 1
            cfg_inout.stim_order = randperm(cfg_inout.N); 
        else
            error('Not tested recently')
        end

        bGenerate_stimuli = 1; % not changed
    else

        crea_extern = [cfg_inout.dir_data_experiment list_other_subjects{bInput} filesep 'Results' filesep];
        files = Get_filenames(crea_extern,['cfgcrea*' cfg_inout.Condition '.mat']);

        if length(files) ~= 1
            error('Multiple create files found for the subject from which the data are being copied from...')
        end
        var = load([crea_extern files{1}]);

        %%%
        if exist(var.cfg_crea.dir_noise,'dir')
            cfg_inout.dir_noise  = var.cfg_crea.dir_noise;
        else
            idx = strfind(cfg_inout.dir_noise,cfg_inout.dir_noise(end));
            dir_noise = cfg_inout.dir_noise(idx(end-1)+1:end-1);
            dir_new   = [cfg_inout.dir_data_experiment list_other_subjects{bInput} filesep];
            cfg_inout.dir_noise = [dir_new dir_noise filesep];
        end
        
        if isfield(var.cfg_crea,'dir_speech')
            warning('dir_speech is deprecated (old cfg_crea is being used)')
            var.cfg_crea.dir_target = var.cfg_crea.dir_speech;
        end
        if exist(var.cfg_crea.dir_target,'dir')
            cfg_inout.dir_target = var.cfg_crea.dir_target;
        else
            idx = strfind(cfg_inout.dir_target,cfg_inout.dir_target(end));
            dir_target = cfg_inout.dir_target(idx(end-1)+1:end-1);
            dir_new   = [cfg_inout.dir_data_experiment list_other_subjects{bInput} filesep];
            cfg_inout.dir_target = [dir_new dir_target filesep];

            cfg_inout.seeds_order        = var.cfg_crea.seeds_order; % to be used sequentially
            cfg_inout.seeds_order_method = var.cfg_crea.seeds_order_method;
            cfg_inout.stim_order         = var.cfg_crea.stim_order; 
        end
        %%% 
        
        if isfield(var.cfg_crea,'bRove_level')
            cfg_inout.bRove_level = var.cfg_crea.bRove_level;
            cfg_inout.Rove_level  = var.cfg_crea.Rove_level;
            cfg_inout.Rove_range  = var.cfg_crea.Rove_range;
        end
        fprintf('Metadata successfully copied from %s\n',crea_extern);

        bGenerate_stimuli = 0; % assuming that the folders from the metadata contain all the files already

    end
end
