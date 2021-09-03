function [cfg_game, data_passation] = fastACI_experiment(experiment, Subject_ID, Condition)
% function [cfg_game, data_passation] = fastACI_experiment(experiment, Subject_ID, Condition)
%
%
% Changes by AO:
%   - cfg_game.resume set to 1 (for oui) or 0 (for non)
%   - cfg_game.ordre_aleatoire renamed to 'randorder_idxs'
%   - init_staircase.m renamed to staircase_init.m
%
% TODO:
%   - bSimulation: convert to some option 'artificial_observer'...
%   - Move init_staircase.m to a predefined folder...
%   - Each experiment should have msg_warmup, msg_instructions, msg_mainexp
%   - istarget = (ListStim(n_stim).N_signstartdateal)==2; % in this line it is assumed a 1-I AFC, change in procedures...
%   - convert displayN to silent or something comparable...
%
% TOASK:
%   - Add an additional variable 'dir_results' (splitting dir_main)
%
% Old name: Script2_Passation_EN.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
if nargin < 3
    Condition = []; % No condition
end
if nargin < 2
    Subject_ID = input('Enter the Subject ID (e.g., ''S01''): ');
end
if nargin == 0
    experiments = {'modulationACI','modulationACI_seeds','speechACI_varnet2013','speechACI_varnet2015'};
    Show_cell(experiments);
    bExp = input('Choose the experiment that you want to run from the list above: ');
    
    experiment = experiments{bExp};
end

Subject_ID_full = Subject_ID;
if sum(Subject_ID=='_')
    idx = find(Subject_ID=='_');
    Subject_ID = Subject_ID(1:idx-1);
end

switch Subject_ID
    case {'dau1997','king2019','osses2021'} % alphabetical order
        bSimulation = 1;
        
    case {'dau1997_preproc'}
        error('Model %s is deprecated...',Subject_ID);
    otherwise
        bSimulation = 0;
end

if sum(Subject_ID=='_') && bSimulation == 0
    error('%s: Subject name (Subject_ID) with the character ''_'' is reserved for the use of auditory models. Please define a new Subject_ID without that character',upper(mfilename))
end

experiment_full = experiment; % it will be used for file naming

%%% Checking whether there are separable conditions (separator='-')
if sum(experiment=='-') % Checking whether it contains an hyphen
    experiment   = strsplit(experiment,'-');
    experiment = experiment{1};
end
%%%

% -------------------------------------------------------------------------
% 1. Loading set-up: 
%    1.1. Loads cfgcrea*.mat
[dir_results, dir_results_completed] = Check_local_dir_data(experiment_full,Subject_ID);

filter2use = [Subject_ID '_' experiment_full];
if ~isempty(Condition)
    filter2use = [filter2use '_' Condition];
end
stored_cfg = Get_filenames(dir_results,['cfgcrea*' filter2use '.mat']);
N_stored_cfg = length(stored_cfg);
if N_stored_cfg==1
    var      = load([dir_results stored_cfg{1}]);
    %%% Check compatibility:
    script_legacy = [experiment '_legacy'];
    if exist([script_legacy '.m'],'file')
        fprintf('%s: Checking forward-compatibility of the stored creation file\n',upper(mfilename))
        exp2eval = sprintf('var.cfg_crea = %s(var.cfg_crea);',script_legacy);
        eval(exp2eval);
    end
    cfg_game = var.cfg_crea;
    % cfg_game.cfg_crea   = var.cfg_crea;
elseif N_stored_cfg > 1
    error('Multiple participants option: has not been validated yet (To do by AO)')
else
    % try
        cfg_game = fastACI_experiment_init(experiment_full,Subject_ID, Condition);
    % catch me
    %     error('%s: fastACI_experiment_init failed\n\t%s',upper(mfilename),me.message);
    % end
end

if ~isfield(cfg_game,'resume')
    cfg_game.resume = []; % 'no' -> new game, 'yes' -> load last saved game, [] -> load last saved game if exists or start a new game
end

% -------------------------------------------------------------------------
%     1.2. Loading parameters: Looks for an existing 'game' (or previous 
%          session for the same participant)
if isempty(cfg_game.resume)
    ListSavegame = Get_filenames(dir_results, ['savegame*' filter2use '.mat']);
    if isempty(ListSavegame)
        cfg_game.resume = 0; % 'non';
    else
        cfg_game.resume = 1; % 'oui';
    end
else
	error('Not validated yet...')
end

if ~isfield(cfg_game,'load_name')
    cfg_game.load_name = [];
else
    error('Not validated yet...')
end

if ~isfield(cfg_game,'Language')
    cfg_game.Language = 'EN';
    warning('Using the default language (''EN''). Please specify in the *_cfg.m file another interface language if required.');
end

switch cfg_game.resume
    case {1,'oui','yes'}
        if isfield(cfg_game,'load_name')
            if isempty(cfg_game.load_name)
                exp2filter = ['savegame*' filter2use '.mat'];
                try
                    load_name = Get_savenames(dir_results, exp2filter, dir_results_completed);
                catch
                    error('%s: No savegame to be loaded',upper(mfilename));
                end
            end
            
            Language = cfg_game.Language;
            
            cfg_game = []; % it will be re-loaded now:
            data_passation = []; % it will be re-loaded now:
            ListStim = [];
            
            load(load_name,'cfg_game'); % loads: cfg_game, data_passation
            script_legacy = [experiment '_legacy'];
            if exist([script_legacy '.m'],'file')
                exp2eval = sprintf('cfg_game = %s(cfg_game,0);',script_legacy);
                eval(exp2eval);
            end
            
            load(load_name,'data_passation');
            cfg_game.load_name = load_name;
            i_current  = data_passation.i_current+1;
            i_savegame = i_current;
            
            if ~isfield(cfg_game,'Language')
                cfg_game.Language = Language;
            end
            
            if cfg_game.is_experiment
                % display welcome message
                msg_welcomeback
            end
             
            cfg_game.resume = 1;
            
            if ~isfield(cfg_game,'dBFS')
                warning('You are loading an ''old'' participant')
                cfg_tmp = [];
                exp2eval = sprintf('cfg_tmp = %s_set(cfg_game);',experiment); % experiment dependent
                eval(exp2eval);
                cfg_game.dBFS = cfg_tmp.dBFS;
                cfg_game.dir_stim = cfg_tmp.dir_stim;
            end
            if ~isfield(cfg_game,'experiment')
                cfg_game.experiment = 'modulationACI';
            end
            
            % cfg_game.is_simulation =  bSimulation;
            % cfg_game.is_experiment = ~bSimulation;
        else
            error('%s: No savegame with the specified name',upper(mfilename))
        end
        
        data_passation.resume_trial = [data_passation.resume_trial, i_savegame];
        clock_str = Get_date_and_time_str;
        data_passation.start_date = {data_passation.date_start, clock_str};
        
    case {0,'non','no'}
            
        % Parameters for targets
        exp2eval = sprintf('cfg_game = %s_set(cfg_game);',experiment); % experiment dependent
        eval(exp2eval);
        
        exp2eval = sprintf('cfg_game = %s_cfg(cfg_game);',experiment);
        eval(exp2eval);
        
        % % Parameters for game
        % cfg_game.is_simulation =  bSimulation;
        % cfg_game.is_experiment = ~bSimulation;
        
        cfg_game.script_name{1} = [mfilename('fullpath') '.m'];

        data_passation.resume_trial = 0;
        data_passation.date_start{1} = Get_date_and_time_str;
        
        i_current = 0; % dummy variable
end

cfg_game.is_simulation =  bSimulation;
cfg_game.is_experiment = ~bSimulation;

%%%
% Simulation parameters
if cfg_game.is_simulation == 1
    
    cfg_game.warmup = 0; % warm up is disabled
    cfg_game.sessionsN = 400; warning('Temporal')
    
    if cfg_game.resume
        % Check wether it already completed the task
        if i_current >= cfg_game.N
            % Then we need to start a new 'run'
            cfg_game.resume = 0; % setting back to 0 to start a new simulation
            data_passation = []; % emptying data_passation
            data_passation.date_start{1} = Get_date_and_time_str;
        end
    end
    
    % Model parameters;
    if ~isfield(cfg_game,'model_cfg_script')
        % First time the model is run, then the configuration file is read 
        %   and backed-up locally:
       
        model_cfg_src = [Subject_ID '_cfg'];
        if exist(model_cfg_src,'file')
            exp2eval = sprintf('def_sim = %s;',model_cfg_src);
            eval(exp2eval);
            
            fprintf('Model configuration found on disk. Check whether the configuration is what you expect:\n');
            def_sim
            
            fprintf('Pausing for 10 s. Press ctr+c to cancel the simulations.\n');
            fprintf('If you want to change the simulation parameters, change manually the script %s and re-run the current script (%s).\n',model_cfg_src,mfilename);
            pause(10);
            
        else
            error('Validate here...')
            def_sim = [];
            def_sim.modelname = Subject_ID;

            %%%
            def_sim = fastACI_set_simulation_config(Subject_ID,def_sim);
            %%%
        end
        
        file_config = cfg_game.experiment;
        if isfield(cfg_game,'Cond_extra_1')
            file_config = [file_config '_' cfg_game.Cond_extra_1];
        end
        if isfield(cfg_game,'Cond_extra_2')
            file_config = [file_config '_' cfg_game.Cond_extra_2];
        end
        file_config = sprintf('%s_%s_cfg_%.0f_%.0f_%.0f_%.0fh_%.0fm_%.0fs.m',file_config,Subject_ID,cfg_game.date);
        cfg_game.model_cfg_script      = file_config;
        cfg_game.model_cfg_script_full = [dir_results file_config];
        
        copyfile([fileparts(which(model_cfg_src)) filesep model_cfg_src '.m'],cfg_game.model_cfg_script_full);
    else
        % The the simulation is resuming, we need to read the configuration file:
        addpath(dir_results);
        exp2eval = sprintf('def_sim = %s;',cfg_game.model_cfg_script(1:end-2));
        eval(exp2eval);
        rmpath(dir_results);
    end
        
    % def_sim.template_every_trial = 0;
    % def_sim.templ_num = 10; % 1;
    % def_sim.det_lev = -6; % NaN of the expvar

    switch cfg_game.experiment
        case 'speechACI_Logatome'
            switch cfg_game.Condition
                case 'bump'
                    warning('det_lev for ''bump'' noises is fixed at the moment...')
                    def_sim.det_lev = 0; % -6 of the expvar
            end
    end

    sim_work = [];
    % sim_work.templ_ref = [];
    % sim_work.templ_tar = [];
    cfg_game.def_sim = def_sim;
    if ~isfield(cfg_game,'sessionsN')
        cfg_game.sessionsN = cfg_game.N;
    end
end

if cfg_game.is_simulation == 1
    % % display welcome message
    % msg_welcome
end
%%%

if ~isfield(cfg_game,'feedback')
    cfg_game.feedback = 0; % feedback is disabled by default
end

%% Load stims, create templates for simulation
if cfg_game.resume == 0
    if isfield(cfg_game,'ListStim')
        %%% The field 'ListStim' must exist for experiments that load files
        %     from disc
        ListStim = cfg_game.ListStim; 
        
        if cfg_game.N ~= length(ListStim)
            error('Number of stimuli does not match.')
        end
    end
        
    [list_signals, list_target_signals] = Get_n_signals(cfg_game);
        
    cfg_game.n_targets_sorted = list_signals;
    cfg_game.n_response_correct_target_sorted = list_target_signals;
end

%% Experiment
if cfg_game.resume == 0
    if ~isfield(cfg_game,'stim_order')
        warning('Assigning the stimulus order, this option has been moved to the _init file and will be removed from script %s soon',upper(mfilename))
        if cfg_game.randorder == 1
            cfg_game.stim_order = randperm(cfg_game.N); 
        else
            cfg_game.stim_order = 1:cfg_game.N; 
        end
    end
    debut_i=1;
    data_passation.next_session_stop = debut_i+cfg_game.sessionsN;
else
    debut_i=i_savegame;
    
    if mod(debut_i,cfg_game.sessionsN) == 1
        % Then we update 'next_session_stop'
        data_passation.next_session_stop = debut_i+cfg_game.sessionsN;
    else
        % Nothing to do: the participant chose to take a break
    end
end

%%% Initialises the staircase
str_inout = [];
str_inout.debut_i = debut_i;

[str_inout, cfg_game] = staircase_init(str_inout,cfg_game);
data_passation.reversal_current = str_inout.reversal_current; % always initialised

response   = str_inout.response;
n_correctinarow = str_inout.n_correctinarow;
expvar     = str_inout.expvar;
i_current  = str_inout.i_current;
stepsize   = str_inout.stepsize;
isbreak    = str_inout.isbreak;
%%% Ends: Initialises

is_warmup = cfg_game.warmup;
if cfg_game.is_experiment == 1
    if is_warmup
        % display instructions warmup
        msg_warmup
        
        data_passation_init = data_passation; % initial data_passation
    else
        % display instructions main exp
        msg_mainexp
    end
end
 
N = cfg_game.N;

% cfg_game.dir_path    = dir_main;
cfg_game.dir_results = dir_results;
cfg_game.dir_results_completed = dir_results_completed;

while i_current <= N && i_current~=data_passation.next_session_stop && isbreak == 0
    
    n_stim = cfg_game.stim_order(i_current);
    
    if is_warmup
        data_passation.i_current = i_current;
        data_passation.n_stim(i_current) = n_stim;
        data_passation.expvar(i_current) = expvar;
    end
        
    % Pre-stores information of the current trial
    if ~is_warmup
        clock_now = clock;
        data_passation.i_current = i_current;
        data_passation.n_stim(i_current) = n_stim;
        data_passation.expvar(i_current) = expvar;
        data_passation.date(i_current,:) = clock_now;
    end
    
    if cfg_game.displayN == 1 || cfg_game.is_simulation
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        if isfield(cfg_game,'expvar_description')
            expvar_description = [', ' cfg_game.expvar_description];
        else
            expvar_description = '';
        end
            
        fprintf('\nDependent variable: expvar = %.4f%s \n',expvar,expvar_description);
    end
    
    str_stim = [];
    str2eval = sprintf('[str_stim,data_passation]=%s_user(cfg_game,data_passation);',cfg_game.experiment);
    eval(str2eval);
    stim_normal = str_stim.tuser;
    
    % stim_normal = il_user(str_inout,cfg_game);
    %%% Create signal: end
    
    if cfg_game.is_experiment
        sil4playing = zeros(0.1*cfg_game.fs,1);
        player = audioplayer([sil4playing; stim_normal],cfg_game.fs);
    end
    if cfg_game.is_simulation
        sil4playing = [];
    end
     
    % Display message, play sound and ask for response
    tic
     
    if is_warmup
        switch cfg_game.Language
            case 'EN'
                fprintf('\n    * WARM-UP PHASE *\n\n');
            case 'FR'
                fprintf('\n    * PHASE D''ECHAUFFEMENT *\n\n');
        end
    else
        switch cfg_game.Language
            case 'EN'
                fprintf('\n    * MAIN EXPERIMENT *\n\n');
                fprintf('    Playing stimulus # %.0f of %.0f (Next session stop in %.0f trials)\n',i_current,cfg_game.N,data_passation.next_session_stop-i_current-1);
            case 'FR'
                fprintf('\n    * EXPERIENCE PRINCIPALE *\n\n');
                fprintf('    Ecoute numero %.0f de %.0f (prochaine pause dans %.0f ecoutes)\n',i_current,cfg_game.N,data_passation.next_session_stop-i_current-1);
        end
    end
    
    if cfg_game.is_experiment
        
        play(player)
        if is_warmup
            
            switch cfg_game.Language
                case 'EN'
                    text2show = {'to play the stim again' ['to play a ' cfg_game.response_names{1}] ['to play a ' cfg_game.response_names{2}] 'to leave the warm-up phase'};
                case 'FR'
                    text2show = {'pour rejouer le son' ['pour ecouter un ' cfg_game.response_names{1}] ['pour ecouter un ' cfg_game.response_names{2}] 'pour quitter l''echauffement'};
            end
            response = Response_keyboard([cfg_game.response_names text2show], cfg_game);
            
        else
            switch cfg_game.Language
                case 'EN'
                    text2show = {'to take a break'};
                case 'FR'
                    text2show = {'pour faire une pause'};
            end
            response = Response_keyboard([cfg_game.response_names text2show], cfg_game, 3.14);
        end
        stop(player)
    elseif cfg_game.is_simulation
        % warning('Temporal')
        switch def_sim.decision_script
            case 'aci_detect' 
                % Default for 'dau1997' and 'osses2021'
                [response,sim_work] = aci_detect(cfg_game,data_passation,def_sim,sim_work);
                if i_current == 1
                    if isfield(sim_work,'thres_for_bias')
                        data_passation.thres_for_bias = sim_work.thres_for_bias;
                    end
                    if isfield(sim_work,'type_decision')
                        data_passation.type_decision = sim_work.type_decision;
                    end
                end
                data_passation.decision_var_mue2choose(i_current,:) = sim_work.decision_var_mue2choose(i_current,:);
                
            case 'king2019_detect'
                % Default for 'king2019'
                [response,sim_work,def_sim] = king2019_detect(cfg_game,data_passation,def_sim,sim_work);
        end
        cfg_game.def_sim = def_sim;
                        
        %%% Template:
        %%% Simulation in itself
        % [response,ir,cfg ] = auditorymodel_TMdetect( in, template, cfg );
        
        % fprintf(['analyse stim # ' num2str(i) ' of ' num2str(cfg_game.N) '\n']);
        % Stim_IR = auditorymodel(stim_normal, cfg_game.fs, cfg_game.model);
        % % redefine the template
        % Signal{2} = create_AM(cfg_game.fc, cfg_game.fm, 10^(m/10), cfg_game.stim_dur, cfg_game.fs)';
        % % Signal{2} = [zeros(cfg_game.cfg_crea.stim_dur*cfg_game.fs-length(Signal{2}),1); Signal{2}];
        % cfg_game.IR{2} = auditorymodel(Signal{2}, cfg_game.fs, cfg_game.model);
        % 
        % % audiowrite(['StimsStims/', ListStim(n_stim).name], stim_normal/max(abs(stim_normal)), cfg_game.fs)
        % % [ response ] = auditorymodel_detection( {cfg_game.IR{1}, cfg_game.IR{2}}, Stim_IR, cfg_game.model );
        % % [ response ] = auditorymodel_PEMO( Stim_IR, {cfg_game.IR{1}, cfg_game.IR{2}}, cfg_game.model );
        % % [ response ] = auditorymodel_PEMO( Stim_IR, cfg_game.IR{2}, cfg_game.model );
        % [ response ] = auditorymodel_TMdetect( Stim_IR, cfg_game.IR{2}, cfg_game.model );
    end
    data_passation.response_time(i_current) = toc;
     
    switch response
        case 3.14 % This is a ''pause''
            clock_str = Get_date_and_time_str;
            data_passation.date_end{length(data_passation.date_start)} = clock_str;
            savename = il_get_savename(experiment_full,Subject_ID,clock_str);
            save([dir_results savename], 'i_current', 'ListStim', 'cfg_game', 'data_passation');
            fprintf('  Saving game to "%s.mat" (folder path: %s)\n',savename,dir_results);
            
        case 3 % play again (if warm-up) or take a break (if main experiment)
            if ~is_warmup
                isbreak = 1;
            end
        case 4 % play pure tone
            str_stim = [];
            data_passation_tmp = data_passation;
            idx = find(cfg_game.n_response_correct_target_sorted == 1); % looks for all '1's
            data_passation_tmp.n_stim(i_current) = idx(round( (length(idx)-1)*random('unif',0,1) )+1); % picks up one randomly
                        
            exp2eval = sprintf('str_stim =  %s_user(cfg_game,data_passation_tmp);',experiment);
            eval(exp2eval);
            stim_normal = str_stim.stim_tone_alone;
            
            player = audioplayer([sil4playing; stim_normal],cfg_game.fs);
            playblocking(player)
            
            switch cfg_game.Language
                case 'EN'
                    fprintf('\n    Press any key\n');
                case 'FR'
                    fprintf('\n    Appuyez sur une touche\n');
            end
            pause;

        case 5 % play modulated tone
            str_stim = [];
            data_passation_tmp = data_passation;
            idx = find(cfg_game.n_response_correct_target_sorted == 2); % looks for all '2's
            data_passation_tmp.n_stim(i_current) = idx(round( (length(idx)-1)*random('unif',0,1) )+1); % picks up one randomly
            exp2eval = sprintf('str_stim =  %s_user(cfg_game,data_passation_tmp);',experiment);
            eval(exp2eval);
            stim_normal = str_stim.stim_tone_alone;
            
            player = audioplayer([sil4playing; stim_normal],cfg_game.fs);
            playblocking(player)
            
            switch cfg_game.Language
                case 'EN'
                    fprintf('\n    Press any key\n');
                case 'FR'
                    fprintf('\n    Appuyez sur une touche\n');
            end
            pause;

        case 6 % escape training
            is_warmup = 0;
            clc
            
            cfg_game.warmup = is_warmup;
            
            str_inout = [];
            str_inout.debut_i = debut_i;
            
            [str_inout,cfg_game] = staircase_init(str_inout,cfg_game); % actual initialisation
            
            response   = str_inout.response;
            n_correctinarow = str_inout.n_correctinarow;
            expvar     = str_inout.expvar;
            i_current  = str_inout.i_current;
            stepsize   = str_inout.stepsize;
            isbreak    = str_inout.isbreak;
            
            data_passation_init.data_training = data_passation;
            data_passation = data_passation_init; % 'empty' data_passation to start the main experiment
            
            data_passation.reversal_current = str_inout.reversal_current;
            %%%% TODO %%%%
            % display instructions main exp

        case {1,2} % responded 1 or 2
            
            data_passation.n_response_correct_target(i_current) = cfg_game.n_response_correct_target_sorted(n_stim);
            resp_num = data_passation.n_response_correct_target(i_current);
            iscorrect = (response == resp_num);
            
            % save trial data
            if ~is_warmup
                data_passation.n_responses(i_current) = response;
                data_passation.n_targets(i_current)   = cfg_game.n_targets_sorted(n_stim);
                data_passation.n_response_correct_target(i_current) = resp_num;
                data_passation.is_correct(i_current) = iscorrect;
            end
            if is_warmup || cfg_game.displayN || cfg_game.feedback
                % ListStim(n_stim).response 
                switch iscorrect
                    case 1
                        switch cfg_game.Language
                            case 'EN'
                                txt_extra = 'You were right';
                            case 'FR'
                                txt_extra = 'Correct';
                        end
                    case 0
                        switch cfg_game.Language
                            case 'EN'
                                txt_extra = 'You were wrong';
                            case 'FR'
                                txt_extra = 'Erreur';
                        end
                end
                
                if isfield(cfg_game,'response_names')
                    resp_name = cfg_game.response_names{cfg_game.n_response_correct_target_sorted(n_stim)};
                else
                    error('Continue validating here...')
                    resp_name = num2str(cfg_game.n_response_correct_target(n_stim));
                end
                % feedback
                switch cfg_game.Language
                    case 'EN'
                        fprintf('\n\t%s => Correct answer was : %.0f (%s)\n\n Press any key to continue.\n',txt_extra,resp_num,resp_name);
                    case 'FR'
                        fprintf('\n\t%s => La bonne reponse etait : %.0f (%s)\n\n Appuyez sur une touche pour continuer.\n',txt_extra,resp_num,resp_name);
                end
                pause(.5); % .5 s pause
            end
            
            if iscorrect
                n_correctinarow = n_correctinarow+1;
            else
                n_correctinarow = 0;
            end
             
            if cfg_game.adapt
                
                switch cfg_game.adapt
                    case {1,'transformed-up-down'}
                        cfg_game = Ensure_field(cfg_game,'rule',[1 2]); % [up down]-rule: [1 2] = 1-up 2-down   
                        cfg_game = Ensure_field(cfg_game,'step_up',1);
                        cfg_game = Ensure_field(cfg_game,'step_down',1);
                        
                    case {2, 'weighted-up-down'}
                        if ~isfield(cfg_game,'rule')
                            error('Check defaults')
                            cfg_game.rule = [1 1]; 
                        end
                        target_score = .707;
                        if ~isfield(cfg_game,'step_up')
                            cfg_game.step_up    = target_score/(1-target_score);
                        end
                        if ~isfield(cfg_game,'step_down')
                            cfg_game.step_down  = 1;
                        end
                    otherwise
                        error('Not validated yet...')
                end
                
                str_inout = [];
                str_inout.iscorrect = iscorrect;
                str_inout.expvar    = expvar;
                str_inout.stepsize  = stepsize;
                str_inout.n_correctinarow = n_correctinarow;
                str_inout.reversal_current = data_passation.reversal_current;
                if isfield(data_passation,'staircase_direction')
                    str_inout.staircase_direction = data_passation.staircase_direction;
                end
                
                [str_inout,cfg_game] = staircase_update(str_inout,cfg_game);
                
                % load updated parameters
                expvar    = str_inout.expvar;
                stepsize  = str_inout.stepsize;
                n_correctinarow = str_inout.n_correctinarow;
                
                data_passation.reversal_current = str_inout.reversal_current;
                data_passation.staircase_direction = str_inout.staircase_direction;
            end
             
            i_current=i_current+1;
        otherwise
            warning('%s: Keyboard response not recognised',upper(mfilename))
    end
    
end
 
%% Save game
clock_str = Get_date_and_time_str;
data_passation.date_end{length(data_passation.date_start)} = clock_str;
savename = il_get_savename(experiment_full, Subject_ID_full, Condition, clock_str);
save([dir_results savename],'cfg_game', 'data_passation');
msg_close

if i_current > N
    % So, the sessions are complete now. 
    
    % 1. Then Get_savenames is run once more and only the last save file will 
    %    be kept in the 'Results' directory:
    Get_savenames(dir_results, exp2filter, dir_results_completed);
    
    % 2. The folder of past sessions will be moved inside the 'Result' folder:
    movefile(dir_results_completed, [dir_results 'Results_past_sessions' filesep]);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fname = il_get_savename(experiment, Subject_ID, Condition, clock_str)

if nargin < 3
    clock_str = Get_date_and_time_str;
end
filter2use = [Subject_ID '_' experiment];

if ~isempty(Condition)
    filter2use = [filter2use '_' Condition];
end

savename = ['savegame_' clock_str];
fname = [savename '_' filter2use];

