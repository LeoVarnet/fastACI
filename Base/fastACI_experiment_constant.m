function [cfg_game, data_passation] = fastACI_experiment_constant(experiment, Subject_ID, Condition, expvar, varargin)
% function [cfg_game, data_passation] = fastACI_experiment_constant(experiment, Subject_ID, Condition, expvar, varargin)
%
%
% Changes by AO:
%   - cfg_game.resume set to 1 (for oui) or 0 (for non)
%   - cfg_game.ordre_aleatoire renamed to 'randorder_idxs'
%   - init_staircase.m renamed to staircase_init.m
%
% Needs to have:
%   - An existing cfg_crea
%   - The waveforms on disk
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
    case {'dau1997','king2019','relanoiborra2019','maxwell2020','osses2021','osses2022a'} % alphabetical order
        bSimulation = 1;
        
    case {'dau1997_preproc'}
        error('Model %s is deprecated...',Subject_ID);
    otherwise
        bSimulation = 0;
end

% From argument function:
definput.import={'fastACI_experiment'}; % arg_fastACI_experiment.m
if bSimulation
    definput.import{end+1} = 'fastACI_simulations';
    definput.import{end+1} = 'fastACI_simulation_detect';
end
% From argument function:
[flags,keyvals]  = ltfatarghelper({},definput,varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    error('%s: No init file will be run, please run a staircase first\n',upper(mfilename));
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
        %%% Disabling the possibility to resume automatically (added 14/04/2022):
        if cfg_game.resume == 1
            cfg_game.resume = 0;
            fprintf('\t%s: By-passing ''resume'', as we are controlling the start and end trial using Ni and Nf\n',upper(mfilename));
        end
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
            
            [cfg_game, data_passation, ListStim] = Convert_ACI_data_type(load_name);
            script_legacy = [experiment '_legacy'];
            if exist([script_legacy '.m'],'file')
                exp2eval = sprintf('cfg_game = %s(cfg_game,0);',script_legacy);
                eval(exp2eval);
            end
            
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

    if ~strcmp(cfg_game.Language,'EN')
        fprintf('%s: Switching the language of the simulations to English\n',upper(mfilename));
        cfg_game.Language = 'EN';
    end
    cfg_game.warmup = 0; % warm up is disabled
    cfg_game.sessionsN = 400; 
    
    if cfg_game.resume
        % Check wether it already completed the task
        if i_current >= cfg_game.N
            % Then we need to start a new 'run'
            cfg_game.resume = 0; % setting back to 0 to start a new simulation
            data_passation = []; % emptying data_passation
            data_passation.date_start{1} = Get_date_and_time_str;
        end
    end
    
    if ~isempty(keyvals.model_cfg_script)
        if isfield(cfg_game,'model_cfg_script')
            warning('Ignoring the field cfg_game.model_cfg_script and using the model_cfg_script indicated in the keyvalues');
        end
        cfg_game.model_cfg_script = keyvals.model_cfg_script;
    end
    
    % Model parameters;
    if ~isfield(cfg_game,'model_cfg_script')
        % First time the model is run, then the configuration file is read 
        %   and backed-up locally:
       
        path_where_supposed = [fastACI_basepath 'Local' filesep];
        model_cfg_src = [Subject_ID '_cfg'];
        
        if exist([path_where_supposed model_cfg_src],'file')
            path = fileparts(which(model_cfg_src));
            path = [path filesep];
            if ~strcmp(path,path_where_supposed)
                error('The script %s is supposed to be in %s,\n(not in %s)',model_cfg_src,path_where_supposed,path);
            end
            exp2eval = sprintf('def_sim = %s(keyvals);',model_cfg_src);
            eval(exp2eval);
            
            fprintf('Model configuration found on disk. Check whether the configuration is what you expect:\n');
            def_sim
            
            fprintf('Pausing for 10 s. Press ctr+c to cancel the simulations.\n');
            fprintf('%s: If you want to change the simulation parameters, change manually the script %s and re-run the current script.\n',upper(mfilename),model_cfg_src);
            pause(10);
            
        else
            def_sim = [];
            def_sim.modelname = Subject_ID;

            %%%
            def_sim = fastACI_set_simulation_config(Subject_ID,def_sim,keyvals);
            %%%
        end
        
        file_config = cfg_game.experiment;
        if isfield(cfg_game,'Cond_extra_1')
            file_config = [file_config '_' cfg_game.Cond_extra_1];
        end
        if isfield(cfg_game,'Cond_extra_2')
            file_config = [file_config '_' cfg_game.Cond_extra_2];
        end
        file_config = sprintf('%s_%s_cfg_%.0f_%.0f_%.0f_%.0fh%.0fm.m',file_config,Subject_ID,cfg_game.date(1:5));
        cfg_game.model_cfg_script      = file_config;
        cfg_game.model_cfg_script_full = [dir_results file_config];
        dir_here = [fileparts(which(model_cfg_src)) filesep];
        copyfile([dir_here model_cfg_src '.m'],cfg_game.model_cfg_script_full);
    elseif exist([dir_results cfg_game.model_cfg_script '.m'],'file')
        % The the simulation is resuming, we need to read the configuration file:
        addpath(dir_results);
        exp2eval = sprintf('def_sim = %s(keyvals);',cfg_game.model_cfg_script(1:end-2));
        eval(exp2eval);
        rmpath(dir_results);
    else
        path_src = [fastACI_basepath 'Simulations' filesep 'Stored_cfg' filesep];
        path_dst = [fastACI_basepath 'Simulations' filesep];
        model_cfg_dst = [Subject_ID '_cfg'];
        
        fsrc = [path_src cfg_game.model_cfg_script '.m'];
        fdst = [path_dst model_cfg_dst '.m'];
        if exist(fsrc,'file')
            copyfile(fsrc,fdst);
        end
        
        exp2eval = sprintf('def_sim = %s(keyvals);',model_cfg_dst);
        eval(exp2eval);  
    end
        
    sim_work = [];
    sim_work.templ_ref = [];
    sim_work.templ_tar = [];
    
    if def_sim.bStore_template == 1
        if exist('fastACI_file_template.m','file')
            file_template = fastACI_file_template(cfg_game.experiment_full, ...
                                    Subject_ID, def_sim.type_decision, keyvals);
            if exist(file_template,'file')
                fprintf('Pausing for 10 s. Press ctr+c to cancel the simulations.\n');
                fprintf('%s: Template found on disk, if this is not what you want, remove/rename the file and re-run the simulations\n (template file: %s).\n',upper(mfilename),file_template);
                pause(10);

                templ_ref = [];
                templ_tar = [];
                load(file_template,'templ_ref');
                load(file_template,'templ_tar');

                fprintf('Stored template found at: %s\n',file_template);
                def_sim.bStore_template = 0;
                
                sim_work.templ_ref = templ_ref;
                sim_work.templ_tar = templ_tar;
            end
        else
           % Nothing to do here: the template will be later assessed    
        end
    end
    
    sim_work.bStore_template = def_sim.bStore_template;
    cfg_game.def_sim = def_sim;
    if ~isfield(cfg_game,'sessionsN')
        cfg_game.sessionsN = cfg_game.N;
    end
end

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
            warning('Number of stimuli does not match.')
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

outs_trial = [];

%%% Initialises the 'staircase' from which the expvar value will be obtained
str_inout = [];
str_inout.debut_i = debut_i;
[str_inout, cfg_game] = staircase_init(str_inout,cfg_game);
if nargin < 4
    warning('No expvar value was specified. The experiment will be run at expvarstart')
    expvar = str_inout.expvar;
end
    
i_current  = str_inout.i_current;
cfg_game.adapt = 0;
isbreak = 0;
%%% Ends: Initialises

data_passation_init = [];
is_warmup = cfg_game.warmup;
if cfg_game.is_experiment == 1
    if is_warmup
        % display instructions warmup
        msg_warmup
        
        data_passation_init = data_passation; % initial data_passation
    else
        % display instructions main exp
        msg_mainexp;
    end
end
 
if isempty(keyvals.Ni)
    Ni = data_passation.i_current;
else
    Ni = keyvals.Ni;
    i_current = Ni;
end
if isempty(keyvals.Nf)
    % N = cfg_game.sessionsN;
else
    Nf = keyvals.Nf;
    % N = Nf-Ni+1;
end
N = cfg_game.N;

data_passation.next_session_stop = Nf;
data_passation.i_current = Ni;

cfg_game.dir_results = dir_results;
cfg_game.dir_results_completed = dir_results_completed;

%%% Checking the calibration:
global global_vars
if cfg_game.is_experiment
    if isfield(global_vars,'dBFS')
        dBFS_playback = global_vars.dBFS;
    else
        dBFS_playback = cfg_game.dBFS; % same dBFS as for the stored sounds...
    end
    cfg_game.dBFS_playback = dBFS_playback;
end
if isfield(cfg_game,'dBFS_playback')
    cfg_game.gain_play = 10^((cfg_game.dBFS_playback-cfg_game.dBFS)/20);
    fprintf('%s: Actual dBFS=%.1f dB SPL\n',upper(mfilename),cfg_game.dBFS_playback);
else
    cfg_game.gain_play = 1;
end
%%% End checking the calibration

while i_current <= N && i_current~=data_passation.next_session_stop && isbreak == 0
     
    ins_trial = [];
    if cfg_game.adapt
        ins_trial.stepsize        = outs_trial.stepsize;
        ins_trial.n_correctinarow = outs_trial.n_correctinarow;
    end
    if cfg_game.is_simulation
        if isfield(cfg_game,'def_sim')
            ins_trial.def_sim = cfg_game.def_sim; % normally this is the case after running the first trial
        else
            ins_trial.def_sim  = def_sim;
        end
        ins_trial.sim_work = sim_work;
    end
    ins_trial.debut_i = debut_i;
    if debut_i == 1
        ins_trial.data_passation_init = data_passation_init;
    end

    data_passation.i_current = i_current;
    %%%
    [cfg_game, data_passation, outs_trial, sim_work] = ...
        fastACI_trial_current(cfg_game, data_passation, expvar, ins_trial, flags, keyvals);
    response  = outs_trial.response;
    i_current = outs_trial.i_current;
    isbreak   = outs_trial.isbreak;
    expvar    = outs_trial.expvar;    
end 
