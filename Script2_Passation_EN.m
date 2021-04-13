function Script2_Passation_EN(experiment, Subject_ID, Condition)
% function Script2_Passation_EN(experiment, Subject_ID, Condition)
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

switch Subject_ID
    case {'dau1997','dau1997_preproc'}
        bSimulation = 1;
    otherwise
        bSimulation = 0;
end

% -------------------------------------------------------------------------
% 1. Loading set-up: 
%    1.1. Loads cfgcrea*.mat
[path,name,ext]=fileparts(which(mfilename)); % path will be the folder where this file is located...
dir_main = [path filesep];    %'C:\Users\Varnet L�o\Dropbox\Professionnel\Matlab\MyScripts\modulationACI\AM';
dir_results = [dir_main 'Interim_results' filesep];
dir_results_subj = [dir_results Subject_ID filesep];

filter2use = [Subject_ID '_' experiment];
if ~isempty(Condition)
    filter2use = [filter2use '_' Condition];
end
stored_cfg = Get_filenames(dir_results_subj,['cfgcrea*' filter2use '.mat']);
N_stored_cfg = length(stored_cfg);
if N_stored_cfg==1
    var      = load([dir_results_subj stored_cfg{1}]);
    cfg_game = var.cfg_crea;
    % cfg_game.cfg_crea   = var.cfg_crea;
elseif N_stored_cfg > 1
    error('Multiple participants option: has not been validated yet (To do by AO)')
else
    try
        cfg_game = Script1_Initialisation_EN(experiment,Subject_ID, Condition);
    catch me
        error('%s: no cfg_crea available\n\t%s',upper(mfilename),me.message);
    end
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

switch cfg_game.resume
    case {1,'oui','yes'}
        if isfield(cfg_game,'load_name')
            if isempty(cfg_game.load_name)
                ListSavegame = dir([dir_results 'savegame*' filter2use '.mat']);
                if isempty(ListSavegame)
                    error('%s: No savegame to be loaded',upper(mfilename));
                end
                index_savegame = 0;
                bytes = 0;
                
                index_all = 1:length(ListSavegame); % index of all the files that were found
                for j=index_all
                    if ListSavegame(j).bytes > bytes
                        index_savegame=j;
                        bytes = ListSavegame(j).bytes; % looks for the largest MAT file
                    end
                end
                load_name = [dir_results ListSavegame(index_savegame).name];
                
                % --- Now we will move the old (completed sessions) to the 
                %     subject's folder.
                index_all(index_savegame) = [];
                if ~isempty(index_all)
                    for j=1:length(index_all)
                        movefile([dir_results      ListSavegame(index_all(j)).name], ... % src
                                 [dir_results_subj ListSavegame(index_all(j)).name]);
                    end
                end
                % ---
            end
            
            cfg_game = []; % it will be re-loaded now:
            ListStim = [];
            
            load(load_name); % loads: cfg_game, data_passation
            cfg_game.load_name = load_name;
            i_current  = data_passation.i_current+1;
            i_savegame = i_current;
            
            % display welcome message
            msg_welcomeback
             
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
            
            cfg_game.is_simulation =  bSimulation;
            cfg_game.is_experiment = ~bSimulation;
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
        
        % Parameters for game
        cfg_game.is_simulation =  bSimulation;
        cfg_game.is_experiment = ~bSimulation;
        % Simulation parameters
        if cfg_game.is_simulation == 1
            % modelparameters;
            % cfg_game.fadein_s           = 0;
            % cfg_game.N_template         = 0;
            cfg_game.warmup = 0; % warm up is disabled
            % cfg_game.sessionsN          = 500;
            % cfg_game.stim_dur           = 0.75;
            def_sim = [];
            def_sim.template_script = 'model_template';
            sim_work = [];
            
            cfg_game.sessionsN = cfg_game.N;
        end

        cfg_game.script_name{1} = [mfilename('fullpath') '.m'];

        data_passation.resume_trial = 0;
        data_passation.date_start{1} = Get_date_and_time_str;

        if cfg_game.is_simulation == 1
            % % display welcome message
            % msg_welcome
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
            error('Number of stimuli does not match.')
        end
    end
        
    [list_signals, list_target_signals] = Get_n_signals(cfg_game);
        
    cfg_game.n_targets_sorted = list_signals;
    cfg_game.n_response_correct_target_sorted = list_target_signals;
end

% Create template
 
if cfg_game.is_simulation == 1
    
    warning('Modelling under construction: Not validated yet...')
    % Signal{1} = create_AM(cfg_game.fc, cfg_game.fm, 0, cfg_game.stim_dur, cfg_game.fs)';
    % Signal{2} = create_AM(cfg_game.fc, cfg_game.fm, 10^(cfg_game.m_start/10), cfg_game.stim_dur, cfg_game.fs)';
    % 
    % if cfg_game.N_template>0 
    %     cfg_game.IR{1} = 0; cfg_game.IR{2} = 0;
    %     for i=1:cfg_game.N_template
    %         fprintf(['Template: stim #' num2str(i) ' of ' num2str(cfg_game.N_template) '\n']);
    % 
    %         WavFile = strcat(cfg_game.folder_name, '/', ListStim(i).name);
    %         [noise,cfg_game.fs] = audioread(WavFile); noise = noise/std(noise);
    %         fadein_samples = cfg_game.fs*cfg_game.fadein_s;
    %         Target{1} = generate_stim( Signal{1}, noise, cfg_game.SNR, fadein_samples);
    %         Target{2} = generate_stim( Signal{2}, noise, cfg_game.SNR, fadein_samples);
    % 
    %         %plot_modep(cfg_game.fc,cfg_game.fmc, auditorymodel(Target{1}, cfg_game.fs, cfg_game.model)/cfg_game.N_template)
    % 
    %         IRind1(:,:,:,i) = auditorymodel(Target{1}, cfg_game.fs, cfg_game.model)/cfg_game.N_template; %cfg_game.IR{1} +
    %         IRind2(:,:,:,i) = auditorymodel(Target{2}, cfg_game.fs, cfg_game.model)/cfg_game.N_template; %cfg_game.IR{2} +
    %     end
    %     cfg_game.IR{1} = mean(IRind1,4);
    %     cfg_game.IR{2} = mean(IRind2,4);
    % else
    %     cfg_game.IR{1} = auditorymodel(Signal{1}, cfg_game.fs, cfg_game.model);
    %     cfg_game.IR{2} = auditorymodel(Signal{2}, cfg_game.fs, cfg_game.model);
    % end
    % %cfg_game.IR{1} = zeros(size(cfg_game.IR{2}));
    % cfg_game.Template = cfg_game.IR{2} - cfg_game.IR{1};     
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

str_inout = staircase_init(str_inout,cfg_game);

% expvar_i   = str_inout.expvar;

response   = str_inout.response;
n_correctinarow = str_inout.n_correctinarow;
expvar     = str_inout.expvar;
i_current  = str_inout.i_current;
stepsize   = str_inout.stepsize;
isbreak    = str_inout.isbreak;
%%% Ends: Initialises

iswarmup = cfg_game.warmup;
if cfg_game.is_experiment == 1
    if iswarmup
        % display instructions warmup
        msg_warmup
        
        data_passation_init = data_passation; % initial data_passation
    else
        % display instructions main exp
        msg_mainexp
    end
end
 
N = cfg_game.N;
% if ~isfield(cfg_game,'ListStim')
%     cfg_game.ListStim = ListStim;
% end
% warning('This is temporary')

cfg_game.dir_path    = dir_main;
cfg_game.dir_results = dir_results;

while i_current <= N && i_current~=data_passation.next_session_stop && isbreak == 0
    
    n_stim = cfg_game.stim_order(i_current);
    
    if iswarmup
        data_passation.i_current = i_current;
        data_passation.n_stim(i_current) = n_stim;
        data_passation.expvar(i_current) = expvar;
    end
        
    % Pre-stores information of the current trial
    if ~iswarmup
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
     
    if iswarmup
        fprintf('\n    * WARM-UP PHASE *\n\n');
    else
        fprintf('\n    * MAIN EXPERIMENT *\n\n');
        fprintf('    Playing stimulus # %.0f of %.0f (Next session stop in %.0f trials)\n',i_current,cfg_game.N,data_passation.next_session_stop-i_current-1);
    end
    
    if cfg_game.is_experiment
        
        play(player)
        if iswarmup
            response = Reponse_clavier([cfg_game.response_names {'to play the stim again' ['to play a ' cfg_game.response_names{1}] ['to play a ' cfg_game.response_names{2}] 'to leave the warm-up phase'}]);
        else
            response = Reponse_clavier([cfg_game.response_names {'to take a break'}], 3.14);
        end
        stop(player)
    elseif cfg_game.is_simulation
        [response,sim_work] = aci_detect(cfg_game,data_passation,def_sim,sim_work);
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
            savename = il_get_savename(experiment,Subject_ID,clock_str);
            save([dir_results savename], 'i_current', 'ListStim', 'cfg_game', 'data_passation');
            fprintf('  Saving game to "%s.mat" (folder path: %s)\n',savename,dir_results);
            
        case 3 % play again (if warm-up) or take a break (if main experiment)
            if ~iswarmup
                isbreak = 1;
            end
        case 4 % play pure tone
            str_stim = [];
            data_passation_tmp = data_passation;
            idx = find(cfg_game.n_response_correct_target == 1); % looks for all '1's
            data_passation_tmp.n_stim(i_current) = idx(round( (length(idx)-1)*random('unif',0,1) )+1); % picks up one randomly
                        
            exp2eval = sprintf('str_stim =  %s_user(cfg_game,data_passation_tmp);',experiment);
            eval(exp2eval);
            stim_normal = str_stim.stim_tone_alone;
            
            player = audioplayer([sil4playing; stim_normal],cfg_game.fs);
            playblocking(player)
            
            fprintf(['\n    Press any key\n']);
            pause;

        case 5 % play modulated tone
            str_stim = [];
            data_passation_tmp = data_passation;
            idx = find(cfg_game.n_response_correct_target == 2); % looks for all '2's
            data_passation_tmp.n_stim(i_current) = idx(round( (length(idx)-1)*random('unif',0,1) )+1); % picks up one randomly
            exp2eval = sprintf('str_stim =  %s_user(cfg_game,data_passation_tmp);',experiment);
            eval(exp2eval);
            stim_normal = str_stim.stim_tone_alone;
            
            player = audioplayer([sil4playing; stim_normal],cfg_game.fs);
            playblocking(player)
            
            fprintf(['\n    Press any key\n']);
            pause;

        case 6 % escape training
            iswarmup = 0;
            clc
            
            cfg_game.warmup = iswarmup;
            
            str_inout = [];
            str_inout.debut_i = debut_i;
            
            str_inout = staircase_init(str_inout,cfg_game); % actual initialisation
            
            response   = str_inout.response;
            n_correctinarow = str_inout.n_correctinarow;
            expvar     = str_inout.expvar;
            i_current  = str_inout.i_current;
            stepsize   = str_inout.stepsize;
            isbreak    = str_inout.isbreak;
            
            data_passation_init.data_training = data_passation;
            data_passation = data_passation_init; % 'empty' data_passation to start the main experiment
            %%%% TODO %%%%
            % display instructions main exp

        case {1,2} % responded 1 or 2
            
            data_passation.n_response_correct_target(i_current) = cfg_game.n_response_correct_target_sorted(n_stim);
            resp_num = data_passation.n_response_correct_target(i_current);
            iscorrect = (response == resp_num);
            
            % save trial data
            if ~iswarmup
                data_passation.n_responses(i_current) = response;
                data_passation.n_targets(i_current)   = cfg_game.n_targets_sorted(n_stim);
                data_passation.n_response_correct_target(i_current) = resp_num;
                data_passation.is_correct(i_current) = iscorrect;
            end
            if iswarmup || cfg_game.displayN || cfg_game.feedback
                % ListStim(n_stim).response 
                switch iscorrect
                    case 1
                        txt_extra = 'You were right';
                    case 0
                        txt_extra = 'You were wrong';
                end
                
                if isfield(cfg_game,'response_names')
                    resp_name = cfg_game.response_names{cfg_game.n_response_correct_target_sorted(n_stim)};
                else
                    error('Continue validating here...')
                    resp_name = num2str(cfg_game.n_response_correct_target(n_stim));
                end
                % feedback
                fprintf('\n\t%s => Correct answer was : %.0f (%s)\n\n Press any key to continue.\n',txt_extra,resp_num,resp_name);
                pause(.5); % .5 s pause
            end
            
            if iscorrect
                n_correctinarow = n_correctinarow+1;
            else
                n_correctinarow = 0;
            end
             
            if cfg_game.adapt
                str_inout = [];
                str_inout.iscorrect = iscorrect;
                str_inout.expvar    = expvar;
                str_inout.stepsize  = stepsize;
                str_inout.n_correctinarow = n_correctinarow;
                
                str_inout = staircase_update(str_inout,cfg_game);
                
                % load updated parameters
                expvar    = str_inout.expvar;
                stepsize  = str_inout.stepsize;
                n_correctinarow = str_inout.n_correctinarow;
            end
             
            i_current=i_current+1;
        otherwise
            warning('%s: Keyboard response not recognised',upper(mfilename))
    end
    
end
 
%% Save game
clock_str = Get_date_and_time_str;
data_passation.date_end{length(data_passation.date_start)} = clock_str;
savename = il_get_savename(experiment, Subject_ID, Condition, clock_str);
save([dir_results savename],'cfg_game', 'data_passation');
msg_close

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

