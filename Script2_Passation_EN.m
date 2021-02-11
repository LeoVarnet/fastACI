function Script2_Passation_EN(experiment, Subject_ID)
% function Script2_Passation_EN(experiment, Subject_ID)
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

% def.expvarunit = 'dB';	% unit of the tracking variable
% def.expvardescription = 'modulation depth';	%'modulation depth' or 'modulation frequency' %description of the tracking variable

if nargin < 2
    Subject_ID = input('Enter the Subject ID (e.g., ''S01''): ');
end
if nargin == 0
    % experiment = 'modulationACI'; 
    experiment = 'speechACI_varnet2015';     
end

bSimulation = 0;
bDebug      = 1;

% -------------------------------------------------------------------------
% 1. Loading set-up: 
%    1.1. Loads cfgcrea*.mat
[path,name,ext]=fileparts(which(mfilename)); % path will be the folder where this file is located...
dir_main = [path filesep];    %'C:\Users\Varnet L�o\Dropbox\Professionnel\Matlab\MyScripts\modulationACI\AM';
dir_results = [dir_main 'Interim_results' filesep];

filter2use = [Subject_ID '_' experiment];
stored_cfg = Get_filenames(dir_results,['cfgcrea*' filter2use '.mat']);
N_stored_cfg = length(stored_cfg);
if N_stored_cfg==1
    var      = load([dir_results stored_cfg{1}]);
    cfg_game = var.cfg_crea;
    % cfg_game.cfg_crea   = var.cfg_crea;
elseif N_stored_cfg > 1
    error('Multiple participants option: has not been validated yet (To do by AO)')
else
    error('%s: no cfg_crea available',upper(mfilename));
end

% TODO: integrate this into Script1
% TODO: change N_noise and N_signal by other more relevant names...
if ~isfield(cfg_game,'N')
    cfg_game.N = cfg_game.N_noise*cfg_game.N_signal;
end

if ~isfield(cfg_game,'resume')
    cfg_game.resume = []; % 'no' -> new game, 'yes' -> load last saved game, [] -> load last saved game if exists or start a new game
end

% -------------------------------------------------------------------------
%     1.2. Loading parameters: Looks for an existing 'game' (or previous 
%          session for the same participant)
if isempty(cfg_game.resume)
    ListSavegame = dir([dir_results 'savegame*' filter2use '.mat']);
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
                for j=1:length(ListSavegame)
                    if ListSavegame(j).bytes > bytes
                        index_savegame=j;
                        bytes = ListSavegame(j).bytes; % looks for the largest MAT file
                    end
                end
                load_name = [dir_results ListSavegame(index_savegame).name];
            end
            
            cfg_game = []; % it will be re-loaded now:
            ListStim = [];
            
            load(load_name);
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
        if ~isfield(cfg_game,'experiment')
            cfg_game.experiment = 'modulationACI';
        end
            
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
            error('Not validated yet...')
            % modelparameters;
            % cfg_game.fadein_s           = 0;
            % cfg_game.N_template         = 0;
            % cfg_game.warmup             = 'no';
            % cfg_game.sessionsN          = 500;
            % cfg_game.stim_dur           = 0.75;
        end

    % cfg_game.current_folder{1} = dir_main; % cd;
    cfg_game.script_name{1} = [mfilename('fullpath') '.m'];

    data_passation.resume_trial = 0;
    data_passation.date_start{1} = Get_date_and_time_str;

	if cfg_game.is_simulation == 1
        error('Not validated yet...')
        % % display welcome message
        % msg_welcome
	end        
end

% TODO: change this name by dir_noise or something like that...
dir_stim = cfg_game.dir_stim;
% clear temp bytes ListSavegame index_savegame clock_now
 
%% Load stims, create templates for simulation

if cfg_game.resume == 0
    if ~isfield(cfg_game,'folder_name')
        disp('    No folder_name is specified...')
        cfg_game.folder_name = '';
    end
    ListStim = dir(strcat([dir_stim cfg_game.folder_name filesep], '*.wav'));
    
    ListStim = rmfield(ListStim,{'date','datenum','bytes', 'isdir'});
    if cfg_game.N ~= length(ListStim)
        switch experiment
            case 'speechACI_varnet2015'
                disp('    No folder_name is specified... / temporal arrangement here...')
            otherwise
                error('Number of stimuli does not match.')
        end
    end
    
    if cfg_game.N_signal > 2
        if ~isfield(cfg_game,'response_correct_target')
            error('A maximum of two alternatives ''1'' and ''2'' can be asked to the participants. Your experiment seems to use %.0f sounds, but they should be mapped to only two responses. Specify the field ''response_correct_target''.',cfg_game.N_signal);
        end
        % cfg_pass.CorrectResponses = [1,2,1,2]; % reponse correcte pour chaque signal (signaux dans l'ordre alphabetique)
    elseif cfg_game.N_signal == 2
        cfg_game.response_correct_target = [1 2];
    end
    
    list_signals = [];
    list_target_signals = [];
    for i=1:cfg_game.N_signal
        % for N_signal == 2 ==> first half equal to one, second half equal to two
        N_conditions = cfg_game.N_noise; % TODO: this is still not the most suitable naming...
        list_signals        = [list_signals i*ones(1,ceil(N_conditions))];%
        list_target_signals = [list_target_signals cfg_game.response_correct_target(i)*ones(1,ceil(N_conditions))];
    end
    
    for i=1:cfg_game.N
        ListStim(i).n_signal = list_signals(i);
    end
    cfg_game.n_signals = list_signals;
    cfg_game.n_response_correct_target = list_target_signals;
end

% Create template
 
if cfg_game.is_simulation == 1
    error('Not validated yet...')
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
 
%cfg_game.N = length(ListStim);
% clear liste_signaux ordre_alea i j IRind1 IRind2
 
%% Experiment
if cfg_game.resume == 0
    if cfg_game.randorder == 1
        cfg_game.randorder_idxs = randperm(cfg_game.N); 
    else
        cfg_game.randorder_idxs = 1:cfg_game.N; 
    end
    debut_i=1;
else
    debut_i=i_savegame;
end

cfg_game.bDebug = bDebug;

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
if ~isfield(cfg_game,'ListStim')
    cfg_game.ListStim = ListStim;
end

while i_current <= N && (cfg_game.is_simulation == 1 || i_current~=debut_i+cfg_game.sessionsN) && isbreak == 0
    
    % cfg_game.i_current = i_current;
    n_stim = cfg_game.randorder_idxs(i_current);
    
    if iswarmup
        
        % i_current = 1;
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
        data_passation.n_signal(i_current) = ListStim(n_stim).n_signal;
        data_passation.date(i_current,:) = clock_now;
    end
    
    if bDebug == 1
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        if isfield(cfg_game,'expvar_description')
            expvar_description = [', ' cfg_game.expvar_description];
        else
            expvar_description = '';
        end
            
        fprintf('\nDependent variable: expvar = %.4f%s \n',expvar,expvar_description);
    end
    
    str_stim = [];
    str2eval = sprintf('str_stim=%s_user(cfg_game,data_passation);',cfg_game.experiment);
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
     
    if cfg_game.is_experiment
        % display
        if iswarmup
            fprintf('\n    * WARM-UP PHASE *\n\n');
        else
            fprintf('\n    * MAIN EXPERIMENT *\n\n');
            if cfg_game.displayN == 1
                fprintf('    Playing stimulus # %.0f of %.0f\n',i_current,cfg_game.N);
            else
                fprintf('    Playing stimulus\n');
            end
        end
        % response = Reponse_clavier([cfg_game.response_names {'r\351ecouter le stim'}], 3.14);
        play(player)
        if iswarmup
            response = Reponse_clavier([cfg_game.response_names {'to play the stim again' ['to play a ' cfg_game.response_names{1}] ['to play a ' cfg_game.response_names{2}] 'to leave the warm-up phase'}]);
        else
            response = Reponse_clavier([cfg_game.response_names {'to take a break'}], 3.14);
        end
        stop(player)
    elseif cfg_game.is_simulation
        error('Not validated yet...')
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
     
    % responsetime = toc;
    % ListStim(n_stim).responsetime  = responsetime; 
    data_passation.responsetime(i_current) = toc;
     
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
            % data_passation_tmp.expvar = expvar;
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
            % iscorrect = (response == cfg_game.ListStim(n_stim).n_signal);
            iscorrect = (response == cfg_game.n_response_correct_target(n_stim));
            
            % save trial data
            if ~iswarmup
                % ListStim(n_stim).n_response = response; 
                data_passation.n_response(i_current) = response;
                % ListStim(n_stim).response = cfg_game.response_names{ListStim(n_stim).n_response};
                % ListStim(n_stim).is_correct = iscorrect; 
                data_passation.is_correct(i_current) = iscorrect;
            end
            if iswarmup || bDebug
                % ListStim(n_stim).response 
                switch iscorrect
                    case 1
                        txt_extra = 'You were right';
                    case 0
                        txt_extra = 'You were wrong';
                end
                
                resp_num = cfg_game.n_response_correct_target(n_stim);
                if isfield(cfg_game,'response_names')
                    resp_name = cfg_game.response_names{cfg_game.n_response_correct_target(n_stim)};
                else
                    error('Continue validating here...')
                    resp_name = num2str(cfg_game.n_response_correct_target(n_stim));
                end
                % feedback
                fprintf('\n %s => Correct answer was : %.1f (%s)\n\n Press any key to continue.\n',txt_extra,resp_num,resp_name);
                  % ListStim(n_stim).n_signal ,cfg_game.response_names{ListStim(n_stim).n_signal}
                pause;
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
savename = il_get_savename(experiment, Subject_ID, clock_str);
save([dir_results savename],'ListStim', 'cfg_game', 'data_passation');
msg_close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fname = il_get_savename(experiment, Subject_ID, clock_str)

if nargin < 3
    clock_str = Get_date_and_time_str;
end

savename = ['savegame_' clock_str];
fname = [savename '_' Subject_ID '_' experiment];