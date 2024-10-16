function [cfg_game, data_passation, outs_trial, sim_work] = fastACI_trial_current(cfg_game, data_passation, expvar, ins_trial, flags, keyvals)
% function [cfg_game, data_passation, outs_trial, sim_work] = fastACI_trial_current(cfg_game, data_passation, expvar, ins_trial, flags, keyvals)
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
sim_work = [];
outs_trial = [];
if nargin < 6
    keyvals = [];
end
if nargin < 5
    flags = Get_idle_flag;
end
if isfield(data_passation,'i_current')
    i_current = data_passation.i_current;
else
    i_current = 1;
end
is_warmup = cfg_game.warmup;
isbreak = 0;

if cfg_game.adapt
    stepsize = ins_trial.stepsize;
    n_correctinarow = ins_trial.n_correctinarow;
end
if cfg_game.is_simulation
    bSimulation = 1;
    
    def_sim  = ins_trial.def_sim;
    sim_work = ins_trial.sim_work;
else
    bSimulation = 0;
end
bExperiment = ~bSimulation;

%%%% Trial start     
n_stim = cfg_game.stim_order(i_current);

if is_warmup
    data_passation.i_current = i_current;
    data_passation.n_stim(i_current) = n_stim;
    data_passation.expvar(i_current) = expvar;
end

% Pre-stores information of the current trial
if ~is_warmup
    % clock_now = clock; % Commented on 17/04/2022
    data_passation.i_current = i_current;
    data_passation.n_stim(i_current) = n_stim;
    data_passation.expvar(i_current) = expvar;
    % data_passation.date(i_current,:) = clock_now; % Commented on 17/04/2022
end

% if cfg_game.displayN == 1 || cfg_game.is_simulation
if cfg_game.is_simulation
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

tic

if bExperiment && isfield(cfg_game,'probe_periodicity') && cfg_game.probe_periodicity > 0
    %%% Defining whether this trial is a probe
    if mod(i_current,cfg_game.probe_periodicity) == 1
        switch cfg_game.Language
            case 'EN'
                clc
                fprintf('\n*** Probe sound ***\n');
            case 'FR'
                clc
                fprintf('\n\t*** Son de r\351f\351rence ***\n');
        end

        % generating a dummy data/cfg
        data_probe = data_passation;
        cfg_probe = cfg_game;
        data_probe.i_current = 1;
        data_probe.n_stim = 1;
        data_probe.expvar = -10;
        cfg_probe.n_targets_sorted = 1;

        str_stim = [];
        str2eval = sprintf('[str_stim,data_probe]=%s_user(cfg_probe,data_probe);',cfg_game.experiment);
        eval(str2eval);
        stim_normal = str_stim.tuser;

        sil4playing = zeros(0.1*cfg_game.fs,size(stim_normal,2));
        player = audioplayer(cfg_game.gain_play*[sil4playing; stim_normal],cfg_game.fs);
        N_samples_stim = size(stim_normal,1) + size(sil4playing,1);

        play(player)
        pause(N_samples_stim/cfg_game.fs);

        switch cfg_game.Language
            case 'EN'
                fprintf('\n    Press any key\n');
            case 'FR'
                fprintf('\n    Appuyez sur une touche\n');
        end
        pause;
        pause(0.5)
    end
end

str_stim = [];
str2eval = sprintf('[str_stim,data_passation]=%s_user(cfg_game,data_passation);',cfg_game.experiment);
eval(str2eval);
stim_normal = str_stim.tuser;

if bExperiment
    sil4playing = zeros(0.1*cfg_game.fs,size(stim_normal,2));
    if ~isfield(cfg_game, 'intervalnum') || cfg_game.intervalnum == 1
        % one single sound is played
        player = audioplayer(cfg_game.gain_play*[sil4playing; stim_normal],cfg_game.fs);
    elseif cfg_game.intervalnum == 2
        % two sounds in random order

        % first determine the random positions 
        randpos = randperm(cfg_game.intervalnum);
        if randpos(1) == 1
            stim_1 = str_stim.tuser;
            stim_2 = str_stim.tref;
        else
            stim_1 = str_stim.tref;
            stim_2 = str_stim.tuser;
        end
        data_passation.randpos(:,i_current) = randpos;

        % then plays the two sounds
        player = audioplayer(cfg_game.gain_play*[sil4playing; stim_1; sil4playing; stim_2],cfg_game.fs);
    else
        error('For the moment the toolbox can only handle experiments with two intervals or less...')
    end

    if cfg_game.intervalnum == 1
        N_samples_stim = size(stim_normal,1) + size(sil4playing,1);
    elseif cfg_game.intervalnum == 2
        N_samples_stim = size(stim_1,1) + size(sil4playing,1) + size(stim_2,1) + size(sil4playing,1);
    end
end
if bSimulation
    sil4playing = [];
end

%%% TODO : change the response names 'target present/absent' to
% 'first/second interval'

% Display message, play sound and ask for response

if is_warmup
    switch cfg_game.Language
        case 'EN'
            fprintf('\n    * WARM-UP PHASE *\n\n');
            if cfg_game.feedback == 1
                fprintf('\tDependent variable: %s = %.2f\n\n',cfg_game.expvar_description,expvar);
            end
        case 'FR'
            fprintf('\n    * PHASE D''\311CHAUFFEMENT *\n\n');
            if cfg_game.feedback == 1
                fprintf('\tNiveau %s = %.2f \n\n',cfg_game.expvar_description,expvar);
            end
    end
else
    N_for_next_stop = min(data_passation.next_session_stop,cfg_game.N_trials+1)-i_current-1;
    msg_currenttrial
end

if bExperiment

    play(player)
    pause(N_samples_stim/cfg_game.fs);
    time_before_response = toc;

    if cfg_game.intervalnum == 1 
        % simple case: each response is the name of an experiment
        if is_warmup
            switch cfg_game.Language
                case 'EN'
                    text2show = {'to play the stim again' ['to play a ' cfg_game.response_names{1}] ['to play a ' cfg_game.response_names{2}] 'to leave the warm-up phase'};
                case 'FR'
                    text2show = {'pour rejouer le son' ['pour \351couter un ' cfg_game.response_names{1}] ['pour \351couter un ' cfg_game.response_names{2}] 'pour quitter l''\351chauffement'};
            end
            response = Response_keyboard([cfg_game.response_names text2show], cfg_game);

        else
            switch cfg_game.Language
                case 'EN'
                    text2show = {'to take a break'};
                case 'FR'
                    text2show = {'pour prendre une pause'};
            end
            response = Response_keyboard([cfg_game.response_names text2show], cfg_game);
        end

    elseif cfg_game.intervalnum == 2
        % in forced choice the responses should be about the interval
        if is_warmup
            switch cfg_game.Language
                case 'EN'
                    FCresponses = {[cfg_game.response_names{1} ' first ' cfg_game.response_names{2} ' second'], [cfg_game.response_names{2} ' first and ' cfg_game.response_names{1} ' second']};
                    text2show = {'to play the stim again' ['to play a ' cfg_game.response_names{1}] ['to play a ' cfg_game.response_names{2}] 'to leave the warm-up phase'};
                case 'FR'
                    FCresponses = {[cfg_game.response_names{1} ' first ' cfg_game.response_names{2} ' second'], [cfg_game.response_names{2} ' first and ' cfg_game.response_names{1} ' second']};
                    text2show = {'pour rejouer le son' ['pour \351couter un ' cfg_game.response_names{1}] ['pour \351couter un ' cfg_game.response_names{2}] 'pour quitter l''\351chauffement'};
            end
            response = Response_keyboard([FCresponses text2show], cfg_game);

        else
            switch cfg_game.Language
                case 'EN'
                    FCresponses = {[cfg_game.response_names{1} ' first ' cfg_game.response_names{2} ' second'], [cfg_game.response_names{2} ' first and ' cfg_game.response_names{1} ' second']};
                    text2show = {'to take a break'};
                case 'FR'
                    FCresponses = {[cfg_game.response_names{1} ' en premier et ' cfg_game.response_names{2} ' en second'], [cfg_game.response_names{2} ' en premier et ' cfg_game.response_names{1} ' en second']};
                    text2show = {'pour prendre une pause'};
            end
            response = Response_keyboard([FCresponses text2show], cfg_game);
        end
        % and now the 
    end
    stop(player)
    
elseif bSimulation
    if cfg_game.intervalnum == 2
        error('Simulations for forced choice are not implemented yet');
        randpos = [1,2];
    end
    time_before_response = toc;
    
    switch def_sim.decision_script
        case 'aci_detect'
            % Default for 'dau1997', 'osses2021', 'osses2022a'
            [response,sim_work,def_sim] = aci_detect(cfg_game,data_passation,def_sim,sim_work, ...
                'argimport',flags,keyvals); 
            data_passation.decision_var_mue2choose(i_current,:) = sim_work.decision_var_mue2choose(i_current,:);

            disp('')
        case 'aci_detect_debug' 
            % This script is not in the fastACI repository but in the 
            %   private fastACI_sim
            [response,sim_work,def_sim] = aci_detect_debug(cfg_game,data_passation,def_sim,sim_work, ...
                'argimport',flags,keyvals); 
            data_passation.decision_var_mue2choose(i_current,:) = sim_work.decision_var_mue2choose(i_current,:);

        case 'king2019_detect'
            % Default for 'king2019'
            [response,sim_work,def_sim] = king2019_detect(cfg_game,data_passation,def_sim,sim_work);
    end
    if i_current == 1
        if isfield(sim_work,'thres_for_bias')
            data_passation.thres_for_bias = sim_work.thres_for_bias;
        end
        if isfield(sim_work,'type_decision')
            data_passation.type_decision = sim_work.type_decision;
        end
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
trial_end_time = toc;
data_passation.response_time(i_current) = trial_end_time-time_before_response;
data_passation.trial_time(i_current) = trial_end_time;

switch response
    case length(cfg_game.response_names)+1 % play again (if warm-up) or take a break (if main experiment)
        if is_warmup
            % It is the training session:
            outs_trial = ins_trial;
        else
            % 'Take a break':
            isbreak = 1;
            i_current = i_current-1; % to start with this same trial when resuming the experiment
            data_passation.i_current = i_current;
        end
    case length(cfg_game.response_names)+2 % play pure tone
        str_stim = [];
        data_passation_tmp = data_passation;
        idx = find(cfg_game.n_response_correct_target_sorted == 1); % looks for all '1's
        data_passation_tmp.n_stim(i_current) = idx(round( (length(idx)-1)*random('unif',0,1) )+1); % picks up one randomly

        exp2eval = sprintf('str_stim =  %s_user(cfg_game,data_passation_tmp);',cfg_game.experiment);
        eval(exp2eval);
        stim_normal = str_stim.stim_tone_alone;

        player = audioplayer(cfg_game.gain_play*[sil4playing; stim_normal],cfg_game.fs);
        playblocking(player)

        switch cfg_game.Language
            case 'EN'
                fprintf('\n    Press any key\n');
            case 'FR'
                fprintf('\n    Appuyez sur une touche\n');
        end
        pause;
        outs_trial = ins_trial;

    case length(cfg_game.response_names)+3 % play modulated tone
        str_stim = [];
        data_passation_tmp = data_passation;
        idx = find(cfg_game.n_response_correct_target_sorted == 2); % looks for all '2's
        data_passation_tmp.n_stim(i_current) = idx(round( (length(idx)-1)*random('unif',0,1) )+1); % picks up one randomly
        exp2eval = sprintf('str_stim =  %s_user(cfg_game,data_passation_tmp);',cfg_game.experiment);
        eval(exp2eval);
        stim_normal = str_stim.stim_tone_alone;

        player = audioplayer(cfg_game.gain_play*[sil4playing; stim_normal],cfg_game.fs);
        playblocking(player)

        switch cfg_game.Language
            case 'EN'
                fprintf('\n    Press any key\n');
            case 'FR'
                fprintf('\n    Appuyez sur une touche\n');
        end
        pause;
        outs_trial = ins_trial;

    case length(cfg_game.response_names)+4 % escape training
        is_warmup = 0;
        clc

        cfg_game.warmup = is_warmup;

        if cfg_game.adapt
            str_inout = [];
            str_inout.debut_i = ins_trial.debut_i;

            [str_inout,cfg_game] = staircase_init(str_inout,cfg_game); % actual initialisation

            response   = str_inout.response;
            outs_trial.n_correctinarow = str_inout.n_correctinarow;
            expvar     = str_inout.expvar;
            i_current  = str_inout.i_current;
            outs_trial.stepsize = str_inout.stepsize;
            isbreak    = str_inout.isbreak;

            ins_trial.data_passation_init.data_training = data_passation;
            data_passation = ins_trial.data_passation_init; % 'empty' data_passation to start the main experiment

            data_passation.reversal_current = str_inout.reversal_current;
        end
        
    case num2cell(1:length(cfg_game.response_names)) % gave a response

        if cfg_game.intervalnum == 1
            % this is the simple case
            resp_num = cfg_game.n_response_correct_target_sorted(n_stim);
            if length(cfg_game.response_names) == cfg_game.N_target % general case: there is one possible response per target
                iscorrect = (response == resp_num);
            else % special case, the number of possible responses is different from the number of targets
                iscorrect = (cfg_game.correctness_matrix(response) == resp_num);
                response_scale = response;
                response = cfg_game.correctness_matrix(response);
            end

            % save trial data
            if ~is_warmup
                data_passation.n_responses(i_current) = response;
                data_passation.n_targets(i_current)   = cfg_game.n_targets_sorted(n_stim);
                data_passation.n_response_correct_target(i_current) = resp_num;
                data_passation.is_correct(i_current) = iscorrect;
                if length(cfg_game.response_names) ~= cfg_game.N_target % special case, the number of possible responses is different from the number of targets
                    data_passation.n_response_scale(i_current) = response_scale;
                end
            end
        elseif cfg_game.intervalnum == 2
            %in the case of a forced choice the responses should be
            %recoded: the datapassation does not store the selected
            %interval but the selected identity of the stim. Index
            %i_current corresponds to the first interval, index
            %i_current+cfg_game.N/2 to the second interval
            resp_num = find(randpos==1);
            iscorrect = (response == resp_num);
            % save trial data for both intervals
            if ~is_warmup
                if response == 1
                    data_passation.n_responses(i_current) = randpos(1);
                    data_passation.n_responses(i_current+cfg_game.N/2) = randpos(2);
                elseif response == 2
                    data_passation.n_responses(i_current) = randpos(2);
                    data_passation.n_responses(i_current+cfg_game.N/2) = randpos(1);
                else
                    error('The toolbox can only handle 2AFC at this time, not 3AFC or more')
                end
                data_passation.n_targets(i_current)   = randpos(1);
                data_passation.n_targets(i_current+cfg_game.N/2)   = randpos(2);
                %%% TODO: at the moment, the combination of options
                %%% forced-choice + complex responses is not allowed
                data_passation.n_response_correct_target(i_current) = randpos(1);
                data_passation.n_response_correct_target(i_current+cfg_game.N/2) = randpos(2);
                data_passation.is_correct(i_current) = iscorrect;
                data_passation.is_correct(i_current+cfg_game.N/2) = iscorrect;
                data_passation.expvar(i_current+cfg_game.N/2) = expvar;
                data_passation.n_stim(i_current+cfg_game.N/2) = data_passation.n_stim(i_current)+cfg_game.N/2;

                if length(cfg_game.response_names) ~= cfg_game.N_target % special case, the number of possible responses is different from the number of targets
                    error('at the moment, the combination of options forced-choice + complex responses is not allowed')
                    data_passation.n_response_scale(i_current) = response_scale;
                    data_passation.n_response_scale(i_current+cfg_game.N/2) = response_scale;
                end
            end
        end

        if is_warmup || cfg_game.feedback % || cfg_game.displayN
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

            if cfg_game.intervalnum == 2

                switch cfg_game.Language
                    case 'EN'
                        resp_name  = [cfg_game.response_names{randpos(1)} ' first ' cfg_game.response_names{randpos(2)} ' second'];
                    case 'FR'
                        resp_name = [cfg_game.response_names{randpos(1)} ' en premier et ' cfg_game.response_names{randpos(2)} ' en second'];
                end
            else
                if isfield(cfg_game,'response_names')
                    resp_name = cfg_game.response_names{cfg_game.n_response_correct_target_sorted(n_stim)};
                else
                    error('Continue validating here...')
                    resp_name = num2str(cfg_game.n_response_correct_target(n_stim));
                end
            end
            % feedback
            switch cfg_game.Language
                case 'EN'
                    fprintf('\n ** %s => Correct answer was : %.0f (%s) ** \n\n',txt_extra,resp_num,resp_name);
                case 'FR'
                    fprintf('\n ** %s => La bonne reponse \351tait : %.0f (%s) **\n\n',txt_extra,resp_num,resp_name);
            end
            pause(1); 
            if ~cfg_game.is_simulation
                % Refreshing only for real participants...
                clc; % pause before refreshing
            end
        end

        if cfg_game.adapt

            if iscorrect
                n_correctinarow = n_correctinarow+1;
            else
                n_correctinarow = 0;
            end
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
            
            outs_trial.stepsize = stepsize;
            outs_trial.n_correctinarow = n_correctinarow;
        else
%             if isfield(cfg_game,'probe_periodicity') && cfg_game.probe_periodicity > 0
%                 %%% Defining whether the next trial is a probe
%                 if mod(i_current+1,cfg_game.probe_periodicity) == 1
%                     expvar = cfg_game.startvar + 10;
%                 else
%                     expvar = cfg_game.startvar;
%                 end
%             end
        end

        i_current = i_current+1; 
        
        if cfg_game.is_simulation
            if i_current == data_passation.next_session_stop
                if isfield(cfg_game.def_sim,'ir_reference')
                    cfg_game.def_sim = Remove_field(cfg_game.def_sim,'ir_reference'); % removing the individual internal representations for the template
                end
                if isfield(def_sim,'ir_signal')
                    cfg_game.def_sim = Remove_field(cfg_game.def_sim,'ir_signal'); % removing the individual internal representations for the template
                end
            end
        end
        
    otherwise
        warning('%s: Keyboard response not recognised',upper(mfilename))
end

outs_trial.response = response;
outs_trial.i_current = i_current; % to be passed
outs_trial.isbreak = isbreak;

outs_trial.expvar = expvar;