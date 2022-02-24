function [cfg_game, data_passation, outs_trial] = fastACI_trial_current(cfg_game, data_passation, expvar, ins_trial, keyvals)
% function [cfg_game, data_passation, outs_trial] = fastACI_trial_current(cfg_game, data_passation, expvar, ins_trial, keyvals)
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

outs_trial = [];
if nargin < 5
    keyvals = [];
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
    def_sim  = ins_trial.def_sim;
    sim_work = ins_trial.sim_work;
end

%%%% Trial start     
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

end

str_stim = [];
str2eval = sprintf('[str_stim,data_passation]=%s_user(cfg_game,data_passation);',cfg_game.experiment);
eval(str2eval);
stim_normal = str_stim.tuser;

% stim_normal = il_user(str_inout,cfg_game);
%%% Create signal: end

if cfg_game.is_experiment
    sil4playing = zeros(0.1*cfg_game.fs,1);
    player = audioplayer(cfg_game.gain_play*[sil4playing; stim_normal],cfg_game.fs);
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
            fprintf('\n    * PHASE D''\311CHAUFFEMENT *\n\n');
    end
else
    N_for_next_stop = min(data_passation.next_session_stop,cfg_game.N+1)-i_current-1;
    msg_currenttrial
end

if cfg_game.is_experiment

    play(player)
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
        response = Response_keyboard([cfg_game.response_names text2show], cfg_game, 3.14);
    end
    stop(player)
elseif cfg_game.is_simulation
    
    switch def_sim.decision_script
        case 'aci_detect' 
            % Default for 'dau1997' and 'osses2021'
            [response,sim_work] = aci_detect(cfg_game,data_passation,def_sim,sim_work,keyvals);
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
data_passation.response_time(i_current) = toc;

switch response
    case 3.14 % Pause: Moved out of fastACI_trial_current
        % Nothing to do
        
    case 3 % play again (if warm-up) or take a break (if main experiment)
        if is_warmup
            % It is the training session:
            outs_trial = ins_trial;
        else
            % 'Take a break':
            isbreak = 1;
            i_current = i_current-1; % to start with this same trial when resuming the experiment
            data_passation.i_current = i_current;
        end
    case 4 % play pure tone
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

    case 5 % play modulated tone
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

    case 6 % escape training
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
                    fprintf('\n ** %s => Correct answer was : %.0f (%s) ** \n\n',txt_extra,resp_num,resp_name);
                case 'FR'
                    fprintf('\n ** %s => La bonne reponse etait : %.0f (%s) **\n\n',txt_extra,resp_num,resp_name);
            end
            pause(.5); 
            if ~cfg_game.is_simulation
                % Refreshing only for real participants...
                clc; % .5 s pause before refreshing
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
        end

        i_current = i_current+1; 
    otherwise
        warning('%s: Keyboard response not recognised',upper(mfilename))
end

outs_trial.response = response;
outs_trial.i_current = i_current; % to be passed
outs_trial.isbreak = isbreak;

outs_trial.expvar = expvar;