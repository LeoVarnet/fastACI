function [data_passation,cfg_game] = Read_ACI_data(file2load,version)
% function [data_passation,cfg_game] = Read_ACI_data(file2load,version)
%
% Variables with these format:
%   2013    
%   2015
%   2015.1 (use g20210226_Script)

% Ascii for e (as in module): 0232
% UTF   for e               = 351 (seems to be)
%
% fprintf('Química, Matemáticas, Español.\n')
% disp('Química, Matemáticas, Español.\n')
% fprintf(native2unicode('Química, Matemáticas, Español.','latin1'))
% encoding = slCharacterEncoding()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    % No version specified
    version = Check_ACI_data_type(file2load);
end

cfg_game = [];
data_passation = [];

var = load(file2load);

% foldername = 'NoiseStims';

switch version
    case 'current'
        data_passation = var.data_passation;
        cfg_game       = var.cfg_game;
        
    otherwise % Needs convertion
        if isfield(var,'i')
            cfg_game.N = var.i;
            var = Remove_field(var,'i');
        end

        switch version
            % case {2015,'2015'}
            %     cfg_game.N = length(var.correct_answer);

            case {2015.1,'2015.1',2021.1,'2021.1'}
                cfg_game.N = cfg_game.N-1;
        end

        if isfield(cfg_game,'N')
            tmp = nan(1,cfg_game.N);
            str = [];
            for i = 1:cfg_game.N
                str = [str sprintf('tmp(%.0f),',i)];
            end
        end

        switch version
            case {2015,'2015'}
                cfg_game.stim_order = var.stim_order;

                data_passation.expvar = var.SNR;

                % n_stim == 1 and 3 is 'da'
                % n_stim == 2 and 4 is 'ga'

                N = length(data_passation.expvar);
                cfg_game.N = N;
                cfg_game.N_target = length(unique(var.n_signal));
                cfg_game.N_presentation = round(cfg_game.N/cfg_game.N_target);

                if ~isfield(data_passation,'n_response')
                    target_response = nan(1,N); % memory allocation
                    target_response_wrong = nan(size(target_response)); % memory allocation
                    n_response      = nan(1,N);

                    idx = find(var.n_signal == 1 | var.n_signal == 3);
                    target_response(idx) = 1;
                    target_response_wrong(idx) = 2;
                    idx = find(var.n_signal == 2 | var.n_signal == 4);
                    target_response(idx) = 2;
                    target_response_wrong(idx) = 1;

                    idx = find(var.correct_answer==1);
                    n_response(idx) = target_response(idx);
                    idx = find(var.correct_answer==0);
                    n_response(idx) = target_response_wrong(idx);
                    data_passation.n_reponse  = n_response;
                else
                    data_passation.n_reponse  = data_pa.n_response;
                end

                data_passation.n_signal   = var.n_signal;

                if ~isfield(data_passation,'is_correct')
                    data_passation.is_correct = double(var.correct_answer);
                else
                    data_passation.is_correct = data_pa.is_correct;
                end

                fn = Get_filenames(cd,'*Bruit*');
                if isempty(fn)
                    fn = Get_filenames(cd,'*Noise*');
                end
                if isempty(fn)
                    error('No noise folder found')
                end
                cfg_game.FolderInBruit   = fn{1}; 

                fn = Get_filenames(cd,'*Signal*');
                if isempty(fn)
                    error('No signal folder found')
                end
                cfg_game.FolderInSignal  = fn{1};

                cfg_game.stim_order = var.stim_order;
                cfg_game.CorrectResponses = [1 2 1 2]; % var.correct_answer;
                cfg_game.NameResponse = {'Da','Ga'}; % manually put
                cfg_game.target_names = cfg_game.NameResponse;
                
                cfg_game.NameSignal = {'Alda','Alga','Arda','Arga'}; % manually put
                
                ListStim = dir(strcat(cfg_game.FolderInBruit, '/*.wav'));
                ListStim = rmfield(ListStim,{'date','datenum','bytes', 'isdir'});

                data_passation.ListStim = ListStim;

                % data_passation.is_correct = double(var.data_passation.is_correct);
                % var.data_passation = Remove_field(var.data_passation,'is_correct');
                % 
                % if isfield(var.ListStim,'n_response')
                %     exp2eval= ['[' str(1:end-1) ']=var.ListStim(:).n_response;'];
                %     eval(exp2eval);
                %     data_passation.n_response = tmp;
                % 
                %     var.ListStim = Remove_field(var.ListStim,'n_response');
                % end
                % 
                % if isfield(var.ListStim,'m')
                %     exp2eval= ['[' str(1:end-1) ']=var.ListStim(:).m;'];
                %     eval(exp2eval);
                %     data_passation.expvar = tmp;
                % 
                %     data_passation.expvar_description = 'modulation depth (dB)';
                % 
                %     var.ListStim = Remove_field(var.ListStim,'m');
                % end
                % 
                % if isfield(var.ListStim,'n_presentation')
                %     exp2eval= ['[' str(1:end-1) ']=var.ListStim(:).n_presentation;'];
                %     eval(exp2eval);
                %     data_passation.n_stim = tmp;
                % 
                %     var.ListStim = Remove_field(var.ListStim,'n_presentation');
                % end
                % 
                % data_passation.ListStim = var.ListStim;
                % 
                % disp('')
                cfg_game.N_response = length(cfg_game.target_names);

            case {2015.1,'2015.1',2021.1,'2021.1'}
                % error('Continue here...')
                cfg_game.stim_order = var.cfg_game.ordre_aleatoire;
                var.cfg_game = Remove_field(var.cfg_game,'ordre_aleatoire');

                data_passation.is_correct = double(var.data_passation.is_correct);
                var.data_passation = Remove_field(var.data_passation,'is_correct');

                if isfield(var.data_passation,'N_signal')
                    data_passation.n_targets = var.data_passation.N_signal;
                    var.data_passation = Remove_field(var.data_passation,'N_signal');
                end
                if isfield(var.cfg_game,'N_signal')
                    cfg_game.N_target = var.cfg_game.N_signal;
                    var.cfg_game = Remove_field(var.cfg_game,'N_signal');
                end
                
                if isfield(var.cfg_game,'N_noise')
                    cfg_game.N_presentation = var.cfg_game.N_noise;
                    var.cfg_game = Remove_field(var.cfg_game,'N_noise');
                else
                    cfg_game.N_presentation = round(cfg_game.N/cfg_game.N_target);
                end
                
                if isfield(var.cfg_game,'response_names')
                    cfg_game.target_names = var.cfg_game.response_names;
                    var.cfg_game = Remove_field(var.cfg_game,'response_names');
                    
                    cfg_game.N_response = length(cfg_game.target_names);
                end
                
                if isfield(var.ListStim,'n_response')
                    exp2eval= ['[' str(1:end-1) ']=var.ListStim(:).n_response;'];
                    eval(exp2eval);
                    data_passation.n_responses = tmp;

                    var.ListStim = Remove_field(var.ListStim,'n_response');
                end

                if isfield(var.ListStim,'m')
                    exp2eval= ['[' str(1:end-1) ']=var.ListStim(:).m;'];
                    eval(exp2eval);
                    data_passation.expvar = tmp;

                    cfg_game.expvar_description = 'modulation depth (dB)';

                    var.ListStim = Remove_field(var.ListStim,'m');
                end

                if isfield(var.ListStim,'n_presentation')
                    exp2eval= ['[' str(1:end-1) ']=var.ListStim(:).n_presentation;'];
                    eval(exp2eval);
                    data_passation.n_stim = tmp;

                    var.ListStim = Remove_field(var.ListStim,'n_presentation');
                end
                
                if isfield(var.data_passation,'responsetime')
                    data_passation.response_time = var.data_passation.responsetime;
                    var.data_passation = Remove_field(var.data_passation,'responsetime');
                end

                if isfield(var.data_passation,'resume_trial')
                    data_passation.resume_trial = var.data_passation.resume_trial;
                    var.data_passation = Remove_field(var.data_passation,'resume_trial');
                end
                
                if isfield(var.cfg_game,'fs')
                    cfg_game.fs = var.cfg_game.fs;
                    var.cfg_game = Remove_field(var.cfg_game,'fs');
                end
                if isfield(var.cfg_game,'m_start')
                    cfg_game.startvar = var.cfg_game.m_start;
                    var.cfg_game = Remove_field(var.cfg_game,'m_start');
                end
                if isfield(var.cfg_game,'min_stepsize')
                    cfg_game.min_stepsize = var.cfg_game.min_stepsize;
                    var.cfg_game = Remove_field(var.cfg_game,'min_stepsize');
                end
                if isfield(var.cfg_game,'adapt_stepsize')
                    cfg_game.adapt_stepsize = var.cfg_game.adapt_stepsize;
                    var.cfg_game = Remove_field(var.cfg_game,'adapt_stepsize');
                end
                
                if isfield(var.cfg_game,'simulation')
                    switch var.cfg_game.simulation
                        case {1,'yes','oui'}
                            cfg_game.is_simulation = 1;
                        case {0,'no','non'}    
                            cfg_game.is_simulation = 0;
                    end
                    cfg_game.is_experiment = ~cfg_game.is_simulation;
                    
                    var.cfg_game = Remove_field(var.cfg_game,'simulation');
                end
                
                cfg_game.Subject_ID = '';
                if isfield(var.cfg_game,'path')
                    try
                        tmp = strsplit(var.cfg_game.path,'/'); % tmp = strsplit(var.cfg_game.path,'\');
                        if length(tmp) == 1
                            tmp = strsplit(var.cfg_game.path,'\'); % assuming windows-formatted path
                        end
                        cfg_game.Subject_ID = tmp{end};
                    catch
                        warning('Subject_ID could not be extracted');
                    end
                end
                                
                cfg_game.ListStim = var.ListStim;
                var = Remove_field(var,'ListStim');

                disp('')

            case {2013,'2013'}
                if isfield(var,'ordre_aleatoire')
                    cfg_game.stim_order = var.ordre_aleatoire{1};
                    var = Remove_field(var,'ordre_aleatoire');
                    stim_order = cfg_game.stim_order;
                end

                if isfield(var,'ListStim')
                    ListStim = var.ListStim;

                    var = Remove_field(var,'ListStim');

                    if isfield(ListStim,'n_presentation')
                        tmp = horzcat(ListStim.n_presentation);
                        data_passation.n_presentation = tmp;
                        data_passation.n_presentation_description = 'trial in which soundXXX was presented';
                        data_passation.n_stim = tmp(stim_order);

                        ListStim = Remove_field(ListStim,'n_presentation');
                    end

                    if isfield(ListStim,'RSB')
                        % data_passation.expvar(1:cfg_game.N) = ListStim(:).RSB;
                        expvar = horzcat(ListStim.RSB);
                        data_passation.expvar = expvar(stim_order);

                        % exp2eval= ['[' str(1:end-1) ']=ListStim(:).RSB;'];
                        % eval(exp2eval);
                        % data_passation.expvar = tmp;

                        cfg_game.expvar_description = 'signal-to-noise ratio';

                        ListStim = Remove_field(ListStim,'RSB');
                    end
                    % expvar = data_passation.expvar(cfg_game.stim_order);

                    if isfield(ListStim,'n_reponse')
                        tmp = horzcat(ListStim.n_reponse);
                        data_passation.n_response = tmp(stim_order);
                        cfg_game.n_response       = tmp;
                        ListStim = Remove_field(ListStim,'n_reponse');
                    end

                    if isfield(ListStim,'n_signal')
                        tmp = horzcat(ListStim.n_signal);
                        data_passation.n_signal = tmp(stim_order);
                        cfg_game.n_signal       = tmp;

                        ListStim = Remove_field(ListStim,'n_signal');
                    end

                    if isfield(ListStim,'is_correct')
                        tmp = horzcat(ListStim.is_correct);
                        data_passation.is_correct = double(tmp(stim_order));

                        ListStim = Remove_field(ListStim,'is_correct');
                    end

                    cfg_game.n_signals(1:cfg_game.N/2)            = 1;
                    cfg_game.n_signals(cfg_game.N/2+1:cfg_game.N) = 2;
                    cfg_game.CorrectResponses = [1 2];
                    cfg_game.NameResponse = {'Aba','Ada'};
                    ListStim = Remove_field(ListStim,'signal');

                    cfg_game.n_response_correct_target       = cfg_game.n_signals;
                    data_passation.n_response_correct_target = cfg_game.n_signals(stim_order);

                    for i = 1:cfg_game.N
                        switch ListStim(i).reponse{1}
                            case cfg_game.NameResponse{1}
                                cfg_game.n_signals(i) = 1;

                            case cfg_game.NameResponse{2}
                                cfg_game.n_signals(i) = 2;
                        end
                    end
                    data_passation.n_signals = cfg_game.n_signals(cfg_game.stim_order);
                    ListStim = Remove_field(ListStim,'reponse');

                    if isfield(ListStim,'score')
                        % score_signal  = horzcat(ListStim.score.score_signal);
                        % score_general = horzcat(ListStim.score.score_general);
                        ListStim = Remove_field(ListStim,'score'); warning('Field score still has to be integrated...')
                    end

                    if ~isfield(data_passation,'is_correct')
                        % If not yet assigned, then 'is_correct' is reconstructed
                        data_passation.is_correct = double(data_passation.n_signals == data_passation.n_response_correct_target);
                        % cfg_game.is_correct = double(cfg_game.n_signals == cfg_game.n_response_correct_target);
                    end
                    data_passation.ListStim = ListStim;
                    disp('')
                end
                
                try
                    cfg_game.N_response       = length(cfg_game.NameResponse);
                catch me
                    warning('cfg_game.N_response not assigned')
                end
        end
end

try
    if isfield(cfg_game,'target_names')
        for i = 1:length(cfg_game.target_names)
            tmp = strsplit(cfg_game.target_names{i},'\'); 
            if length(tmp) > 1
                % if there is this separator, then there is an unrecognised format:
                
                disp('')
                switch tmp{end}
                    case '351' % this is an e'
                        % Mapping between 232 (ASCII - I believe) and 351 was done manually
                        str = [];
                        for j = 1:length(tmp)-1
                            str = [str tmp{j}];
                        end
                        str = [str native2unicode(232,'ISO-8859-1')];
                        cfg_game.target_names{i} = str; 
                        
                    otherwise
                        warning('Unrecognised character was not fixed')
                end
            end
        end
    end
end

disp('')