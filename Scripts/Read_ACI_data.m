function [data_passation,cfg_game] = Read_ACI_data(file2load,version)
% function [data_passation,cfg_game] = Read_ACI_data(file2load,version)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    % No version specified
    version = Check_ACI_data_type(file2load);
end

cfg_game = [];
data_passation = [];

var = load(file2load);

% n_signal = data_passation.N_signal(1:length(n_response));
% foldername = 'NoiseStims';

if isfield(var,'i')
    cfg_game.N = var.i;
    var = Remove_field(var,'i');
end

switch version
    case {2015,'2015'}
        cfg_game.N = cfg_game.N-1;
end
   
tmp = nan(1,cfg_game.N);
str = [];
for i = 1:cfg_game.N
    str = [str sprintf('tmp(%.0f),',i)];
end

switch version
    case {2015,'2015'}
        cfg_game.stim_order = var.cfg_game.ordre_aleatoire;
        var.cfg_game = Remove_field(var.cfg_game,'ordre_aleatoire');
        
        data_passation.is_correct = double(var.data_passation.is_correct);
        var.data_passation = Remove_field(var.data_passation,'is_correct');
        
        if isfield(var.ListStim,'n_response')
            exp2eval= ['[' str(1:end-1) ']=var.ListStim(:).n_response;'];
            eval(exp2eval);
            data_passation.n_response = tmp;

            var.ListStim = Remove_field(var.ListStim,'n_response');
        end
            
        if isfield(var.ListStim,'m')
            exp2eval= ['[' str(1:end-1) ']=var.ListStim(:).m;'];
            eval(exp2eval);
            data_passation.expvar = tmp;

            data_passation.expvar_description = 'modulation depth (dB)';

            var.ListStim = Remove_field(var.ListStim,'m');
        end
        
        if isfield(var.ListStim,'n_presentation')
            exp2eval= ['[' str(1:end-1) ']=var.ListStim(:).n_presentation;'];
            eval(exp2eval);
            data_passation.n_stim = tmp;

            var.ListStim = Remove_field(var.ListStim,'n_presentation');
        end
            
        data_passation.ListStim = var.ListStim;
        
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
                
                data_passation.expvar_description = 'signal-to-noise ratio';
                
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
            cfg_game.ListStim = ListStim;
            disp('')
        end
end

disp('')