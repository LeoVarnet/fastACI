function [cfg_pass, data_passation, ListStim] = Convert_ACI_data_type(file_savegame,opts,version)
% function [cfg_pass, data_passation, ListStim] = Convert_ACI_data_type(file_savegame,opts,version)
%
% cfg_pass.expvar            data_passation.RSB
%                            data_passation.n_reponse
%                            data_passation.n_signal
% cfg_pass.ordre_aleatoire   cfg_pass.stim_order
% cfg_pass.response_correct_target  cfg_pass.CorrectResponses
% cfg_pass.response_names      cfg_pass.response_names
% ListStim                   cfg_pass.ListStim

% cfg_pass.dir_noise
% cfg_pass.dir_target
% 
% cfg_pass.N_response = length(cfg_pass.response_names);
% cfg_pass.target_names = {'Alda','Alga','Arda','Arga'}; % manually put
% cfg_pass.N_signal = length(cfg_pass.target_names);
% cfg_pass.N = length(data_passation.n_signal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_pass = []; % will be filled in...
data_passation = [];
var = load(file_savegame);

if nargin < 3
    [version,version_description] = Check_ACI_data_type(file_savegame);
    if ~strcmp(version,'current')
        fprintf('%s: The file was automatically recognised as: ''%s''\n',upper(mfilename),version_description);
    end
end

if nargin < 2
    opts = [];
end
[data_pa,cfg_pa] = Read_ACI_data(file_savegame,version);

switch version
    case {'current',2015.1,'2015.1'}
        cfg_pass = cfg_pa;
        data_passation = data_pa;
        
        if isfield(cfg_pass,'dir_speech')
            cfg_pass.dir_target = cfg_pass.dir_speech;
            cfg_pass = Remove_field(cfg_pass,'dir_speech');
        end
        if isfield(opts,'dir_target') 
            if ~isempty(opts.dir_target)
                cfg_pass.dir_target = opts.dir_target; % Assigning the keyval if given in opts
            end
        end
        
        if isfield(opts,'dir_noise') 
            if ~isempty(opts.dir_noise)
                cfg_pass.dir_noise = opts.dir_noise; % Assigning the keyval if given in opts
            end
        end
        if isfield(cfg_pass,'dir_noise')
            % Nothing to do...
        else
            fn = Get_filenames(cd,'*Bruit*');
            if isempty(fn)
                fn = Get_filenames(cd,'*Noise*');
            end
            if ~isempty(fn)
                cfg_pass.dir_noise = fn{1}; 
            else
                cfg_pass.dir_noise = '';
                warning('dir_noise: No noise folder found')
            end
            
        end
        
        if isfield(cfg_pass,'ListStim')
            ListStim = cfg_pass.ListStim;
        elseif isfield(data_passation,'ListStim')
            ListStim = data_passation.ListStim;
            cfg_pass.ListStim = data_passation.ListStim;
        else
            ListStim = dir(strcat(cfg_pass.dir_noise, '/*.wav'));
            ListStim = rmfield(ListStim,{'date','datenum','bytes', 'isdir'});
            
            if isfield(ListStim,'folder')
                ListStim = rmfield(ListStim,'folder');
            end
            cfg_pass.ListStim = ListStim; 
        end
                
        if isfield(data_passation,'ListStim')
            data_passation = Remove_field(data_passation,'ListStim');
        end
        
        if length(cfg_pass.response_names) == cfg_pass.N_target
            cfg_pass.target_names = cfg_pa.response_names;
        else
            cfg_pass.target_names = cfg_pa.response_names;
            % TODO
            %error('Continue here...') % Automate this to be read from experiment scripts...
            % cfg_pass.NameTarget = {'Alda','Alga','Arda','Arga'}; % manually put
        end
        cfg_pass.N_response = length(cfg_pa.response_names);
        
        disp('')
    case {2021.1,'2021.1'}    
        cfg_pass = cfg_pa; 
        data_passation = data_pa;
        
        cfg_pass.experiment = 'modulationACI';
        cfg_pass.response_correct_target = [1,2];
        cfg_pass.response_names = cfg_pass.target_names;
        data_passation.i_current = length(data_passation.is_correct);
        cfg_pass.n_targets_sorted = sort(data_passation.n_targets); %ones(1,data_passation.i_current);
        
        if isfield(opts,'dir_noise')
            if ~isempty(opts.dir_noise)
                cfg_pass.dir_noise = opts.dir_noise; % Assigning the keyval if given in opts
            end
        end
        if isfield(cfg_pass,'dir_noise')
            % Nothing to do...
        else % if isempty(cfg_pass,'dir_noise') || ~isfield(cfg_pass,'dir_noise')
            dir_savegame = fileparts(file_savegame);
            if isempty(dir_savegame)
                dir_savegame = [cd filesep]; % assuming that the sound folders can be in the current folder
            else
                dir_savegame = [dir_savegame filesep];
            end
            
            fn = Get_filenames(dir_savegame,'*Bruit*');
            if isempty(fn)
                fn = Get_filenames(dir_savegame,'*Noise*');
            end
            if ~isempty(fn)
                cfg_pass.dir_noise = [dir_savegame fn{1} filesep]; 
            else
                cfg_pass.dir_noise = '';
                warning('dir_noise: No noise folder found')
            end            
        end
        
        %%%
        if isfield(opts,'dir_target')
            if ~isempty(opts.dir_target)
                cfg_pass.dir_target = opts.dir_target; % Assigning the keyval if given in opts
            end
        end
        if isfield(cfg_pass,'dir_target')
            % Nothing to do...
        else % if isempty(cfg_pass,'dir_target') || ~isfield(cfg_pass,'dir_target')
            dir_savegame = fileparts(file_savegame);
            if isempty(dir_savegame)
                dir_savegame = cd; % assuming that the sound folders can be in the current folder
            else
                dir_savegame = [dir_savegame filesep];
            end
            
            fn = Get_filenames(dir_savegame,'*Target*');
            if ~isempty(fn)
                cfg_pass.dir_target = [dir_savegame fn{1} filesep]; 
            else
                cfg_pass.dir_target = '';
                warning('dir_target: No target folder found')
            end            
        end
        
        if isfield(data_passation,'ListStim')
            ListStim = data_passation.ListStim;
            cfg_pass.ListStim = data_passation.ListStim;
            data_passation = Remove_field(data_passation,'ListStim');
        else
            ListStim = dir(strcat(cfg_pass.dir_noise, '/*.wav'));
            ListStim = rmfield(ListStim,{'date','datenum','bytes', 'isdir'});
            
            if isfield(ListStim,'folder')
                ListStim = rmfield(ListStim,'folder');
            end
            cfg_pass.ListStim = ListStim;
        end
                
        % if length(cfg_pass.response_names) == cfg_pass.N_target
        %     cfg_pass.target_names = cfg_pa.response_names;
        % else
        %     error('Continue here...') % Automate this to be read from experiment scripts...
        %     % cfg_pass.NameTarget = {'Alda','Alga','Arda','Arga'}; % manually put
        % end
        if isfield(cfg_pa,'target_names') && ~isfield(cfg_pass,'N_response')
            cfg_pass.N_response = length(cfg_pa.target_names);
        end
        
    case {2013,'2013'}
        data_passation.expvar = data_pa.expvar;
        
        if isfield(opts,'dir_noise')
            cfg_pass.dir_noise = opts.dir_noise;
        else
            fn = Get_filenames(cd,'*Bruit*');
            if isempty(fn)
                fn = Get_filenames(cd,'*Noise*');
            end
            if isempty(fn)
                error('dir_noise: No noise folder found')
            end
            cfg_pass.dir_noise   = fn{1}; 
        end
        
        if isfield(opts,'dir_signal')
            cfg_pass.dir_target = opts.dir_signal;
        else
            if isfield(opts,'dir_target')
                cfg_pass.dir_target = opts.dir_target;
            else
                fn = Get_filenames(cd,'*Signal*');
                if isempty(fn)
                    error('dir_target: No signal folder found')
                end
                cfg_pass.dir_target  = fn{1};
            end
        end
        
        stim_order = cfg_pa.stim_order;
        cfg_pass.n_response_correct_target_sorted = cfg_pa.n_response_correct_target;
        cfg_pass.stim_order  = stim_order;
        cfg_pass.response_correct_target = cfg_pa.CorrectResponses;
        cfg_pass.response_names   = cfg_pa.NameResponse;
        
        data_passation.n_response_correct_target = cfg_pa.n_response_correct_target(stim_order);
        data_passation.i_current = length(stim_order);
        
        N_target = length(cfg_pass.response_correct_target);
        
        if length(cfg_pass.response_names) == N_target
            cfg_pass.target_names = {'aba','ada'}; % cfg_pass.response_correct_target;
        else
            error('Continue here...') % Automate this to be read from experiment scripts...
            % cfg_pass.NameTarget = {'Alda','Alga','Arda','Arga'}; % manually put
        end
        
        if ~isfield(data_pa,'n_signal')
            data_passation.n_targets = cfg_pa.n_signals(stim_order);
        else
            data_passation.n_targets = data_pa.n_signal;
        end
        
        if isfield(data_pa,'n_responses')
            data_passation.n_responses  = data_pa.n_responses;
        elseif isfield(data_pa,'n_response')
            data_passation.n_responses  = data_pa.n_response;
        else
            error('Continue here')
        end
        if ~isfield(data_pa,'is_correct')
            error('Continue here')
        else
            data_passation.is_correct  = data_pa.is_correct;
        end
        
        if ~isfield(cfg_pa,'Subject_ID')
            disp('No Subject_ID found...')
            cfg_pass.Subject_ID  = '';
        else
            cfg_pass.Subject_ID  = cfg_pa.Subject_ID;
        end
        
        cfg_pass.N_target = N_target;
        cfg_pass.N = length(data_pa.expvar);
        
        disp('')
        
        ListStim = data_pa.ListStim;
        cfg_pass.ListStim = data_pa.ListStim;
        cfg_pass.N_response = length(cfg_pass.response_names);
        
    case {2015,'2015'}
        data_passation = data_pa;
        cfg_pass       = cfg_pa;
        
        if isfield(cfg_pass,'FolderInBruit')
            cfg_pass.dir_noise = cfg_pass.FolderInBruit;
            cfg_pass = rmfield(cfg_pass,'FolderInBruit');
        end
        if isfield(opts,'dir_noise')
            cfg_pass.dir_noise = opts.dir_noise; % Assigning the keyval if given in opts
        end
        
        if isfield(cfg_pass,'FolderInSignal')
            cfg_pass.dir_target = cfg_pass.FolderInSignal;
            cfg_pass = rmfield(cfg_pass,'FolderInSignal');
        end
        if isfield(opts,'dir_target')
            cfg_pass.dir_target = opts.dir_target; % Assigning the keyval if given in opts
        end
        
        if isfield(cfg_pass,'N_signal')
            cfg_pass.N_target = cfg_pass.N_signal;
            cfg_pass = rmfield(cfg_pass,'N_signal');
        end
        
        if isfield(cfg_pass,'CorrectResponses')
            cfg_pass.response_correct_target = cfg_pass.CorrectResponses;
            cfg_pass = rmfield(cfg_pass,'CorrectResponses');
        end
        
        if isfield(cfg_pass,'NameResponse')
            cfg_pass.response_names = cfg_pass.NameResponse;
            cfg_pass = rmfield(cfg_pass,'NameResponse');
        end
                
        if isfield(cfg_pass,'NameSignal')
            cfg_pass.target_names = cfg_pass.NameSignal;
            cfg_pass = rmfield(cfg_pass,'NameSignal');
        end
        
        %%%
        if isfield(data_passation,'n_reponse')
            data_passation.n_responses = data_passation.n_reponse;
            data_passation = rmfield(data_passation,'n_reponse');
        end
        
        if isfield(data_passation,'n_signal')
            data_passation.n_targets = data_passation.n_signal;
            data_passation = rmfield(data_passation,'n_signal');
        end
        
        % data_passation.n_signal   = var.n_signal;
        % 
        % if ~isfield(data_pa,'is_correct')
        %     data_passation.is_correct = double(var.correct_answer);
        % else
        %     data_passation.is_correct = data_pa.is_correct;
        % end
        % 
        % fn = Get_filenames(cd,'*Bruit*');
        % if isempty(fn)
        %     fn = Get_filenames(cd,'*Noise*');
        % end
        % if isempty(fn)
        %     error('No noise folder found')
        % end
        % cfg_pass.dir_noise   = fn{1}; 
        % 
        % fn = Get_filenames(cd,'*Signal*');
        % if isempty(fn)
        %     error('No signal folder found')
        % end
        % cfg_pass.dir_target  = fn{1};
        % 
        % cfg_pass.stim_order = var.stim_order;
        % cfg_pass.response_correct_target = [1 2 1 2]; % var.correct_answer;
        % cfg_pass.NameResponse = {'Da','Ga'}; % manually put
        % 
        % cfg_pass.NameTarget = {'Alda','Alga','Arda','Arga'}; % manually put
        % 
        % ListStim = dir(strcat(cfg_pass.dir_noise, '/*.wav'));
        % ListStim = rmfield(ListStim,{'date','datenum','bytes', 'isdir'});
        
        ListStim = data_passation.ListStim;
        cfg_pass.ListStim = data_passation.ListStim;
        data_passation = Remove_field(data_passation,'ListStim');
        
end

if isfield(cfg_pa,'experiment')
    cfg_pass.experiment = cfg_pa.experiment;
    switch cfg_pass.experiment
        case 'AM_revcorr'
            % Updating to the current name
            cfg_pass.experiment = 'modulationACI';
    end
end

if ~isfield(cfg_pass,'experiment_full')
    cfg_pass.experiment_full = cfg_pass.experiment;
end

if ~isfield(data_passation,'resume_trial')
    data_passation.resume_trial = 1;
end

%%%% New 28/03/2025 to make it compatible with older versions
if ~isfield(data_passation,'n_stim')
    data_passation.n_stim = cfg_pass.stim_order;
end