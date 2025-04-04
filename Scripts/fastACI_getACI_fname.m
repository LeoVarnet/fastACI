function [fnameACI, cfg_game, data_passation, ListStim, flags, keyvals, bCalculation] = fastACI_getACI_fname(savegame_file,varargin)
% function [fnameACI, cfg_game, data_passation, ListStim, flags, keyvals, bCalculation] = fastACI_getACI_fname(savegame_file,varargin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Condition can be renamed to NameCond

if nargin == 0
    error('%s: Please spefify the identifier of the subject from whom you want to process the data',upper(mfilename));
end

% From argument function:
definput.import={'fastACI_getACI'}; % arg_fastACI_getACI.m

idx = find(strcmp(varargin,'l1glm'));
if ~isempty(idx)
    definput.import{end+1} = 'lasso'; % requesting to load from arg_lasso.m (if needed)
end
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

%% 1. Reading the experimental data (*.mat file):
[cfg_game, data_passation, ListStim] = Convert_ACI_data_type(savegame_file,keyvals);
N = cfg_game.N;

% default parameters -- should already be there at this point
% this is only to ensure compatibility
if ~isfield(cfg_game, 'intervalnum')
    cfg_game.intervalnum = 1; % default experiment is a yes-no
end

if ~isfield(cfg_game, 'N_trials')
    cfg_game.N_trials = cfg_game.N/cfg_game.intervalnum;
end

if isempty(keyvals.dir_out)
    % curr_dir = [pwd filesep]; % current directory
    [path,name,ext]=fileparts(which(savegame_file));
    if isempty(path)
        [path,name,ext]=fileparts( savegame_file );
    end
    if ~isempty(path)
        path = [path filesep 'Results_ACI' filesep];
    else
        path = ['Results_ACI' filesep];
    end
    if ~exist(path,'dir')
        mkdir(path);
    end
    dir_out = path;
else
    dir_out = keyvals.dir_out;
end

%% 2. Reading or setting options for calculation:
% if isfield(opts_ACI,'IdxTrialsLoad')
%     error('%s: Please redefine the field ''IdxTrialsLoad'' (deprecated name) to ''idx_trialselect'' (new name)',upper(mfilename));
% end

%%% 2.1 Options for the ACI calculation: ----------------------------------

% General parameters
TF_type = flags.TF_type;
glmfct  = flags.glmfct;
% END: From argument function:

if isempty(keyvals.idx_trialselect) && (data_passation.i_current == cfg_game.N_trials)
    keyvals.idx_trialselect = 1:N;
    str_last_trial = '';
else
    if ~isempty(keyvals.idx_trialselect)
        N_here = length(keyvals.idx_trialselect);
    else
        N_here = data_passation.i_current;
    end
    if N_here == cfg_game.N
        str_last_trial = '';
    else
        str_last_trial = ['-' Get_abbreviation_for_filename('up-to-trial') num2str(N_here)];
    end
end

Condition = '';
if isfield(cfg_game,'Condition')
    Condition = ['-' cfg_game.Condition];
end

if isempty(cfg_game.Subject_ID)
    try
        % Trying to find a Subject ID if empty
        di = cfg_game.dir_noise;
        if strcmp(di(end),filesep)
            di = di(1:end-1);
        end
        di = strsplit(di,filesep);
        Subject_ID = di{end-1};

        cfg_game.Subject_ID = Subject_ID;
    catch
        warning('No Subject_ID was found...')
        cfg_game.Subject_ID = '';
    end
end

%%%
str_trialtype_analysis = [];
if ~isempty(keyvals.trialtype_analysis)
    switch keyvals.trialtype_analysis
        case 'total'
            %%% Nothing to do: just an empty name
        otherwise
            str_trialtype_analysis = ['-' keyvals.trialtype_analysis];
    end
end

if flags.do_no_bias
    % Number of trials for targets 1 or 2 will be 'equalised'. This 
    %     processing is introduced in the _preprocessing script.
    str_trialtype_analysis = [str_trialtype_analysis '-' Get_abbreviation_for_filename('no_bias')];
end

if keyvals.add_signal
    str_add_signal = ['-' Get_abbreviation_for_filename('addsignal')];
else
    str_add_signal = '';
end

if isempty(keyvals.perc)
    bMaybe_do_expvar_limits = 1;
    
elseif isnan(keyvals.perc(1)) || isnan(keyvals.perc(2)) 
    bMaybe_do_expvar_limits = 1;
    
else
    bMaybe_do_expvar_limits = 0;
    if isempty(keyvals.expvar_limits)
        % Nothing to do here...
    else
        error('keyvals.perc AND keyvals.expvar_limits are non empty arrays. Only one of the two options is allowed...' );
    end
end

if bMaybe_do_expvar_limits == 0
    keyvals.expvar_limits  = [prctile(data_passation.expvar,keyvals.perc(1)) prctile(data_passation.expvar,keyvals.perc(2))];
    str_out = Get_abbreviation_for_filename('perc');
    str_expvar_limits = sprintf('-%s%.0f_%.0f',str_out,keyvals.perc);
else
    if ~isempty(keyvals.expvar_limits)
        str_out = Get_abbreviation_for_filename('expvar');
        sign_str1 = num2str(abs(keyvals.expvar_limits(1)));
        if keyvals.expvar_limits(1) < 0
            sign_str1 = ['m' sign_str1];
        end
        sign_str2 = num2str(abs(keyvals.expvar_limits(2)));
        if keyvals.expvar_limits(2) < 0
            sign_str2 = ['m' sign_str2];
        end
        str_expvar_limits = sprintf('-%s_%s_%s',str_out,sign_str1,sign_str2);
    else
        str_expvar_limits = [];
    end
end

if keyvals.expvar_after_reversal ~= 0
    str_expvar_rev = sprintf('-rev%.0f',keyvals.expvar_after_reversal);
else
    str_expvar_rev = [];
end

try
    if contains(savegame_file,'swap-tar.mat')
        str_target_order = '-swap';
    else
        str_target_order = '';
    end
catch
    str_target_order = '';
end

str_TF_type = Get_abbreviation_for_filename(TF_type);
str_glmfct  = Get_abbreviation_for_filename(glmfct);  
if strcmp(keyvals.pyramid_script,'imgaussfilt')
    str_glmfct = [str_glmfct '+pyrga'];
end
    
fnameACI = [dir_out 'ACI-' cfg_game.Subject_ID '-' cfg_game.experiment Condition ...
            str_trialtype_analysis '-' str_TF_type '-' str_glmfct str_last_trial ...
            str_add_signal str_expvar_limits str_expvar_rev str_target_order '.mat'];

bCalculation = ~exist(fnameACI,'file');