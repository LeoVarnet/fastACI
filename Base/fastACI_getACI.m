function [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(savegame_file,varargin)
% function [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(savegame_file,varargin)
%
% 1. Description (FR):
%       Se placer dans le dossier contenant le dossier contenant les donnees 
%       du participant a analyser
%
%    For shortening the calculation process in this script, the following 
%       simplifications can be done (please run with NO simplification when
%       preparing data for publication):
%         - IdxTrialsLoad set to a lower range (e.g. IdxTrialsLoad = 1:1000)
%         - cfg_ACI.lambda0   = 85; % Valeur initiale de lambda
%
% To set:
%       DimCI: ('tf' or 'lyon') - Choice of auditory model: Lyon or something else (see Varnet2015: 'Cochleograms')
%       opts_ACI.glmfct (default: glm)
%
% 1. Reading data from a MAT file and check compatibility (il_convert_tb_ACI_data)
% 2. Reading/setting options for calculation
% 3. Reading the sound waveforms, getting T-F representations (data_load)
% 4. Preprocessing before the ACI assessment
%
% Scripts from where 'Script4_Calcul_ACI_debug' is called from:
%   1. g20210301_recreating_varnet2013.m
%   2. g20210413_SAO5000_varnet2013.m
%
% Changing the parameter names:
% New name        Old name             Changed on:
% -- (removed)    withX                7/05/2021
%                 CI_glmqpoptim_fct
% TF_type         DimCI
%
% Old name: Script4_Calcul_ACI.m (changed on 7 July 2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Condition can be renamed to NameCond

if nargin == 0
    error('%s: Please spefify the identifier of the subject from whom you want to process the data',upper(mfilename));
end

Data_matrix = [];

% From argument function:
definput.import={'fastACI_getACI'}; % arg_fastACI_getACI.m
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

%% 1. Reading the experimental data (*.mat file):
[cfg_game, data_passation, ListStim] = Convert_ACI_data_type(savegame_file,keyvals);
N = cfg_game.N;

if isempty(keyvals.dir_out)
    warning('No output directory (opts_ACI.dir_out) has been specified, the same folder where the MAT file is will be used...');
    % curr_dir = [pwd filesep]; % current directory
    [path,name,ext]=fileparts(which(savegame_file));
    if isempty(path)
        [path,name,ext]=fileparts( savegame_file );
    end
    path = [path filesep];
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
% opts_ACI = Ensure_field(opts_ACI,'Analysis_condition','total'); % tableau de cellules contenant les noms de la ou des conditions à calculer
%Analysis_condition = ''; % '_total'; warning('Remove soon...')

% do_permutation = flags.do_permutation; % By default the permutation test is 'on'
do_recreate_validation = flags.do_recreate_validation;
TF_type = flags.TF_type;
glmfct  = flags.glmfct;
% END: From argument function:

if isempty(keyvals.idx_trialselect)
    keyvals.idx_trialselect = 1:N;
    str_last_trial = '';
else
    str_last_trial = ['-up-to-trial-' num2str(length(keyvals.idx_trialselect))];
end

Condition = '';
if isfield(cfg_game,'Condition')
    Condition = ['-' cfg_game.Condition];
end

if isempty(cfg_game.Subject_ID)
    % Trying to find a Subject ID if empty
    di = cfg_game.dir_noise;
    if strcmp(di(end),filesep)
        di = di(1:end-1);
    end
    di = strsplit(di,filesep);
    Subject_ID = di{end-1};
    
    cfg_game.Subject_ID = Subject_ID;
end

%%%
trialtype_analysis = [];
if ~isempty(keyvals.trialtype_analysis)
    switch keyvals.trialtype_analysis
        case 'total'
            %%% Nothing to do: just an empty name
        otherwise
            trialtype_analysis = ['-' keyvals.trialtype_analysis];
    end
end

if flags.do_no_bias
    % Number of trials for targets 1 or 2 will be 'equalised'. This 
    %     processing is introduced in the _preprocessing script.
    trialtype_analysis = [trialtype_analysis '+nobias'];
end

if keyvals.add_signal
    label_add_signal = '+addsignal';
else
    label_add_signal = '';
end

if ~isempty(keyvals.expvar_limits)
    label_expvar_limits = sprintf('-expvar-%.0f-to-%.0f',keyvals.expvar_limits);
else
    label_expvar_limits = [];
end
fnameACI = [dir_out cfg_game.Subject_ID '_' cfg_game.experiment Condition '-ACI' trialtype_analysis '-' TF_type '-' glmfct str_last_trial label_add_signal label_expvar_limits '.mat'];

%%%

bCalculation = ~exist(fnameACI,'file');

if bCalculation == 0
    if keyvals.skip_if_on_disk == 0
        disp('ACI found on disk:')
        bCalculation = input('  Enter 1 to re-calculate (overwrite) or 0 to read stored results:');
    end
    
    if bCalculation == 1
        fname_old = [fnameACI(1:end-4) '-old.mat'];
        if exist(fname_old,'file')
            error('%s: Trying to back up your old-exsiting results, but there is already a file named %s',upper(mfilename),fname_old);
        end
        movefile(fnameACI,fname_old);
        fprintf('%s: Old mat file successfully backed up as: %s\n',upper(mfilename),fname_old);
    end
end

%%% 2.2 Creating structure cfg_ACI: ---------------------------------------
% cfg_game = Ensure_field(cfg_game,'response_names',{'Signal#1','Signal#2'});
cfg_ACI = import_cfg(cfg_game, 'dir_noise', 'dir_target', 'N', 'N_target', ... % 'dir_target', 'N_response'
    'stim_order', 'target_names', 'response_correct_target','response_names');

if isfield(cfg_game,'Subject_ID')
    cfg_ACI.Subject_ID = cfg_game.Subject_ID;
end
if isfield(cfg_game,'experiment')
    cfg_ACI.experiment = cfg_game.experiment;
end
if isfield(cfg_game,'Condition')
    cfg_ACI.Condition = cfg_game.Condition;
end
%

cfg_ACI = arg_TF_type(cfg_ACI, flags, keyvals);
cfg_ACI = arg_glmfct(cfg_ACI, flags);

cfg_ACI.flags = flags;
cfg_ACI.keyvals = keyvals;
cfg_ACI.fnameACI = fnameACI;

cfg_ACI.idx_trialselect   = cfg_ACI.keyvals.idx_trialselect; % numeros des essais utilises pour le calcul (defaut = 1:cfg_ACI.N), si possible les essais sont pris dans l'ordre de presentation
cfg_ACI.withU             = 1; % 'yes'; % Ajouter deux paramètres U au modèle
 
switch cfg_ACI.glmfct
    case 'glmfitqp'
        check_cfg(cfg_ACI, 'prior','lambda0', 'stepsize', 'maxiter', 'nobreak', 'minDiffSecondRound');
        cfg_ACI.N_folds   = cfg_ACI.keyvals.N_folds;
        
    case {'lassoglm','lasso','lassoslow','lassoglmslow'}
        cfg_ACI.lambda0   = [];
        cfg_ACI.N_folds    = cfg_ACI.keyvals.N_folds; 
        
    case 'classic_revcorr'
        
end
if cfg_ACI.zscore == 0
    error('%s: the glmfct options require that cfg_ACI.zscore is 1',upper(mfilename))
end
 
cfg_ACI.f_limits = cfg_ACI.keyvals.f_limits; % bande de frequences pour analyse
cfg_ACI.t_limits = cfg_ACI.keyvals.t_limits; % bande de temps pour analyse
 
if flags.do_lyon
    check_cfg(cfg_ACI, 'freq_analysis', 'time_analysis', 'decimation', 'earQ', 'stepfactor');
end

cfg_ACI = set_default_cfg(cfg_ACI, 'N_trialselect', length(cfg_ACI.idx_trialselect));

%% 3. Loading the data: Reading the sound waveforms

% ADD HERE: RECONSTRUCTION OF NOISE WAVEFORMS FOR SEEDS EXPERIMENTS
if bCalculation || do_recreate_validation
    
    if ~exist(cfg_ACI.dir_noise,'dir') && isempty(cfg_ACI.keyvals.dir_noise)
        % If it does not exist
        if isfield(cfg_game,'seeds_order')
            % Nothing to do, the waveforms will be generated inside data_load
            cfg_ACI.cfg_game = cfg_game;
            error('Under development...')            
        else
            error('%s: No valid noise directory (''dir_noise''). cfg_game contains a folder that was not found, please enter a valid dir_noise.',upper(mfilename))
        end
        
        if isfield(cfg_game,'Rove_level')
            cfg_ACI.Rove_level = cfg_game.Rove_level;
        end
        
    elseif ~isempty(cfg_ACI.keyvals.dir_noise)
        fprintf('%s: Using dir_noise specified as input parameter by the user\n',upper(mfilename));
        fprintf('\tNew cfg_game.dir_noise=%s\n',cfg_ACI.keyvals.dir_noise);
        if isfield(cfg_ACI,'dir_noise')
            fprintf('\t The directory used during the experiments has been stored as dir_noise_original=%s\n',cfg_ACI.dir_noise);
            cfg_ACI.dir_noise_original = cfg_ACI.dir_noise;
        end
        cfg_ACI.dir_noise = cfg_ACI.keyvals.dir_noise;
        
        if ~isempty(cfg_ACI.keyvals.dir_target)
            cfg_ACI.dir_target_original = cfg_ACI.dir_target;
            cfg_ACI.dir_target = cfg_ACI.keyvals.dir_target;
        end
    end
    
    if ~strcmp(fileparts(fileparts(cfg_ACI.dir_noise)), fileparts(fileparts(cfg_ACI.dir_target))) 
        % Extra check: if someone enters this part of the code, maybe he/she
        % does not have compatible dir_noise and dir_target directories, and 
        % therefore we throw a warning that appears for 10 s.
        warning('dir_noise and dir_target were found to be located in different root folders, please check that this is correct (ignore this message otherwise)');
        pause(10)
    end
    
    if data_passation.i_current ~= cfg_game.N
        fprintf('%s: Less trials have been tested by this participant than the expected cfg_game.N=%.0f trials\n',upper(mfilename),cfg_game.N);
        fprintf('\tPress ctrl+C to cancel the current ACI calculation, otherwise, the ACI will be obtained for less trials...\n');
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        pause(10);
        
        fprintf('\tToo late: ACI will be assessed using %.0f trials only\n',data_passation.i_current);
        
        % N = 150; warning('temporal')% data_passation.i_current;
        N = data_passation.i_current;
        cfg_ACI.N = N;
        cfg_ACI.stim_order = cfg_ACI.stim_order(1:N);
        cfg_ACI.idx_trialselect = 1:N;
        cfg_ACI.N_trialselect = N;
    end
    
    if isempty(keyvals.Data_matrix)
        % Loading the data regularly:
        [Data_matrix,cfg_ACI] = fastACI_getACI_dataload(cfg_ACI, ListStim, cfg_game);
    else
        Data_matrix = keyvals.Data_matrix;
        
        N_ref = 10;
        N = cfg_ACI.N;
        cfg_ACI_ref = cfg_ACI;
        cfg_ACI_ref.N = N_ref;
        [Data_matrix_ref,cfg_ACI] = fastACI_getACI_dataload(cfg_ACI_ref, ListStim, cfg_game);
        cfg_ACI.N = N; % restoring the initial N
        
        for ii = 1:N_ref
            diffe(ii) = sum(Data_matrix(ii,:)-Data_matrix_ref(ii,:));
        end
        if sum(diffe) ~= 0
            error('The input Data_matrix seems to be different from the expected Data_matrix loaded from ListStim');
        end
        
    end
end
 
%% 4. Preprocessing of the data, before the ACI calculation
if bCalculation || do_recreate_validation        
    
    [y, y_correct, X, U, cfg_ACI] = fastACI_getACI_preprocess(cfg_ACI, data_passation, Data_matrix);
    
end

%% 5. Calculation of the ACI
if bCalculation
    
    [ACI, results, cfg_ACI] = fastACI_getACI_calculate(cfg_ACI, y, y_correct, X, U);
    
    if isfield(cfg_ACI.keyvals,'Data_matrix')
        cfg_ACI.keyvals = Remove_field(cfg_ACI.keyvals,'Data_matrix');
    end
    
    save(fnameACI, 'ACI', 'cfg_ACI', 'results')
    results.fnameACI = fnameACI;
    results.fnameACI_description = 'File name where the fastACI results were stored...';
    
else
    ACI = [];
    cfg_ACI = [];
    results = [];
    load(fnameACI,'ACI','cfg_ACI','results');
 
    % if isfield(results,'lambdas')
    %     % This data display is just to make the user aware what for data are being loaded
    %     [round(results.lambdas') round(10*results.cvgofs')/10]
    % end
end
 
results.fnameACI = fnameACI;
results.fnameACI_description = 'File name where the fastACI results were stored...';

%% 6. Recreate validation, if requested:
if do_recreate_validation
    %%% Need to apply normalisation to X (check that it is the same)
    
    if isfield(cfg_ACI,'folds')
        folds = cfg_ACI.folds;
    else
        warning('The exact folds used during the ACI assessment were not found, using new folds...')
        folds = getcvfolds(length(y),cfg_ACI.N_folds);
    end
    
    if isfield(results,'finalfit')
        w = results.finalfit.w;
    end
    
    if cfg_ACI.withU
        Xin = [X U];
    else
        Xin = X;
    end
    
    if isfield(results,'finalfit')
        eval_opts = results.evaluation.opts;
        
        for i = 1:size(folds,2)
            idxs = find(folds(:,i) == 1);
            [ll(i),~,~,~,lleach,y_est] = evalGlmLikelihood(y(idxs),Xin(idxs,:),w,eval_opts.baseline(idxs),eval_opts.family,eval_opts.familyextra,eval_opts.weights(idxs));
        end
        ll = ll/length(idxs);

        [ll_tot,~,~,~,lleach,y_est] = evalGlmLikelihood(y,Xin,w,eval_opts.baseline,eval_opts.family,eval_opts.familyextra,eval_opts.weights);
        ll_tot = ll_tot/length(y);

        figure; 
        plot(y_est,'bs-'); hold on, grid on
        plot(y(idxs),'ro')
    else
        warning('Skipping the re-validation (bRecreate_validation) because ''results.finalfit'' was not found...')
    end
    
end

%% 7. Preparing final output and plotting:
 
% affichage de l'ACI

Max_here = max(max(ACI));
Min_here = min(min(ACI));
ACI_norm = 2*(ACI-Min_here)/(Max_here-Min_here)-1;

if flags.do_plot || nargout == 0
    if isfield(results,'ACI_perm')
        disp('Plotting permuation test results...')
        %%% Plotting permutation test:
        idxs         = find(ACI(:)>results.ACI_perm_CI_high(:) |  ACI(:)< results.ACI_perm_CI_low(:));
        % idxs_exclude = find(ACI(:)<=results.ACI_perm_CI_high(:) & ACI(:)>=results.ACI_perm_CI_low(:));

        ACI_ex = zeros(size(ACI));
        ACI_ex(idxs) = ACI(idxs);
        figure;
        subplot(1,2,1)
        % affichage_tf(results.ACI_norm, 'CI', 'cfg',cfg_ACI); hold on
        affichage_tf(ACI, 'CI', 'cfg',cfg_ACI);
        title(glmfct)

        subplot(1,2,2)
        % affichage_tf(results.ACI_norm, 'CI', 'cfg',cfg_ACI); hold on
        affichage_tf(ACI_ex, 'CI', 'cfg',cfg_ACI);
        title([glmfct ' (only within CI)'])
        %%% End plotting permutation test
    else
        
        figure;
        affichage_tf(ACI, 'CI', 'cfg',cfg_ACI);
        title(glmfct)
        
        % figure
        % affichage_tf(ACI_norm, 'CI', 'cfg',cfg_ACI); affichage_tf(ACI, 'CInorm', 'cfg', cfg_ACI)
        % affichage_tf(ACI, 'tvalue', 'cfg', cfg_ACI); affichage_tf(ACI, 'prob', 'cfg', cfg_ACI)
    end
end

results.ACI_norm       = ACI_norm;
results.ACI_norm_description = 'normalised ACI image using the optimal lambda, with weights between -1 and 1';

results.cfg_game       = cfg_game;
results.data_passation = data_passation;
