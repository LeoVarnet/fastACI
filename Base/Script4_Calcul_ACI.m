function [ACI,cfg_ACI,results] = Script4_Calcul_ACI(savegame_file,varargin)
% function [ACI,cfg_ACI,results] = Script4_Calcul_ACI(savegame_file,varargin)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Condition can be renamed to NameCond

if nargin == 0
    error('%s: Please spefify the identifier of the subject from whom you want to process the data',upper(mfilename));
end

opts_ACI = []; warning('Temporal')
%% 1. Reading the experimental data (*.mat file):
[cfg_game, data_passation, ListStim] = Convert_ACI_data_type(savegame_file,opts_ACI);
N = cfg_game.N;

if ~isfield(opts_ACI,'dir_out')
    warning('No output directory (opts_ACI.dir_out) has been specified, the same folder where the MAT file is will be used...');
    % dir_subject = [dir_where Subject_ID filesep];
    % 
    % curr_dir = [pwd filesep]; % current directory
    [path,name,ext]=fileparts(savegame_file);
    path = [path filesep];
    dir_out = path;
else
    dir_out = opts.dir_out;
end

%% 2. Reading or setting options for calculation:
% if isfield(opts_ACI,'IdxTrialsLoad')
%     error('%s: Please redefine the field ''IdxTrialsLoad'' (deprecated name) to ''idx_trialselect'' (new name)',upper(mfilename));
% end

%%% 2.1 Options for the ACI calculation: ----------------------------------

% General parameters
opts_ACI = Ensure_field(opts_ACI,'Analysis_condition','total'); % tableau de cellules contenant les noms de la ou des conditions à calculer
Analysis_condition = opts_ACI.Analysis_condition;

% From argument function:
definput.import={'Script4_Calcul_ACI'}; % arg_Script4_Calcul_ACI
[flags,keyvals]  = ltfatarghelper([],definput,varargin);

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

fnameACI = [dir_out cfg_game.Subject_ID '_' cfg_game.experiment Condition '-ACI_' Analysis_condition '-' TF_type '-' glmfct str_last_trial '.mat'];

bCalculation = ~exist(fnameACI,'file');

if bCalculation == 0
    disp('ACI found on disk:')
    bCalculation = input('  Enter 1 to re-calculate (overwrite) or 0 to read stored results:');
    
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
cfg_ACI = import_cfg(cfg_game, 'dir_noise', 'N', 'N_target', ... % 'dir_target', 'N_response'
    'stim_order', 'target_names', 'response_correct_target','response_names');
cfg_ACI = arg_TF_type(cfg_ACI, flags);
cfg_ACI = arg_glmfct(cfg_ACI, flags);

cfg_ACI.flags = flags;
cfg_ACI.keyvals = keyvals;
cfg_ACI.fnameACI = fnameACI;

cfg_ACI.idx_trialselect   = cfg_ACI.keyvals.idx_trialselect; % numeros des essais utilises pour le calcul (defaut = 1:cfg_ACI.N), si possible les essais sont pris dans l'ordre de presentation

cfg_ACI.WithSNR           = 'no'; % baser le calcul sur les échantillons de bruit uniquement, au niveau de bruit de l'essai (default 'non')
cfg_ACI.WithSignal        = 'no'; % baser le calcul sur le stimulus complet (signal + bruit) (default 'non'). Si WithSNR et WithSignal sont négatifs, le calcul se base sur le bruit uniquement, sans prise en compte du SNR
cfg_ACI.dB                = 'no'; % Utiliser des prédicteurs en dB plutôt qu'en amplitude

cfg_ACI.withU             = 1; % 'yes'; % Ajouter deux paramètres U au modèle
 
switch cfg_ACI.glmfct
    case 'glmfitqp'
        check_cfg(cfg_ACI, 'prior','lambda0', 'stepsize', 'maxiter', 'nobreak', 'minDiffSecondRound');
        cfg_ACI.N_folds   = cfg_ACI.keyvals.N_folds;
        
    case 'lassoglm'
        cfg_ACI.lambda0   = [];
        cfg_ACI.N_folds    = cfg_ACI.keyvals.N_folds; 
        
    case 'classic_revcorr'
        
end
if cfg_ACI.zscore == 0
    error('%s: the glmfct options require that cfg_ACI.zscore is 1',upper(mfilename))
end
 
% cfg_ACI.filtrage          = 'non'; % Effectuer un filtrage de l'ACI
% if isyes(cfg_ACI.filtrage)   
%     % caracteristiques du filtre :
%     cfg_ACI.flt_freqcoup  = 1/0.050; % frequence de coupure en Hz
%     cfg_ACI.flt_quefrcoup = 1/500;   % quefrence de coupure en s
% end

cfg_ACI.freq_analysis = cfg_ACI.keyvals.f_limits; % bande de frequences pour analyse
cfg_ACI.time_analysis = cfg_ACI.keyvals.t_limits; % bande de temps pour analyse
 
if flags.do_tf
    check_cfg(cfg_ACI, 'freq_analysis', 'time_analysis', 'overlap', 'Nwindow', 'NFFT');
end
   
if flags.do_lyon
    check_cfg(cfg_ACI, 'freq_analysis', 'time_analysis', 'decimation', 'earQ', 'stepfactor');
end

% cfg_ACI = set_default_cfg(cfg_ACI, 'dir_noise', 'ListeBruit', 'dir_target', ... 
%     'ListeSignal', 'N', N, 'N_target', N_target, 'N_response', N_response, 'target_names', {'signal#1', 'signal#2'}, 'response_names', {'signal#1', 'signal#2'}, ...
%     'IdxTrialsLoad', 1:N, 'WithSNR', 'non'); %'Basis', 'laplacian', 'folds', 5
cfg_ACI = set_default_cfg(cfg_ACI, 'N_trialselect', length(cfg_ACI.idx_trialselect));

%% 3. Loading the data: Reading the sound waveforms

% ADD HERE: RECONSTRUCTION OF NOISE WAVEFORMS FOR SEEDS EXPERIMENTS
if bCalculation || do_recreate_validation
    
    if ~isempty(cfg_ACI.keyvals.dir_noise)
        fprintf('%s: Using dir_noise specified as input parameter by the user\n',upper(mfilename));
        fprintf('\tNew cfg_game.dir_noise=%s\n',cfg_ACI.keyvals.dir_noise);
        cfg_ACI.dir_noise = cfg_ACI.keyvals.dir_noise;
    end
    
    if ~exist(cfg_ACI.dir_noise,'dir')
        % If it does not exist
        if isfield(opts_ACI,'dir_noise')
            dir_noise = opts_ACI.dir_noise;
            if ~exist(dir_noise,'dir')
                error('%s: No valid noise directory (''dir_noise''). cfg_game contains a folder that was not found, please a valid opts_ACI.dir_noise.',upper(mfilename))
            end
            cfg_ACI.dir_noise = dir_noise;
        end
    end
        
    [Data_matrix,cfg_ACI] = data_load(cfg_ACI, ListStim, 'WithSNR', cfg_ACI.WithSNR, 'WithSignal', cfg_ACI.WithSignal, 'dB', cfg_ACI.dB);
end
 
%% 4. Preprocessing of the data, before the ACI calculation
if bCalculation || do_recreate_validation        
    
    [y, y_correct, X, U, cfg_ACI] = Script4_Calcul_ACI_preprocess(cfg_ACI, data_passation, Data_matrix);
    
end

%% 5. Calculation of the ACI
if bCalculation
    
    [ACI, results, cfg_ACI] = Script4_Calcul_ACI_calculate(cfg_ACI, y, y_correct, X, U);
    
    save(fnameACI, 'ACI', 'cfg_ACI', 'results')
    results.fnameACI = fnameACI;
    results.fnameACI_description = 'File name where the fastACI results were stored...';
    
else
    ACI = [];
    cfg_ACI = [];
    results = [];
    load(fnameACI,'ACI','cfg_ACI','results');
 
    if isfield(results,'lambdas')
        % This data display is just to make the user aware what for data are being loaded
        [round(results.lambdas') round(10*results.cvgofs')/10]
    end
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

%%% TODO: check that ACI is already reshaped after the ACI calculation...
N_f = length(cfg_ACI.f);
N_t = length(cfg_ACI.t);
switch cfg_ACI.glmfct
    case {'glmfitqp','CI_glmqp_fct'}
        % %%% results.fits{1,idx_lambdaoptim}.w(1:end-2) is incorrect, as that
        % % is the fitting for the first fold using the 'lambdaoptim', however
        % % with lambdaoptim there is a very final fit, after finishing with
        % % the folds and those results are stored in 'results.finalfit':
        % %
        % % % ACI=reshape(results.fits{1,idx_lambdaoptim}.w(1:end-2),length(cfg_ACI.f),length(cfg_ACI.t));
        % 
        % ACI=reshape(results.finalfit.w(1:end-2),length(cfg_ACI.f),length(cfg_ACI.t));

    case 'lassoglm'
        w = ACI;
        if size(ACI,1)~=N_f || size(ACI,2)~=N_t
            ACI=reshape(w,N_f,N_t);
        end
        
    case 'classic_revcorr'
        % Nothing to do, ACI should already be reshaped
end

Max_here = max(max(ACI));
Min_here = min(min(ACI));
ACI_norm = 2*(ACI-Min_here)/(Max_here-Min_here)-1;

% if isfield(results,'cvgofs')
%     % determination du lambda optimal
%     [mincvgof, idx_lambdaoptim] = min(mean(results.cvgofs,1));
% 
%     figure
%     plot(results.lambdas,results.cvgofs)
%     hold on;
%     plot(results.lambdas,mean(results.cvgofs,1),'k','LineWidth',2.5);
%     plot(results.lambdas(idx_lambdaoptim),mincvgof,'k*','MarkerSize',10)
%     set(gca, 'XScale', 'log'); set(gca,'YGrid','on') ;
%     title('CV deviance and optimal lambda'); 
%     ylabel(['CV deviance (' num2str(cfg_ACI.N_folds) '-fold)']); xlabel('lambda')
% 
%     fprintf(['Valeur du lambda optimal : ' num2str(results.lambdas(idx_lambdaoptim)) '\n'])
% end

bPlot = 1;
if bPlot || nargout == 0
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
        % affichage_tf(ACI_norm, 'CI', 'cfg',cfg_ACI)
        
        % figure
        % affichage_tf(ACI, 'CInorm', 'cfg', cfg_ACI)
        % 
        % figure
        % affichage_tf(ACI, 'tvalue', 'cfg', cfg_ACI)
        % 
        % figure
        % affichage_tf(ACI, 'prob', 'cfg', cfg_ACI)
    end
end

results.ACI_norm       = ACI_norm;
results.ACI_norm_description = 'normalised ACI image using the optimal lambda, with weights between -1 and 1';

results.cfg_game       = cfg_game;
results.data_passation = data_passation;