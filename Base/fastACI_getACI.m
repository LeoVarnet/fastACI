function [ACI,cfg_ACI,results,Data_matrix,extra_outs] = fastACI_getACI(savegame_file,varargin)
% function [ACI,cfg_ACI,results,Data_matrix,extra_outs] = fastACI_getACI(savegame_file,varargin)
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
% % Example:
%   dir_where = [fastACI_paths('dir_data') 'speechACI_Logatome-abda-S43M' filesep 'SLV' filesep 'Results' filesep];
%   savefile = [dir_where 'savegame_2021_11_19_12_50_SLV_speechACI_Logatome-abda-S43M_bumpv1p2_10dB.mat'];
%   [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(savefile);
% 
%   [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(savefile,'force_dataload');
%
% Old name: Script4_Calcul_ACI.m (changed on 7 July 2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Condition can be renamed to NameCond

if nargin == 0
    error('%s: Please spefify the identifier of the subject from whom you want to process the data',upper(mfilename));
end

Data_matrix = [];

%%%
[fnameACI, cfg_game, data_passation, ListStim, flags, keyvals] = fastACI_getACI_fname(savegame_file,varargin{:});
bCalculation = ~exist(fnameACI,'file');
if isempty(keyvals.dir_noise)
    if ~exist(cfg_game.dir_noise,'dir') && isempty(keyvals.Data_matrix)
        % Prepare ACI analysis
        cfg_game = Check_cfg_crea_dirs(cfg_game);
        
        if bCalculation == 1
            files = Get_filenames(cfg_game.dir_noise,'*.wav');
            if isempty(files) % then dir_noise does not exist yet...
                if isfield(cfg_game,'Condition')
                    fprintf('Please indicate the folder where the stimuli can be found (dir_noise)\n');
                    str_here = sprintf('Experiment=%s, subject=%s, condition=%s',cfg_game.experiment_full,cfg_game.Subject_ID, cfg_game.Condition);
                    cfg_game.dir_noise = uigetdir([pwd filesep],str_here);
                end
                cfg_game.dir_noise = [cfg_game.dir_noise filesep];
            end
        end
    end
else
    disp('')
    % cfg_game.dir_noise = keyvals.dir_noise;
end

% General parameters
do_recreate_validation = flags.do_recreate_validation;
glmfct  = flags.glmfct;
N = cfg_game.N;
% END: From argument function:

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
cfg_ACI = import_cfg(cfg_game, 'dir_noise', 'N', 'N_target', ... % 'dir_target', 'N_response'
    'stim_order', 'target_names', 'response_correct_target','response_names');

if isfield(cfg_game,'dir_target')
    cfg_ACI.dir_target = cfg_game.dir_target;
end
if isfield(cfg_game,'sessionsN')
    cfg_ACI.L_session = cfg_game.sessionsN; % used in fastACI_getACI_preprocess.m
end
if isfield(cfg_game,'Subject_ID')
    cfg_ACI.Subject_ID = cfg_game.Subject_ID;
end
if isfield(cfg_game,'experiment')
    cfg_ACI.experiment = cfg_game.experiment;
end
if isfield(cfg_game,'Condition')
    cfg_ACI.Condition = cfg_game.Condition;
end

extra_outs.bCalculation = bCalculation;

cfg_ACI = arg_TF_type(cfg_ACI, flags, keyvals);
cfg_ACI = arg_glmfct( cfg_ACI, flags, keyvals);

cfg_ACI.flags = flags;
cfg_ACI.keyvals = keyvals;
cfg_ACI.fnameACI = fnameACI;

cfg_ACI.idx_trialselect   = cfg_ACI.keyvals.idx_trialselect; % numeros des essais utilises pour le calcul (defaut = 1:cfg_ACI.N), si possible les essais sont pris dans l'ordre de presentation
cfg_ACI.withU             = 1; % 'yes'; % Ajouter deux paramètres U au modèle
 
switch cfg_ACI.glmfct
    case 'glmfitqp'
        check_cfg(cfg_ACI, 'prior','lambda0', 'stepsize', 'maxiter', 'nobreak', 'minDiffSecondRound');
        cfg_ACI.N_folds   = cfg_ACI.keyvals.N_folds;
        
    case {'lassoglm','lasso','l1lm','l1glm'}
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

if ~isempty(cfg_ACI.keyvals.ACI_crosspred)
    
    if ischar(cfg_ACI.keyvals.ACI_crosspred)
        % then it is a char
        ACI_crosspred = {cfg_ACI.keyvals.ACI_crosspred}; % it is converted into a cell array
    else
        ACI_crosspred = cfg_ACI.keyvals.ACI_crosspred;
    end
    N_crosspred = length(ACI_crosspred);
    
    fname_cross = Save_crosspred(fnameACI,[],mfilename,keyvals.suffix_crosspred);
    fname_cross = [fname_cross '.mat'];
    if ~exist(fname_cross,'file')
        bCrossPred = 1;
    else
        var = load(fname_cross);
        N_crosspred_here = length(var.crosspred);
        count = 0;
        for i = 1:N_crosspred_here
            if var.crosspred(i).bProcessed ~= 0
                count = count+1;
            end
        end
        
        fprintf('Cross prediction found on disk\n')
        
        if N_crosspred_here ~= N_crosspred
            error('Mismatch between the requested cross predictions and the number of expected cross predictions in the existing crosspred.mat file\n To proceed, check your local files and remove the previously existing crosspred.mat file. You can then re-run the script')
        end
        if count ~= N_crosspred_here
            % This means that the cross prediction is incomplete
            bCrossPred = 1;
        else
            % This means that the cross prediction is complete
            bCrossPred = 0;
        end
    end
else
    bCrossPred = 0;
end
%% 3. Loading the data: Reading the sound waveforms

% ADD HERE: RECONSTRUCTION OF NOISE WAVEFORMS FOR SEEDS EXPERIMENTS
if bCalculation || do_recreate_validation || flags.do_force_dataload || bCrossPred
    
    if ~exist(cfg_ACI.dir_noise,'dir') && isempty(cfg_ACI.keyvals.dir_noise) && isempty(keyvals.Data_matrix)
    
        % If it does not exist
        if isfield(cfg_game,'seeds_order')
            % Nothing to do, the waveforms will be generated inside data_load
            cfg_ACI.cfg_game = cfg_game;
            % error('Under development...')            
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
    
    if isfield(cfg_ACI,'dir_target') && ~strcmp(fileparts(fileparts(cfg_ACI.dir_noise)), fileparts(fileparts(cfg_ACI.dir_target))) 
        % Extra check: if someone enters this part of the code, maybe he/she
        % does not have compatible dir_noise and dir_target directories, and 
        % therefore we throw a warning that appears for 10 s.
        warning('dir_noise and dir_target were found to be located in different root folders, please check that this is correct (ignore this message otherwise)');
        pause(10)
    end
    
    if (data_passation.i_current ~= cfg_game.N && bCalculation) || (data_passation.i_current ~= cfg_game.N && bCrossPred) || flags.do_force_dataload
        if data_passation.i_current ~= cfg_game.N
            fprintf('%s: Less trials have been tested by this participant than the expected cfg_game.N=%.0f trials\n',upper(mfilename),cfg_game.N);
            fprintf('\tPress ctrl+C to cancel the current ACI calculation, otherwise, the ACI will be obtained for less trials...\n');
            fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            pause(10);
            
            fprintf('\tToo late: ACI will be assessed using %.0f trials only\n',data_passation.i_current);
        end
        
        N = data_passation.i_current;
        cfg_ACI.N = N;
        cfg_ACI.stim_order = cfg_ACI.stim_order(1:N);
        cfg_ACI.idx_trialselect = 1:N;
        cfg_ACI.N_trialselect = N;
    end
    
    if isempty(keyvals.Data_matrix)
        % Loading the data regularly:
        if isempty(keyvals.script_dataload)
            % If empty: then it looks for candidates for script to load the data
            script4dataload = [cfg_game.experiment '_dataload'];
            if exist([script4dataload '.m'],'file')
                if bCalculation
                    fprintf('%s: dataload script found for this experiment (%s.m)\n',upper(mfilename),script4dataload);
                    fprintf('\t If you want to use the default fastACI_getACI_dataload.m file instead, abort this\n');
                    fprintf('\t processing now (press ctrl+c) and enter a new keyval called ''force_default_dataload'' to 1 \n');
                    pause(10); % Leo: I put this pause for a reason, don't remove or modify this please
                end
                exp2eval = sprintf('[Data_matrix,cfg_ACI] = %s(cfg_ACI, ListStim, cfg_game, data_passation);',script4dataload);
                eval(exp2eval);
            else
                % Default:
                [Data_matrix,cfg_ACI] = fastACI_getACI_dataload(cfg_ACI, ListStim, cfg_game, data_passation);
            end
        else
            script4dataload = keyvals.script_dataload;
            if strcmp(script4dataload(end-1:end),'.m') % in case the extension is added
                script4dataload = script4dataload(1:end-2); % removing the extension
            end
                
            exp2eval = sprintf('[Data_matrix,cfg_ACI] = %s(cfg_ACI, ListStim, cfg_game, data_passation);',script4dataload);
            eval(exp2eval);
        end
    else
        Data_matrix = keyvals.Data_matrix;
        
        if ~exist(cfg_ACI.dir_noise,'dir')
            if keyvals.consistency_check == 1
                warning('The consistency check will be skipped because the sounds were not found on disk...');
                fprintf('Pausing for 5 second...')
                pause(5);
                keyvals.consistency_check = 0;
            end
        end
        switch flags.TF_type
            case {'spect','gammatone'} % Then Data_matrix is checked for consistency
                if keyvals.consistency_check
                    N_ref = 10;
                    N = cfg_ACI.N;
                    cfg_ACI_ref = cfg_ACI;
                    cfg_ACI_ref.N = N_ref;
                    [Data_matrix_ref,cfg_ACI] = fastACI_getACI_dataload(cfg_ACI_ref, ListStim, cfg_game, data_passation);
                    cfg_ACI.N = N; % restoring the initial N

                    for ii = 1:N_ref
                        diffe(ii) = sum(Data_matrix(ii,:)-Data_matrix_ref(ii,:));
                    end
                    if sum(diffe) ~= 0
                        error('The input Data_matrix seems to be different from the expected Data_matrix loaded from ListStim');
                    end
                end
                bConsistency_check = keyvals.consistency_check;
            otherwise
                bConsistency_check = 0;
        end
        
        if bConsistency_check == 0
            if ~isfield(cfg_ACI,'f')
                cfg_ACI.f = transpose(1:size(Data_matrix,2));
                cfg_ACI.f_description = 'Frequency bin';
            end
            if ~isfield(cfg_ACI,'t')
                cfg_ACI.t = 1:size(Data_matrix,3);
                cfg_ACI.t_description = 'Time bin';
            end
            % cfg_ACI.N_t = lengt;
            % N_f = cfg_ACI.N_f;
            warning('flags.TF_type=%s seems to be a custom T-F configuration. No Data_matrix consistency will be checked...',flags.TF_type);
        end
    end
end
 
%% 4. Preprocessing of the data, before the ACI calculation
if bCalculation || bCrossPred || do_recreate_validation        
    
    [y, y_correct, X, U, cfg_ACI] = fastACI_getACI_preprocess(cfg_ACI, data_passation, Data_matrix);
    
end

%% 5. Calculation of the ACI
if bCalculation
    
    [ACI, results, cfg_ACI] = fastACI_getACI_calculate(cfg_ACI, y, y_correct, X, U);
    
    if isfield(cfg_ACI.keyvals,'Data_matrix')
        cfg_ACI.keyvals = Remove_field(cfg_ACI.keyvals,'Data_matrix');
    end
    info_toolbox = Get_toolbox_info(mfilename);
    try
        save(fnameACI, 'ACI', 'cfg_ACI', 'results','info_toolbox');
    catch
        [dir_where,file,ext] = fileparts(fnameACI);
        mkdir(dir_where);
        % warning('No destination folder could be created...')
        save(fnameACI, 'ACI', 'cfg_ACI', 'results','info_toolbox');
    end
    results.fnameACI = fnameACI;
    results.fnameACI_description = 'File name where the fastACI results were stored...';
    
else
    ACI = [];
    cfg_ACI = [];
    results = [];
    
    bCheck_target = 0;
    if isfield(cfg_ACI,'dir_target')
        dir_target = cfg_ACI.dir_target;
        bCheck_target = 1;
    end
    bCheck_noise = 0;
    if isfield(cfg_ACI,'dir_noise')
        dir_noise = cfg_ACI.dir_noise;
        bCheck_noise = 1;
    end
    
    load(fnameACI,'ACI','cfg_ACI','results');
 
    if isfield(cfg_ACI,'fnameACI')
        if ~strcmp(cfg_ACI.fnameACI,fnameACI)
            % In this case fnameACI is probably a location where the ACI MAT 
            %   file was moved to, after the creation of the MAT file:
            cfg_ACI.fnameACI = fnameACI;
        end
    end
    if bCheck_noise
        if ~strcmp(cfg_ACI.dir_noise,dir_noise)
            cfg_ACI.dir_noise = dir_noise;
        end
    end
    if bCheck_target
        if ~strcmp(cfg_ACI.dir_target,dir_target)
            cfg_ACI.dir_target = dir_target;
        end
    end
    % if isfield(results,'lambdas')
    %     % This data display is just to make the user aware what for data are being loaded
    %     [round(results.lambdas') round(10*results.cvgofs')/10]
    % end
end
 
results.fnameACI = fnameACI;
results.fnameACI_description = 'File name where the fastACI results were stored...';

%% 6. Validation, if requested:

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

%cross validation 

if bCrossPred
    
    if exist(fname_cross,'file')
        load(fname_cross);
    end
    
    for i = 1:N_crosspred
        
        if exist(fname_cross,'file')
            bProcessed_already = crosspred(i).bProcessed;
        else
            bProcessed_already = 0;
        end
        
        if exist(ACI_crosspred{i},'file')
            bProceed = 1;
        else
            bProceed = 0;
        end
        
        crosspred(i).Subject_test = cfg_ACI.Subject_ID;
        crosspred(i).fnameACI_ref = ACI_crosspred{i};
        crosspred(i).bProcessed = bProceed;
        
        if bProceed && ~bProcessed_already
            fprintf('%s: Cross validating participant %s''s data using %s\n', ...
                upper(mfilename),crosspred(i).Subject_test,crosspred(i).fnameACI_ref);
            var = load(ACI_crosspred{i},'cfg_ACI');
            cfg_crosspred = var.cfg_ACI;

            var = load(ACI_crosspred{i},'results');
            results_crosspred = var.results;

            crosspred(i).ACI_crosspred = ACI_crosspred{i};
            %%%%TODO%%%%
            switch glmfct
                case 'l1glm'
                    bFolds_from_ref_ACI = 1; 
                    if bFolds_from_ref_ACI
                        % This is the default
                        CV=results.FitInfo.CV;
                    else
                        %  Tested on 12/07/2022
                        error('Using the folds from %s',ACI_crosspred{i});
                        CV=var.results.FitInfo.CV; % reads the folds from the ACI to be adopted
                    end
                        
                    % To do: change CV.training into CV.train

                    %[~, idx_lambda_optim] = min(mean(results_crosspred.FitInfo.Devtest,2));
                    crosspred(i).lambdas = results_crosspred.lambdas;
                    for i_lambda = 1:length(results_crosspred.lambdas)
                        for i_fold = 1:cfg_crosspred.N_folds
                            %the following lines are copied from function lassoglmslow.
                            %Maybe we should use a separate function
                            idx_training = CV.training(i_fold); % idxs of the training set in this fold
                            idx_test     = CV.test(i_fold); % idxs for the test (validation) in this fold

                            coef = [results_crosspred.FitInfo.Intercept(i_lambda,i_fold); results_crosspred.FitInfo.B(:,i_lambda,i_fold)];

                            %%% Training:
                            yhat_train = glmval(coef,X(idx_training,:),'logit');
                            [PC,MSE,Dev,MSE_rounded] = Get_prediction_metrics(yhat_train,y,idx_training);

                            crosspred(i).MSE_train(i_lambda,i_fold) = MSE; % changed from MSEtrain_crosspred to MSE_train
                            crosspred(i).PC_train(i_lambda,i_fold) = PC;
                            crosspred(i).yhat_train(i_lambda,i_fold,1:length(yhat_train)) = yhat_train;

                            %%% Test (or validation):
                            yhat_test = glmval(coef,X(idx_test,:),'logit'); % X(CV.test(i_fold),:)*B_temp + FitInfo_temp.Intercept;
                            [PC,MSE,Dev,~, yhat_test_rounded,PC_t,MSE_t,Dev_t] = Get_prediction_metrics(yhat_test,y,idx_test);

                            crosspred(i).Dev_test(i_lambda,i_fold) = Dev; % -2*(sum(log(binopdf(y(idx_test),1,yhat_test))) - sum(log(binopdf(y(idx_test),1,y(idx_test)))));
                            crosspred(i).MSE_test(i_lambda,i_fold) = MSE;
                            crosspred(i).PC_test(i_lambda,i_fold)  = PC;
                            crosspred(i).yhat_test(i_lambda,i_fold,1:length(yhat_test)) = yhat_test;

                            crosspred(i).PC_test_t(i_lambda,i_fold,1:length(yhat_test))  = PC_t; % yhat_test_rounded==y(idx_test);
                            crosspred(i).Dev_test_t(i_lambda,i_fold,1:length(yhat_test)) = Dev_t; % -2*((log(binopdf(y(idx_test),1,yhat_test))) - (log(binopdf(y(idx_test),1,y(idx_test)))));
                        end
                    end

                otherwise
                    error('cross validation not implemented yet for this glmfct');
            end
        end
    end
    Save_crosspred(fnameACI,crosspred,mfilename,keyvals.suffix_crosspred);
    
    results.crosspred = crosspred;
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
        out_affichage = affichage_tf(ACI, 'CI', 'cfg',cfg_ACI);
        title(glmfct)

        subplot(1,2,2)
        % affichage_tf(results.ACI_norm, 'CI', 'cfg',cfg_ACI); hold on
        affichage_tf(ACI_ex, 'CI', 'cfg',cfg_ACI);
        title([glmfct ' (only within CI)'])
        %%% End plotting permutation test
    else
        
        figure;
        out_affichage = affichage_tf(ACI, 'CI', 'cfg',cfg_ACI);
        title(glmfct)
        
        if isfield(cfg_ACI,'t_description')
            xlabel(cfg_ACI.t_description);
        end
        if isfield(cfg_ACI,'f_description')
            ylabel(cfg_ACI.f_description);
        end
        
        % figure
        % affichage_tf(ACI_norm, 'CI', 'cfg',cfg_ACI); affichage_tf(ACI, 'CInorm', 'cfg', cfg_ACI)
        % affichage_tf(ACI, 'tvalue', 'cfg', cfg_ACI); affichage_tf(ACI, 'prob', 'cfg', cfg_ACI)
    end
    extra_outs.out_affichage = out_affichage;
    
end

results.ACI_norm       = ACI_norm;
results.ACI_norm_description = 'normalised ACI image using the optimal lambda, with weights between -1 and 1';

results.cfg_game       = cfg_game;
results.data_passation = data_passation;
