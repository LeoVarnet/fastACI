function [ACI,results,cfg_ACI] = fastACI_getACI_calculate(cfg_ACI, y, y_correct, X, U)
% function [ACI,results,cfg_ACI] = fastACI_getACI_calculate(cfg_ACI)
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
% New name      Old name        Changed on:
%
% Old name: Script4_Calcul_ACI_calculate.m (changed on 21/05/2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

glmfct = cfg_ACI.glmfct;

flags = cfg_ACI.flags; 

%% 5. Calculation of the ACI

N_f = length(cfg_ACI.f);
N_t = length(cfg_ACI.t);

fprintf('Starting the ACI assessment\n');% ,Analysis_condition);
tic

do_permutation = flags.do_permutation;
if do_permutation
    cfg_perm = cfg_ACI.cfg_perm;
    N_perm = cfg_perm.N_perm;
end

switch glmfct
    case {'glmfitqp','glm_L2'}
        [ACI, results, cfg_ACI] = CI_glmqpoptim_fct(cfg_ACI, y, y_correct, X, U); 

        if do_permutation
            cfg_perm = Merge_structs(cfg_perm,cfg_ACI);
            cfg_perm.maxiter = 1;
            cfg_perm.lambda0 = results.finallambda * sqrt(cfg_ACI.stepsize); % crossValidate.m: stepsize will be 'compensated'

            for i = 1:N_perm
                fprintf('\t Assessing permuted ACI: %.0f of %.0f\n',i,N_perm);

                U_perm_here = squeeze(cfg_perm.U_perm(:,:,i));
                [ACI_perm_here, results_perm, cfg_perm] = CI_glmqpoptim_fct(cfg_perm, cfg_perm.y_perm(:,i), [], X, U_perm_here); 
                ACI_perm(:,:,i) = reshape(ACI_perm_here, [N_f,N_t]);
            end

            ACI_perm_CI_low  = prctile(ACI_perm,5,3);
            ACI_perm_CI_high = prctile(ACI_perm,95,3);

            results.idxs_perm = cfg_perm.idxs_perm;
            results.ACI_perm  = ACI_perm;
            results.ACI_perm_CI_low   = ACI_perm_CI_low;
            results.ACI_perm_CI_high  = ACI_perm_CI_high;
        end

        % The following warning appears when using the Optimisation toolbox 
        %   of MATLAB in version R2020a and it seems to have originated as 
        %   of version R2018a (according to some forums on the internet).
        [msg,warnID] = lastwarn;
        msg_ref='The quasi-newton algorithm does not use analytic Hessian. Hessian flag in options will be ignored (supplied Hessian will not be used).';
        warnID_ref = 'optim:fminunc:HessIgnored';
        if strcmp(msg,msg_ref) || strcmp(warnID,warnID_ref)
            warning('Turn on manually the use of the old optimisation toolbox of MATLAB (go to glmfitqp.m->function il_irls)')
        end

        % results.evaluation.opts

    case {'classic_revcorr','correlation'}
        % The following fitting uses all elements of y and X (instead of 
        %     splitting the elements into training and test trials):
        [ACI, pval_ACI] = corr(y,X);
        ACI = reshape(ACI, [N_f,N_t]);
        pval_ACI = reshape(pval_ACI, [N_f,N_t]);
        
        for i_folds = 1:10
            % Added by Alejandro on 20/04/2023, see pres_osses2023_04_WASdag('fig3a') % and 'fig3b'
            % Because the fitting stored in ACI was obtained from all trials, 
            %     the fitting is performed again using 90% of the data for 
            %     training and 10% of the data for test. These fittings are
            %     stored in ACI_subset and are purely used to obtain an 
            %     appropriate goodness-of-fit metric (NOT overestimated).
            %     ACI_subset is overwritten in every 'i_folds'-cycle:
            idx90 = round(.9*size(X,1));
            idx_random = randperm(size(X,1));
            idx_test = idx_random(idx90+1:end);
            idx_random = idx_random(1:idx90);
            ACI_subset = corr(y(idx_random),X(idx_random,:));
        
            FitInfo_local = Get_PA_for_classic_revcorr(ACI_subset,X(idx_test,:),y(idx_test)); % FitInfo from trials that were not used to obtain ACI_subset
            FitInfo.PC(i_folds,1) = FitInfo_local.PC;
            FitInfo.idx_test(i_folds,:) = idx_test;
        end
                
        results.ACI = ACI;
        results.pval_ACI = pval_ACI;
        results.FitInfo = FitInfo;
        
        if do_permutation
            for i = 1:N_perm
                fprintf('\t Assessing permuted ACI: %.0f of %.0f\n',i,N_perm);

                U_perm_here = squeeze(cfg_perm.U_perm(:,:,i));
                ACI_perm_here = corr(cfg_perm.y_perm(:,i), X); 
                ACI_perm(:,:,i) = reshape(ACI_perm_here, [N_f,N_t]);
            end

            ACI_perm_CI_low  = prctile(ACI_perm,5,3);
            ACI_perm_CI_high = prctile(ACI_perm,95,3);

            results.idxs_perm = cfg_perm.idxs_perm;
            results.ACI_perm  = ACI_perm;
            results.ACI_perm_CI_low   = ACI_perm_CI_low;
            results.ACI_perm_CI_high  = ACI_perm_CI_high;
        end
    case 'weighted_sum'
        X_11 = mean(X(y==1 & y_correct == 1,:));
        X_10 = mean(X(y==1 & y_correct == 0,:));
        X_01 = mean(X(y==0 & y_correct == 1,:));
        X_00 = mean(X(y==0 & y_correct == 0,:));

        ACI = X_11 - X_01 + X_10 -X_00;
        ACI = reshape(ACI, [N_f,N_t]);
        results.ACI     = ACI;
    case 'glm'
        % Implementation by Leo Varnet on 15/11/2023. It requires the 
        %   statistical toolbox. This fitting was used for our 'segmentation'
        %   experiment (segmentation_user.m, 2022-2023).
        tic
        [B,Dev,Stat] = glmfit(X,y,'binomial','link','logit');
        ACI_size = [length(cfg_ACI.f),length(cfg_ACI.t)];
        ACI = reshape(B(2:end),ACI_size); % AO: Why is B(1) excluded?
        results.B = B;
        results.Dev = Dev;
        results.dfe = Stat.dfe;
        results.covb = Stat.covb;
        results.se = Stat.se;
        results.coeffcorr = Stat.coeffcorr;
        results.t = Stat.t;
        results.p = Stat.p;
        
        for i_folds = 1:10
            % Added by Alejandro on 9/08/2023
            % Because the fitting stored in ACI was obtained from all trials, 
            %     the fitting is performed again using 90% of the data for 
            %     training and 10% of the data for test. These fittings are
            %     stored in ACI_subset and are purely used to obtain an 
            %     appropriate goodness-of-fit metric (NOT overestimated).
            %     ACI_subset is overwritten in every 'i_folds'-cycle:
            warning_id = 'stats:glmfit:IllConditioned';
            warning('off',warning_id);

            FitInfo_local = Get_PA_for_glm(X,y); % FitInfo from trials that were not used to obtain ACI_subset
            FitInfo.PC(i_folds,1) = FitInfo_local.PC;
            FitInfo.idx_test(i_folds,:) = FitInfo_local.idx_test;
            
            warning('on',warning_id);
        end
        results.FitInfo = FitInfo;
        
    case 'lassoglm_original'
        % Implementation by Leo Varnet on 13/01/2023
        error('The GLM fitting using %s has not been validated yet using fastACI',glmfct);
        tic
        [B,Fitinfo] = lassoglm(X,y,'binomial','link','logit');
        ACI = reshape(B(2:end),[length(cfg_ACI.f),length(cfg_ACI.t)]);
        % results.B = B;
        % results.Dev = Dev;
        % results.dfe = Stat.dfe;
        % results.covb = Stat.covb;
        % results.se = Stat.se;
        % results.coeffcorr = Stat.coeffcorr;
        % results.t = Stat.t;
        % results.p = Stat.p;

    case {'lassoglm','lasso','l1lm','l1glm','glm_L1_GB','lm_L1_GB'}
        N_folds = cfg_ACI.N_folds; % k_folds validation
        % cfg_ACI = Ensure_field(cfg_ACI,'lambda0',[]);
        lambda0 = cfg_ACI.lambda0;

        tic
        switch glmfct
            case 'lassoglm'
                [B,FitInfo] = lassoglm(X,y,'binomial','CV',N_folds,'lambda',lambda0,'link','logit');
                results.cvgofs  = FitInfo.Deviance; % this is already the mean value
                idxlambda = FitInfo.IndexMinDeviance;
                
            case 'lasso'
                [B,FitInfo] = lasso(X,y,'CV',N_folds);
                idxlambda = FitInfo.IndexMinMSE;
                
            case {'l1lm','lm_L1_GB'}
                lambda0 = cfg_ACI.keyvals.lambda;
                if isempty(lambda0)
                    [B,FitInfo] = lassoslow(X,y,N_folds); % default lambdas will be tested
                else
                    [B,FitInfo] = lassoslow(X,y,N_folds,lambda0); 
                end
                [~, idxlambda] = min(mean(FitInfo.MSE_test,2));
                
            case {'l1glm','glm_L1_GB'}
                lambda0 = cfg_ACI.keyvals.lambda;
                if isempty(lambda0)
                    [B,FitInfo] = lassoglmslow(X,y,N_folds); % default lambdas will be tested
                else
                    [B,FitInfo] = lassoglmslow(X,y,N_folds,lambda0); 
                end
                [~, idxlambda] = min(mean(FitInfo.Dev_test,2));
        end
        toc
        
        results.idxlambda = idxlambda;
        [ACI, cfg_ACI, sumReWeight] = Convert_lasso_B2ACI(B, cfg_ACI, idxlambda, cfg_ACI.keyvals);
        
        %%%
        if do_permutation
            for i = 1:N_perm
                fprintf('\t Assessing permuted ACI: %.0f of %.0f\n',i,N_perm);
                 
                % U_perm_here = squeeze(cfg_perm.U_perm(:,:,i));
                % [ACI_perm_here, results_perm, cfg_perm] = CI_glmqpoptim_fct(cfg_perm, cfg_perm.y_perm(:,i), [], X, U_perm_here); 
                % ACI_perm(:,:,i) = reshape(ACI_perm_here, [N_f,N_t]);
                
                Lambda = FitInfo.Lambda(idxlambda);
                switch glmfct
                    case 'lassoglm'
                        [B_perm,FitInfo_perm] = lassoglm(X,cfg_perm.y_perm(:,i),'binomial','CV',N_folds,'lambda',lambda0,'link','logit');
                        
                    case 'lasso'
                        [B_perm,FitInfo_perm] = lasso(X,cfg_perm.y_perm(:,i),'CV',N_folds,'Lambda',Lambda);
                    otherwise
                        error('permutation test not available yet for this option')
                end
                ACI_perm(:,:,i) = Convert_lasso_B2ACI(B_perm, cfg_ACI, cfg_ACI.keyvals);
            end
            
            ACI_perm_CI_low  = prctile(ACI_perm,5,3);
            ACI_perm_CI_high = prctile(ACI_perm,95,3);
            
            results.idxs_perm = cfg_perm.idxs_perm;
            results.ACI_perm  = ACI_perm;
            results.ACI_perm_CI_low   = ACI_perm_CI_low;
            results.ACI_perm_CI_high  = ACI_perm_CI_high;
        end
        %%%

        results.lambdas = FitInfo.Lambda; 
        % results.cvgofs  = FitInfo.Deviance; % this is already the mean value
        results.FitInfo = FitInfo;
        results.B       = B;
        results.ACI     = ACI;
        results.ACIs    = sumReWeight;

    otherwise
        error('Function %s not debugged yet...',glmfct);
        % [ACI, results, cfg_ACI] = feval(cfg_ACI.glmfct, cfg_ACI, y, y_correct, X, U );
end
fprintf('Completed condition \n');
toc