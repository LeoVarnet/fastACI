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
    case 'glmfitqp'
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

        % The following warning appears when using the Optimization
        %   toolbox of MATLAB in version R2020a and it seems to
        %   have originated as of version R2018a (according to
        %   some forums on the internet).
        [msg,warnID] = lastwarn;
        msg_ref='The quasi-newton algorithm does not use analytic Hessian. Hessian flag in options will be ignored (supplied Hessian will not be used).';
        warnID_ref = 'optim:fminunc:HessIgnored';
        if strcmp(msg,msg_ref) || strcmp(warnID,warnID_ref)
            warning('Turn on manually the use of the old optimisation toolbox of MATLAB (go to glmfitqp.m->function il_irls)')
        end

        % results.evaluation.opts

    case 'classic_revcorr'
        [ACI, pval_ACI] = corr(y,X);
        ACI = reshape(ACI, [N_f,N_t]);
        pval_ACI = reshape(pval_ACI, [N_f,N_t]);
        results.ACI = ACI;
        results.pval_ACI = pval_ACI;

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

    case {'lassoglm','lasso','lassoslow','lassoglmslow'}
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
                
            case 'lassoslow'
                lambda0 = cfg_ACI.keyvals.lambda;
                if isempty(lambda0)
                    [B,FitInfo] = lassoslow(X,y,N_folds); % default lambdas will be tested
                else
                    [B,FitInfo] = lassoslow(X,y,N_folds,lambda0); 
                end
                [~, idxlambda] = min(mean(FitInfo.MSEtest,2));
                
            case 'lassoglmslow'
                lambda0 = cfg_ACI.keyvals.lambda;
                if isempty(lambda0)
                    [B,FitInfo] = lassoglmslow(X,y,N_folds); % default lambdas will be tested
                else
                    [B,FitInfo] = lassoglmslow(X,y,N_folds,lambda0); 
                end
                [~, idxlambda] = min(mean(FitInfo.Devtest,2));
        end
        toc
        
        results.idxlambda = idxlambda;
        [ACI, cfg_ACI, sumReWeight] = Convert_lasso_B2ACI(B, cfg_ACI, results);
        
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
                end
                ACI_perm(:,:,i) = Convert_lasso_B2ACI(B_perm, cfg_ACI);
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
        % results.w       = B(:,FitInfo.IndexMinDeviance);
        % results.finalfit.w = results.w;
        
        results.ACI = ACI;
        results.ACIs = sumReWeight;
        
        if flags.do_plot
            figure; 
            % h=pcolor(t_X, f_X, -sumReWeight(:,:,idxlambda)); set(h, 'EdgeColor', 'none'); xlabel('time (s)'); ylabel('freq (Hz)'); colorbar; title('betaSmooth: ACI obtained with a lasso regression on smooth basis'); colorbar; caxis([-1 1]*max(abs(caxis)));
            h=pcolor(cfg_ACI.t, cfg_ACI.f, ACI); 
            set(h, 'EdgeColor', 'none'); 
            xlabel('time (s)'); 
            ylabel('freq (Hz)'); colorbar; 
            title('betaSmooth: ACI obtained with a lasso regression on smooth basis'); 
            colorbar; % caxis([-1 1]*max(abs(caxis)));
        end

    otherwise
        error('Function %s not debugged yet...',glmfct);
        % [ACI, results, cfg_ACI] = feval(cfg_ACI.glmfct, cfg_ACI, y, y_correct, X, U );
end
fprintf('Completed condition \n');
toc