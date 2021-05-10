function [ACI,results,cfg_ACI] = Script4_Calcul_ACI_calculate(cfg_ACI, y, y_correct, X, U)
% function [ACI,results,cfg_ACI] = Script4_Calcul_ACI_calculate(cfg_ACI)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

glmfct = cfg_ACI.glmfct;

%% 5. Calculation of the ACI

N_f = length(cfg_ACI.f);
N_t = length(cfg_ACI.t);

fprintf('Starting the ACI assessment\n');% ,Analysis_condition);
tic

do_permutation = cfg_ACI.flags.do_permutation;
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
        ACI = corr(y,X);
        ACI = reshape(ACI, [N_f,N_t]);
        results.ACI = ACI;

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

    case 'lassoglm'
        N_folds = cfg_ACI.N_folds; % k_folds validation
        % cfg_ACI = Ensure_field(cfg_ACI,'lambda0',[]);
        lambda0 = cfg_ACI.lambda0;

        tic
        [B,FitInfo] = lassoglm(X,y,'binomial','CV',N_folds,'lambda',lambda0,'link','logit');
        toc

        results.lambdas = FitInfo.Lambda; 
        results.cvgofs  = FitInfo.Deviance; % this is already the mean value
        results.FitInfo = FitInfo;
        results.B       = B;
        results.w       = B(:,FitInfo.IndexMinDeviance);
        % results.finalfit.w = results.w;
        ACI = B(:,FitInfo.IndexMinDeviance);

    otherwise
        error('Function %s not debugged yet...',glmfct);
        % [ACI, results, cfg_ACI] = feval(cfg_ACI.glmfct, cfg_ACI, y, y_correct, X, U );
end
fprintf('Completed condition \n');
toc
