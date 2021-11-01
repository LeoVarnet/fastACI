function cfg_inout = arg_glmfct(cfg_inout,flags)
% function cfginout = arg_glmfct(cfg_inout,flags)
%
% Enabled so far:
%   'permutation' or 'no_permutation': to perform or not the permutation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(cfg_inout,'glmfct')
    cfg_inout.glmfct = flags.glmfct;
end
glmfct = cfg_inout.glmfct;

zscore = 1;
if ~isfield(cfg_inout,'zscore')
    cfg_inout.zscore = zscore; % Parametres de calcul du cochleogramme
else
    if zscore ~= cfg_inout.zscore
        fprintf('\t%s: non-default value for field ''zscore'' is being used\n',glmfct);
    end
end

% Default number of permutations N_perm (if do_permutation == 1):
switch glmfct
    case {'glmfitqp','CI_glmqpoptim_fct'}
        
        if strcmp(glmfct,'CI_glmqpoptim_fct')
            warning('glmfct name ''%s'' will be soon deprecated, please use ''glmfitqp'' instead...',glmfct)
        end
        
        prior = 'smoothness';
        lambda0 = 5;
        stepsize = 1.5; % Progression step for lambda
        maxiter  = 30; % Number of iterations
        nobreak  = 1;
        minDiffSecondRound = 10;
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'prior',prior); 
        if bAssigned == 0
            if ~strcmp(cfg_inout.prior,prior)
                fprintf('\t%s: non-default value for field ''prior'' is being used\n',glmfct);
            end
        end
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'lambda0',lambda0); 
        if bAssigned == 0
            if ~strcmp(cfg_inout.lambda0,lambda0)
                fprintf('\t%s: non-default value for field ''lambda0'' is being used\n',glmfct);
            end
        end
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'stepsize',stepsize); 
        if bAssigned == 0
            if ~strcmp(cfg_inout.stepsize,stepsize)
                fprintf('\t%s: non-default value for field ''stepsize'' is being used\n',glmfct);
            end
        end
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'maxiter',maxiter); 
        if bAssigned == 0
            if ~strcmp(cfg_inout.maxiter,maxiter)
                fprintf('\t%s: non-default value for field ''maxiter'' is being used\n',glmfct);
            end
        end
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'nobreak',nobreak); 
        if bAssigned == 0
            if ~strcmp(cfg_inout.nobreak,nobreak)
                fprintf('\t%s: non-default value for field ''nobreak'' is being used\n',glmfct);
            end
        end
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'minDiffSecondRound',minDiffSecondRound);
        if bAssigned == 0
            if ~strcmp(cfg_inout.minDiffSecondRound,minDiffSecondRound)
                fprintf('\t%s: non-default value for field ''minDiffSecondRound'' is being used\n',glmfct);
            end
        end
        
    case {'lasso','lassoglm','lassoslow','lassoglmslow'}
        
        Nlevel    = 5; % number of levels (= degrees of filtering) in the Gaussian pyramid
        Nlevelmin = 2; % minimum level considered in the analysis (default=2, 
                       % Leo suggests 1)
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'lasso_Nlevel',Nlevel);
        if bAssigned == 0
            if ~strcmp(cfg_inout.lasso_Nlevel,Nlevel)
                fprintf('\t%s: non-default value for field ''lasso_Nlevel'' is being used\n',glmfct);
            end
        end
                    
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'lasso_Nlevelmin',Nlevelmin);
        if bAssigned == 0
            if ~strcmp(cfg_inout.lasso_Nlevelmin,Nlevelmin)
                fprintf('\t%s: non-default value for field ''lasso_Nlevelmin'' is being used\n',glmfct);
            end
        end
        
    case 'classic_revcorr'
end