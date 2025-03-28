function cfg_inout = arg_glmfct(cfg_inout,flags,keyvals)
% function cfginout = arg_glmfct(cfg_inout,flags,keyvals)
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
    case {'glmfitqp','glm_L2'}
        
        % extract relevant parameters from keyval, or default values
        
        prior = 'smoothness';
        % if ~isfield(keyvals,'lambda0')
        %     lambda0 = 5;
        % else
             lambda0 = keyvals.lambda0;
        % end
        % if ~isfield(keyvals,'stepsize')
        %     stepsize = 1.5; % Progression step for lambda
        % else
             stepsize = keyvals.stepsize;
        % end
        % if ~isfield(keyvals,'maxiter')
        %     maxiter  = 30; % Number of iterations
        % else
             maxiter = keyvals.maxiter;
        % end
        % if ~isfield(keyvals,'nobreak')
        %     nobreak  = 1;
        % else
             nobreak = keyvals.nobreak;
        % end
        % if ~isfield(keyvals,'precision')
        %     %%%TODO : default precision
        % else
             precision = keyvals.precision;
        % end
        minDiffSecondRound = 3;
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'prior',prior); 
        if bAssigned == 0
            if ~strcmp(cfg_inout.prior,prior)
                fprintf('\t%s: non-default value for field ''prior'' is being used\n',glmfct);
            end
        end
        
         [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'lambda0',lambda0); 
        % if bAssigned == 0
        %     if ~strcmp(cfg_inout.lambda0,lambda0)
        %         fprintf('\t%s: non-default value for field ''lambda0'' is being used\n',glmfct);
        %     end
        % end
        % 
         [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'stepsize',stepsize); 
        % if bAssigned == 0
        %     if ~strcmp(cfg_inout.stepsize,stepsize)
        %         fprintf('\t%s: non-default value for field ''stepsize'' is being used\n',glmfct);
        %     end
        % end
        % 
         [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'maxiter',maxiter); 
        % if bAssigned == 0
        %     if ~strcmp(cfg_inout.maxiter,maxiter)
        %         fprintf('\t%s: non-default value for field ''maxiter'' is being used\n',glmfct);
        %     end
        % end
        % 
         [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'nobreak',nobreak); 
        % if bAssigned == 0
        %     if ~strcmp(cfg_inout.nobreak,nobreak)
        %         fprintf('\t%s: non-default value for field ''nobreak'' is being used\n',glmfct);
        %     end
        % end
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'minDiffSecondRound',minDiffSecondRound);
        if bAssigned == 0
            if ~strcmp(cfg_inout.minDiffSecondRound,minDiffSecondRound)
                fprintf('\t%s: non-default value for field ''minDiffSecondRound'' is being used\n',glmfct);
            end
        end

         [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'precision',precision);
        % if bAssigned == 0
        %     if ~strcmp(cfg_inout.precision,precision)
        %         fprintf('\t%s: non-default value for field ''precision'' is being used\n',glmfct);
        %     end
        % end
        
    case {'lasso','lassoglm','l1lm','l1glm','glm_L1_GB','lm_L1_GB'}
        
        Nlevel    = keyvals.lasso_Nlevel; % 5; % number of levels (= degrees of filtering) in the Gaussian pyramid
        Nlevelmin = 2; % minimum level considered in the analysis 
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'lasso_Nlevel',Nlevel);
        if bAssigned == 0
            if cfg_inout.lasso_Nlevel~=Nlevel
                fprintf('\t%s: non-default value for field ''lasso_Nlevel'' is being used\n',glmfct);
            end
        end
                    
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'lasso_Nlevelmin',Nlevelmin);
        if bAssigned == 0
            if cfg_inout.lasso_Nlevelmin~=Nlevelmin
                fprintf('\t%s: non-default value for field ''lasso_Nlevelmin'' is being used\n',glmfct);
            end
        end
        
    case 'glm'
        % no additional parameter
    case {'classic_revcorr','correlation'}
        % no additional parameter
    case 'weighted_sum'
        % no additional parameter
end