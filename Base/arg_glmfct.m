function cfg_inout = arg_glmfct(cfg_inout,flags)
% function cfginout = arg_glmfct(cfg_inout,flags)
%
% Enabled so far:
%   'permutation' or 'no_permutation': to perform or not the permutation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_inout.glmfct = flags.glmfct;

zscore = 1;
if ~isfield(cfg_inout,'zscore')
    cfg_inout.zscore = zscore; % Parametres de calcul du cochleogramme
else
    if zscore ~= cfg_inout.zscore
        fprintf('\t%s: non-default value for field ''zscore'' is being used\n',flags.glmfct);
    end
end

% Default number of permutations N_perm (if do_permutation == 1):
switch flags.glmfct
    case 'CI_glmqpoptim_fct'
        
        prior = 'smoothness';
        lambda0 = 5;
        stepsize = 1.5; % Progression step for lambda
        maxiter  = 30; % Number of iterations
        nobreak  = 1;
        minDiffSecondRound = 10;
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'prior',prior); 
        if bAssigned == 0
            if ~strcmp(cfg_inout.prior,prior)
                fprintf('\t%s: non-default value for field ''prior'' is being used\n',flags.glmfct);
            end
        end
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'lambda0',lambda0); 
        if bAssigned == 0
            if ~strcmp(cfg_inout.lambda0,lambda0)
                fprintf('\t%s: non-default value for field ''lambda0'' is being used\n',flags.glmfct);
            end
        end
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'stepsize',stepsize); 
        if bAssigned == 0
            if ~strcmp(cfg_inout.stepsize,stepsize)
                fprintf('\t%s: non-default value for field ''stepsize'' is being used\n',flags.glmfct);
            end
        end
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'maxiter',maxiter); 
        if bAssigned == 0
            if ~strcmp(cfg_inout.maxiter,maxiter)
                fprintf('\t%s: non-default value for field ''maxiter'' is being used\n',flags.glmfct);
            end
        end
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'nobreak',nobreak); 
        if bAssigned == 0
            if ~strcmp(cfg_inout.nobreak,nobreak)
                fprintf('\t%s: non-default value for field ''nobreak'' is being used\n',flags.glmfct);
            end
        end
        
        [cfg_inout,bAssigned] = Ensure_field(cfg_inout,'minDiffSecondRound',minDiffSecondRound);
        if bAssigned == 0
            if ~strcmp(cfg_inout.minDiffSecondRound,minDiffSecondRound)
                fprintf('\t%s: non-default value for field ''minDiffSecondRound'' is being used\n',flags.glmfct);
            end
        end
        
    case 'lassoglm'
        % decimation = 350;
        % earQ = 8;
        % stepfactor = 0.4;
        % 
        % if ~isfield(cfg_inout,'decimation')
        %     cfg_inout.decimation = decimation; % Parametres de calcul du cochleogramme
        % else
        %     if decimation ~= cfg_inout.decimation
        %         fprintf('\t%s: non-default value for field ''decimation'' is being used\n',flags.DimCI);
        %     end
        % end
        % 
        % if ~isfield(cfg_inout,'earQ')
        %     cfg_inout.earQ = earQ; 
        % else
        %     if earQ ~= cfg_inout.earQ
        %         fprintf('\t%s: non-default value for field ''earQ'' is being used\n',flags.DimCI);
        %     end
        % end
        % 
        % if ~isfield(cfg_inout,'stepfactor')
        %     cfg_inout.stepfactor = stepfactor; 
        % else
        %     if stepfactor ~= cfg_inout.stepfactor
        %         fprintf('\t%s: non-default value for field ''stepfactor'' is being used\n',flags.DimCI);
        %     end
        % end
    case 'classic_revcorr'
end