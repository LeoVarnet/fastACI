function cfg_inout = arg_TF_type(cfg_inout,flags)
% function cfginout = arg_DimCI(cfg_inout,flags)
%
% Enabled so far:
%   'permutation' or 'no_permutation': to perform or not the permutation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% definput.flags.permutation = {'permutation','no_permutation'};
% definput.flags.recreate_validation = {'no_recreate_validation','recreate_validation'};


% Default number of permutations N_perm (if do_permutation == 1):
switch flags.TF_type
    case 'tf'
        % definput.keyvals.prior = 'smoothness'; % default number of ACI assessments for the permutation test
        % definput.keyvals.lambda0 = 5;
        % definput.keyvals.stepsize = 1.5; % Progression step for lambda
        % definput.keyvals.maxiter  = 30; % Number of iterations
        % definput.keyvals.nobreak  = 1;
        % definput.keyvals.minDiffSecondRound = 10;
        % definput.keyvals.nobreak  = 1;
        % %definput.groups.adt_dau = {'tau',[0.005 0.050 0.129 0.253 0.500]};
        % Parametres pour analyse spectrogramme       
        
        overlap = 0;    % Parametres de calcul du spectrogramme
        Nwindow = 512;  %
        NFFT    = 512;  %
        
        if ~isfield(cfg_inout,'overlap')
            cfg_inout.overlap = overlap; % Parametres de calcul du cochleogramme
        else
            if overlap ~= cfg_inout.overlap
                fprintf('\t%s: non-default value for field ''overlap'' is being used\n',flags.DimCI);
            end
        end
        
        if ~isfield(cfg_inout,'Nwindow')
            cfg_inout.Nwindow = Nwindow; % Parametres de calcul du cochleogramme
        else
            if Nwindow ~= cfg_inout.Nwindow
                fprintf('\t%s: non-default value for field ''Nwindow'' is being used\n',flags.DimCI);
            end
        end
        
        if ~isfield(cfg_inout,'NFFT')
            cfg_inout.NFFT = NFFT; % Parametres de calcul du cochleogramme
        else
            if NFFT ~= cfg_inout.NFFT
                fprintf('\t%s: non-default value for field ''NFFT'' is being used\n',flags.DimCI);
            end
        end
        
    case 'lyon'
        decimation = 350;
        earQ = 8;
        stepfactor = 0.4;
        
        if ~isfield(cfg_inout,'decimation')
            cfg_inout.decimation = decimation; % Parametres de calcul du cochleogramme
        else
            if decimation ~= cfg_inout.decimation
                fprintf('\t%s: non-default value for field ''decimation'' is being used\n',flags.DimCI);
            end
        end
        
        if ~isfield(cfg_inout,'earQ')
            cfg_inout.earQ = earQ; 
        else
            if earQ ~= cfg_inout.earQ
                fprintf('\t%s: non-default value for field ''earQ'' is being used\n',flags.DimCI);
            end
        end
        
        if ~isfield(cfg_inout,'stepfactor')
            cfg_inout.stepfactor = stepfactor; 
        else
            if stepfactor ~= cfg_inout.stepfactor
                fprintf('\t%s: non-default value for field ''stepfactor'' is being used\n',flags.DimCI);
            end
        end
end

definput.keyvals.zscore = 1;