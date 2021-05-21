function cfg_inout = arg_TF_type(cfg_inout,flags,keyvals)
% function cfginout = arg_DimCI(cfg_inout,flags,keyvals)
%
% Enabled so far:
%   'permutation' or 'no_permutation': to perform or not the permutation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% definput.flags.permutation = {'permutation','no_permutation'};
% definput.flags.recreate_validation = {'no_recreate_validation','recreate_validation'};

% Default number of permutations N_perm (if do_permutation == 1):
switch flags.TF_type
    case {'spect','tf'}
        % definput.keyvals.prior = 'smoothness'; % default number of ACI assessments for the permutation test
        % definput.keyvals.lambda0 = 5;
        % definput.keyvals.stepsize = 1.5; % Progression step for lambda
        % definput.keyvals.maxiter  = 30; % Number of iterations
        % definput.keyvals.nobreak  = 1;
        % definput.keyvals.minDiffSecondRound = 10;
        % definput.keyvals.nobreak  = 1;
        
        % Parametres pour analyse spectrogramme       
        if strcmp(flags.TF_type,'tf')
            warning('%s: flag ''tf'' will be soon deprecated, use ''spect'' instead...',upper(filename));
        end
        
        cfg_inout.spect_overlap = keyvals.spect_overlap; % default = 0
        cfg_inout.spect_Nwindow = keyvals.spect_Nwindow; % default = 512
        cfg_inout.spect_NFFT    = keyvals.spect_NFFT;    % default = 512
        cfg_inout.spect_unit    = keyvals.spect_unit;    % default = 'dB'
                
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
        
    case 'noise_logspect'
        cfg_inout.logspect_unit = keyvals.logspect_unit; % default = 'dB'
end

% definput.keyvals.zscore = 1;