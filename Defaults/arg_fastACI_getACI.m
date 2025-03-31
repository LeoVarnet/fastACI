function definput = arg_fastACI_getACI(definput)
% function definput = arg_fastACI_getACI(definput)
%
% Enabled so far:
%   'permutation' or 'no_permutation': to perform or not the permutation test
% Old name: arg_Script4_Calcul_ACI.m (changed on 21/05/2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.flags.TF_type = {'gammatone','spect','lyon','noise_logspect','adapt','tf'}; 
definput.flags.glmfct = {'correlation','glm','glm_L1_GB','glm_L2','lm_L1_GB','weighted_sum',...
    'glmfitqp','lassoglm_original','lassoglm','lasso','lassoslow','l1glm','l1lm','classic_revcorr',...
};
% new names (26/09/2024): 'l1glm' is now 'glm_L1_GB'; 'glmfitqp' is now
% 'glm_L2', 'classic_revcorr' is now 'correlation'
definput.flags.permutation = {'no_permutation','permutation'};
definput.flags.force_dataload = {'no_force_dataload','force_dataload'};
definput.flags.bias   = {'bias'  ,'no_bias'};
definput.flags.recreate_validation = {'no_recreate_validation','recreate_validation'};
definput.flags.plot = {'plot','no_plot'};

%%% Parameters for lasso fitting:
%       The following parameters are only relevant if glmfct is one of the 
%       lasso functions:
definput.keyvals.pyramid_script = 'imgaussfilt'; % script used to reduce/expand the Gaussian pyramid
% definput.keyvals.pyramid_script = 'imresize';
definput.keyvals.pyramid_shape = 0; % default, other possible values: -1
definput.keyvals.pyramid_padding = 'replicate'; % 'Padding',0
%%% End: parameters for lasso fitting

definput.keyvals.apply_SNR  = 0; % old: WithSNR
definput.keyvals.add_signal = 0; % old: WithSignal 

definput.keyvals.skip_if_on_disk = 1;
definput.keyvals.consistency_check = 1;
definput.keyvals.Data_matrix = [];

definput.keyvals.zscore = 1;

definput.keyvals.ACI_crosspred = [];
definput.keyvals.suffix_crosspred = ''; % if empty (default) generates Crosspred.mat
                                        % if non empty generates Crosspred-noise.mat (if suffix_crosspred == 'noise') 
% Default number of permutations N_perm (if do_permutation == 1):
definput.keyvals.N_perm = 20; % default number of ACI assessments for the permutation test
definput.keyvals.N_folds = 10; 

definput.keyvals.dir_noise = []; 
definput.keyvals.idx_trialselect = []; 

definput.keyvals.trialtype_analysis = 'total';  %'total', 'incorrect', 'correct', 't1', 't2'

definput.keyvals.f_limits = [1 10000]; % Hz, arbitrary frequencies to be used as limits 
definput.keyvals.t_limits = [0 1]; 

%%% Used in fastACI_getACI_preprocess.m:
definput.keyvals.expvar_limits = [];
definput.keyvals.expvar_after_reversal = 0; % 0 = all trials are used; 4 = all trials 
                                    % after reversal 4 are kept. New option as of 24/01/2022
definput.keyvals.perc = [];
%%%
definput.keyvals.dir_target = [];
definput.keyvals.dir_noise  = [];
definput.keyvals.dir_out    = [];

definput.groups.varnet2013 = {'spect','glmfitqp', ...
    'pyramid_script','imresize', ... % this keyval is irrelevant for processing
    'f_limits',[0 4050],'t_limits',[0 0.3425]};

definput.groups.glmfitqp         = {'pyramid_script',[]};
definput.groups.glm              = {'pyramid_script',[]};
definput.groups.classic_revcorr  = {'pyramid_script',[]};

%%% Only used in 'spect' is used:
definput.keyvals.spect_overlap = 0;    % Parametres de calcul du spectrogramme
definput.keyvals.spect_Nwindow = 512;  %
definput.keyvals.spect_NFFT    = 512;  %
definput.keyvals.spect_unit    = 'dB'; % 'linear'

%%% Only used in 'noise_logspect' is used:
definput.keyvals.logspect_unit = 'dB'; % 'linear'

%%% Only used in 'gammatone'
definput.keyvals.bwmul    = .5; % recommended either 0.5 or 1 ERB
definput.keyvals.binwidth = 10e-3; % 10 ms, time resolution for ACI

%%% Only used in 'l1glm'
definput.keyvals.lambda = [];

%%% Only used in 'glmfitqp'
definput.keyvals.lambda0 = 5;
definput.keyvals.precision = [];
definput.keyvals.stepsize = 1.5;

definput.keyvals.script_dataload = ''; % fastACI_getACI_dataload