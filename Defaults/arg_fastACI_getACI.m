function definput = arg_fastACI_getACI(definput)
% function definput = arg_fastACI_getACI(definput)
%
% Enabled so far:
%   'permutation' or 'no_permutation': to perform or not the permutation test
% Old name: arg_Script4_Calcul_ACI.m (changed on 21/05/2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.flags.TF_type = {'spect','lyon','noise_logspect','gammatone','adapt','tf'};
definput.flags.glmfct = {'glmfitqp'         ,'lassoglm','lasso','lassoslow','lassoglmslow','classic_revcorr'};
% Old names:   glmfct = {'CI_glmqpoptim_fct','lassoglm','lasso',                           'classic_revcorr'};
definput.flags.permutation = {'permutation','no_permutation'};
definput.flags.force_dataload = {'no_force_dataload','force_dataload'};
definput.flags.bias   = {'bias'  ,'no_bias'};
definput.flags.recreate_validation = {'no_recreate_validation','recreate_validation'};
definput.flags.plot = {'plot','no_plot'};

definput.keyvals.apply_SNR  = 0; % old: WithSNR
definput.keyvals.add_signal = 0; % old: WithSignal 

definput.keyvals.skip_if_on_disk = 1;
definput.keyvals.Data_matrix = [];
% 
definput.keyvals.cfg_crosspred = {};
definput.keyvals.results_crosspred = {};

% Default number of permutations N_perm (if do_permutation == 1):
definput.keyvals.N_perm = 20; % default number of ACI assessments for the permutation test
definput.keyvals.N_folds = 10; 

definput.keyvals.dir_noise = []; 
definput.keyvals.idx_trialselect = []; 

definput.keyvals.trialtype_analysis = 'total'; 

definput.keyvals.f_limits = [1 10000]; % Hz, arbitrary frequencies to be used as limits 
definput.keyvals.t_limits = [0 1]; 
definput.keyvals.expvar_limits = [];
definput.keyvals.perc = [];

definput.keyvals.dir_target = [];
definput.keyvals.dir_noise  = [];
definput.keyvals.dir_out    = [];

definput.groups.varnet2013 = {'spect','glmfitqp','f_limits',[0 4050],'t_limits',[0 0.3425]};

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

definput.keyvals.lambda = [];

% % Managing some default values: time and frequency limits
% if flags.do_tf
%     opts_ACI = Ensure_field(opts_ACI,'freq_analysis',[1 10000]);
%     opts_ACI = Ensure_field(opts_ACI,'time_analysis',[0 0.6]);
% end
% 
% if flags.do_lyon
%     opts_ACI = Ensure_field(opts_ACI,'freq_analysis',[1  8000]);
%     opts_ACI = Ensure_field(opts_ACI,'time_analysis',[0 1]); 
% end