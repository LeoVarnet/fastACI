function definput = arg_Script4_Calcul_ACI(definput)
% function definput = arg_Script4_Calcul_ACI(definput)
%
% Enabled so far:
%   'permutation' or 'no_permutation': to perform or not the permutation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.flags.TF_type = {'spect','lyon','noise_logspect', ...
                          'tf'};
definput.flags.glmfct = {'glmfitqp'         ,'lassoglm','lasso','classic_revcorr'};
% Old names:   glmfct = {'CI_glmqpoptim_fct','lassoglm','lasso','classic_revcorr'};
definput.flags.permutation = {'permutation','no_permutation'};
definput.flags.recreate_validation = {'no_recreate_validation','recreate_validation'};
definput.flags.plot = {'plot','no_plot'};

definput.keyvals.apply_SNR  = 0; % old: WithSNR
definput.keyvals.add_signal = 0; % old: WithSignal

% Default number of permutations N_perm (if do_permutation == 1):
definput.keyvals.N_perm = 20; % default number of ACI assessments for the permutation test
definput.keyvals.N_folds = 10; 

definput.keyvals.dir_noise = []; 
definput.keyvals.idx_trialselect = []; 

definput.keyvals.f_limits = [1 10000]; % Hz, arbitrary frequencies to be used as limits 
definput.keyvals.t_limits = [0 1]; 
definput.keyvals.expvar_limits = [];

definput.keyvals.dir_target = [];
definput.keyvals.dir_noise  = [];
definput.keyvals.dir_out    = [];

definput.groups.varnet2013 = {'spect','glmfitqp','f_limits',[0 4050],'t_limits',[0 0.3425]};

definput.keyvals.spect_overlap = 0;    % Parametres de calcul du spectrogramme
definput.keyvals.spect_Nwindow = 512;  %
definput.keyvals.spect_NFFT    = 512;  %
definput.keyvals.spect_unit    = 'dB'; % 'linear'

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