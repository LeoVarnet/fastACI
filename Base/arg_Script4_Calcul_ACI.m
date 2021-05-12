function definput = arg_Script4_Calcul_ACI(definput)
% function definput = arg_Script4_Calcul_ACI(definput)
%
% Enabled so far:
%   'permutation' or 'no_permutation': to perform or not the permutation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.flags.TF_type = {'tf','lyon','noise_logspect'};
definput.flags.glmfct = {'glmfitqp'         ,'lassoglm','classic_revcorr'};
% Old names:   glmfct = {'CI_glmqpoptim_fct','lassoglm','classic_revcorr'};
definput.flags.permutation = {'permutation','no_permutation'};
definput.flags.recreate_validation = {'no_recreate_validation','recreate_validation'};

% Default number of permutations N_perm (if do_permutation == 1):
definput.keyvals.N_perm = 20; % default number of ACI assessments for the permutation test
definput.keyvals.N_folds = 10; 

definput.keyvals.dir_noise = []; 
definput.keyvals.idx_trialselect = []; 

definput.keyvals.f_limits = [1 10000]; % Hz, arbitrary frequencies to be used as limits 
definput.keyvals.t_limits = [0 1]; 
definput.keyvals.expvar_limits = [];

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