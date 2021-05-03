function definput = arg_Script4_Calcul_ACI(definput)
% function definput = arg_Script4_Calcul_ACI(definput)
%
% Enabled so far:
%   'permutation' or 'no_permutation': to perform or not the permutation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.flags.DimCI  = {'tf','lyon','noise_logspect'};
definput.flags.glmfct = {'CI_glmqpoptim_fct','lassoglm','classic_revcorr'};

definput.flags.permutation = {'permutation','no_permutation'};
definput.flags.recreate_validation = {'no_recreate_validation','recreate_validation'};

% Default number of permutations N_perm (if do_permutation == 1):
definput.keyvals.N_perm = 20; % default number of ACI assessments for the permutation test

% %definput.groups.adt_dau = {'tau',[0.005 0.050 0.129 0.253 0.500]};