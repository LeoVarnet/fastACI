function [y, y_correct, X, U, cfg_ACI] = Script4_Calcul_ACI_preprocess(cfg_ACI, data_passation, Data_matrix)
% function [y, y_correct, X, U, cfg_ACI] = Script4_Calcul_ACI_preprocess(cfg_ACI, data_passation, Data_matrix)
%
% Changing the parameter names:
% New name      Old name        Changed on:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_response_here = data_passation.n_responses(cfg_ACI.IdxTrialsLoad);
n_targets_here  = data_passation.n_targets(cfg_ACI.IdxTrialsLoad);
expvar          = data_passation.expvar(cfg_ACI.IdxTrialsLoad);
    
% Selection of the condition
% cfg_ACI.NameCond = Condition;
cfg_ACI.n_signal_analysis = []; % perform analysis only on trials corresponding to specific signals numbers
cfg_ACI.SNR_analysis      = []; % select a range of SNR
cfg_ACI.n_trials_analysis = []; % select a range of trial numbers

do_permutation = cfg_ACI.flags.do_permutation; % By default the permutation test is 'on'
if do_permutation
    N_perm = cfg_ACI.keyvals.N_perm; 
end 

% switch Analysis_condition % cfg_ACI.NameCond
%     case 'total'
        % nothing to do
%     case 'Al'
%         cfg_ACI.n_signal_analysis = [1 2];
%     case 'Ar'
%         cfg_ACI.n_signal_analysis = [3 4];
%     case 'Alda'
%         cfg_ACI.n_signal_analysis = [1];
%     case 'Alga'
%         cfg_ACI.n_signal_analysis = [2];
%     case 'Arda'
%         cfg_ACI.n_signal_analysis = [3];
%     case 'Arga'
%         cfg_ACI.n_signal_analysis = [4];
%     case 'lowSNR'
%         cfg_ACI.SNR_analysis = [min(expvar) median(expvar)];
%     case 'highSNR'
%         cfg_ACI.SNR_analysis = [median(expvar) max(expvar)];
%     case 'firsttrials'
%         cfg_ACI.n_trials_analysis = [1 N/2];
%     case 'lasttrials'
%         cfg_ACI.n_trials_analysis = [N/2+1 N];
%     otherwise
%         error(['Condition inconnue : ' cfg_ACI.NameCond])
% end

cfg_ACI = set_default_cfg(cfg_ACI, 'SNR_analysis', [min(expvar) max(expvar)], 'n_signal_analysis', 1:cfg_ACI.N_target, 'n_trials_analysis', [1 cfg_ACI.N_TrialsLoad]);

N_TrialsLoad = cfg_ACI.N_TrialsLoad;
select_n_trials = (1:N_TrialsLoad>=cfg_ACI.n_trials_analysis(1) & 1:N_TrialsLoad<=cfg_ACI.n_trials_analysis(2));
 
select_SNR = (expvar>=cfg_ACI.SNR_analysis(1) & expvar<=cfg_ACI.SNR_analysis(2));
select_n_signal = zeros(1,N_TrialsLoad);

for i=cfg_ACI.n_signal_analysis
    select_n_signal = ((n_targets_here==i) | select_n_signal);
end
idx_analysis = find(select_n_trials & select_SNR & select_n_signal);

% % Egalisation du nbr de reponses , si demandé
% if isyes(cfg_ACI.egalisation)
%     error('Not validated yet...')
%     % idx_1=find(n_reponse_here==1);
%     % idx_2=find(n_reponse_here==2);
%     % min_nidx=min(length(idx_1),length(idx_2));
%     % idx_analysis=[idx_1(1:min_nidx),idx_2(1:min_nidx)];
% end

cfg_ACI.N_trials = length(idx_analysis);
y         = double((n_response_here(idx_analysis)==1)');
y_correct = double((cfg_ACI.response_correct_target(n_targets_here(idx_analysis))==1)');
% y_correct indicates the trials where target '1' has been expected to 
%     be chosen or 'target 1 present' trials
    
%%% Start: Testing Alejandro on 27/04/2021
if do_permutation
    cfg_perm = [];
    cfg_perm.N_perm = N_perm;
    for i = 1:N_perm
        cfg_perm.idxs_perm(:,i) = transpose(randperm(cfg_ACI.N_trials));
        cfg_perm.y_perm(:,i)    = y(cfg_perm.idxs_perm(:,i)); % random answers
        n_targets_perm = n_targets_here(cfg_perm.idxs_perm(:,i));
        cfg_perm.y_correct_perm(:,i) = double((cfg_ACI.response_correct_target(n_targets_perm(idx_analysis))==1)');

        n_responses_perm = data_passation.n_responses(idx_analysis);
        n_responses_perm = n_responses_perm(cfg_perm.idxs_perm(:,i));
        cfg_perm.perc_correct_if_perm(i) = sum(n_responses_perm == data_passation.n_response_correct_target(idx_analysis))/length(idx_analysis);
        if i == 1
            cfg_perm.perc_correct = sum(data_passation.is_correct(idx_analysis))/length(idx_analysis);
        end
    end
end
    
% y_correct_perm(:,i), X, U_perm_here

%%% End: Testing Alejandro on 27/04/2021
    
switch cfg_ACI.withX
    case {1,'yes','oui'}
        X = Data_matrix(idx_analysis,:);
    otherwise
        X = [];
end
    
switch cfg_ACI.withU
    case {1,'yes','oui'}
        U = [y_correct ones(cfg_ACI.N_trials, 1)];
        if do_permutation
            for i = 1:N_perm
                cfg_perm.U_perm(:,:,i) = [cfg_perm.y_correct_perm(:,i) ones(cfg_ACI.N_trials, 1)]; %%% Testing Alejandro
            end
        end

    otherwise
        U = [];
        if do_permutation
            for i = 1:N_perm
                cfg_perm.U_perm(:,:,i) = []; %%% Testing Alejandro
            end
        end
end
cfg_ACI.N_trials_true  = sum(y==1);
cfg_ACI.N_trials_false = sum(y==0);

if do_permutation
    cfg_ACI.cfg_perm = cfg_perm;
end
    
% warning('Temporal')
%%% Commented by Leo on 28/04/2021
% if cfg_ACI.zscore
%     X=(X-mean(X(:)))/std(X(:)); % Again a normalisation...
% end
%%% End commented by Leo on 28/04/2021

if cfg_ACI.zscore
    %%%% Added by Leo on 28/04/2021:
    % Here each row of X (along dimension 1, each trial) is 
    %     normalised to have a mean of 0 and a std of 1.
    X = zscore(X,[],1);
    % U = zscore(U,[],1);
end